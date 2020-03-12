# neigh_delaunay.py
# A script to find the nearest neighbours of a cell using Delaunay 
# triangulation. 

import sys

# This code must be run with python3!
if (sys.version_info < (3, 5)):
    print("This code mus tbe run with Python version 3.5 or higher")
    sys.exit(1)

import numpy as np
from scipy.spatial import Delaunay

args = sys.argv
if (len(args) != 11):
    print("usage: neigh_delaunay.py npoints lx ly xbuff ybuff tstart tend tinc pos_file out_file")
    sys.exit(1)

npoints = int(args.pop(1))
lx = int(args.pop(1))
ly = int(args.pop(1))
xbuff = float(args.pop(1))
ybuff = float(args.pop(1))
tstart = int(args.pop(1))
tend = int(args.pop(1))
tinc = int(args.pop(1))
pos_file = args.pop(1)
out_file = args.pop(1)
xbuff *= lx
ybuff *= ly
nframes = int((tend-tstart)//tinc+1)

periodic_loc = [(lx,-ly),(lx,0),(lx,ly),(0,-ly),(0,ly),(-lx,-ly),
                (-lx,0),(-lx,ly)]

# Read position data
nlines = npoints+2
reader = open(pos_file, 'r')
writer = open(out_file, 'w')

while True:
    # Read header section (including time info)
    for i in range(2):
        line = reader.readline()
    if (not line): break
    data = line.split()
    time = int(data[1])
    if (time > tend):
        break
    elif (time < tstart or (time-tstart)%tinc != 0):
        for i in range(npoints):
            line = reader.readline()
    else:
        pos = []
        index_map = [n for n in range(npoints)]
        frame = int((time-tstart)//tinc)
        for n in range(npoints):
            line = reader.readline()
            data = line.split()
            x = float(data[0])
            y = float(data[1])
            if (x < lx+xbuff and x > -xbuff and y < ly+ybuff and y > -ybuff):
                pos.append((x,y))
        for n in range(npoints):
            x,y = pos[n]
            for pt in periodic_loc:
                xp = x+pt[0]
                yp = y+pt[1]
                if (xp < lx+xbuff and xp > -xbuff and yp < ly+ybuff and 
                    yp > -ybuff):
                    pos.append((xp,yp))
                    index_map.append(n)
        pos = np.array(pos)
        
        # Perform delaunay triangulation
        tri = Delaunay(pos)

        # Create neighbour lists
        neighbours = [set() for i in range(len(pos))]
        for s in tri.simplices:
            for i in s:
                for j in s:
                    if (i != j):
                        neighbours[i].add(index_map[j])
        
        # Output result
        writer.write("Cells: {:d}\n".format(npoints))
        writer.write("Timestep: {:d}\n".format(time))
        for i in range(npoints):
            for n in sorted(neighbours[i]):
                writer.write("{:d} ".format(n))
            writer.write("\n")

reader.close()
writer.close()
