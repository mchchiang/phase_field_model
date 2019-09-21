# overlap.py

import sys
import numpy as np
from skimage import measure
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

# This code must be run with python3!
if (sys.version_info < (3, 5)):
    print("This code must be run with Python version 3.5 or higher")

args = sys.argv
if (len(args) != 4):
    print("usage: overlap.py time fileroot pos_file")
    sys.exit(1)

npoints = 100
clx = 41
cly = 41
lx = 160
ly = 138
time = int(args.pop(1))
fileroot = args.pop(1)
pos_file = args.pop(1)
xbuff = 0.2*lx
ybuff = 0.2*ly

polygons = []
points = []
pos = []

def iwrap(x):
    global lx
    remainder = x % lx
    if (remainder >= 0):
        return remainder
    return lx + remainder

def jwrap(y):
    global ly
    remainder = y % ly
    if (remainder >= 0):
        return remainder
    return ly + remainder

def add_polygon(cxcm, cycm, xcm, ycm, dx, dy, contour):
    global polygons, lx, ly, xbuff, ybuff
    xcm += dx
    ycm += dy
    if (xcm < lx+xbuff and xcm > -xbuff and ycm < ly+ybuff and ycm > -ybuff):
        poly = np.copy(contour)
        poly[:,0] += (xcm-cxcm) + 0.5
        poly[:,1] += (ycm-cycm) + 0.5
        points.append((xcm,ycm))
        polygons.append(Polygon(poly))

periodic_loc = [(lx,-ly),(lx,0),(lx,ly),(0,-ly),(0,0),(0,ly),(-lx,-ly),
                (-lx,0),(-lx,ly)]

# Read position data
nlines = npoints + 2
reader = open(pos_file, 'r')
while True:
    # Read header section (including time info)
    for i in range(2):
        line = reader.readline()
    if (not line): break
    data = line.split()
    t = int(data[1])
    if (t > time):
        break
    elif (t < time):
        for i in range(npoints):
            line = reader.readline()
    else:
        for i in range(npoints):
            line = reader.readline()
            data = line.split()
            # Read wrapped position data
            x = float(data[0])
            y = float(data[1])
            pos.append((x,y))
reader.close()
    
for n in range(npoints):
    filename = fileroot + ("cell_{:d}.dat.{:d}".format(n,time))
    print("Reading data from cell {:d}".format(n))
    with open(filename,'r') as reader:
        # For computing the local field centre of mass
        xavg = 0.0
        yavg = 0.0
        mass = 0.0
        local_field = np.zeros((clx,cly))
        for l, line in enumerate(reader):
            data = line.split()
            if (len(data) != 3): continue
            i = int(data[0])
            j = int(data[1])
            phi = float(data[2])
            local_field[i,j] = phi
            xavg += (phi*(i+0.5)) # Uset the centre of a lattice element
            yavg += (phi*(j+0.5))
            mass += phi
        if (mass > 0.0):
            xavg /= mass
            yavg /= mass
        else:
            xavg = 0.0
            yavg = 0.0
        contours = measure.find_contours(local_field, 1.0)
        contour = contours[0]
        for pt in periodic_loc:
            add_polygon(xavg, yavg, pos[n][0], pos[n][1], pt[0], pt[1], 
                        contour)

fig, ax = plt.subplots()
ax.set_xlim([0,lx])
ax.set_ylim([0,ly])

# Draw borders but no axes and ticks
plt.tick_params(axis="both", which="both", bottom=False, top=False,
                labelbottom=False, right=False, left=False, labelleft=False)

# Set plot margins
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

# For drawing each cell's interior
patches1 = PatchCollection(polygons, linewidth=0)
plt_polygons1 = ax.add_collection(patches1)
plt_polygons1.set_facecolor((1.0,0.5,0.0,0.7)) # Orange

# For drawing each cell's boundary
patches2 = PatchCollection(polygons, linewidth=1)
plt_polygons2 = ax.add_collection(patches2)
plt_polygons2.set_facecolor((0.0,0.0,0.0,0.0))
plt_polygons2.set_edgecolor("black")

# For drawing each cell's centre of mass
points = np.array(points)
plt_pts, = ax.plot(points[:,0], points[:,1], '.', markersize=5, color="black")

plt.show()
