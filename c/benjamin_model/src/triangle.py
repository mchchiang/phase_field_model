# triangle.py
# Genreate the coordinates of the cells' centre of mass which are
# arranged in a triangular lattice

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("nx", type=int, help="number of cells in each row")
parser.add_argument("ny", type=int, help="number of cells in each col")
parser.add_argument("radius", type=float, help="radius of each cell")
parser.add_argument("outfile", help="name of output file")
args = parser.parse_args()

nx = args.nx
ny = args.ny
r = args.radius
outfile = args.outfile

ncells = nx*ny

lx = int(2*r*nx)
ly = int(3**0.5*r*ny)
#lx = int(3**0.5*r*nx)
#ly = int(2*r*ny)

print "lx = {:d} ly = {:d}".format(lx, ly)

with open (outfile, "w") as writer:
    for i in xrange(0,ncells):
        x = r*((i/nx)%2) + 2*r*(i%nx)
        y = 3**0.5*r*(i/nx)
#        x = 3**0.5*r*(i/nx)
#        y = r*((i/nx)%2) + 2*r*(i%nx)
        writer.write("{:d} {:.5f} {:.5f}\n".format(i,x,y))
