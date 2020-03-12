# triangle.py
# Genreate the coordinates of the cells' centre of mass which are
# arranged in a triangular lattice (with noise)

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("nx", type=int, help="number of cells in each row")
parser.add_argument("ny", type=int, help="number of cells in each col")
parser.add_argument("radius", type=float, help="radius of each cell")
parser.add_argument("cutoff", type=float, help="cutoff radius")
parser.add_argument("seed", type=int, help="seed for random generator")
parser.add_argument("outfile", help="name of output file")
args = parser.parse_args()

nx = args.nx
ny = args.ny
r = args.radius
rc = args.cutoff
seed = args.seed
outfile = args.outfile

rng = np.random.RandomState(seed)
pi = np.pi

ncells = nx*ny

lx = int(2*r*nx)
ly = int(np.ceil(3**0.5*r*ny))
size = int((lx*ly)**0.5)
lx = size
ly = size
dx = lx/float(nx)
dy = ly/float(ny)

print "lx = {:d} ly = {:d}".format(lx, ly)

with open (outfile, "w") as writer:
    for i in xrange(0,ncells):
        x = dx*0.5*((i/nx)%2) + dx*(i%nx)
        y = dy*(i/nx)
        if (rc > 0.0):
            rx = rng.normal(0.0,rc)
            ry = rng.normal(0.0,rc)
            x = x+rx
            y = y+ry
        writer.write("{:d} {:.5f} {:.5f}\n".format(i,x,y))
