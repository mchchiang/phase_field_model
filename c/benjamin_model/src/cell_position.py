# cell_position.py
# Genreate the coordinates of the cells' centre of mass, which are arranged in
# a triangular lattice, and also their velocities

import argparse
import random
import math

parser = argparse.ArgumentParser()
parser.add_argument("nx", type=int, help="Number of cells in each row")
parser.add_argument("ny", type=int, help="Number of cells in each col")
parser.add_argument("radius", type=float, help="Radius of each cell")
parser.add_argument("v", type=float, help="Speed of each cell")
parser.add_argument("seed", type=int, help="Seed for random generator")
parser.add_argument("outfile", help="name of output file")
args = parser.parse_args()

nx = args.nx
ny = args.ny
r = args.radius
v = args.v
seed = args.seed
outfile = args.outfile

ncells = nx*ny

lx = int(2*r*nx)
ly = int(3**0.5*r*ny)
pi = math.pi
random.seed(seed)

# Output how big is the actual system
print "lx = {:d} ly = {:d}".format(lx, ly)

with open (outfile, "w") as writer:
    for i in xrange(0,ncells):
        x = r*((i/nx)%2) + 2*r*(i%nx)
        y = 3**0.5*r*(i/nx)
        # Generate a random 2D direction [0,2*pi)
        theta = random.random()*2.0*pi
        vx = v*math.cos(theta)
        vy = v*math.sin(theta)
        writer.write("{:d} {:.5f} {:.5f} {:.5f} {:.5}\n".format(i,x,y,vx,vy))
