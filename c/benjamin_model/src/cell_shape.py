# cell_shape.py
# Genreate a simple shape of a cell on its own lattice

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("nx", type=int, help="Width of lattice")
parser.add_argument("ny", type=int, help="Height of lattice")
parser.add_argument("phi0", type=float, help="Maximum value of phi")
parser.add_argument("outfile", help="Name of output file")
subparsers = parser.add_subparsers(title="shapes",dest="shape")
parser_circle = subparsers.add_parser("circle")
parser_square = subparsers.add_parser("square")
parser_circle.add_argument("radius", type=float, help="Radius of the cell")
parser_square.add_argument("length", type=float, help="Length of the cell")
args = parser.parse_args()

nx = args.nx
ny = args.ny
phi0 = args.phi0
outfile = args.outfile

def circle(x, y, r):
    x0 = x/2
    y0 = y/2
    phi = [[phi0 if ((i-x0)**2+(j-y0)**2)**0.5 < r else 0.0 \
            for j in xrange(ny)] for i in xrange(nx)]
    return phi

def square(x, y, l):
    x0 = x/2
    y0 = y/2
    lh = l/2
    phi = [[phi0 if (abs(i-x0)<lh and abs(j-y0)<lh) else 0.0 \
            for j in xrange(ny)] for i in xrange(nx)]
    return phi

def zeros(x, y):
    phi = [[0.0 for j in xrange(ny)] for i in xrange(nx)]
    return phi

if args.shape == "circle":
    phi = circle(nx, ny, args.radius)
elif args.shape == "square":
    phi = square(nx, ny, args.length)
else:
    phi = zeros(nx, ny)

with open(outfile, "w") as writer:
    for i in xrange(nx):
        for j in xrange(ny):
            writer.write("{:d} {:d} {:.5f}\n".format(i, j, phi[i][j]))
