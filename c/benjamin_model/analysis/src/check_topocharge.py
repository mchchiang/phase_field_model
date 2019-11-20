# check_topocharge.py

import sys

args = sys.argv

if (len(args) != 3):
    print "usage: check_topocharge.py npoints in_file"
    sys.exit(1)

npoints = int(args.pop(1))
in_file = args.pop(1)

reader = open(in_file, 'r')
while True:
    for i in xrange(2):
        line = reader.readline()
    if (not line): break
    data = line.split()
    time = int(data[1])
    
    charge = 0
    for i in xrange(npoints):
        line = reader.readline()
        data = line.split()
        charge += (len(data)-6)
    if (charge != 0):
        print "Topological charge not conserved at time = {:d}, s = {:d}".format(time, charge)

reader.close()
