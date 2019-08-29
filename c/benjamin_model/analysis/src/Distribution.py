# Distribution.py
#
# Return the probability density distribution of a set of data
#

import sys
import math

def read_array(filename, array, value_col):
    index = 0
    with open(filename, 'r') as f:
        for line in f:
            # Ignore lines starting with '#' or have no data
            if (line.startswith("#")): continue
            data = line.strip().split()
            if (data == []): continue
            value = float(data[value_col])
            # Ignore any nan and inf values
            if (math.isnan(value) or math.isinf(value)): continue
            array.append(value)

            
args = sys.argv
args.pop(0) # Ignore self

if (len(args) < 6):
    print "Usage: [value_col] [min] [max] [bin_size] " \
        "[data_file] [out_file]"
    sys.exit(1)

value_col = int(args.pop(0))
min_val = float(args.pop(0))
max_val = float(args.pop(0))
bin_size = float(args.pop(0))
data_file = args.pop(0)
out_file = args.pop(0)

samples = []

# Read in arrays 
read_array(data_file, samples, value_col)

# Create distribution
nbins = int(math.ceil((max_val-min_val)/bin_size))
distrb = [0.0 for i in xrange(nbins)]
cdf = [0.0 for i in xrange(nbins)]
count = 0

for val in samples:
    if (val < min_val or val > max_val): continue # Skip out of range values
    bin_index = int(math.floor((val-min_val)/bin_size))
    distrb[bin_index] += 1.0
    count += 1

# Normalise distribution
distrb = map(lambda x : x / float(count) if count != 0 else 0.0 , distrb)

# Compute cumulative distribution
count = 0.0
for i,v in enumerate(distrb):
    count += v
    cdf[i] = count

# Output distribution
with open(out_file, 'w') as writer:
    for i, val in enumerate(zip(distrb,cdf)):
        left = i*bin_size+min_val
        mid = (i+0.5)*bin_size+min_val
        right = (i+1.0)*bin_size+min_val
        writer.write("{:.5g} {:.5g} {:.5g} {:.5g} {:.5g}\n".format(
            left, mid, right, val[0], val[1]))


