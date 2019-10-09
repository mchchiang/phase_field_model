# PearsonCorrelation.py
#
# Compute the Pearson correlation coefficient of two arrays
#

import sys
import scipy
import scipy.stats
import math

def read_array(filename, array, index_col, value_col):
    with open(filename, 'r') as f:
        for line in f:
            if (line.startswith("#")): continue
            data = line.strip().split()
            if (data == []): continue
            index = int(data[index_col])
            value = float(data[value_col])
            array[index] = value


args = sys.argv
args.pop(0) # Ignore self

if (len(args) < 3):
    print "Usage: [size] [index_col] [value_col] " \
        "[array1_file] [array2_file] [out_file]"
    sys.exit(1)

size = int(args.pop(0))
index_col = int(args.pop(0))
value_col = int(args.pop(0))
array1_file = args.pop(0)
array2_file = args.pop(0)
out_file = args.pop(0)

array1 = [0.0 for i in xrange(size)]
array2 = [0.0 for i in xrange(size)]

# Read in arrays 
read_array(array1_file, array1, index_col, value_col)
read_array(array2_file, array2, index_col, value_col)

# Remove nan or inf entries
arrays = [[a,b] for a,b in zip(array1,array2) if not 
          (math.isnan(a) or math.isinf(a) or math.isnan(b) or math.isinf(b))]
array1 = [x[0] for x in arrays]
array2 = [x[1] for x in arrays]

array1 = scipy.array(array1)
array2 = scipy.array(array2)

# Compute correlation
corr = scipy.stats.pearsonr(array1, array2)

with open(out_file, 'w') as writer:
    writer.write("{:.5g} {:.5g}\n".format(*corr))
