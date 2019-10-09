# KS2SampleTest.py
#
# Perform the Kolmogorov-Smmirnov two sample test
#

import sys
import scipy
import scipy.stats
import math

def read_array(filename, array, index_col, value_col):
    with open(filename, 'r') as f:
        index = 0
        for line in f:
            if (line.startswith("#")): continue
            data = line.strip().split()
            if (data == []): continue
            value = float(data[value_col])
            if (index_col != -1):
                index = int(data[index_col])                
                array[index] = value
            else:
                array[index] = value
                index += 1 

args = sys.argv
args.pop(0) # Ignore self

if (len(args) != 7):
    print "Usage: [size1] [size2] [index_col] [value_col] " \
        "[array1_file] [array2_file] [out_file]"
    sys.exit(1)

size1 = int(args.pop(0))
size2 = int(args.pop(0))
index_col = int(args.pop(0))
value_col = int(args.pop(0))
array1_file = args.pop(0)
array2_file = args.pop(0)
out_file = args.pop(0)

array1 = [0.0 for i in xrange(size1)]
array2 = [0.0 for i in xrange(size2)]

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

# Compute KS 2-Sample Test
result = scipy.stats.ks_2samp(array1, array2)

with open(out_file, 'w') as writer:
    writer.write("{:.5g} {:.5g}\n".format(*result))


