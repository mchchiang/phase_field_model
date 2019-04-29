import sys
from itertools import izip
import math

args = sys.argv

if (len(args) < 7):
    print "Usage: AverageMultiFiles.py [ref_col] [avg_col] [err_col] " \
        "[startpt] [output_file] [data_files]"
    sys.exit(1)

args.pop(0) # Ignore self
ref_col = int(args.pop(0))
avg_col = int(args.pop(0))
err_col = int(args.pop(0))
startpt = int(args.pop(0))
output_file = args.pop(0)

files = [open(i, "r") for i in args]
writer = open(output_file, "w")

for rows in izip(*files):
    ref = 0.0
    avg = 0.0
    avgSq = 0.0
    error = 0.0
    hasData = False

    for line in rows:
        if (not line.startswith("#")):
            data = line.strip().split()
            if (data == []): #ignore any lines start with \n
                break
            
            if (ref_col != -1):
                ref = float(data[ref_col])

            if (ref_col == -1 or ref >= startpt):
                value = float(data[avg_col])
                avg += value
                avgSq += value*value

                if (err_col != -1):
                    sigma = float(data[err_col])
                    error += sigma*sigma

                hasData = True

    if (hasData == True):
        n = float(len(rows))
        avg /= n
        avgSq /= n

        # use un-biased estimate of variance
        if (n > 1):
            var = n / (n-1) * (avgSq - avg*avg)
            sigma = math.sqrt(var) 
        else:
            var = 0
            sigma = 0
        
        if (err_col != -1):
            error = math.sqrt(error) / n
        else:
            error = sigma / math.sqrt(n)
        
        if (ref_col == -1):
            output = "%.5f %.5f %.5f\n" % (avg, sigma,  error)
        else:
            output = "%.5f %.5f %.5f %.5f\n" % (ref, avg, sigma,  error)
        writer.write(output)

for f in files:
    f.close()
writer.close()
    
        
    
