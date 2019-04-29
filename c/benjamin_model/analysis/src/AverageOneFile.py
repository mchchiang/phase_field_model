import sys
import math

args = sys.argv

if (len(args) != 7):
    print "Usage: Average [avg_col] [err_col] [start_row] [end_row] " \
        "[input_file] [output_file]"
    sys.exit(1)

args.pop(0) # Ignore self
avg_col = int(args.pop(0))
err_col = int(args.pop(0))
start_row = int(args.pop(0))
end_row = int(args.pop(0))
input_file = args.pop(0)
output_file = args.pop(0)

files = [open(i, "r") for i in args]
writer = open(output_file, "w")

avg = 0.0
avgSq = 0.0
error = 0.0
n = 0

with open(input_file, 'r') as f:
    for row, line in enumerate(f):
        if (row > end_row):
            break
        if (row < start_row):
            continue
        if (line.startswith("#")): # Ignore lines with comments
            continue
        data = line.strip().split()
        
        if (data == []): # Ignore any lines start with \n
            continue
        
        value = float(data[avg_col])

        if (math.isnan(value) or math.isinf(value)):
            continue

        avg += value
        avgSq += value*value
        
        if (err_col != -1):
            sigma = float(data[err_col])
            error += sigma*sigma

        n += 1
    
    n = float(n)
    avg /= max(n,1)
    avgSq /= max(n,1)

    # Use un-biased estimate of variance
    if (n > 1):
        var = n / (n-1) * (avgSq - avg*avg)
        sigma = math.sqrt(var) 
    else:
        var = 0.0
        sigma = 0.0
    
    if (err_col != -1):
        error = math.sqrt(error) / max(n,1)
    else:
        error = sigma / math.sqrt(n) if n > 0.0 else 0.0
        
    output = "{:.5f} {:.5f} {:.5f}\n".format(avg, sigma,  error)
    writer.write(output)

writer.close()
