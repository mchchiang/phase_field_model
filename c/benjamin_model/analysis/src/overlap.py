# overlap.py

import sys

# This code must be run with python3!
if (sys.version_info < (3, 5)):
    print("This code must be run with Python version 3.5 or higher")
    sys.exit(1)

import numpy as np
from skimage import measure
import matplotlib as mpl
import matplotlib.cm as mplcm
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

args = sys.argv
if (len(args) != 4):
    print("usage: overlap.py time fileroot pos_file")
    sys.exit(1)

npoints = 100
clx = 41
cly = 41
lx = 160
ly = 138
time = int(args.pop(1))
fileroot = args.pop(1)
pos_file = args.pop(1)
xbuff = 0.2*lx
ybuff = 0.2*ly

polygons = []
points = []
cm = []
local_cm = []
local_field = np.zeros((npoints,clx,cly))
olap_field = np.zeros((npoints,clx,cly))
olap = []
index_map = []

def sgn(val):
    return (0.0 < val) - (val < 0.0)

def iwrap(x):
    global lx
    remainder = x % lx
    if (remainder >= 0):
        return remainder
    return lx + remainder

def jwrap(y):
    global ly
    remainder = y % ly
    if (remainder >= 0):
        return remainder
    return ly + remainder

def idiff(x1, x2):
    global lx
    dx1 = x1-x2
    dx2 = -sgn(dx1)*(lx-abs(dx1))
    if (abs(dx1) < abs(dx2)):
        return dx1
    else:
        return dx2

def jdiff(y1, y2):
    global ly
    dy1 = y1-y2
    dy2 = -sgn(dy1)*(ly-abs(dy1))
    if (abs(dy1) < abs(dy2)):
        return dy1
    else:
        return dy2

def add_polygon(index, cxcm, cycm, xcm, ycm, dx, dy, contour):
    global polygons, lx, ly, xbuff, ybuff
    xcm += dx
    ycm += dy
    if (xcm < lx+xbuff and xcm > -xbuff and ycm < ly+ybuff and ycm > -ybuff):
        poly = np.copy(contour)
        poly[:,0] += (xcm-cxcm) + 0.5
        poly[:,1] += (ycm-cycm) + 0.5
        points.append((xcm,ycm))
        polygons.append(Polygon(poly))
        index_map.append(index)

periodic_loc = [(lx,-ly),(lx,0),(lx,ly),(0,-ly),(0,0),(0,ly),(-lx,-ly),
                (-lx,0),(-lx,ly)]

# Read position data
nlines = npoints + 2
reader = open(pos_file, 'r')
while True:
    # Read header section (including time info)
    for i in range(2):
        line = reader.readline()
    if (not line): break
    data = line.split()
    t = int(data[1])
    if (t > time):
        break
    elif (t < time):
        for i in range(npoints):
            line = reader.readline()
    else:
        for i in range(npoints):
            line = reader.readline()
            data = line.split()
            # Read wrapped position data
            x = float(data[0])
            y = float(data[1])
            cm.append((x,y))
reader.close()
    
for n in range(npoints):
    filename = fileroot + ("cell_{:d}.dat.{:d}".format(n,time))
    with open(filename,'r') as reader:
        # For computing the local field centre of mass
        xavg = 0.0
        yavg = 0.0
        mass = 0.0
        for l, line in enumerate(reader):
            data = line.split()
            if (len(data) != 3): continue
            i = int(data[0])
            j = int(data[1])
            phi = float(data[2])
            local_field[n,i,j] = phi
            xavg += (phi*(i+0.5)) # Uset the centre of a lattice element
            yavg += (phi*(j+0.5))
            mass += phi
        if (mass > 0.0):
            xavg /= mass
            yavg /= mass
        else:
            xavg = 0.0
            yavg = 0.0
        local_cm.append((xavg,yavg))
        contours = measure.find_contours(local_field[n], 1.0)
        contour = contours[0]
        for pt in periodic_loc:
            add_polygon(n, xavg, yavg, cm[n][0], cm[n][1], pt[0], pt[1], 
                        contour)

# Compute degree of overlap
for m in range(npoints):
    olap_avg = 0.0
    area = 0.0
    x0m = iwrap(cm[m][0]-local_cm[m][0])
    y0m = jwrap(cm[m][1]-local_cm[m][1])
    for i in range(clx):
        for j in range(cly):
            if (local_field[m,i,j] < 1.0): continue
            area += 1.0
    for n in range(npoints):
        if (m == n): continue
        x0n = iwrap(cm[n][0]-local_cm[n][0])
        y0n = jwrap(cm[n][1]-local_cm[n][1])
        dx0mn = idiff(x0m,x0n)
        dy0mn = jdiff(y0m,y0n)
        if (abs(dx0mn) > lx or abs(dy0mn) > ly): continue
        for i in range(clx):
            xn = int(round(dx0mn+i))
            if (xn < 0 or xn >= clx): continue
            for j in range(cly):
                phim = local_field[m,i,j]
                if (phim < 1.0): continue
                yn = int(round(dy0mn+j))
                if (yn < 0 or yn >= cly): continue
                phin = local_field[n,xn,yn]
                if (phin < 1.0): continue
                phim = 0.0 if phim < 1.0 else 1.0
                phin = 0.0 if phin < 1.0 else 1.0
                prod = phim*phin
                olap_avg += prod
                olap_field[m,i,j] += prod
    olap_avg /= area
    olap.append(olap_avg)

# Plot settings
data_min = 0.59
data_max = 0.65
#data_min = min(olap)
#data_max = max(olap)
norm = mpl.colors.Normalize(vmin=data_min, vmax=data_max, clip=True)
mapper = mplcm.ScalarMappable(norm=norm, cmap=mplcm.RdYlBu_r)
mapper.set_array([])

fig, ax = plt.subplots()
"""
plt.imshow(olap_field[14])
plt.show()
"""
ax.set_xlim([0,lx])
ax.set_ylim([0,ly])

cbar = plt.colorbar(mapper)

# Draw borders but no axes and ticks
plt.tick_params(axis="both", which="both", bottom=False, top=False,
                labelbottom=False, right=False, left=False, labelleft=False)

# Set plot margins
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

# For drawing each cell's interior
colours = []
for pt in range(len(index_map)):
    colours.append(mapper.to_rgba(olap[index_map[pt]],alpha=0.8))
patches1 = PatchCollection(polygons, linewidth=0)
plt_polygons1 = ax.add_collection(patches1)
plt_polygons1.set_facecolor(colours)

# For drawing each cell's boundary
patches2 = PatchCollection(polygons, linewidth=1)
plt_polygons2 = ax.add_collection(patches2)
plt_polygons2.set_facecolor((0.0,0.0,0.0,0.0))
plt_polygons2.set_edgecolor("black")

# For drawing each cell's centre of mass
points = np.array(points)
plt_pts, = ax.plot(points[:,0], points[:,1], '.', markersize=5, color="black")

plt.show()

