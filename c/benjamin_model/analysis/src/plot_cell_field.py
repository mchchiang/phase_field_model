# plot_cell_field.py

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
import matplotlib.animation as animation
from matplotlib import rc
from matplotlib.animation import FuncAnimation
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

args = sys.argv
if (len(args) != 16):
    print("usage: plot_cell_field.py npoints lx ly clx cly Dr dt tstart tend tinc make_movie print_to_screen fileroot pos_file out_file")
    sys.exit(1)

npoints = int(args.pop(1))
lx = int(args.pop(1))
ly = int(args.pop(1))
clx = int(args.pop(1))
cly = int(args.pop(1))
Dr = float(args.pop(1))
dt = float(args.pop(1))
tstart = int(args.pop(1))
tend = int(args.pop(1))
tinc = int(args.pop(1))
make_movie = bool(int(args.pop(1)))
print_to_screen = bool(int(args.pop(1)))
fileroot = args.pop(1)
pos_file = args.pop(1)
out_file = args.pop(1)
xbuff = 0.2*lx
ybuff = 0.2*ly
nframes = (tend-tstart)//tinc+1

use_label=0 # 1

if (not make_movie):
    tend = tstart

# Data arrays
polygons = [[] for i in range(nframes)]
points = [[] for i in range(nframes)]
pos = [[(0.,0.) for j in range(npoints)] for i in range(nframes)]
index_map = [[] for i in range(nframes)]
time_map = [i*tinc+tstart for i in range(nframes)]

# Useful functions for plotting cell fields
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

def add_polygon(index, cxcm, cycm, xcm, ycm, dx, dy, contour, frame):
    global points, polygons, index_map, lx, ly, xbuff, ybuff
    xcm += dx
    ycm += dy
    if (xcm < lx+xbuff and xcm > -xbuff and ycm < ly+ybuff and ycm > -ybuff):
        poly = np.copy(contour)
        poly[:,0] += (xcm-cxcm) + 0.5
        poly[:,1] += (ycm-cycm) + 0.5
        points[frame].append((xcm,ycm))
        polygons[frame].append(Polygon(poly))
        index_map[frame].append(index)

periodic_loc = [(lx,-ly),(lx,0),(lx,ly),(0,-ly),(0,0),(0,ly),(-lx,-ly),
                (-lx,0),(-lx,ly)]

# Read cell field data
def read_field(index, frame, time):
    global fileroot
    filename = fileroot + ("cell_{:d}.dat.{:d}".format(index, time))
    with open(filename,'r') as field_reader:
        # For computing the local field centre of mass
        xavg = 0.0
        yavg = 0.0
        mass = 0.0
        local_field = np.zeros((clx,cly))
        for l, line in enumerate(field_reader):
            data = line.split()
            if (len(data) != 3): continue
            i = int(data[0])
            j = int(data[1])
            phi = float(data[2])
            local_field[i,j] = phi
            xavg += (phi*(i+0.5)) # Uset the centre of a lattice element
            yavg += (phi*(j+0.5))
            mass += phi
        if (mass > 0.0):
            xavg /= mass
            yavg /= mass
        else:
            xavg = 0.0
            yavg = 0.0
        contours = measure.find_contours(local_field, 1.0)
        contour = contours[0]
        for pt in periodic_loc:
            add_polygon(index, xavg, yavg, pos[frame][index][0], 
                        pos[frame][index][1], pt[0], pt[1], contour, frame)

# Read position data
nlines = npoints + 2
reader = open(pos_file, 'r')
while True:
    # Read header section (including time info)
    for i in range(2):
        line = reader.readline()
    if (not line): break
    data = line.split()
    time = int(data[1])
    if (time > tend):
        break
    elif (time < tstart or (time-tstart)%tinc != 0):
        for i in range(npoints):
            line = reader.readline()
    else:
        print("Reading data at timestep = {:d}".format(time))
        frame = (time-tstart)//tinc
        for n in range(npoints):
            line = reader.readline()
            data = line.split()
            
            # Read wrapped position data
            x = float(data[0])
            y = float(data[1])
            pos[frame][n] = (x,y)
            
            # Read cell field data
            read_field(n, frame, time)
            
reader.close()
    
# Make animation
# Plot settings
data_min = 0
data_max = npoints-1
norm = mpl.colors.Normalize(vmin=data_min, vmax=data_max, clip=True)
mapper = mplcm.ScalarMappable(norm=norm, cmap=mplcm.RdYlBu_r)
mapper.set_array([])

fig, ax = plt.subplots()
ax.set_xlim([0,lx])
ax.set_ylim([0,ly])

if (not make_movie):
# Use Latex typesetting when not making movies
    mpl.rcParams["text.latex.preamble"] = [
        r'\usepackage{amsmath}',
        r'\usepackage{amssymb}',
        r'\usepackage[scaled=1]{helvet}',
        r'\usepackage{sansmath}',
        r'\sansmath']
    plt.rc("text", usetex=True)

# Draw borders but no axes and ticks
plt.tick_params(axis="both", which="both", bottom=False, top=False,
                labelbottom=False, right=False, left=False, labelleft=False)

# Set plot margins
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

# Get the artist for plotting centre of Voronoi cells
for i in range(nframes):
    points[i] = np.array(points[i])
plt_pts, = ax.plot([],[], '.', markersize=5, color="black") # Empty data

# Get the artist for plotting the time label
if (use_label):
    plt_time_txt = ax.text(0.5,0.005, "", fontsize=14,
                           horizontalalignment="center",
                           transform=plt.gcf().transFigure)

# Get the artist for plotting the polygons
patches_int = PatchCollection([], linewidth=0)
plt_polygons_int = ax.add_collection(patches_int)
patches_out = PatchCollection([], linewidth=1)
plt_polygons_out = ax.add_collection(patches_out)
plt_polygons_out.set_facecolor((0.0,0.0,0.0,0.0))
plt_polygons_out.set_edgecolor("black")

def plot_data(frame):
    global pos, index_map, time_map
    
    # Plot centres of cells
    plt_pts.set_xdata(points[frame][:,0])
    plt_pts.set_ydata(points[frame][:,1])
    
    # Plot the cell polygons
    colors = []
    for pt in range(len(polygons[frame])):
        colors.append(mapper.to_rgba(index_map[frame][pt], alpha=0.8))
    plt_polygons_int.set_paths(polygons[frame])
    plt_polygons_int.set_facecolor(colors)
    plt_polygons_out.set_paths(polygons[frame])
    
    # Plot the time label
    if (use_label):
        plt_time_txt.set_text(r"$D_rt = {:.1f}$".format(time_map[frame]*Dr*dt))
    
    return plt_pts, plt_polygons_int, plt_polygons_out

if (make_movie):
    if (print_to_screen):
        ani = FuncAnimation(fig, plot_data, np.arange(nframes), 
                            fargs=[], interval=1)
        plt.show()
    else:
        ani = FuncAnimation(fig, plot_data, np.arange(nframes), 
                            fargs=[], interval=1)
        Writer = animation.writers["ffmpeg"]
        writer = Writer(fps=15, bitrate=1500)
        ani.save(out_file,writer=writer)
else:
    plot_data(0)
    if (print_to_screen):
        plt.show()
    else:
        plt.savefig(out_file, transparent=True)
