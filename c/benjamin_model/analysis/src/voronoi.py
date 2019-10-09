# voronoi.py
# A script to plot a Voronoi representation of the system and to visualise
# a specific data set

import sys
import numpy as np
import scipy.spatial as spatial
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
from matplotlib.animation import FuncAnimation
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

args = sys.argv
if (len(args) != 20):
    print "usage: voronoi.py npoints lx ly Dr dt data_col data_min data_max tic_start tic_end tic_inc tstart tend tinc make_movie print_to_screen pos_file data_file out_file"
    sys.exit(1)

npoints = int(args.pop(1))
lx = int(args.pop(1))
ly = int(args.pop(1))
Dr = float(args.pop(1))
dt = float(args.pop(1))
data_col = int(args.pop(1))
data_min = float(args.pop(1))
data_max = float(args.pop(1))
tic_start = float(args.pop(1))
tic_end = float(args.pop(1))
tic_inc = float(args.pop(1))
tstart = int(args.pop(1))
tend = int(args.pop(1))
tinc = int(args.pop(1))
make_movie = bool(int(args.pop(1)))
print_to_screen = bool(int(args.pop(1)))
pos_file = args.pop(1)
data_file = args.pop(1)
out_file = args.pop(1)
xbuff = 0.2*lx
ybuff = 0.2*ly
nframes = (tend-tstart)/tinc+1

use_label = 0 # 1

if (not make_movie):
    tend = tstart

# Data arrays
pos = [[] for i in xrange(nframes)]
data_val = [[0.0 for j in xrange(npoints)] for i in xrange(nframes)]
index_map = [[] for i in xrange(nframes)]
time_map = [i*tinc+tstart for i in xrange(nframes)]

# Useful functions for plotting
def add_point(index, x, y, frame):
    global pos, data_map, lx, ly, xbuff, ybuff
    if (x < lx+xbuff and x > -xbuff and y < ly+ybuff and y > -ybuff):
        pos[frame].append((x,y))
        index_map[frame].append(index)
#        data_map[frame][(x,y)] = val

periodic_loc = [(lx,-ly),(lx,0),(lx,ly),(0,-ly),(0,0),(0,ly),(-lx,-ly),
                (-lx,0),(-lx,ly)]


# Read position data
nlines = npoints + 2
reader = open(pos_file, 'r')

while True:
    # Read header section (including time info)
    for i in xrange(2):
        line = reader.readline()
    if (not line): break
    data = line.split()
    time = int(data[1])
    if (time > tend):
        break
    elif (time < tstart or (time-tstart) % tinc != 0):
        for i in xrange(npoints):
            line = reader.readline()
    else:
        frame = (time-tstart)/tinc
        for n in xrange(npoints):
            line = reader.readline()
            data = line.split()
            x = float(data[0])
            y = float(data[1])
            for pt in periodic_loc:
                add_point(n, x+pt[0], y+pt[1], frame)
    
reader.close()

if (data_col >= 0):
    data_reader = open(data_file, 'r')
    while True:
        # Read header section (including time info)
        for i in xrange(2):
            line = data_reader.readline()
        if (not line): break
        data = line.split()
        time = int(data[1])
        if (time > tend):
            break
        elif (time < tstart or (time-tstart) % tinc != 0):
            for n in xrange(npoints):
                data_reader.readline()
        else:
            frame = (time-tstart)/tinc
            for n in xrange(npoints):
                line = data_reader.readline()
                data = line.split()
                data_val[frame][n] = float(data[data_col])
                
    data_reader.close()

# Make animation
# Plot settings
norm = mpl.colors.Normalize(vmin=data_min, vmax=data_max, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.RdYlBu_r)
mapper.set_array([])
fig, ax = plt.subplots()
ax.set_xlim([0,lx])
ax.set_ylim([0,ly])
cbar = plt.colorbar(mapper)
cbar.set_ticks(np.arange(tic_start, tic_end+tic_inc/2.0, tic_inc))

if (not make_movie):
# Use Latex typesetting when not making movies
    mpl.rcParams["text.latex.unicode"] = True
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
plt.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05)

# Get the artist for plotting centre of Voronoi cells
plt_pts, = ax.plot([],[], '.', markersize=5, color="black") # Empty data

# Get the artist for plotting the time label
if (use_label):
    plt_time_txt = ax.text(0.45,0.005,"",fontsize=14,
                           horizontalalignment="center",
                           transform=plt.gcf().transFigure)

# Get the artist for plotting the polygons
patches = PatchCollection([], linewidth=1.0)
plt_polygons = ax.add_collection(patches)

def plot_data(frame):
    global pos, data_map, time_map

    # Plot centres of Voronoi cells
    vor = spatial.Voronoi(pos[frame])
    plt_pts.set_xdata(vor.points[:,0])
    plt_pts.set_ydata(vor.points[:,1])
    
    # Plot the Voronoi polygons
    colors = []
    polygons = []
    for r in xrange(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if -1 in region: continue
#        pt = tuple(vor.points[r])
        poly = [vor.vertices[i] for i in region]
        polygons.append(Polygon(poly))
        if (data_col >= 0):
            colors.append(mapper.to_rgba(
                data_val[frame][index_map[frame][r]]))
        else:
            colours.append(mapper.to_rgba(index_map[frame][r]))
    plt_polygons.set_paths(polygons)
    plt_polygons.set_facecolor(colors)
    plt_polygons.set_edgecolor("black")

    # Plot the time label
    if (use_label):
        plt_time_txt.set_text(r"$D_rt = {:.1f}$".format(time_map[frame]*Dr*dt))
    
    return plt_pts, plt_polygons,

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
