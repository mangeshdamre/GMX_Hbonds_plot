#!/usr/bin/env python3.8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy.stats
import math
import string
import sys
import argparse
from tabulate import tabulate
#####################################################################

parser = argparse.ArgumentParser(description='Usage of HBONDS script!')
parser.add_argument("-x",
                    "--x",
                    type=int,
                    required=False,
                    default=0,
                    help="This is the 'xmin' variable : default = %(default)s")
parser.add_argument("-X",
                    "--X",
                    type=int,
                    required=True,
                    default=25,
                    help="This is the 'xmax' variable : default = %(default)s")
parser.add_argument("-y",
                    "--y",
                    type=int,
                    required=False,
                    default=0,
                    help="This is the 'ymin' variable : default = %(default)s")
parser.add_argument("-Y",
                    "--Y",
                    type=int,
                    required=True,
                    default=10,
                    help="This is the 'Ymax' variable : default = %(default)s")
parser.add_argument(
    "-t",
    "--time",
    type=int,
    required=True,
    default=25,
    help="This is the 'Simulation time in ns' variable : default = %(default)s"
)
args = parser.parse_args()
#####################################################################
xmin = args.x
xmax = args.X
ymin = args.y
ymax = args.Y
time = args.time * 1000

print(f'\n')
#####################################################################
#xmax = int(sys.argv[1])	# x-axis limit
#time = time * 1000

inp_nameA = "A/hbnumA.xvg"
inp_nameB = "B/hbnumB.xvg"
inp_nameC = "C/hbnumC.xvg"
inp_nameD = "D/hbnumD.xvg"
inp_nameE = "E/hbnumE.xvg"
inp_nameF = "F/hbnumF.xvg"

inp_dataA = np.loadtxt(inp_nameA, comments=['#', '@', '&'])[0:time]
inp_dataB = np.loadtxt(inp_nameB, comments=['#', '@', '&'])[0:time]
inp_dataC = np.loadtxt(inp_nameC, comments=['#', '@', '&'])[0:time]
inp_dataD = np.loadtxt(inp_nameD, comments=['#', '@', '&'])[0:time]
inp_dataE = np.loadtxt(inp_nameE, comments=['#', '@', '&'])[0:time]
inp_dataF = np.loadtxt(inp_nameF, comments=['#', '@', '&'])[0:time]

framesA = inp_dataA[:, 0] / 1000
framesB = inp_dataB[:, 0] / 1000
framesC = inp_dataC[:, 0] / 1000
framesD = inp_dataD[:, 0] / 1000
framesE = inp_dataE[:, 0] / 1000
framesF = inp_dataF[:, 0] / 1000

hbondsA = inp_dataA[:, 1]
hbondsB = inp_dataB[:, 1]
hbondsC = inp_dataC[:, 1]
hbondsD = inp_dataD[:, 1]
hbondsE = inp_dataE[:, 1]
hbondsF = inp_dataF[:, 1]

print(
    tabulate(
        [['Total time (ns)', time / 1000], ['xmin', xmin], ['xmax', xmax],
         ['ymin', ymin], ['ymax', ymax],
         ['Total number of frames', len(framesA)]],
        tablefmt='fancy_grid'))
print(
    tabulate([['1', framesA[0]], ['2', framesA[1]], ['Last', framesA[-1]]],
             headers=['Frame', 'Time (ns)'],
             tablefmt='fancy_grid',
             floatfmt=".3f"))
#print(len(hbondsA))
#print(hbondsA)
print(f'\n')


def movingaverage(interval, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(interval, window, 'same')


############### Global setting ##############
xmin = xmin
xmax = xmax
xinterval = 10
ymin = ymin
ymax = ymax
yinterval = 2
alpha1 = 0.3
alpha2 = 1.0
x_ticks = np.arange(xmin, xmax, xinterval)
y_ticks = np.arange(ymin, ymax, yinterval)
plt.rcParams['axes.grid'] = True
plt.rcParams["legend.loc"] = 'upper right'
#############################################
# plot with various axes scales
plt.figure()

plt.subplot(321)
plt.plot(framesA,
         hbondsA,
         label='A',
         color='red',
         marker='',
         linewidth=1,
         alpha=alpha1)
y_av = movingaverage(hbondsA, 50)
plt.plot(framesA,
         y_av,
         label='',
         color='red',
         marker='',
         linewidth=1,
         alpha=alpha2)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.title('')
plt.legend()

plt.subplot(322)
plt.plot(framesB,
         hbondsB,
         label='B',
         color='blue',
         marker='',
         linewidth=1,
         alpha=alpha1)
y_av = movingaverage(hbondsB, 50)
plt.plot(framesB,
         y_av,
         label='',
         color='blue',
         marker='',
         linewidth=1,
         alpha=alpha2)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
ax = plt.gca()
ax.axes.yaxis.set_ticklabels([])
plt.legend()

plt.subplot(323)
plt.plot(framesC,
         hbondsC,
         label='C',
         color='green',
         marker='',
         linewidth=1,
         alpha=alpha1)
y_av = movingaverage(hbondsC, 50)
plt.plot(framesC,
         y_av,
         label='',
         color='green',
         marker='',
         linewidth=1,
         alpha=alpha2)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.ylabel('Number')
plt.legend()

plt.subplot(324)
plt.plot(framesD,
         hbondsD,
         label='D',
         color='purple',
         marker='',
         linewidth=1,
         alpha=alpha1)
y_av = movingaverage(hbondsD, 50)
plt.plot(framesD,
         y_av,
         label='',
         color='purple',
         marker='',
         linewidth=1,
         alpha=alpha2)
plt.xticks(x_ticks)
plt.yticks(y_ticks)
ax = plt.gca()
ax.axes.yaxis.set_ticklabels([])
plt.legend()

plt.subplot(325)
plt.plot(framesE,
         hbondsE,
         label='E',
         color='brown',
         marker='',
         linewidth=1,
         alpha=alpha1)
y_av = movingaverage(hbondsE, 50)
plt.plot(framesE,
         y_av,
         label='',
         color='brown',
         marker='',
         linewidth=1,
         alpha=alpha2)
plt.xlabel('time(ns)')
plt.xticks(x_ticks)
plt.yticks(y_ticks)
plt.legend()

plt.subplot(326)
plt.plot(framesF,
         hbondsF,
         label='F',
         color='teal',
         marker='',
         linewidth=1,
         alpha=alpha1)
y_av = movingaverage(hbondsF, 50)
plt.plot(framesF,
         y_av,
         label='',
         color='teal',
         marker='',
         linewidth=1,
         alpha=alpha2)
plt.xlabel('time(ns)')
plt.xticks(x_ticks)
plt.yticks(y_ticks)
ax = plt.gca()
ax.axes.yaxis.set_ticklabels([])
plt.legend()

plt.suptitle('Hydrogen Bonds')

plt.subplots_adjust(top=0.90,
                    bottom=0.10,
                    left=0.10,
                    right=0.95,
                    hspace=0.07,
                    wspace=0.05)
plt.savefig('hbnum-python.png')
