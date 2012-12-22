#!/usr/bin/env python

import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.texmanager import TexManager
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

mp.rc('text', usetex=True)
mp.rcParams['font.sans-serif'] = 'computer modern bright'
TexManager.font_info['computern modern bright'] = (
    ('cmbr', r'\usepackage{cmbright}'))
mp.rcParams['text.latex.preamble'] = (
    r'\PassOptionsToPackage{dvips}{graphicx}',
    r'\usepackage{cmbright}',
    r'\usepackage{units}',
    r'\usepackage{rotating}',
    )
mp.rc('font', family='sans-serif', size=10)
mp.rc('axes', titlesize='small', labelsize='x-small')
mp.rc('legend', fontsize='x-small')
mp.rc('xtick', labelsize='x-small')
mp.rc('ytick', labelsize='x-small')

CMRMap = mp.colors.LinearSegmentedColormap(
    'CMR',
    {
      'blue': ((0.0, 0.0, 0.0),
        (0.125, 0.5, 0.5),
        (0.25, 0.75, 0.75),
        (0.375, 0.5, 0.5),
        (0.5, 0.15, 0.15),
        (0.625, 0.0, 0.0),
        (0.75, 0.1, 0.1),
        (0.875, 0.5, 0.5),
        (1.0, 1.0, 1.0)),
      'green': ((0.0, 0.0, 0.0),
        (0.125, 0.15, 0.15),
        (0.25, 0.15, 0.15),
        (0.375, 0.2, 0.2),
        (0.5, 0.25, 0.25),
        (0.625, 0.5, 0.5),
        (0.75, 0.75, 0.75),
        (0.875, 0.9, 0.9),
        (1.0, 1.0, 1.0)),
      'red': ((0.0, 0.0, 0.0),
        (0.125, 0.15, 0.15),
        (0.25, 0.3, 0.3),
        (0.375, 0.6, 0.6),
        (0.5, 1.0, 1.0),
        (0.625, 0.9, 0.9),
        (0.75, 0.9, 0.9),
        (0.875, 0.9, 0.9),
        (1.0, 1.0, 1.0))}
      )

data = np.loadtxt('tmp.dat')
npv = sum(data[:,2]==data[0,2])
data = data.reshape((npv, -1, data.shape[1]))
Z = data[:, :, 0]
PV = data[:, :, 2]
Ydot = data[:, :, 4]

fig = plt.gcf()
fig.set_size_inches(4, 3)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Z, PV, Ydot, rstride=1, cstride=1, cmap=cm.jet, #CMRMap,
    linewidth=0, antialiased=False, shade='interp')
ax.set_xlabel(r'$z\ \left[-\right]$')
ax.set_ylabel(r'PV\ $\left[-\right]$')
ax.set_zlabel(r'\turnbox{180}{$\dot{\omega}_c\ \left[\unitfrac{1}{s}\right]$}')
#ax.set_zlabel(r'$\dot{\omega}_c\ \left[\unitfrac{1}{s}\right]$')
ax.view_init(azim=220)

#fig.savefig('Ydot.eps')
fig.savefig('Ydot.png', dpi=600)

#plt.show(block=True)
