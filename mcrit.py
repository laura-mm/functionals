# code to plot mcrit
import matplotlib
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import gridspec
plt.rc('text', usetex=True)
plt.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
lett = ['$(a)$', '$(b)$', '$(c)$']

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('axes', titlesize=20)

ggrid = 200
m2grid = 1409
m5grid = 1409

m2 = np.linspace(((2000 - m2grid)/100.0) - 10.0, 10, m2grid + 1)
g2 = np.linspace((1.0/ggrid) - 1, 1, ggrid)
m5 = np.linspace(((2000 - m5grid)/100.0) -10.0, 10, m5grid + 1)
g5 = np.linspace((1.0/ggrid) - 1, 1, ggrid)
g2, m2 = np.meshgrid(g2, m2)
g5, m5 = np.meshgrid(g5, m5)

mcrit2 = np.loadtxt('m2.txt', delimiter = ",")
mcrit2 = mcrit2.reshape(m2grid + 1, ggrid)
mcrit5 = np.loadtxt('m5l.txt', delimiter = ",")
mcrit5 = mcrit5.reshape(m5grid + 1, ggrid)

low = min(mcrit2.min(), mcrit5.min())
high = max(mcrit2.max(), mcrit5.max())

lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']

plt.subplots(1, 2, figsize=(18,9))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.25]) 

plt.subplot(gs[0]) # a = 2
plt.pcolormesh(g2, m2, mcrit2, norm = colors.LogNorm(vmin = low, vmax = high), cmap = cm.gist_ncar)
plt.axhline(y = 0, linestyle = '--', color = 'r')
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\mu$')
plt.xlim([-1, 1])
plt.ylim([-4, 10])
plt.text(0.85, 9.2, lett[0], fontsize=20)
plt.title(r'$a = 2$')

ax1 = plt.subplot(gs[1]) # a = 0.5
plt.pcolormesh(g5, m5, mcrit5, norm = colors.LogNorm(vmin = low, vmax = high), cmap = cm.gist_ncar)
plt.axhline(y = 0, linestyle = '--', color = 'r')
plt.xlabel(r'$\gamma$')
plt.colorbar()
plt.xlim([-1, 1])
plt.ylim([-4, 10])
plt.text(0.85, 9.2, lett[1], fontsize=20)
plt.title(r'$a = 0.5$')
#plt.title(r'$\sigma_c$ in $\mu - \gamma$ plane')

plt.subplots_adjust(left = 0.05, bottom = 0.08, right = 1.0, top = 0.96, wspace = 0.07)
plt.savefig('mcritV2.jpg')

plt.show()
