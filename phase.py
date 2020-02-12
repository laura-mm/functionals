# heat maps and phase diagrams

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
plt.rc('text', usetex=True)
plt.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
lett = ['$(a)$', '$(b)$', '$(c)$']

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('axes', titlesize=20)

grid = 100 # for heat map
gridl = 1000 # lines on heat map
gridg = 200 # 3d plots

# could also look at interpolation - to make it more smooth and less pixelly


# fp solution plot acrit only
################################################

ast = 3 #start for a
gst = 3 # start for gamma

plt.figure(1)
ag = np.linspace((gst*(2.0/gridg)) - 1, 1, gridg + 1 - gst)
a = np.linspace(ast*(2.0/gridg), 2, gridg + 1 - ast)
ag, a = np.meshgrid(ag, a)
acrit = np.loadtxt('acrit_0.txt', delimiter=",")
acrit = acrit.reshape(gridg + 1 - gst, gridg + 1 - ast)
acrit = np.transpose(acrit)
plt.pcolormesh(ag, a, acrit, norm = colors.LogNorm(vmin = acrit.min(), vmax = acrit.max()), cmap = cm.gist_ncar)
plt.axhline(y = 1, linestyle = '--', color = 'r')
plt.colorbar()
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$a$')
plt.xlim([-1, 1])
plt.ylim([0, 2])
plt.subplots_adjust(left = 0.1, right = 1.0, top = 0.96, bottom = 0.1)
plt.savefig('acrit.pdf')

###################################

#heat maps
##############################################

g = np.linspace(-1, 1, grid + 1)
x = np.linspace(-1, 3, grid + 1)
s = np.power(10, x)

la5 = np.loadtxt('phline_5_0.txt', delimiter=",")
la2 = np.loadtxt('phline_20_0.txt', delimiter=",")

gl = np.linspace((2.0/gridl) - 1, 1, gridl)
l1 = math.sqrt(2)/(1+gl);

# holling with 3 values for a
##########################################

plt.subplots(1, 3, figsize=(18,6))

plt.subplot(1, 3, 1) # linear
linear = np.loadtxt('ltran.txt', delimiter=",")
linear = linear.reshape(grid + 1, grid + 1, 3)
linear = np.transpose(linear, (1, 0, 2))
plt.imshow(linear, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(gl, l1, 'k')
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\sigma$')
plt.text(0.8, 500, lett[0], fontsize=20)
plt.title(r'$a \to \infty$')

plt.subplot(1, 3, 2) # holling a = 2
h2 = np.loadtxt('h2.txt', delimiter=",")
h2 = h2.reshape(grid + 1, grid + 1, 3)
h2 = np.transpose(h2, (1, 0, 2))
plt.imshow(h2, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.xlabel(r'$\gamma$')
plt.text(0.8, 500, lett[1], fontsize=20)
plt.title(r'$a = 2$')

plt.subplot(1, 3, 3) # holling a = 5
h5 = np.loadtxt('h5.txt', delimiter=",")
h5 = h5.reshape(grid + 1, grid + 1, 3)
h5 = np.transpose(h5, (1, 0, 2))
plt.imshow(h5, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.xlabel(r'$\gamma$')
plt.ylim([0.1, 1000])
plt.text(0.8, 500, lett[2], fontsize=20)
plt.title(r'$a = 0.5$')

plt.subplots_adjust(left = 0.05, right = 0.98, wspace = 0.17)
plt.savefig('HolheatV2.pdf')

# piecewise heatmaps, 2 values for a
####################################################

plt.subplots(1, 2, figsize=(12,6))

plt.subplot(1, 2, 1) # p2
p2 = np.loadtxt('p2.txt', delimiter=",")
p2 = p2.reshape(grid + 1, grid + 1, 3)
p2 = np.transpose(p2, (1, 0, 2))
plt.imshow(p2, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(gl, l1, 'k')
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\sigma$')
plt.text(0.8, 500, lett[0], fontsize=20)
plt.title(r'$a = 2$')


plt.subplot(1, 2, 2) # p5
p5 = np.loadtxt('p5.txt', delimiter=",")
p5 = p5.reshape(grid + 1, grid + 1, 3)
p5 = np.transpose(p5, (1, 0, 2))
plt.imshow(p5, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(gl, la5, 'k')
plt.xlabel(r'$\gamma$')
plt.ylim([0.1, 1000])
plt.text(0.8, 500, lett[1], fontsize=20)
plt.title(r'$a = 0.5$')


plt.subplots_adjust(left = 0.07, right = 0.98, wspace = 0.17)
plt.savefig('pheatV2.pdf')

# peicewise a = 2 with different values of N
#######################################################

plt.subplots(1, 3, figsize=(18,6)) 

plt.subplot(1, 3, 1) # N = 200
p2 = np.loadtxt('p2.txt', delimiter=",")
p2 = p2.reshape(grid + 1, grid + 1, 3)
p2 = np.transpose(p2, (1, 0, 2))
plt.imshow(p2, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(gl, l1, 'k')
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\sigma$')
plt.text(0.8, 500, lett[0], fontsize=20)
plt.title(r'$N = 200$')

plt.subplot(1, 3, 2) # N = 50
p2_50 = np.loadtxt('p2_50.txt', delimiter=",")
p2_50 = p2_50.reshape(grid + 1, grid + 1, 3)
p2_50 = np.transpose(p2_50, (1, 0, 2))
plt.imshow(p2_50, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(gl, l1, 'k')
plt.xlabel(r'$\gamma$')
plt.text(0.8, 500, lett[1], fontsize=20)
plt.title(r'$N = 50$')

plt.subplot(1, 3, 3) # N = 20
p2_20 = np.loadtxt('p2_20.txt', delimiter=",")
p2_20 = p2_20.reshape(grid + 1, grid + 1, 3)
p2_20 = np.transpose(p2_20, (1, 0, 2))
plt.imshow(p2_20, extent=[-1, 1, 0.1, 1000], origin = 'lower')
plt.yscale('log')
plt.semilogy(gl, l1, 'k')
plt.xlabel(r'$\gamma$')
plt.text(0.8, 500, lett[2], fontsize=20)
plt.title(r'$N = 20$')


plt.subplots_adjust(left = 0.05, right = 0.98, wspace = 0.17)
plt.savefig('Nvar.pdf')


###############################
# end of heat maps


plt.show()




