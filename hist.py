# histograms

import numpy as np
import matplotlib.pyplot as plt
import math
plt.rc('text', usetex=True)
plt.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.rc('axes', labelsize=20)
plt.rc('axes', titlesize=20)


mu = 0
a = 2
N = 200
runs = 100

##### for a = 2
g0 = np.loadtxt('meashist_20_0_0.txt', delimiter=",")
g0 = g0.reshape(len(g0)/9, 9)
g0 = np.transpose(g0)
g2 = np.loadtxt('meashist_20_0_2.txt', delimiter=",")
g2 = g2.reshape(len(g2)/9, 9)
g2 = np.transpose(g2)
g4 = np.loadtxt('meashist_20_0_4.txt', delimiter=",")
g4 = g4.reshape(len(g4)/9, 9)
g4 = np.transpose(g4)

h = np.loadtxt('hist2_20_0.txt', delimiter=",")
h = h.reshape(3, 3, runs*N)

g = [g0, g2, g4]
######


###### for a = 0.5
j0 = np.loadtxt('meashist_5_0_0.txt', delimiter=",")
j0 = j0.reshape(len(j0)/9, 9)
j0 = np.transpose(j0)
j2 = np.loadtxt('meashist_5_0_2.txt', delimiter=",")
j2 = j2.reshape(len(j2)/9, 9)
j2 = np.transpose(j2)
j4 = np.loadtxt('meashist_5_0_4.txt', delimiter=",")
j4 = j4.reshape(len(j4)/9, 9)
j4 = np.transpose(j4)

k = np.loadtxt('hist2_5_0.txt', delimiter=",")
k = k.reshape(3, 3, runs*N)

j = [j0, j2, j4]
######

col = ['r', 'b', 'g']
lab = ['$\gamma = -1$', '$\gamma = 0$', '$\gamma = 1$']
tit = ['$\sigma = 10^{-1}$', '$\sigma = 10^{-0.5}$', '$\sigma = 10^0$']
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$', '$(g)$', '$(h)$', '$(i)$', '$(j)$', '$(k)$', '$(l)$',]

# sigma, z1, z2, top, middle, bottom, M, q, help for g

plt.figure(figsize=(12,8)) # a = 2

for gi in range (3):

	i0 = np.where(abs(np.log10(g[gi][0]) + 1) == np.min(abs(np.log10(g[gi][0]) + 1)))
	i1 = np.where(abs(np.log10(g[gi][0]) + 0.5) == np.min(abs(np.log10(g[gi][0]) + 0.5)))
	i2 = np.where(abs(np.log10(g[gi][0]) + 0) == np.min(abs(np.log10(g[gi][0]) + 0)))
	i = [i0, i1, i2]

	plt.subplot(3, 4, (4*gi) + 1)
	plt.semilogx(g[gi][0], g[gi][1], col[gi])
	plt.semilogx(g[gi][0], -g[gi][2], col[gi])
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlim(10**(-1), 10**1)
	plt.ylim(-10, 20)
	plt.xlabel(r'$\sigma$')
	plt.ylabel(r'$z_1, z_2$')
	plt.yticks(np.arange(-10, 30, 10))
	plt.text(4, 14, lett[4*gi], fontsize=20)

	
	if gi == 2:
		plt.semilogx(0, 0, col[0], label = lab[0])
		plt.semilogx(0, 0, col[1], label = lab[1])
		plt.semilogx(0, 0, col[2], label = lab[2])

	for ii in range (3):

		plt.subplot(3, 4, (4*gi) + 2 + ii)
		gaus_m = (1 + (mu*g[gi][6][i[ii]]))/g[gi][8][i[ii]]
		gaus_sd = (g[gi][0][i[ii]]*np.sqrt(g[gi][7][i[ii]]))/g[gi][8][i[ii]]
		x1 = (g[gi][1][i[ii]]*gaus_sd) + gaus_m
		x2 = (-g[gi][2][i[ii]]*gaus_sd) + gaus_m
		x = np.linspace(x2, x1, 1001)
		plt.hist(h[gi][ii], 50, normed = 1, edgecolor = col[gi], color = col[gi])
		plt.plot(x, np.exp((-(x-gaus_m)**2)/(2*gaus_sd**2))/(gaus_sd*np.sqrt(2*math.pi)), 'k')
		plt.xlim(0, 3)
		plt.ylim(0, 4.5)
		plt.xlabel(r'$x$')
		plt.xticks(np.arange(0, 4, 1))
		plt.yticks(np.arange(0, 5, 1))
		plt.text(2.25, 3.75, lett[4*gi + 1 + ii], fontsize=20)
		if gi == 0:
			plt.title(tit[ii])

plt.subplots_adjust(left = 0.07, right = 0.98, bottom = 0.1, top = 0.94, wspace = 0.2, hspace = 0.3)
plt.savefig('hist_2.pdf')

###################################################

plt.figure(figsize=(12,8)) # a = 0.5

for gi in range (3):

	i0 = np.where(abs(np.log10(j[gi][0]) + 1) == np.min(abs(np.log10(j[gi][0]) + 1)))
	i1 = np.where(abs(np.log10(j[gi][0]) + 0.5) == np.min(abs(np.log10(j[gi][0]) + 0.5)))
	i2 = np.where(abs(np.log10(j[gi][0]) + 0) == np.min(abs(np.log10(j[gi][0]) + 0)))

	if gi == 2:
		i2 = (i2[0][9]) # because there were many places where sigma = 1, we take the one in the middle

	i = [i0, i1, i2]

	plt.subplot(3, 4, (4*gi) + 1)
	plt.semilogx(j[gi][0], j[gi][1], col[gi])
	plt.semilogx(j[gi][0], -j[gi][2], col[gi])
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlim(10**(-1), 10**1)
	plt.ylim(-5, 5)
	plt.xlabel(r'$\sigma$')
	plt.ylabel(r'$z_1, z_2$')
	plt.yticks(np.arange(-4, 6, 2))
	plt.text(4, 3, lett[4*gi], fontsize=20)


	if gi == 2:
		plt.semilogx(0, 0, col[0], label = lab[0])
		plt.semilogx(0, 0, col[1], label = lab[1])
		plt.semilogx(0, 0, col[2], label = lab[2])

	for ii in range (3):

		plt.subplot(3, 4, (4*gi) + 2 + ii)
		gaus_m = (1 + (mu*j[gi][6][i[ii]]))/j[gi][8][i[ii]]
		gaus_sd = (j[gi][0][i[ii]]*np.sqrt(j[gi][7][i[ii]]))/j[gi][8][i[ii]]
		x1 = (j[gi][1][i[ii]]*gaus_sd) + gaus_m
		x2 = (-j[gi][2][i[ii]]*gaus_sd) + gaus_m
		x = np.linspace(x2, x1, 1001)
		plt.hist(k[gi][ii], 50, normed = 1, edgecolor = col[gi], color = col[gi])
		plt.plot(x, np.exp((-(x-gaus_m)**2)/(2*gaus_sd**2))/(gaus_sd*np.sqrt(2*math.pi)), 'k')
		plt.xlim(0.47, 1.53)
		plt.ylim(0, 4.5)
		plt.xlabel(r'$x$')
		plt.xticks(np.arange(0.5, 2, 0.5))
		plt.yticks(np.arange(0, 5, 1))
		plt.text(1.25, 3.75, lett[4*gi + 1 + ii], fontsize=20)
		if gi == 0:
			plt.title(tit[ii])

plt.subplots_adjust(left = 0.07, right = 0.98, bottom = 0.1, top = 0.94, wspace = 0.2, hspace = 0.3)
plt.savefig('hist_5.pdf')

#######################################################

plt.show()



