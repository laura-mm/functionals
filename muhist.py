# histograms for varying mu

import numpy as np
import matplotlib.pyplot as plt
import math
plt.rc('text', usetex=True)
plt.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('axes', titlesize=20)
plt.rc('legend', fontsize=20)

#mu = [-4, -2, 0] # can change these
mgrid = 2000
mstart = 591
m = np.linspace((0.01*mstart) - 10, 10, mgrid + 1 - mstart)


h2 = np.loadtxt('mhist2.txt', delimiter=",")
h2 = h2.reshape(mgrid + 1 - mstart, 5, 8)
h2 = np.transpose(h2, (1, 2, 0))
h5 = np.loadtxt('mhist5.txt', delimiter=",")
h5 = h5.reshape(mgrid + 1 - mstart, 5, 8)
h5 = np.transpose(h5, (1, 2, 0))

minlist = []

ev = 100

lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$', '$(g)$', '$(h)$', '$(i)$', '$(j)$', '$(k)$', '$(l)$',]

colo = ['r', 'm', 'b', 'c', 'g']
mark = ['o', 's', '*', 'x', '+']
lab = ['$\gamma = -1$', '$\gamma = -0.5$', '$\gamma = 0$', '$\gamma = 0.5$', '$\gamma = 1$']

plt.subplots(2, 2, figsize=(12,6)) # sigma and phi at critical point

plt.subplot(2, 2, 1) # sigma, a = 2
for gi in range (1, 5):
	plt.semilogy(m, h2[gi][1], colo[gi], marker = mark[gi], markevery = ev, label = lab[gi]) #sigma
plt.ylabel(r'$\sigma_c$')
plt.xlim([-4, 10])
plt.ylim([10**-0.5, 10**1.5])
plt.text(8.5, 10**1.2, lett[0], fontsize=20)
plt.title(r'$a = 2$')

plt.subplot(2, 2, 3) # phi, a = 2
for gi in range (1, 5):
	plt.plot(m, h2[gi][7], colo[gi], marker = mark[gi], markevery = ev, label = lab[gi]) #phi
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\phi$')
plt.xlim([-4, 10])
plt.ylim([0, 0.5])
plt.text(8.5, 0.42, lett[2], fontsize=20)
plt.legend(bbox_to_anchor=(0, -0.6, 2.14, 0.2), loc='lower left', ncol=4, mode="expand", borderaxespad=0)

plt.subplot(2, 2, 2) # sigma, a = 0.5
for gi in range (1, 5):
	plt.semilogy(m, h5[gi][1], colo[gi], marker = mark[gi], markevery = ev, label = lab[gi]) #sigma
plt.xlim([-4, 10])
plt.ylim([10**-0.5, 10**1.5])
plt.text(8.5, 10**1.2, lett[1], fontsize=20)
plt.title(r'$a = 0.5$')

plt.subplot(2, 2, 4) # phi, a = 0.5
for gi in range (1, 5):
	plt.plot(m, h5[gi][7], colo[gi], marker = mark[gi], markevery = ev, label = lab[gi]) #phi
plt.xlabel(r'$\mu$')
plt.xlim([-4, 10])
plt.ylim([0, 0.5])
plt.text(8.5, 0.42, lett[3], fontsize=20)

plt.subplots_adjust(left = 0.07, right = 0.98, bottom = 0.25, top = 0.94, wspace = 0.14, hspace = 0.27)
plt.savefig('muhistV2.pdf')

####################################################

plt.figure(figsize=(12,12)) # z1, z2, phi, sigma, 3 histograms at various mu, a = 0.5
plt.title('a=0.5')

# gamma, sigma, z2, z1, M, q, help, phi

for gi in range (1, 5):

	im = np.where(np.log10(h5[gi][1]) == np.min(np.log10(h5[gi][1])))
	minlist.append(im)
	i2 = np.where(m == 2.0)
	i4 = np.where(m == -2.0)
	i1 = np.where(m == -1.0)
	i0 = np.where(m == 0.0)

	i = [i4, i0, i2]

	plt.subplot(4, 5, (5*gi) - 4)
	plt.plot(m, h5[gi][2], colo[gi]) #z2
	plt.plot(m, h5[gi][3], colo[gi]) #z1
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlabel('mu')
	plt.ylabel('z1, z2')

	plt.subplot(4, 5, (5*gi) - 3)
	plt.plot(m, np.log10(h5[gi][1]), colo[gi]) #sigma
	plt.plot(m, h5[gi][7], colo[gi]) #phi
	plt.xlabel('mu')
	plt.ylabel('phi, sigma')


	for ii in range (3):

		mu = m[i[ii]]
		plt.subplot(4, 5, (5*gi) - 2 + ii)
		gaus_m = (1 + (mu*h5[gi][4][i[ii]]))/h5[gi][6][i[ii]]
		gaus_sd = (h5[gi][1][i[ii]]*np.sqrt(h5[gi][5][i[ii]]))/h5[gi][6][i[ii]]
		x1 = (h5[gi][3][i[ii]]*gaus_sd) + gaus_m
		x2 = (h5[gi][2][i[ii]]*gaus_sd) + gaus_m
		x = np.linspace(x1, x2, 1001)
		plt.plot(x, np.exp((-(x-gaus_m)**2)/(2*gaus_sd**2))/(gaus_sd*np.sqrt(2*math.pi)), 'k')
		plt.xlim(0.5, 1.5)
		plt.ylim(0, 0.35)
		plt.xlabel('x')
		plt.ylabel('p(x)')

plt.savefig('muhist5.jpg')

###################################################

plt.figure(figsize=(12,12)) # z1, z2, phi, sigma, 3 histograms at various mu, a = 0.5
plt.title('a=2')

# gamma, sigma, z2, z1, M, q, help, phi

for gi in range (1, 5):

	#im = np.where(np.log10(h2[gi][1]) == np.min(np.log10(h2[gi][1])))
	#minlist.append(im)
	i2 = np.where(m == 2.0)
	i4 = np.where(m == -2.0)
	i1 = np.where(m == -1.0)
	i0 = np.where(m == 0.0)

	i = [i4, i0, i2]

	plt.subplot(4, 5, (5*gi) - 4)
	plt.plot(m, h2[gi][2], colo[gi]) #z2
	plt.plot(m, h2[gi][3], colo[gi]) #z1
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlabel('mu')
	plt.ylabel('z1, z2')

	plt.subplot(4, 5, (5*gi) - 3)
	plt.plot(m, np.log10(h2[gi][1]), colo[gi]) #sigma
	plt.plot(m, h2[gi][7], colo[gi]) #phi
	plt.xlabel('mu')
	plt.ylabel('phi, sigma')
	


	for ii in range (3):

		mu = m[i[ii]]
		plt.subplot(4, 5, (5*gi) - 2 + ii)
		gaus_m = (1 + (mu*h2[gi][4][i[ii]]))/h2[gi][6][i[ii]]
		gaus_sd = (h2[gi][1][i[ii]]*np.sqrt(h2[gi][5][i[ii]]))/h2[gi][6][i[ii]]
		x1 = (h2[gi][3][i[ii]]*gaus_sd) + gaus_m
		x2 = (h2[gi][2][i[ii]]*gaus_sd) + gaus_m
		x = np.linspace(x1, x2, 1001)
		plt.plot(x, np.exp((-(x-gaus_m)**2)/(2*gaus_sd**2))/(gaus_sd*np.sqrt(2*math.pi)), 'k')
		plt.xlim(0, 3)
		plt.ylim(0, 0.35)
		plt.xlabel('x')
		plt.ylabel('p(x)')

plt.savefig('muhist2.jpg')

#######################################################

plt.figure(4) # sigma, phi, z1, z2 for a = 2 and a = 0.5
plt.subplot(4, 2, 1)
for gi in range (1, 5):
	plt.plot(m, np.log10(h2[gi][1]), colo[gi]) #sigma
	plt.xlabel('mu')
	plt.ylabel('sigma')
	plt.title('a=2')
plt.subplot(4, 2, 3)
for gi in range (1, 5):
	plt.plot(m, h2[gi][7], colo[gi]) #phi
	plt.xlabel('mu')
	plt.ylabel('phi')
plt.subplot(4, 2, 2)
for gi in range (1, 5):
	plt.plot(m, np.log10(h5[gi][1]), colo[gi]) #sigma
	plt.xlabel('mu')
	plt.ylabel('sigma')
	plt.title('a=0.5')
plt.subplot(4, 2, 4)
for gi in range (1, 5):
	plt.plot(m, h5[gi][7], colo[gi]) #phi
	plt.xlabel('mu')
	plt.ylabel('phi')
plt.subplot(4, 2, 5)
for gi in range (1, 5):
	plt.plot(m, h2[gi][3], colo[gi]) #z1
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlabel('mu')
	plt.ylabel('z1')
plt.subplot(4, 2, 7)
for gi in range (1, 5):
	plt.plot(m, h2[gi][2], colo[gi]) #z2
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlabel('mu')
	plt.ylabel('z2')
plt.subplot(4, 2, 6)
for gi in range (1, 5):
	plt.plot(m, h5[gi][3], colo[gi]) #z1
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlabel('mu')
	plt.ylabel('z1')
plt.subplot(4, 2, 8)
for gi in range (1, 5):
	plt.plot(m, h5[gi][2], colo[gi]) #z2
	plt.axhline(y = 0, linestyle = '--', color = 'k')
	plt.xlabel('mu')
	plt.ylabel('z2')
	
plt.savefig('muhist.jpg')

######################################################################


plt.show()



