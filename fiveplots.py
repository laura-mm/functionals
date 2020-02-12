# plots for phi, M, diversity, d, h
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
grid = 20

crit = np.loadtxt('crit_20_0.txt', delimiter=",")
crih = np.loadtxt('crit_5_0.txt', delimiter=",")

# g is for a = 2
g0 = np.loadtxt('measfp_20_0_0.txt', delimiter=",")
g0 = g0.reshape(len(g0)/9, 9)
g0 = np.transpose(g0)
g1 = np.loadtxt('measfp_20_0_1.txt', delimiter=",")
g1 = g1.reshape(len(g1)/9, 9)
g1 = np.transpose(g1)
g2 = np.loadtxt('measfp_20_0_2.txt', delimiter=",")
g2 = g2.reshape(len(g2)/9, 9)
g2 = np.transpose(g2)
g3 = np.loadtxt('measfp_20_0_3.txt', delimiter=",")
g3 = g3.reshape(len(g3)/9, 9)
g3 = np.transpose(g3)
g4 = np.loadtxt('measfp_20_0_4.txt', delimiter=",")
g4 = g4.reshape(len(g4)/9, 9)
g4 = np.transpose(g4)

# h is for a = 0.5
h0 = np.loadtxt('measfp_5_0_0.txt', delimiter=",")
h0 = h0.reshape(len(h0)/9, 9)
h0 = np.transpose(h0)
h1 = np.loadtxt('measfp_5_0_1.txt', delimiter=",")
h1 = h1.reshape(len(h1)/9, 9)
h1 = np.transpose(h1)
h2 = np.loadtxt('measfp_5_0_2.txt', delimiter=",")
h2 = h2.reshape(len(h2)/9, 9)
h2 = np.transpose(h2)
h3 = np.loadtxt('measfp_5_0_3.txt', delimiter=",")
h3 = h3.reshape(len(h3)/9, 9)
h3 = np.transpose(h3)
h4 = np.loadtxt('measfp_5_0_4.txt', delimiter=",")
h4 = h4.reshape(len(h4)/9, 9)
h4 = np.transpose(h4)

"""
measp2 = np.loadtxt('meas22_20_0.txt', delimiter=",")
measp2 = measp2.reshape(5, grid + 1, 8)
measp2 = np.transpose(measp2, (0, 2, 1)) # gam, meas, sig
measp5 = np.loadtxt('meas22_5_0.txt', delimiter=",")
measp5 = measp5.reshape(5, grid + 1, 8)
measp5 = np.transpose(measp5, (0, 2, 1)) # gam, meas, sig
meash2 = np.loadtxt('meas33_20_0.txt', delimiter=",")
meash2 = meash2.reshape(5, grid + 1, 8)
meash2 = np.transpose(meash2, (0, 2, 1)) # gam, meas, sig
meash5 = np.loadtxt('meas33_5_0.txt', delimiter=",")
meash5 = meash5.reshape(5, grid + 1, 8)
meash5 = np.transpose(meash5, (0, 2, 1)) # gam, meas, sig
"""

# now for version 2:
measp2 = np.loadtxt('measpV2_20_0.txt', delimiter=",")
measp2 = measp2.reshape(5, grid + 1, 8)
measp2 = np.transpose(measp2, (0, 2, 1)) # gam, meas, sig
measp5 = np.loadtxt('measpV2_5_0.txt', delimiter=",")
measp5 = measp5.reshape(5, grid + 1, 8)
measp5 = np.transpose(measp5, (0, 2, 1)) # gam, meas, sig
meash2 = np.loadtxt('meashV2_20_0.txt', delimiter=",")
meash2 = meash2.reshape(5, grid + 1, 8)
meash2 = np.transpose(meash2, (0, 2, 1)) # gam, meas, sig
meash5 = np.loadtxt('meashV2_5_0.txt', delimiter=",")
meash5 = meash5.reshape(5, grid + 1, 8)
meash5 = np.transpose(meash5, (0, 2, 1)) # gam, meas, sig


g = [g0, g1, g2, g3, g4]
h = [h0, h1, h2, h3, h4]
col = ['r', 'm', 'b', 'c', 'g']
colo = ['r.', 'm.', 'b.', 'c.', 'g.']
lab = ['$\gamma = -1$', '$\gamma = -0.5$', '$\gamma = 0$', '$\gamma = 0.5$', '$\gamma = 1$']
mark = ['o', 's', '*', 'x', '+']
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']
x = np.linspace(-1, 1, grid + 1)
s = np.power(10, x)

# sigma, z1, z2, top, middle, bottom, M, q, help for g

# top, middle, bottom, M, q, diversity, dsq, h, for meas

# piecewise phi, M, diversity
###############################################################

plt.subplots(3, 2, figsize=(12,12))

plt.subplot(3, 2, 1) # phi, a = 2
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][4], col[gi])
	plt.semilogx(s, measp2[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$a = 2$')

plt.subplot(3, 2, 3) # M, a = 2
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6], col[gi])
	plt.semilogx(s, measp2[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 2, 5) # diversity, a = 2
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6]*g[gi][6]/g[gi][7], col[gi])
	plt.semilogx(s, measp2[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
plt.text(6, 0.9, lett[4], fontsize=20)
plt.legend(bbox_to_anchor=(0, -0.4, 2.17, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

plt.subplot(3, 2, 2) # phi, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][4], col[gi])
	plt.semilogx(s, measp5[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylim([-0.2, 1.01])
plt.text(6, 0.85, lett[1], fontsize=20)
plt.title(r'$a = 0.5$')

plt.subplot(3, 2, 4) # M, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6], col[gi])
	plt.semilogx(s, measp5[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylim([0.2, 2])
plt.text(6, 1.8, lett[3], fontsize=20)

plt.subplot(3, 2, 6) # diversity, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6]*h[gi][6]/h[gi][7], col[gi])
	plt.semilogx(s, measp5[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylim([0.3, 1])
plt.text(6, 0.9, lett[5], fontsize=20)

plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.97, wspace = 0.17)
plt.savefig('five2V2.pdf')

#end of piecewise m phi diversity plots
##################################################################

#section for piecewise, d, h
#################################################################

plt.subplots(2, 2, figsize=(12,9))

plt.subplot(2, 2, 1) # d, a = 2
for gi in range(5):
	plt.semilogx(s, measp2[gi][6], col[gi], marker = mark[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel('$d$')
plt.text(6, 0.55, lett[0], fontsize=20)
#plt.ylim([0, 0.5])
plt.title(r'$a = 2$')

plt.subplot(2, 2, 3) # h, a = 2
for gi in range(5):
	plt.semilogx(s, measp2[gi][7], col[gi], label = lab[gi], marker = mark[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel('$h$')
plt.ylim([0, 0.03])
plt.text(6, 0.0275, lett[2], fontsize=20)
plt.legend(bbox_to_anchor=(0, -0.3, 2.21, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

plt.subplot(2, 2, 2) # d, a = 0.5
for gi in range(5):
	plt.semilogx(s, measp5[gi][6], col[gi], marker = mark[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.text(6, 0.055, lett[1], fontsize=20)
#plt.ylim([0, 0.05])
plt.title(r'$a = 0.5$')

plt.subplot(2, 2, 4) # h, a = 0.5
for gi in range(5):
	plt.semilogx(s, measp5[gi][7], col[gi], marker = mark[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylim([0, 0.0005])
plt.text(6, 0.00045, lett[3], fontsize=20)

plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.96, wspace = 0.21, hspace = 0.16)
plt.savefig('five2dhV2.pdf')

# end of section peicewise d h
##############################################################################

#section for holling, M diversity
############################################################################

plt.subplots(2, 2, figsize=(12,9))

plt.subplot(2, 2, 1) # M, a = 2
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6], col[gi])
	plt.semilogx(s, meash2[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
plt.text(6, 1.8, lett[0], fontsize=20)
plt.title(r'$a = 2$')

plt.subplot(2, 2, 3) # diversity, a = 2
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6]*g[gi][6]/g[gi][7], col[gi])
	plt.semilogx(s, meash2[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
plt.text(6, 0.9, lett[2], fontsize=20)
plt.legend(bbox_to_anchor=(0, -0.3, 2.17, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

plt.subplot(2, 2, 2) # M, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6], col[gi])
	plt.semilogx(s, meash5[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
plt.xlim([10**(-1), 10**1])
plt.ylim([0.2, 2])
plt.text(6, 1.8, lett[1], fontsize=20)
plt.title(r'$a = 0.5$')

plt.subplot(2, 2, 4) # diversity, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6]*h[gi][6]/h[gi][7], col[gi])
	plt.semilogx(s, meash5[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylim([0.3, 1])
plt.text(6, 0.9, lett[3], fontsize=20)

plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.96, wspace = 0.17)
plt.savefig('five3V2.pdf')

#end of holling M diversity
####################################################################

#section for holling, d h
#####################################################################

plt.subplots(2, 2, figsize=(12,9))

plt.subplot(2, 2, 1) # d, a = 2
for gi in range(5):
	plt.semilogx(s, meash2[gi][6], col[gi], marker = mark[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel('$d$')
plt.text(6, 0.35, lett[0], fontsize=20)
#plt.ylim([0, 0.5])
plt.title(r'$a = 2$')

plt.subplot(2, 2, 3) # h, a = 2
for gi in range(5):
	plt.semilogx(s, meash2[gi][7], col[gi], label = lab[gi], marker = mark[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel('$h$')
plt.ylim([0, 0.01])
plt.text(6, 0.009, lett[2], fontsize=20)
plt.legend(bbox_to_anchor=(0, -0.3, 2.25, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

plt.subplot(2, 2, 2) # d, a = 0.5
for gi in range(5):
	plt.semilogx(s, meash5[gi][6], col[gi], marker = mark[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.text(6, 0.035, lett[1], fontsize=20)
#plt.ylim([0, 0.05])
plt.title(r'$a = 0.5$')

plt.subplot(2, 2, 4) # h, a = 0.5
for gi in range(5):
	plt.semilogx(s, meash5[gi][7], col[gi], marker = mark[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylim([0, 0.00007])
plt.text(6, 0.00006, lett[3], fontsize=20)

plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.96, wspace = 0.25, hspace = 0.18)
plt.savefig('five3dhV2.pdf')

# end of section holling d h
##############################################################

# end of plots in paper
#############################################################

#correction bit
###########################################################

meas2p = np.loadtxt('measp_20_0.txt', delimiter=",")
meas2p = meas2p.reshape(5, grid + 1, 10)
meas2p = np.transpose(meas2p, (0, 2, 1)) # gam, meas, sig
meas5p = np.loadtxt('measp_5_0.txt', delimiter=",")
meas5p = meas5p.reshape(5, grid + 1, 10)
meas5p = np.transpose(meas5p, (0, 2, 1)) # gam, meas, sig
meas2h = np.loadtxt('meash_20_0.txt', delimiter=",")
meas2h = meas2h.reshape(5, grid + 1, 10)
meas2h = np.transpose(meas2h, (0, 2, 1)) # gam, meas, sig
meas5h = np.loadtxt('meash_5_0.txt', delimiter=",")
meas5h = meas5h.reshape(5, grid + 1, 10)
meas5h = np.transpose(meas5h, (0, 2, 1)) # gam, meas, sig


#section for piecewise, d h, after correction
#############################################################

plt.subplots(2, 2, figsize=(12,9))

plt.subplot(2, 2, 1) # d, a = 2
for gi in range(5):
	#plt.semilogx(s, meas2p[gi][6], col[gi], marker = mark[gi])
	#plt.semilogx(s, meas2p[gi][8], col[gi], label = lab[gi], marker = '^')
	plt.semilogx(s, meas2p[gi][9], col[gi], label = lab[gi], marker = '>')
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel('$d$')
plt.text(6, 0.55, lett[0], fontsize=20)
plt.ylim([0, 0.5])
plt.title(r'$a = 2$')

plt.subplot(2, 2, 3) # h, a = 2
for gi in range(5):
	plt.semilogx(s, meas2p[gi][7], col[gi], label = lab[gi], marker = mark[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel('$h$')
plt.ylim([0, 0.025])
plt.text(6, 0.0225, lett[2], fontsize=20)

plt.subplot(2, 2, 2) # d, a = 0.5
for gi in range(5):
	#plt.semilogx(s, meas5p[gi][6], col[gi], marker = mark[gi])
	#plt.semilogx(s, meas5p[gi][8], col[gi], marker = '^')
	plt.semilogx(s, meas5p[gi][9], col[gi], marker = '>')
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.text(6, 0.055, lett[1], fontsize=20)
plt.ylim([0, 0.05])
plt.title(r'$a = 0.5$')

plt.subplot(2, 2, 4) # h, a = 0.5
for gi in range(5):
	plt.semilogx(s, meas5p[gi][7], col[gi], marker = mark[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
#plt.ylabel('$h$')
#plt.ylim([0, 0.0005])
plt.text(6, 0.00045, lett[3], fontsize=20)

plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.07, top = 0.96, wspace = 0.21, hspace = 0.16)
plt.savefig('pp9.jpg')

# end of piecewise, d h, after correction
#################################################################

#section for holling, dh, after correction
#####################################################################

plt.subplots(2, 2, figsize=(12,9))

plt.subplot(2, 2, 1) # d, a = 2
for gi in range(5):
	#plt.semilogx(s, meas2h[gi][6], col[gi], marker = mark[gi])
	#plt.semilogx(s, meas2h[gi][8], col[gi], label = lab[gi], marker = '^')
	plt.semilogx(s, meas2h[gi][9], col[gi], label = lab[gi], marker = '>')
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel('$d$')
plt.text(6, 0.45, lett[0], fontsize=20)
plt.ylim([0, 0.5])
plt.title(r'$a = 2$')

plt.subplot(2, 2, 3) # h, a = 2
for gi in range(5):
	plt.semilogx(s, meas2h[gi][7], col[gi], label = lab[gi], marker = mark[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
#plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel('$h$')
#plt.ylim([0, 0.012])
plt.text(6, 0.011, lett[2], fontsize=20)

plt.subplot(2, 2, 2) # d, a = 0.5
for gi in range(5):
	#plt.semilogx(s, meas5h[gi][6], col[gi], marker = mark[gi])
	#plt.semilogx(s, meas5h[gi][8], col[gi], marker = '^')
	plt.semilogx(s, meas5h[gi][9], col[gi], marker = '>')
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.text(6, 0.0325, lett[1], fontsize=20)
plt.ylim([0, 0.05])
plt.title(r'$a = 0.5$')

plt.subplot(2, 2, 4) # h, a = 0.5
for gi in range(5):
	plt.semilogx(s, meas5h[gi][7], col[gi], marker = mark[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
#plt.ylabel('$h$')
plt.ylim([0, 0.00014])
plt.text(6, 0.000125, lett[3], fontsize=20)

plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.07, top = 0.96, wspace = 0.25, hspace = 0.18)
plt.savefig('hh9.jpg')

# end of holling, dh, after correction
#########################################################################



# varying T with piecewise function
##################################################################

#meas(a)T(Tpower)

meas2T5 = np.loadtxt('measpT_20_200000.txt', delimiter=",")
meas2T5 = meas2T5.reshape(5, grid + 1, 8)
meas2T5 = np.transpose(meas2T5, (0, 2, 1)) # gam, meas, sig

meas2T4 = np.loadtxt('measpT_20_20000.txt', delimiter=",")
meas2T4 = meas2T4.reshape(5, grid + 1, 8)
meas2T4 = np.transpose(meas2T4, (0, 2, 1)) # gam, meas, sig

meas2T3 = np.loadtxt('measpT_20_2000.txt', delimiter=",")
meas2T3 = meas2T3.reshape(5, grid + 1, 8)
meas2T3 = np.transpose(meas2T3, (0, 2, 1)) # gam, meas, sig

meas2T2 = np.loadtxt('measpT_20_200.txt', delimiter=",")
meas2T2 = meas2T2.reshape(5, grid + 1, 8)
meas2T2 = np.transpose(meas2T2, (0, 2, 1)) # gam, meas, sig

meas5T5 = np.loadtxt('measpT_5_200000.txt', delimiter=",")
meas5T5 = meas5T5.reshape(5, grid + 1, 8)
meas5T5 = np.transpose(meas5T5, (0, 2, 1)) # gam, meas, sig

meas5T4 = np.loadtxt('measpT_5_20000.txt', delimiter=",")
meas5T4 = meas5T4.reshape(5, grid + 1, 8)
meas5T4 = np.transpose(meas5T4, (0, 2, 1)) # gam, meas, sig

meas5T3 = np.loadtxt('measpT_5_2000.txt', delimiter=",")
meas5T3 = meas5T3.reshape(5, grid + 1, 8)
meas5T3 = np.transpose(meas5T3, (0, 2, 1)) # gam, meas, sig

meas5T2 = np.loadtxt('measpT_5_200.txt', delimiter=",")
meas5T2 = meas5T2.reshape(5, grid + 1, 8)
meas5T2 = np.transpose(meas5T2, (0, 2, 1)) # gam, meas, sig


# a = 2
###################################################


plt.subplots(3, 4, figsize=(20,20))
plt.title(r'$a = 2$')

plt.subplot(3, 4, 1) # phi, a = 2, N = 200000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][4], col[gi])
	plt.semilogx(s, meas2T5[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^5$')

plt.subplot(3, 4, 5) # M, a = 2, N = 200000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6], col[gi])
	plt.semilogx(s, meas2T5[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 9) # diversity, a = 2, N = 200000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6]*g[gi][6]/g[gi][7], col[gi])
	plt.semilogx(s, meas2T5[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)
plt.legend(bbox_to_anchor=(0, -0.4, 4.5, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

plt.subplot(3, 4, 2) # phi, a = 2, N = 20000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][4], col[gi])
	plt.semilogx(s, meas2T4[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^4$')

plt.subplot(3, 4, 6) # M, a = 2, N = 20000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6], col[gi])
	plt.semilogx(s, meas2T4[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 10) # diversity, a = 2, N = 20000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6]*g[gi][6]/g[gi][7], col[gi])
	plt.semilogx(s, meas2T4[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)

plt.subplot(3, 4, 3) # phi, a = 2, N = 2000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][4], col[gi])
	plt.semilogx(s, meas2T3[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^3$')

plt.subplot(3, 4, 7) # M, a = 2, N = 2000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6], col[gi])
	plt.semilogx(s, meas2T3[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 11) # diversity, a = 2, N = 2000
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6]*g[gi][6]/g[gi][7], col[gi])
	plt.semilogx(s, meas2T3[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)

plt.subplot(3, 4, 4) # phi, a = 2, N = 200
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][4], col[gi])
	plt.semilogx(s, meas2T2[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^2$')

plt.subplot(3, 4, 8) # M, a = 2, N = 200
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6], col[gi])
	plt.semilogx(s, meas2T2[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 12) # diversity, a = 2, N = 200
for gi in range(5):
	plt.semilogx(g[gi][0], g[gi][6]*g[gi][6]/g[gi][7], col[gi])
	plt.semilogx(s, meas2T2[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crit[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)

#plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.97, wspace = 0.17)
plt.savefig('fiveTa2.pdf')
plt.savefig('fiveTa2.jpg')


# a = 0.5
###################################################


plt.subplots(3, 4, figsize=(20,20))
plt.title(r'$a = 0.5$')

plt.subplot(3, 4, 1) # phi, a = 0.5, N = 200000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][4], col[gi])
	plt.semilogx(s, meas5T5[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^5$')

plt.subplot(3, 4, 5) # M, a = 0.5, N = 200000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6], col[gi])
	plt.semilogx(s, meas5T5[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 9) # diversity, a = 0.5, N = 200000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6]*h[gi][6]/h[gi][7], col[gi])
	plt.semilogx(s, meas5T5[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)
plt.legend(bbox_to_anchor=(0, -0.4, 4.5, 0.2), loc='lower left', ncol=5, mode="expand", borderaxespad=0)

plt.subplot(3, 4, 2) # phi, a = 0.5, N = 20000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][4], col[gi])
	plt.semilogx(s, meas5T4[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^4$')

plt.subplot(3, 4, 6) # M, a = 0.5, N = 20000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6], col[gi])
	plt.semilogx(s, meas5T4[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 10) # diversity, a = 0.5, N = 20000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6]*h[gi][6]/h[gi][7], col[gi])
	plt.semilogx(s, meas5T4[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)

plt.subplot(3, 4, 3) # phi, a = 0.5, N = 2000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][4], col[gi])
	plt.semilogx(s, meas5T3[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^3$')

plt.subplot(3, 4, 7) # M, a = 0.5, N = 2000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6], col[gi])
	plt.semilogx(s, meas5T3[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 11) # diversity, a = 0.5, N = 2000
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6]*h[gi][6]/h[gi][7], col[gi])
	plt.semilogx(s, meas5T3[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)

plt.subplot(3, 4, 4) # phi, a = 0.5, N = 200
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][4], col[gi])
	plt.semilogx(s, meas5T2[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$\phi$')
plt.ylim([-0.2, 1.01])
#plt.text(6, 0.85, lett[0], fontsize=20)
plt.title(r'$T = 2 \times 10^2$')

plt.subplot(3, 4, 8) # M, a = 0.5, N = 200
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6], col[gi])
	plt.semilogx(s, meas5T2[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylabel(r'$M^*$')
plt.ylim([0.2, 2])
#plt.text(6, 1.8, lett[2], fontsize=20)

plt.subplot(3, 4, 12) # diversity, a = 0.5, N = 200
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6]*h[gi][6]/h[gi][7], col[gi])
	plt.semilogx(s, meas5T2[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\text{diversity }\frac{M^2}{q}$')
plt.ylim([0.3, 1])
#plt.text(6, 0.9, lett[4], fontsize=20)

#plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.97, wspace = 0.17)
plt.savefig('fiveTa5.pdf')
plt.savefig('fiveTa5.jpg')

# end of varying T
#############################################







"""



plt.subplot(3, 2, 2) # phi, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][4], col[gi])
	plt.semilogx(s, measp5[gi][1], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylim([-0.2, 1.01])
plt.text(6, 0.85, lett[1], fontsize=20)
plt.title(r'$a = 0.5$')

plt.subplot(3, 2, 4) # M, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6], col[gi])
	plt.semilogx(s, measp5[gi][3], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.ylim([0.2, 2])
plt.text(6, 1.8, lett[3], fontsize=20)

plt.subplot(3, 2, 6) # diversity, a = 0.5
for gi in range(5):
	plt.semilogx(h[gi][0], h[gi][6]*h[gi][6]/h[gi][7], col[gi])
	plt.semilogx(s, measp5[gi][5], colo[gi], marker = mark[gi], label = lab[gi])
	plt.axvline(x = crih[gi], linestyle = '--', color = col[gi])
plt.xlim([10**(-1), 10**1])
plt.xlabel(r'$\sigma$')
plt.ylim([0.3, 1])
plt.text(6, 0.9, lett[5], fontsize=20)

plt.subplots_adjust(left=0.09, right = 0.97, bottom = 0.13, top = 0.97, wspace = 0.17)
plt.savefig('five2V2.pdf')
"""



plt.show()

