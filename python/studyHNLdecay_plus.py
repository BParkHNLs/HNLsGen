'''
Plotting script to plot N lifetime and N->lpi branching ratio, figure 6 of AN
'''


from decays import *
from common import *

import matplotlib.pyplot as plt
plt.style.use('./mystyle.mplstyle')
#masses = [0.5+i*0.1 for i in range(0,60) ] # GeV

masses = [float('{:.2f}'.format(0.5+i*0.1)) for i in range(0,55) ] # GeV


# for historical reasons first couplings is muonic!
#decays_muonic = [HNLDecays(m,1.,0,0) for m in masses]
#ctaus_bondarenko_tot_muonic = [2*ctau_from_gamma(decays_muonic[im].decay_rate['tot']    ) for im,m in enumerate(masses)]

#decays_electronic = [HNLDecays(m,0.,1.,0) for m in masses]
#ctaus_bondarenko_tot_electronic = [2*ctau_from_gamma(decays_electronic[im].decay_rate['tot']    ) for im,m in enumerate(masses)]

decays_mixed = [HNLDecays(m,0.5,0.5,0) for m in masses]
ctaus_bondarenko_tot_mixed = [2*ctau_from_gamma(decays_mixed[im].decay_rate['tot']    ) for im,m in enumerate(masses)]
print ctaus_bondarenko_tot_mixed





####
## Total Ctau of N, as a function of mass
###
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
#ax.plot(masses, ctaus_bondarenko_tot_muonic,     'r', label='$V_e^2=0,  V_{\\mu}^2=1,  V_{\\tau}^2=0$')
#ax.plot(masses, ctaus_bondarenko_tot_electronic, 'g', label='$V_e^2=1,  V_{\\mu}^2=0,  V_{\\tau}^2=0$')
ax.plot(masses, ctaus_bondarenko_tot_mixed,      'b') # label='$V_e^2=0.5, V_{\\mu}^2=0.5, V_{\\tau}^2=0$')
ax.set_ylabel('$c\\tau$ (mm) ')
ax.set_xlabel('HNL mass (GeV)')
ax.text(0.22, 1.05, '$V_e^2=0.5, V_{\\mu}^2=0.5, V_{\\tau}^2=0$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
#ax.set_ylim(1e-7,1e-5)
ax.set_xlim(0.5,6.)
#ax.set_yticks([0.1,0.2,0.4,0.6,0.8,1.0], minor=True)
ax.legend(loc='upper right', frameon=True, fontsize=18) 
ax.grid(which='both', axis='both')
ax.set_yscale('log')
#ax.set_xscale('log')

fig.savefig('HNL_ctau.pdf')
fig.savefig('HNL_ctau.png')


####
## Branching ratio of N to lpi, as a function of mass
####
gammas_bondarenko_tot =   [decays_mixed[im].decay_rate['tot']  for im,m in enumerate(masses)]
gammas_bondarenko_Nmupi = [decays_mixed[im].decay_rate['mupi'] for im,m in enumerate(masses)]
gammas_bondarenko_Nelpi = [decays_mixed[im].decay_rate['elpi']  for im,m in enumerate(masses)]
brs_bondarenko_Nmupi = [gp / gtot for gp,gtot in zip(gammas_bondarenko_Nmupi, gammas_bondarenko_tot)]
brs_bondarenko_Nelpi = [gp / gtot for gp,gtot in zip(gammas_bondarenko_Nelpi, gammas_bondarenko_tot)]

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.plot(masses, brs_bondarenko_Nmupi, 'b', label='N $\\to \mu\pi$')
ax.plot(masses, brs_bondarenko_Nelpi, 'r', label='N $\\to e \pi$')
ax.set_ylabel('BR$(N\\rightarrow\\ell\\pi)$')
ax.set_xlabel('HNL mass (GeV)')
ax.text(0.22, 1.05, '$V_e^2=0.5, V_{\\mu}^2=0.5, V_{\\tau}^2=0$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
ax.grid(which='both', axis='both')
ax.set_yscale('log')
ax.set_xlim(0.5,6.)
ax.legend(loc='upper right', frameon=True, fontsize=18) 

fig.savefig('HNL_BR_comp.pdf')
fig.savefig('HNL_BR_comp.png')






