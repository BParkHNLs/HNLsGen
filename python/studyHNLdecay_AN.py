'''
Script to plot N lifetime and N->lpi branching ratios, figures 5 and 6 of AN
'''


from decays import *
from common import *

import matplotlib.pyplot as plt
import pandas as pd

#plt.style.use('./mystyle.mplstyle')
#masses = [0.5+i*0.1 for i in range(0,60) ] # GeV

masses = [float('{:.2f}'.format(0.5+i*0.1)) for i in range(0,55) ] # GeV


# for historical reasons first couplings is muonic!
#decays_muonic = [HNLDecays(m,1.,0,0) for m in masses]
#ctaus_bondarenko_tot_muonic = [2*ctau_from_gamma(decays_muonic[im].decay_rate['tot']    ) for im,m in enumerate(masses)]

#decays_electronic = [HNLDecays(m,0.,1.,0) for m in masses]
#ctaus_bondarenko_tot_electronic = [2*ctau_from_gamma(decays_electronic[im].decay_rate['tot']    ) for im,m in enumerate(masses)]

decays_mixed = [HNLDecays(m,0.5,0.5,0) for m in masses]
#ctaus_bondarenko_tot_mixed = [2*ctau_from_gamma(decays_mixed[im].decay_rate['tot']    ) for im,m in enumerate(masses)]
ctaus_bondarenko_tot_mixed = [ctau_from_gamma(2*decays_mixed[im].decay_rate['tot']    ) for im,m in enumerate(masses)]
print ctaus_bondarenko_tot_mixed

# get data extracted from Bondarenko paper
file_bondarenko_quark = pd.read_csv('./bondarenko_data/bondarenko_quark.csv')
mass_bondarenko_quark = file_bondarenko_quark['mass']
BR_bondarenko_quark   = file_bondarenko_quark['BR']

file_bondarenko_lepton = pd.read_csv('./bondarenko_data/bondarenko_lepton.csv')
mass_bondarenko_lepton = file_bondarenko_lepton['mass']
BR_bondarenko_lepton   = file_bondarenko_lepton['BR']

file_bondarenko_invisible = pd.read_csv('./bondarenko_data/bondarenko_invisible.csv')
mass_bondarenko_invisible = file_bondarenko_invisible['mass']
BR_bondarenko_invisible   = file_bondarenko_invisible['BR']


####
## Total Ctau of N, as a function of mass
## Figure 6 left
###
#fig, ax = plt.subplots(1, 1, figsize=(8, 8))
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

fig.savefig('./plots/HNL_ctau.pdf')
fig.savefig('./plots/HNL_ctau.png')


####
## Branching ratio of N to lpi, as a function of mass
## Figure 6 right
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

fig.savefig('./plots/HNL_BR_comp.pdf')
fig.savefig('./plots/HNL_BR_comp.png')


#### 
## Branching ratio composition for different processes, as a function of mass 
## Figure 5 
####

masses = [float('{:.2f}'.format(1.+i*0.1)) for i in range(0,40) ] # GeV
decays_onethird = [HNLDecays(m,0.3333,0.3333,0.3333) for m in masses]

gammas_bondarenko_lep = [2*decays_onethird[im].decay_rate['tot_lep'] for im,m in enumerate(masses)]
gammas_bondarenko_had = [2*decays_onethird[im].decay_rate['tot_had'] for im,m in enumerate(masses)]
gammas_bondarenko_neu = [2*decays_onethird[im].decay_rate['tot_neu'] for im,m in enumerate(masses)]
gammas_bondarenko_tot = [2*decays_onethird[im].decay_rate['tot']     for im,m in enumerate(masses)]

gammas_bondarenko_lepFrac = [x_lep / x_tot for x_lep,x_tot in zip(gammas_bondarenko_lep,gammas_bondarenko_tot)]
gammas_bondarenko_hadFrac = [x_had / x_tot for x_had,x_tot in zip(gammas_bondarenko_had,gammas_bondarenko_tot)]
gammas_bondarenko_neuFrac = [x_neu / x_tot for x_neu,x_tot in zip(gammas_bondarenko_neu,gammas_bondarenko_tot)]

#fig, ax = plt.subplots(1, 1, figsize=(8,8))
fig, ax = plt.subplots(1, 1, figsize=(7, 5))
p1, = ax.plot(masses, gammas_bondarenko_lepFrac, 'r', label='leptons')
p2, = ax.plot(masses, gammas_bondarenko_hadFrac, 'b',  label='quarks')
p3, = ax.plot(masses, gammas_bondarenko_neuFrac, 'g',  label='neutrinos')

p4, = ax.plot(mass_bondarenko_quark, BR_bondarenko_quark, linewidth=0,  marker='*', markersize=7, color='black', label='Bondarenko et al.')
p5, = ax.plot(mass_bondarenko_lepton, BR_bondarenko_lepton, linewidth=0,  marker='*', markersize=7, color='black')#, label='Bondarenko et al.')
p6, = ax.plot(mass_bondarenko_invisible, BR_bondarenko_invisible, linewidth=0,  marker='*', markersize=7, color='black')#, label='Bondarenko et al.')
ax.legend(handles=[p1, p2, p3, p4],numpoints =1)

ax.set_ylabel('BR')
ax.set_xlabel('HNL mass (GeV)')
ax.set_ylim(0.1,1.0)
ax.set_xlim(1,5)
ax.set_yticks([0.1,0.2,0.4,0.6,0.8,1.0], minor=True)
#ax.legend(loc='center right', frameon=True) 
ax.grid(which='both', axis='both')
ax.set_yscale('log')
ax.set_xscale('log')

fig.savefig('./plots/HNL_BR_bon_log.pdf')
fig.savefig('./plots/HNL_BR_bon_log.png')


