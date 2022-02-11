'''
This script is used to produce several plots regarding the HNL decay (lifetime, BRs)
In particular it was used to perform the Bondarenko vs Peskin studies shown in 
https://indico.cern.ch/event/978280/contributions/4120549/attachments/2148757/3653052/2020_11_26_Bpark_HNL.pdf
'''

from decays import *
from old_common import *

import matplotlib.pyplot as plt



#masses = [0.5+i*0.1 for i in range(0,60) ] # GeV
masses = [float('{:.2f}'.format(1.+i*0.1)) for i in range(0,40) ] # GeV

ctaus_peskin_lep = [ctau_from_gamma(gamma_lep(mass=m,vv=1)) for m in masses] 
ctaus_peskin_had = [ctau_from_gamma(gamma_had(mass=m,vv=1)) for m in masses]   
ctaus_peskin_neu = [ctau_from_gamma(gamma_neu(mass=m,vv=1)) for m in masses]   
ctaus_peskin_tot = [ctau_from_gamma(gamma_total(mass=m,vv=1)) for m in masses]    

ctaus_peskin = [ctaus_peskin_lep, ctaus_peskin_had, ctaus_peskin_neu, ctaus_peskin_tot]

decays = [HNLDecays(m,1.) for m in masses]

# factor 2 is because Peskin considers Majorana
ctaus_bondarenko_lep = [ctau_from_gamma(2*decays[im].decay_rate['tot_lep']) for im,m in enumerate(masses)]
ctaus_bondarenko_had = [ctau_from_gamma(2*decays[im].decay_rate['tot_had']) for im,m in enumerate(masses)]
ctaus_bondarenko_neu = [ctau_from_gamma(2*decays[im].decay_rate['tot_neu']) for im,m in enumerate(masses)]
ctaus_bondarenko_tot = [ctau_from_gamma(2*decays[im].decay_rate['tot']    ) for im,m in enumerate(masses)]

ctaus_bondarenko = [ctaus_bondarenko_lep, ctaus_bondarenko_had, ctaus_bondarenko_neu, ctaus_bondarenko_tot]

labels = ['Leptonic', 'Hadronic', 'Neutrinos', 'Total']


####
## Partial widths, comparison between Peskin and Bondarenko, as a function of HNL mass
## slide 6 
####

### Log Plot
fig, axes = plt.subplots(1, 4, figsize=(5.5*4, 5))
for i,ax in enumerate(axes):
  ax.plot(masses, ctaus_peskin[i], 'b', label='{} (Peskin)'.format(labels[i]))
  ax.plot(masses, ctaus_bondarenko[i], 'r', label='{} (Bondarenko)'.format(labels[i]))
  ax.set_ylabel('c$\\tau$ (mm)')
  ax.set_xlabel('HNL mass (GeV)')
  ax.legend(loc='upper right', frameon=True) 
  ax.grid(which='both', axis='both')
  ax.set_yscale('log')
 

### Lin Plot
fig.savefig('./plots/HNLlifetime_check_log.pdf')
fig.savefig('./plots/HNLlifetime_check_log.png')

fig, axes = plt.subplots(1, 4, figsize=(5.5*4, 5))
for i,ax in enumerate(axes):
  ax.plot(masses, ctaus_peskin[i], 'b', label='{} (Peskin)'.format(labels[i]))
  ax.plot(masses, ctaus_bondarenko[i], 'r', label='{} (Bondarenko)'.format(labels[i]))
  ax.set_ylabel('c$\\tau$ (mm)')
  ax.set_xlabel('HNL mass (GeV)')
  ax.legend(loc='upper right', frameon=True) 
  ax.grid(which='both', axis='both')
  #ax.set_yscale('log')
 
fig.savefig('./plots/HNLlifetime_check_lin.pdf')
fig.savefig('./plots/HNLlifetime_check_lin.png')


### Ratio plot
ctaus_ratios = []
for pes,bon in zip(ctaus_peskin,ctaus_bondarenko):
  ctaus_ratio = []
  for im,mass in enumerate(masses):
    ctaus_ratio.append(pes[im]/bon[im])
  ctaus_ratios.append(ctaus_ratio)

fig, axes = plt.subplots(1, 4, figsize=(5.5*4, 5))
for i,ax in enumerate(axes):
  ax.plot(masses, ctaus_ratios[i], 'r', label='{}'.format(labels[i]))
  ax.set_ylabel('Ratio of c$\\tau$ (mm)')
  ax.set_xlabel('HNL mass (GeV)')
  ax.set_ylim(0.4,1.4)
  ax.legend(loc='upper right', frameon=True) 
  ax.grid(which='both', axis='both')
  #ax.set_yscale('log')
fig.suptitle('Ratio Peskin/Bondarenko', fontsize=14) 
fig.savefig('./plots/HNLlifetime_ratio_lin.pdf')
fig.savefig('./plots/HNLlifetime_ratio_lin.png')


####
## Branching ratio composition for different processes, as a function of HNL mass
## Figure 13 right of Bondarenko
## Slide 7,
####
gammas_bondarenko_lep = [2*decays[im].decay_rate['tot_lep'] for im,m in enumerate(masses)]
gammas_bondarenko_had = [2*decays[im].decay_rate['tot_had'] for im,m in enumerate(masses)]
gammas_bondarenko_neu = [2*decays[im].decay_rate['tot_neu'] for im,m in enumerate(masses)]
gammas_bondarenko_tot = [2*decays[im].decay_rate['tot']     for im,m in enumerate(masses)]

gammas_bondarenko_lepFrac = [x_lep / x_tot for x_lep,x_tot in zip(gammas_bondarenko_lep,gammas_bondarenko_tot)]
gammas_bondarenko_hadFrac = [x_had / x_tot for x_had,x_tot in zip(gammas_bondarenko_had,gammas_bondarenko_tot)]
gammas_bondarenko_neuFrac = [x_neu / x_tot for x_neu,x_tot in zip(gammas_bondarenko_neu,gammas_bondarenko_tot)]

fig, ax = plt.subplots(1, 1, figsize=(5.5*1, 5))
ax.plot(masses, gammas_bondarenko_lepFrac, 'r', label='leptons')
ax.plot(masses, gammas_bondarenko_hadFrac, 'b',  label='quarks')
ax.plot(masses, gammas_bondarenko_neuFrac, 'g',  label='neutrinos')
ax.set_ylabel('BR')
ax.set_xlabel('HNL mass (GeV)')
ax.set_ylim(0.1,1.0)
ax.set_xlim(1,5)
ax.set_yticks([0.1,0.2,0.4,0.6,0.8,1.0], minor=True)
ax.legend(loc='center right', frameon=True) 
ax.grid(which='both', axis='both')
ax.set_yscale('log')
ax.set_xscale('log')

fig.savefig('./plots/HNL_BR_bon_log.pdf')
fig.savefig('./plots/HNL_BR_bon_log.png')


####
## Lifetime of N as a function of HNL mass
## slide 7 right
####

taus_bondarenko_tot = [ctau / const_c * 0.001 for ctau in ctaus_bondarenko_tot ]

fig, ax = plt.subplots(1, 1, figsize=(5.5*1, 5))
ax.plot(masses, taus_bondarenko_tot, 'r', label='$V_e=V_{\\tau}=0$')
ax.set_ylabel('$\\tau$ (s) ')
ax.set_xlabel('HNL mass (GeV)')
ax.set_ylim(1e-16,1e-04)
ax.set_xlim(0.05,5)
#ax.set_yticks([0.1,0.2,0.4,0.6,0.8,1.0], minor=True)
ax.legend(loc='center right', frameon=False) 
ax.grid(which='both', axis='both')
ax.set_yscale('log')
ax.set_xscale('log')

fig.savefig('./plots/HNL_tau_bon_log.pdf')
fig.savefig('./plots/HNL_tau_bon_log.png')


####
## Width of N / Mass, as a function of mass
## To check narrow width approximation
####
fig, ax = plt.subplots(1, 1, figsize=(5.5*1, 5))
ax.plot(masses, gammas_bondarenko_tot_Over_mass, 'r', label='bondarenko')
ax.set_ylabel('$\\Gamma/m $')
ax.set_xlabel('HNL mass (GeV)')
ax.grid(which='both', axis='both')
ax.set_yscale('log')
fig.savefig('./plots/HNL_gamma_over_mass_log.pdf')
fig.savefig('./plots/HNL_gamma_over_mass_log.png')


####
## BR comparison between Peskin and Bondarenko, as a function of HNL mass 
## see slide 8
####
# compare BRs (N->mupi)
gammas_bondarenko_Nmupi = [2*decays[im].decay_rate['mupi'] for im,m in enumerate(masses)]
brs_bondarenko_Nmupi = [gp / gtot for gp,gtot in zip(gammas_bondarenko_Nmupi, gammas_bondarenko_tot)]
brs_peskin_Nmupi = [BR_HNLmupion(mass=m) for im,m in enumerate(masses)]
brs_peskin_Nmupi_corrected = [2*BR_HNLmupion(mass=m) for im,m in enumerate(masses)]
gammas_peskin_Nmupi_corrected = [2*gamma_partial(mass=m,vv=1) for im,m in enumerate(masses)]
ratios = [p/b for p,b in zip(brs_peskin_Nmupi,brs_bondarenko_Nmupi)]
ratios_corrected = [p/b for p,b in zip(brs_peskin_Nmupi_corrected,brs_bondarenko_Nmupi)]
ratios_gamma_corrected = [p/b for p,b in zip(gammas_peskin_Nmupi_corrected,gammas_bondarenko_Nmupi)]
fig, ax = plt.subplots(2, 1, figsize=(5.5*1, 5*2))
ax[0].plot(masses, brs_peskin_Nmupi_corrected, 'b', label='peskin')
#ax[0].plot(masses, brs_peskin_Nmupi, 'black', label='peskin (old)')
ax[0].plot(masses, brs_bondarenko_Nmupi, 'r', label='bondarenko')
#ax[0].plot(masses, gammas_peskin_Nmupi_corrected, 'b', label='peskin')
#ax[0].plot(masses, gammas_bondarenko_Nmupi, 'r', label='bondarenko')
ax[0].set_ylabel('BR$(N\\rightarrow\\mu\\pi)$')
ax[0].set_xlabel('HNL mass (GeV)')
ax[0].grid(which='both', axis='both')
ax[0].set_yscale('log')
ax[0].legend(loc='upper right', frameon=False) 

ax[1].plot(masses, ratios_corrected, 'b', label='P/B')
#ax[1].plot(masses, ratios          , 'black', label='P(old)/B')
#ax[1].plot(masses, ratios_gamma_corrected, 'b', label='P/B')
ax[1].set_ylabel('Ratio')
ax[1].set_xlabel('HNL mass (GeV)')
ax[1].grid(which='both', axis='both')
#ax[1].set_yscale('log')
ax[1].legend(loc='upper right', frameon=False) 
ax[1].set_ylim(0.2,1.4)

fig.savefig('./plots/HNL_BR.pdf')
fig.savefig('./plots/HNL_BR.png')




