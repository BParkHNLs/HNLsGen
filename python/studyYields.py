'''
This script was used to estimate the relative importance between two contributions: (B -> mu N X) and (B -> X D -> mu N X)
The calculations are purely based on BRs.
See the presentation for more details:
https://indico.cern.ch/event/1045239/contributions/4390951/attachments/2270279/3855507/2021_06_24_Bpark_BtoD.pdf
'''

import os
import sys
import numpy as np
from glob import glob
from array import array

from python.common import *
from python.decays import Decays


#from python.old_common import BR__B_l_N, BR__B_D0_l_N, BR__B_l_nu, BR__B_D0_l_nu, # BR_HNLmupion

# NOTE: whenever an upper limit is presented in the PDG, the maximal value is taken
# NOTE: For inclusive branching fractions, e.g., B => D anything, the values usually are multiplicities, not branching fractions. They can be greater than one.
#       what does this mean ?? And are these estimates impaired by this?

def getExpNevtsD(mass,which='BtoDX'): # 'BtoDX' , 'BtoDuvu'
  dec    = Decays(mass=mass, mixing_angle_square=1)


  tot_D0_to_HNL = dec.D0_to_pieHNL.BR + dec.D0_to_KeHNL.BR + dec.D0_to_piuHNL.BR + dec.D0_to_KuHNL.BR + dec.D0_to_KstareHNL.BR + dec.D0_to_KstaruHNL.BR
  tot_D0_to_uHNL = dec.D0_to_piuHNL.BR + dec.D0_to_KuHNL.BR + dec.D0_to_KstaruHNL.BR
  
  # note, contributions with a tau should be null for HNL mass = 1, 1.5 GeV (mass of the tau 1.78 GeV, mass of the D~1.9 GeV)
  tot_D_to_HNL  = dec.D_to_eHNL.BR +  dec.D_to_uHNL.BR +  dec.D_to_tHNL.BR + \
                  dec.D_to_K0eHNL.BR + dec.D_to_pi0eHNL.BR + dec.D_to_K0uHNL.BR + dec.D_to_pi0uHNL.BR + \
                  dec.D_to_K0stareHNL.BR + dec.D_to_K0staruHNL.BR
  tot_D_to_uHNL = dec.D0_to_KstaruHNL.BR + dec.D_to_K0uHNL.BR + dec.D_to_pi0uHNL.BR + dec.D_to_K0staruHNL.BR  

  tot_Ds_to_HNL = dec.Ds_to_eHNL.BR + dec.Ds_to_uHNL.BR + dec.Ds_to_tHNL.BR 
  tot_Ds_to_uHNL = dec.Ds_to_uHNL.BR

  ### B+/- 
  # https://pdglive.lbl.gov/Particle.action?init=0&node=S041&home=MXXX045#decayclump_J  
  # inclusive modes
  gamma38,gamma39 = 0.086,0.79     
  gamma40,gamma41 = 0.025,0.099     
  gamma42,gamma43 = 0.079,0.011

  # semileptonic modes
  gamma4,gamma5,gamma6,gamma7,gamma8,gamma11,gamma18 = 0.0235,0.077,0.0566,0.0188,0.044,0.0188,0.00061

  if which=='BtoDX':
    B_w =  (gamma38+gamma39) * tot_D0_to_uHNL + \
           (gamma40+gamma41) * tot_D_to_uHNL  + \
           (gamma42+gamma43) * tot_Ds_to_uHNL 

  elif which=='BtoDuvu':
    B_w =  (gamma4*2 + gamma5 + gamma6*2 + gamma7 ) * tot_D0_to_HNL + \
           (gamma8*2 + gamma11*2 )                  * tot_D_to_HNL + \
           (gamma18*2)                              * tot_Ds_to_HNL
  else:
    B_w = 0.

  ### B0  
  #https://pdglive.lbl.gov/Particle.action?init=0&node=S042&home=MXXX045#decayclump_I
  # inclusive
  gamma23,gamma24,gamma25,gamma26,gamma27,gamma28 = 0.081,0.474,0.039,0.369,0.103,0.026
  # semi-leptonic
  gamma4,gamma5,gamma6,gamma7,gamma8,gamma11 = 0.0231,0.0108,0.0506,0.0157,0.041,0.023  

  if which=='BtoDX':
    B0_w = (gamma23+gamma24)  * tot_D0_to_uHNL + \
           (gamma25+gamma26)  * tot_D_to_uHNL + \
           (gamma27+gamma28)  * tot_Ds_to_uHNL

  elif which=='BtoDuvu':
    B0_w = (gamma4*2 + gamma5 + gamma6*2 + gamma7) * tot_D_to_HNL + \
           (gamma8*2 + gamma11)                    * tot_D0_to_HNL 
           # no d_S contribution
  else:
    B0_w = 0.

  ### Bs
  # https://pdglive.lbl.gov/Particle.action?init=0&node=S086&home=MXXX046
  gamma1 = 0.93  
  gamma5,gamma6 = 0.0252,0.054 # Ds(*) u vu

  if which=='BtoDX':
    Bs_w = gamma1 * tot_Ds_to_uHNL
  elif which=='BtoDuvu':
    Bs_w = (gamma5 + gamma6) * tot_Ds_to_HNL # taking everything into account well is a mess
  else: 
    Bs_w = 0.

  ### Bc
  # https://pdglive.lbl.gov/Particle.action?init=0&node=S091&home=MXXX049#decayclump_A
  Bc_w = 0  # Bc->D negligible, and Bc contribution negligible at HNL masses: 1-1.8 GeV

  Bweights = [B_w, B0_w, Bs_w, Bc_w]
  Bfracs =   [0.4, 0.4, 0.1, 0.001]
  Bspecies = ['B', 'B0', 'Bs', 'Bc']
  expNevts = 0
  
  for ib,b in enumerate(Bspecies):
    expNevts        += Bfracs[ib] * Bweights[ib]

  return expNevts


def getExpNevtsB(mass):
 
  dec    = Decays(mass=mass, mixing_angle_square=1)

  B_w =     (dec.B_to_uHNL.BR+       \
             dec.B_to_D0uHNL.BR+     \
             dec.B_to_D0staruHNL.BR+ \
             dec.B_to_pi0uHNL.BR+    \
             dec.B_to_rho0uHNL.BR)   

  B0_w =    (dec.B0_to_DuHNL.BR+     \
             dec.B0_to_DstaruHNL.BR+ \
             dec.B0_to_piuHNL.BR+    \
             dec.B0_to_rhouHNL.BR)   
   
  Bs_w =    (dec.Bs_to_DsuHNL.BR+ \
             dec.Bs_to_DsstaruHNL.BR+ \
             dec.Bs_to_KuHNL.BR+ \
             dec.Bs_to_KstaruHNL.BR)  
  
  Bc_w =    (dec.Bc_to_uHNL.BR) 

  Bweights = [B_w, B0_w, Bs_w, Bc_w]
  Bfracs =   [0.4, 0.4, 0.1, 0.001]
  Bspecies = ['B', 'B0', 'Bs', 'Bc']
  expNevts = 0
  
  for ib,b in enumerate(Bspecies):
    expNevts        += Bfracs[ib] * Bweights[ib]

  return expNevts



if __name__ == "__main__":

  
  print 'Mass=0.5GeV'
  print 'B                => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsB(mass=0.5))
  print 'B => D X         => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=0.5, which='BtoDX'))
  print 'B => mu v_mu D X => l  HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=0.5, which='BtoDuvu'))
  print ''

  print 'Mass=1.0GeV'
  print 'B                => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsB(mass=1.0))
  print 'B => D X         => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=1.0, which='BtoDX'))
  print 'B => mu v_mu D X => l  HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=1.0, which='BtoDuvu'))
  print ''
  
  print 'Mass=1.5GeV'
  print 'B                => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsB(mass=1.5))
  print 'B => D X         => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=1.5, which='BtoDX'))
  print 'B => mu v_mu D X => l  HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=1.5, which='BtoDuvu'))
  print ''
  
  print 'Mass=1.8GeV'
  print 'B                => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsB(mass=1.8))
  print 'B => D X         => mu HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=1.8, which='BtoDX'))
  print 'B => mu v_mu D X => l  HNL X:   Expected = {:.3f}'.format(getExpNevtsD(mass=1.8, which='BtoDuvu'))
  print ''

#  # use the B->DX and subtract B->DuX
#  B_w_inc =  (gamma38+gamma39  -  (gamma4*2+gamma5+gamma6*2+gamma7)) * tot_D0_to_uHNL + \
#             (gamma40+gamma41  -  (gamma8*2 + gamma11*2           )) * tot_D_to_uHNL  + \
#             (gamma42+gamma43  -  (gamma18*2                      )) * tot_Ds_to_uHNL 

#  B0_w_inc = (gamma23+gamma24  - (gamma4*2+gamma5+gamma6*2+gamma7)) * tot_D0_to_uHNL + \
#             (gamma25+gamma26  - (gamma11*2))                       * tot_D_to_uHNL + \
#             (gamma27+gamma28)                                      * tot_Ds_to_uHNL

