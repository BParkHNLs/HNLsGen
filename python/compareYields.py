import os
import sys
import numpy as np
from glob import glob
from array import array

from python.common import *
from python.decays import Decays

#from python.old_common import BR__B_l_N, BR__B_D0_l_N, BR__B_l_nu, BR__B_D0_l_nu, # BR_HNLmupion



def getExpNevtsD(mass):
  dec    = Decays(mass=mass, mixing_angle_square=1)


  tot_D0_to_HNL = dec.D0_to_pieHNL.BR + dec.D0_to_KeHNL.BR + dec.D0_to_piuHNL.BR + dec.D0_to_KuHNL.BR + dec.D0_to_KstareHNL.BR + dec.D0_to_KstaruHNL.BR
  tot_D0_to_uHNL = dec.D0_to_piuHNL.BR + dec.D0_to_KuHNL.BR + dec.D0_to_KstaruHNL.BR
  
  tot_D_to_HNL  = dec.D_to_eHNL.BR +  dec.D_to_uHNL.BR +  dec.D_to_tHNL.BR + \
                  dec.D_to_K0eHNL.BR + dec.D_to_pi0eHNL.BR + dec.D_to_K0uHNL.BR + dec.D_to_pi0uHNL.BR + \
                  dec.D_to_K0stareHNL.BR + dec.D_to_K0staruHNL.BR
  tot_D_to_uHNL = dec.D0_to_KstaruHNL.BR + dec.D_to_K0uHNL.BR + dec.D_to_pi0uHNL.BR + dec.D_to_K0staruHNL.BR  

  tot_Ds_to_HNL = dec.Ds_to_eHNL.BR + dec.Ds_to_uHNL.BR + dec.Ds_to_tHNL.BR 
  tot_Ds_to_uHNL = dec.Ds_to_uHNL.BR

  # B+/-
  SM_B_to_D0X_BR =   0.086 + 0.79     # incluive modes numbers taken from https://pdglive.lbl.gov/Particle.action?init=0&node=S041&home=MXXX045#decayclump_J
  SM_B_to_DX_BR =    0.025 + 0.099    # 
  SM_B_to_DsX_BR =   0.079 + 0.011

  B_w_inc =  SM_B_to_D0X_BR * tot_D0_to_uHNL + \
             SM_B_to_DX_BR  * tot_D_to_uHNL  + \
             SM_B_to_DsX_BR * tot_Ds_to_uHNL 

  gamma4,gamma5,gamma6,gamma7,gamma8,gamma11,gamma18 = 0.0235,0.077,0.0566,0.0188,0.044,0.0188,0.00061

  B_w_lep =  (gamma4*2 + gamma5 + gamma6*2 + gamma7 ) * tot_D0_to_HNL + \
             (gamma8*2 + gamma11*2 )                  * tot_D_to_HNL + \
             (gamma18*2)                              * tot_Ds_to_HNL

  B_w = B_w_inc # + B_w_lep 

  # B0  https://pdglive.lbl.gov/Particle.action?init=0&node=S042&home=MXXX045#decayclump_I
  gamma23,gamma24,gamma25,gamma26,gamma27,gamma28 = 0.081,0.474,0.039,0.369,0.103,0.026
  B0_w_inc = (gamma23 + gamma24) * tot_D0_to_uHNL + \
             (gamma25 + gamma26) * tot_D_to_uHNL + \
             (gamma27 + gamma28) * tot_Ds_to_uHNL

  #B0_w_lep = (gamma4*2 + gamma5 + gamma6*2 + gamma7) * tot_D_to_HNL + \
  #           (gamma8*2 + gamma11)                    * tot_D0_to_HNL 
             # no d_S contribution

  B0_w = B0_w_inc

  # Bs
  gamma1 = 0.93   # https://pdgprod.lbl.gov/pdgprod/pdgLive/Particle.action?init=0&node=S086&home=MXXX046
  Bs_w = gamma1 * tot_Ds_to_uHNL

  # Bc
  Bc_w = 0

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

  print 'Mass=1GeV'
  print 'B => mu HNL X:        Expected = {}'.format(getExpNevtsB(mass=1))
  print 'B => X D => mu HNL:   Expected = {}'.format(getExpNevtsD(mass=1))
  print ''
  print 'Mass=1.5GeV'
  print 'B => mu HNL X:        Expected = {}'.format(getExpNevtsB(mass=1.5))
  print 'B => X D => mu HNL:   Expected = {}'.format(getExpNevtsD(mass=1.5))
  print ''
  print 'Mass=1.8GeV'
  print 'B => mu HNL X:        Expected = {}'.format(getExpNevtsB(mass=1.8))
  print 'B => X D => mu HNL:   Expected = {}'.format(getExpNevtsD(mass=1.8))



