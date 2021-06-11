from python.decays import HNLDecays

import numpy as np
import math
from scipy.stats import expon  

# constants 
const_pi = math.pi
const_hbar = 6.582119569e-22 * 1e-03 # GeV s  # from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
const_c = 299792458. # m / s                  # from http://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf


def ctau_from_gamma(gamma):
    tau_natural = 1. / gamma                  # 1/GeV
    tau = tau_natural * const_hbar            # s
    ctau = tau * const_c * 1000               # mm
    return ctau

def gamma_total(mass,vv):
    '''
    Total width for N (Dirac)
    '''
    gamma_total =  HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['tot']   # GeV
    return gamma_total

def gamma_partial(mass,vv):
    '''
    Partial width for N->mupi (Dirac)
    '''
    gamma_partial = HNLDecays(mass=mass,mixing_angle_square=vv).decay_rate['mupi'] # GeV
    return gamma_partial

def BR_HNLmupion(mass): # vv is irrelevant, as it cancels out in the ratio
    return gamma_partial(mass=mass,vv=1.)/gamma_total(mass=mass,vv=1.)

def getCtau(mass=-99,vv=-99,ismaj=True):
    '''
    Helper function to go from vv,m -> ctau
    '''
    mult = 2. if ismaj else 1.
    return ctau_from_gamma(mult*gamma_total(mass=mass,vv=vv))

def getVV(mass=-99.,ctau=-99.,ismaj=True):
    '''
    Helper function to go from ctau,m -> vv
    '''
    mult = 2. if ismaj else 1.
    ref_m = 1. # GeV
    ref_vv = 1. 
    ref_ctau = ctau_from_gamma(mult*gamma_total(mass=ref_m,vv=ref_vv))

    k = ref_ctau * np.power(ref_m,5) * ref_vv

    return k/(np.power(mass, 5) * ctau)
   

def getCtauEff(ctau,cut=1000.):
  '''
  Returns an approximate estimate of the displacement filter efficiency (assuming beta*gamma=1) !!!
  '''
  rv = expon(scale=ctau)
  return rv.cdf(cut) # evaluate the cdf at ct=1000mm=1m

class Point(object):
  '''
  Class that contains information on mass,ctau,vv of a given signal point
  '''
  def __init__(self,mass,ctau=None,vv=None,is_ctau_rw=False,orig_vv=None,is_reco_rw=False,ismaj=True,myfilterEff=None):
    self.mass = mass
    self.is_ctau_rw = is_ctau_rw
    self.ismaj = ismaj
    self.myfilterEff = myfilterEff
    self.cfg = None

    if not vv: 
      self.ctau=ctau 
      self.vv=getVV(mass=self.mass, ctau=self.ctau, ismaj=self.ismaj)
    elif not ctau:
      self.vv = vv
      self.ctau=getCtau(mass=self.mass, vv=self.vv, ismaj=self.ismaj)
    else:
      self.vv = vv
      self.ctau = ctau

    if self.is_ctau_rw:
      self.orig_vv = orig_vv
      self.orig_ctau = getCtau(mass=self.mass,vv=orig_vv, ismaj=self.ismaj)
    else:
      self.orig_vv = self.vv 
      self.orig_ctau = self.ctau
   
    self.is_reco_rw = is_reco_rw
    self.name = '{:.1f}_{:.1f}'.format(self.mass,self.vv)
  #def getExpMedian():
    rv = expon(scale=self.ctau) 
    self.median = rv.median()
    #return rv.mean(),rv.median()

  def stamp(self):
    attrs=[]
    for k,v in self.__dict__.items():
      attrs.append(' {}={} '.format(k,v))
    attrs=' '.join(attrs)
    print(attrs)

  def stamp_simpli(self):
    attrs=[]
    attrs.append('m={}'.format(self.mass))
    attrs.append('vv={:.1e}'.format(self.vv))
    attrs.append('ctau={:.1e}'.format(self.ctau))
    #attrs.append('{}'.format(self.median))
    attrs=' '.join(attrs)
    print(attrs)
    return attrs

  def setConfig(self,cfg):
    if not self.cfg: 
      self.cfg = cfg

class Config(object):
  def __init__(self,nevtseff=10,muoneff=1.0,displeff=1.0,timeevt=80,timejob=23,contingency=1.5):
    self.nevtseff = int(nevtseff)     # number of wanted effective events 
    self.muoneff = float(muoneff)     #
    self.displeff = float(displeff)   #
    self.filtereff = float(muoneff*displeff) # efficiency of the generator filter
    self.timeevt = float(timeevt)     # time needed for each event (all steps),   in seconds    # take the one of the worst case scenario
    self.timejob = int(timejob)       # time for each job, including contingency, in hours
    self.contingency = float(contingency)   # factor 

    self.nevts       = int(float(self.nevtseff) / self.filtereff)                  # not needed for prodHelper.py  
    self.nevtseffjob = int( (self.timejob*3600./self.contingency) / self.timeevt )
    self.njobs       = int(float(self.nevtseff) / float(self.nevtseffjob)) + 1    
    self.nevtsjob    = int(float(self.nevts) / float(self.njobs) )                      # this will regulate maxEvents in cmsRun 


  def stamp(self):
    attrs=[]
    for k,v in self.__dict__.items():
      attrs.append(' {}={} '.format(k,v))
    attrs=' '.join(attrs)
    print(attrs)


if __name__ == "__main__":
  import matplotlib.pyplot as plt

  # test old vs new relation 
  mass = 1.0

  vvs = [5e-03, 1e-03, 5e-04, 1e-04, 5e-05, 1e-05, 5e-06, 1e-06, 5e-07]

