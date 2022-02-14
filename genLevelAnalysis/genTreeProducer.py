'''
Script to produce gen-level ntuples out of GEN-SIM files
'''

import sys
import os
import ROOT
import numpy as np
import glob
import re
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from scipy.constants import c as speed_of_light

from python.common import Point,getVV,getCtau 

##############
# Globals
#############
branches = [
    'run',  
    'lumi', 
    'event',
     
    # the mother
    'b_pt',
    'b_eta',
    'b_y',
    'b_phi',
    'b_mass',
    'b_q',
    'b_pdgid',
    'b_ct_reco',
    
    # daughters of the B
    # # the HNL
    'hnl_pt',
    'hnl_eta',
    'hnl_phi',
    'hnl_mass',
    'hnl_q',
    'hnl_qprop',
    'hnl_pdgid',
    #'hnl_ct_lhe', # from LHE information
    'hnl_ct_reco', # from Lxyz = ct\gamma\beta --> ct = Lxyz/(\gamma\beta)
    'hnl_beta',  # Lorentz
    'hnl_gamma', # Lorentz

    # # the meson
    'd_pt',
    'd_eta',
    'd_phi',
    'd_mass',
    'd_q',
    'd_pdgid',

    # # the trigger lepton
    'l0_pt',
    'l0_eta',
    'l0_phi',
    'l0_mass',
    'l0_q',
    'l0_pdgid',

    # daughters of the D meson
#    # # the pion
#    'pi_pt',
#    'pi_eta',
#    'pi_phi',
#    'pi_mass',
#    'pi_q',
#    'pi_pdgid',
#    
#    # # the kaon
#    'k_pt',
#    'k_eta',
#    'k_phi',
#    'k_mass',
#    'k_q',
#    'k_pdgid',
    
    # daughters of the HNL
    # # the lepton
    'l1_pt',
    'l1_eta',
    'l1_phi',
    'l1_mass',
    'l1_q',
    'l1_pdgid',

    # # the pion
    'pi1_pt',
    'pi1_eta',
    'pi1_phi',
    'pi1_mass',
    'pi1_q',
    'pi1_pdgid',
   
    # invariant masses
    'lep_pi_invmass',
    'lep_mu_invmass',
    'k_pi_invmass',
    #'b_invmass',
    'bpartial_invmass',

    'Lz',  # 1D longitudinal displacement for the HNL
    'Lxy', # 2D transverse displacement for the HNL
    'Lxyz',  # 3D displacement for the HNL
    'Lxy_cos', # cosine of the pointing angle in the transverse plane

    'Lxyz_b', #3D displacement of the B wrt to primary vertex
    'Lxyz_l0', #3D displacement of the trigger lepton wrt to B vertex

    'weight_reco', # reconstruction efficiency weight, as obtained by A-M's study
]

# couplings to be tested, for which the reweight is run
new_vvs = [
    1e-10, 
    5e-10, 
    1e-09, 
    5e-09, 
    1e-08, 
    5e-08, 
    1e-07, 
    2e-07, 
    5e-07, 
    8e-07, 
    1e-06, 
    2e-06, 
    5e-06, 
    8e-06, 
    1e-05, 
    2e-05, 
    3e-05, 
    4e-05, 
    5e-05, 
    6e-05, 
    7e-05, 
    8e-05, 
    9e-05, 
    0.0001, 
    0.0002, 
    0.0003, 
    0.0004, 
    0.0005, 
    0.0006, 
    0.0007, 
    0.0008, 
    0.0009,
    0.0010,
]

# add weight branches
for vv in new_vvs:
    branches.append('weight_%s'      %(str(vv).replace('-', 'm')))
    branches.append('ctau_%s'        %(str(vv).replace('-', 'm')))
    branches.append('xs_scale_to_%s' %(str(vv).replace('-', 'm')))

import pandas as pd
bins_pt  = [(0,7), (7,10), (10,15),(15,30), (30,1000)] #rows
bins_lxy = [(0,10),(10,30),(30,50),(50,100),(100,150),(150,300),(300,500),(500,1e10)] #columns, everything in mm
reco_weights = pd.read_csv('reco_weights_updated.csv', sep=',', comment='#') # version of 07/07 

###################
# Functions and classes
##################
def stampDaugthers(particle, indent=0):
    for idau in range(particle.numberOfDaughters()):
        daughter = particle.daughter(idau)
        print ' '*4*indent, '|---->', 'pdgid', daughter.pdgId(), '\t pt {0:>5}  eta {1:>5}  phi {2:>5}'.format('%.1f'%daughter.pt(), '%.2f'%daughter.eta(), '%.2f'%daughter.phi())
        stampDaugthers(daughter, indent+1)

def stampParticle(particle, indent=0):
   print ' '*4*indent, '', 'pdgid', particle.pdgId(), '\t pt {0:>5}  eta {1:>5}  phi {2:>5}'.format('%.1f'%particle.pt(), '%.2f'%particle.eta(), '%.2f'%particle.phi())

def isAncestor(a, p):
    if a == p :
        return True
    for i in xrange(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
            return True
    return False

def weight_to_new_ctau(old_ctau, old_v2, new_v2, ct):
    '''
    Returns an event weight based on the ratio of the normalised lifetime distributions.
    old_ctau: reference ctau
    old_v2  : reference coupling squared
    new_v2  : target coupling squared
    ct      : heavy neutrino lifetime in the specific event
    '''
    new_ctau = old_ctau * old_v2 / new_v2
    weight = old_ctau/new_ctau * np.exp( (1./old_ctau - 1./new_ctau) * ct )
    return weight, new_ctau

def scale_to_new_xs(old_v2, new_v2):
    '''
    At equal mass, the xs is simply proportional to v2 
    '''
    return new_v2 / old_v2


class Vertex(object):
   def __init__(self, index=0, x=0, y=0, z=0, pT=0):    
      self.index = index
      self.x     = x
      self.y     = y
      self.z     = z
      self.pT    = pT

def getpT(input):
   return input.pT


def runGenTreeProducer(infiles='./step*root',outfilename='out.root',this_mass=1,this_ctau=500,this_vv=0.0013,doFromMini=False):
  # input and output
  files = glob.glob(infiles)

  if len(files)==0: raise RuntimeError('No files to be run!, glob expression = {}'.format(infiles))
  outfile = ROOT.TFile.Open(outfilename, 'recreate')

  handles = OrderedDict()
  if doFromMini:
    handles['genP'] = ('packedGenParticles' , Handle('std::vector<pat::PackedGenParticle>')) # does not seem to be working quite well
  else:
    handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>'))
  #handles['lhe']         = ('externalLHEProducer', Handle('LHEEventProduct'))
  
  # output file and tree gymnastics
  ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
  tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       
    
  # get to the real thing
  print 'loading the file ...'
  print files
  events = Events(files)
  print '... done!'

  for i, event in enumerate(events):
    if opt.maxEvts != -1 and i+1 > opt.maxEvts: break    

    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())
  
    if i%1000==0:
    #if i%1==0:
        percentage = float(i)/events.size()*100.
        print '\t===> processing %d / %d event \t completed %.1f%s' %(i, events.size(), percentage, '%')
  
    #print '\n Event {a}'.format(a=i)
        
    # get the heavy neutrino
    the_hns = [ip for ip in event.genP if abs(ip.pdgId())==9900015 and ip.isLastCopy()] 
    if len(the_hns)!=0:
       event.the_hn = the_hns[0] # one per event
    else:
       event.the_hn = None
   
    # find the B; if the event is a BHNL, take the mother of the HNL, else, take the genParticles which are Bs
    if event.the_hn:
      event.the_hn.mothers = [event.the_hn.mother(jj) for jj in range(event.the_hn.numberOfMothers())]
      the_b_mothers = sorted([ii for ii in event.the_hn.mothers if (abs(ii.pdgId())==521 or abs(ii.pdgId())==511 or abs(ii.pdgId())==531 \
                                                                 or abs(ii.pdgId())==541)], key = lambda x : x.pt(), reverse=True)
      if len(the_b_mothers)!=0:
        event.the_b_mother = the_b_mothers[0]
    else:
      #the_b_mothers = sorted([ii for ii in event.genP if abs(ii.pdgId())==521 and ()  and ii.isLastCopy()], key = lambda x : x.pt(), reverse=True)
      for ip in event.genP:
        if abs(ip.pdgId())==521  and ip.isLastCopy():
          hasJPsiDaughter = False
          hasKDaughter = False
          for idau in range(ip.numberOfDaughters()):
            if abs(ip.daughter(idau).pdgId())==443:
              hasJPsiDaughter = True
            if abs(ip.daughter(idau).pdgId())==321:
              hasKDaughter = True
          if hasJPsiDaughter and hasKDaughter:
            event.the_b_mother = ip
      
    # get the other daughters of the B meson
    event.the_b_mother.daughters = [event.the_b_mother.daughter(jj) for jj in range(event.the_b_mother.numberOfDaughters())]

    if opt.doDebug: 
      stampParticle(event.the_b_mother)
      stampDaugthers(event.the_b_mother)
  
    # # first the daughter meson if it exists (we will call it D, but it could be something else)
    the_ds = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId()) in [111,211,113,213,321,323,411,421,413,423,431,433]], 
                     key = lambda x : x.pt(), reverse=True)
    if len(the_ds)!=0:
      event.the_d = the_ds[0]
    else:
      event.the_d = None

    # # then the trigger lepton
    the_pls = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt(), reverse=True)
    if len(the_pls)!=0:
      event.the_pl = the_pls[0]
    else:
      event.the_pl = None
     
    if event.the_hn: 
      # HNL daughters
      event.the_hn.initialdaus = [event.the_hn.daughter(jj) for jj in range(event.the_hn.numberOfDaughters())]
    
      # # the lepton
      the_lep_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt(), reverse=True)
      if len(the_lep_daughters)!=0:
         event.the_hn.lep = the_lep_daughters[0]
      else: 
         event.the_hn.lep = None
     
      # # the pion
      the_pi_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId())==211], key = lambda x : x.pt(), reverse=True)
      if len(the_pi_daughters)!=0:
         event.the_hn.pi = the_pi_daughters[0]
      else:
         event.the_hn.pi = None
     
      # invariant masses   
      # # to get the invariant mass of the HNL daughters
      if len(the_pi_daughters)!=0 and len(the_lep_daughters)!=0:   
        event.the_hnldaughters = event.the_hn.lep.p4() + event.the_hn.pi.p4()
    
      # # to get the full or partial invariant mass of the B 
      # # # partial
      if len(the_pls)!=0 and len(the_hns)!=0:
        event.the_bdaughters_partial = event.the_hn.p4() + event.the_pl.p4()
        event.the_bdaughters_mu_lep  = event.the_pl.p4() + event.the_hn.lep.p4()
      # # # total
      #event.the_bdaughters_all = sum([ii.p4() for ii in event.the_b_mother.daughters])
      
      # identify the primary vertex
      # for that, needs the trigger lepton
      if len(the_pls)!=0:
        event.the_hn.the_pv = event.the_pl.vertex()
      
      # identify the secondary vertex
      if len(the_lep_daughters)!=0:
        event.the_hn.the_sv = event.the_hn.lep.vertex()
      
      # 2D transverse and 3D displacement, Pythagoras
      if len(the_pls)!=0 and len(the_lep_daughters)!=0:
        event.Lz  =  np.sqrt((event.the_hn.the_pv.z() - event.the_hn.the_sv.z())**2) * 10 # everything in mm
        
        event.Lxy  = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                             (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2) * 10 # everything in mm
  
        event.Lxyz = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                             (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2 + \
                             (event.the_hn.the_pv.z() - event.the_hn.the_sv.z())**2) * 10 # everything in mm
   
      # per event ct, as derived from the flight distance and Lorentz boost  
      event.the_hn.beta  = event.the_hn.p4().Beta()
      event.the_hn.gamma = event.the_hn.p4().Gamma()
      
      # we get the lifetime from the kinematics 
      if len(the_pls)!=0 and len(the_lep_daughters)!=0:
        event.the_hn.ct_reco = event.Lxyz / (event.the_hn.beta * event.the_hn.gamma)
    
      # pointing angle    
      #if len(the_pls) and len(the_lep_daughters):
      #  hn_pt_vect = ROOT.math.XYZVector(event.the_hnldaughters.px(),
      #                                   event.the_hnldaughters.py(),
      #                                   0.)
    
      #  Lxy_vect = ROOT.GlobalPoint(-1*((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())), 
      #                              -1*((event.the_hn.the_pv.y() - event.the_hn.the_sv.y())),
      #                               0)
    
      #  vperptau = ROOT.math.XYZVector(Lxy_vect.x(), Lxy_vect.y(), 0.)
     
        #event.cos_pointing = vperptau.Dot(hn_pt_vect)/(vperptau.R()*hn_pt_vect.R())
    
    if len(the_pls)!=0:
      event.Lxyz_pl = np.sqrt((event.the_b_mother.vx() - event.the_pl.vx())**2 + \
                              (event.the_b_mother.vy() - event.the_pl.vy())**2 + \
                              (event.the_b_mother.vz() - event.the_pl.vz())**2) * 10 # everything in mm
  
    # set the reconstruction efficiency weigh
    event.weight_reco = 1.
    if opt.doWeightReco:

    #  while event.weight_reco == 1.:
        for i,bin_pt in enumerate(bins_pt):
            for j,bin_lxy in enumerate(bins_lxy):
              #print i,bin_pt,j,bin_lxy,reco_weights.iloc[i][j]
              if event.Lxy > bin_lxy[0] and event.Lxy < bin_lxy[1] and event.the_hn.pt() > bin_pt[0] and event.the_hn.pt() < bin_pt[1]:
                event.weight_reco = reco_weights.iloc[i][j] 
                #print i,bin_pt,j,bin_lxy,event.weight_reco

    # get the lifetime of the B
    event.the_b_mother.beta  = event.the_b_mother.p4().Beta()
    event.the_b_mother.gamma = event.the_b_mother.p4().Gamma()   
    b_p4 = event.the_b_mother.p4()
    event.the_b_mother.y     = 0.5 * ROOT.TMath.Log( (b_p4.E()+b_p4.Pz()) / (b_p4.E()-b_p4.Pz()) ) if (b_p4.E()-b_p4.Pz())!=0 else -99
    #event.the_b_mother.ct_reco = event.Lxyz_pl / (event.the_b_mother.beta * event.the_b_mother.gamma)
  
    tofill['run'        ] = event.eventAuxiliary().run()
    tofill['lumi'       ] = event.eventAuxiliary().luminosityBlock()
    tofill['event'      ] = event.eventAuxiliary().event()
     
    if event.the_b_mother:
        tofill['b_pt'   ] = event.the_b_mother.pt()     
        tofill['b_eta'  ] = event.the_b_mother.eta()    
        tofill['b_y'    ] = event.the_b_mother.y 
        tofill['b_phi'  ] = event.the_b_mother.phi()    
        tofill['b_mass' ] = event.the_b_mother.mass()   
        tofill['b_q'    ] = event.the_b_mother.charge()   
        tofill['b_pdgid'] = event.the_b_mother.pdgId()   
        #tofill['b_ct_reco'] = event.the_b_mother.ct_reco   
    
    if event.the_hn: 
      tofill['hnl_pt'     ] = event.the_hn.pt()     
      #tofill['hnl_p'      ] = event.the_hn.pt() * ROOT.TMath.CosH(event.the_hn.eta())
      tofill['hnl_eta'    ] = event.the_hn.eta()    
      tofill['hnl_phi'    ] = event.the_hn.phi()    
      tofill['hnl_mass'   ] = event.the_hn.mass()   
      #tofill['hnl_ct_lhe' ] = event.hnl_ct_lhe 
      if len(the_lep_daughters)!=0:
        tofill['hnl_ct_reco'] = event.the_hn.ct_reco 
      tofill['hnl_beta'   ] = event.the_hn.beta  
      tofill['hnl_gamma'  ] = event.the_hn.gamma  
      tofill['hnl_pdgid'  ] = event.the_hn.pdgId()  
  
    if event.the_d:
      tofill['d_pt'   ] = event.the_d.pt()     
      tofill['d_eta'  ] = event.the_d.eta()    
      tofill['d_phi'  ] = event.the_d.phi()    
      tofill['d_mass' ] = event.the_d.mass()   
      tofill['d_q'    ] = event.the_d.charge()   
      tofill['d_pdgid'] = event.the_d.pdgId()   
   
    if event.the_pl:
       tofill['l0_pt'      ] = event.the_pl.pt()
       tofill['l0_eta'     ] = event.the_pl.eta()
       tofill['l0_phi'     ] = event.the_pl.phi()
       tofill['l0_mass'    ] = event.the_pl.mass()
       tofill['l0_q'       ] = event.the_pl.charge()
       tofill['l0_pdgid'   ] = event.the_pl.pdgId()
 
    if event.the_hn: 
      if event.the_hn.lep:
         tofill['l1_pt'      ] = event.the_hn.lep.pt()
         tofill['l1_eta'     ] = event.the_hn.lep.eta()
         tofill['l1_phi'     ] = event.the_hn.lep.phi()
         tofill['l1_mass'    ] = event.the_hn.lep.mass()
         tofill['l1_q'       ] = event.the_hn.lep.charge()
         tofill['l1_pdgid'   ] = event.the_hn.lep.pdgId()
    
      if event.the_hn.pi:   
         tofill['pi1_pt'      ] = event.the_hn.pi.pt()
         tofill['pi1_eta'     ] = event.the_hn.pi.eta()
         tofill['pi1_phi'     ] = event.the_hn.pi.phi()
         tofill['pi1_mass'    ] = event.the_hn.pi.mass()
         tofill['pi1_q'       ] = event.the_hn.pi.charge()
         tofill['pi1_pdgid'   ] = event.the_hn.pi.pdgId()
  
      tofill['lep_pi_invmass' ] = event.the_hnldaughters.mass()
      tofill['bpartial_invmass'] = event.the_bdaughters_partial.mass()
      tofill['lep_mu_invmass'] = event.the_bdaughters_mu_lep.mass()
    #tofill['b_invmass'] = event.the_bdaughters_all.mass()
    
      # hnl charge
      if event.the_hn.lep and event.the_hn.pi:
        tofill['hnl_q'      ] = event.the_hn.lep.charge() + event.the_hn.pi.charge() 
      tofill['hnl_qprop']   = event.the_hn.charge()
    
    if len(the_pls)!=0 and len(the_lep_daughters)!=0:
      tofill['Lz'        ] = event.Lz
      tofill['Lxy'        ] = event.Lxy
      tofill['Lxyz'       ] = event.Lxyz
      #tofill['Lxy_cos'    ] = event.cos_pointing
  
    if event.the_pl and event.the_b_mother:
      tofill['Lxyz_l0'] = event.Lxyz_pl

    # weights for ctau reweighting 
    if event.the_hn:
      for vv in new_vvs:
          tofill['weight_%s'      %(str(vv).replace('-', 'm'))] = weight_to_new_ctau(old_ctau=this_ctau, old_v2=this_vv, new_v2=vv, ct=event.the_hn.ct_reco)[0]
          tofill['ctau_%s'        %(str(vv).replace('-', 'm'))] = weight_to_new_ctau(old_ctau=this_ctau, old_v2=this_vv, new_v2=vv, ct=event.the_hn.ct_reco)[1]
          tofill['xs_scale_to_%s' %(str(vv).replace('-', 'm'))] = vv / this_vv
  
    tofill['weight_reco'] = event.weight_reco
         
    ntuple.Fill(array('f',tofill.values()))
  
  outfile.cd()
  ntuple.Write()
  outfile.Close()



def getOptions():
   from argparse import ArgumentParser
   parser = ArgumentParser(description='options for Gen tree producer', add_help=True)
   parser.add_argument('--pl', type=str, dest='pl', help='production label', default='V02_muFromB_pt5_eta1p6_njt30')
   #parser.add_argument('--expr', type=str, dest='expr', help='file regular expression', default='step1*root')
   parser.add_argument('--points', type=str, dest='pointFile', help='name of file contaning information on scan to be run', default='points.py')
   parser.add_argument('--maxEvts', type=int, dest='maxEvts', help='max number of events to be run', default=-1)
   parser.add_argument('--doDebug', dest='doDebug', help='', action='store_true', default=False)
   parser.add_argument('--doFromMini', dest='doFromMini', help='run on the miniAOD, as opposed to the GEN-SIM', action='store_true', default=False)
   parser.add_argument('--doSaveScratch', dest='doSaveScratch', help='save on scratch, not on work', action='store_true', default=False)
   parser.add_argument('--doWeightReco', dest='doWeightReco', help='calculate weight for reconstruction efficiency', action='store_true', default=False)
   return parser.parse_args()
  

if __name__ == "__main__":

  opt = getOptions()
  # point objects are imported
  sys.path.append('../slurm/.')
  ps = __import__(opt.pointFile.split('.py')[0]) 

  for p in ps.points:

    user = os.environ["USER"] 
    fileExpr = 'step1*nj*root' if not opt.doFromMini else 'step4*root'
    expr = '/pnfs/psi.ch/cms/trivcat/store/user/{usr}/BHNLsGen/{pl}/mass{m}_ctau{ctau}/{ex}'.format(usr=user,pl=opt.pl,m=p.mass,ctau=p.ctau,ex=fileExpr)
    #expr = '/work/mratti/GEN_HNL/CMSSW_10_2_15/src/HNLsGen/slurm/pilot_V15_control_BfilterNoMufilter/BPH-step1_numEvent357.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/fragments_test/CMSSW_10_2_15/src/ok_BToNMuX_NToEMuPi_test.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/fragments_test/CMSSW_10_2_15/src/BToNMuX_NToEMuPi_test.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/fragments_test/CMSSW_10_2_27/src/BToNMuX_NToEMuPi_test.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/fragments_test/CMSSW_10_2_15/src/BcToNMuX_NToEMuPi_test.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/fragments_test/CMSSW_10_2_27/src/BcToNMuX_NToEMuPi_test.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/fragments_test/CMSSW_10_2_15/src/BHNL_test*.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/CMSSW_10_2_3/src/HNLsGen/slurm/testManyChan_n100_njt1/BPH-step1_numEvent1000.root'
    #expr = '/work/mratti/GEN_HNL_newPythia/CMSSW_10_2_3/src/HNLsGen/slurm/testBc_n10_njt1/BPH-step1_numEvent1000.root'
    opath = './' if not opt.doSaveScratch else '/scratch/mratti/BHNLsGenDump/'
    outfilename = '{opath}/outputfiles/{pl}/mass{m}_ctau{ctau}_miniGenTree{suffix}.root'.format(opath=opath, pl=opt.pl,m=p.mass,ctau=p.ctau,suffix='_fromMini' if opt.doFromMini else '')
    os.system('mkdir -p {opath}/outputfiles/{pl}'.format(opath=opath ,pl=opt.pl))    

    print('************************************************************************') 
    print('\nGoing to run gen tree producer over')
    print('   pl  : {}'.format(opt.pl))
    print('   mass: {} GeV'.format(p.mass))
    print('   ctau: {} mm'.format(p.ctau))
    print('   VV  : {} \n'.format(p.vv))

    runGenTreeProducer(infiles=expr,outfilename=outfilename,this_mass=p.mass,this_ctau=p.ctau,this_vv=p.vv,doFromMini=opt.doFromMini)

