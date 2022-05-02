# http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf

import os
import sys
import glob
import ROOT
import numpy as np
import pandas as pd
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from scipy.constants import c as speed_of_light
from os import path


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

branches = [
    'run',  
    'lumi', 
    'event',
     
    # the mother
    'b_pt',
    'b_eta',
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
    'hnl_pdgid',
    #'hnl_ct_lhe', # from LHE information
    'hnl_ct_reco', # from Lxyz = ct\gamma\beta --> ct = Lxyz/(\gamma\beta)
    'hnl_beta',  # Lorentz
    'hnl_gamma', # Lorentz

    # # the accompagnying meson
    'X_pt',
    'X_eta',
    'X_phi',
    'X_mass',
    'X_q',
    'X_pdgid',

    # # the B lepton
    'mu_fromB_pt',
    'mu_fromB_eta',
    'mu_fromB_phi',
    'mu_fromB_mass',
    'mu_fromB_q',
    'mu_fromB_pdgid',
    'mu_fromB_satisfies_BParkHLT_cond',
    
    # daughters of the HNL
    # # the lepton
    'mu_fromHNL_pt',
    'mu_fromHNL_eta',
    'mu_fromHNL_phi',
    'mu_fromHNL_mass',
    'mu_fromHNL_q',
    'mu_fromHNL_pdgid',
    'mu_fromHNL_satisfies_BParkHLT_cond',

    # # the pion
    'pi_fromHNL_pt',
    'pi_fromHNL_eta',
    'pi_fromHNL_phi',
    'pi_fromHNL_mass',
    'pi_fromHNL_q',
    'pi_fromHNL_pdgid',

    # event topology
    'hnl_tag_side',
    'number_triggering_muon_hnlside',
    'atleast_1triggeringmuon_otherside',
    'mu_fromHNL_only_satisfies_BParkHLT_cond_tagside',
    'onetriggering_mu_fromB_tagside',
    'onetriggering_mu_fromHNL_tagside',
   
    # invariant masses
    'mu_pi_invmass',
    'mu_HNL_X_invmass',

    # displacement
    'Lxy',   # 2D transverse displacement for the HNL
    'Lxyz',  # 3D displacement for the HNL
    #'Lxy_cos', # cosine of the pointing angle in the transverse plane

    #'Lxyz_b', #3D displacement of the B wrt to primary vertex
    'Lxyz_l0', #3D displacement of the triggered lepton wrt to B vertex

     # reconstruction efficiency weight
    'weight_reco',

    # ct weight (used for lifetime reweighting study)
    #'weight_ct',
]

# couplings to be tested, for which the reweight is run
new_v2s = [
    1e-10, 
    5e-10, 
    1e-9, 
    5e-9, 
    1e-8, 
    5e-8, 
    1e-7, 
    5e-7, 
    1e-6, 
    5e-6, 
    6e-06, 
    8e-06, 
    1e-5, 
    2e-5, 
    3e-5, 
    4e-5, 
    5e-5, 
    7e-05, 
    0.0001, 
    0.0002, 
    0.00025, 
    0.0003, 
    0.0005, 
    0.0012,
]

# ctau weights, each coupling gets its own weight
weights = OrderedDict(zip(new_v2s, np.ones(len(new_v2s))))

# add weight branches
#for vv in new_v2s:
#    branches.append('weight_%s'      %(str(vv).replace('-', 'm')))
#    branches.append('ctau_%s'        %(str(vv).replace('-', 'm')))
#    branches.append('xs_scale_to_%s' %(str(vv).replace('-', 'm')))


#handles = OrderedDict()
#handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>'))
#
## output file and tree gymnastics
##outfile = ROOT.TFile.Open('outputfiles/test_modfilter_v4_m3_ctau1000.root', 'recreate')
#outfile = ROOT.TFile.Open('outputfiles/test.root', 'recreate')
#ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
#tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       
#
## get to the real thing
#print 'loading the file ...'
#infiles = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/test_modfilter_v3_n10000000_njt500/mass4.5_ctau1.0/step1*.root'
##infiles = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/mass4.5_ctau1.0/step1*root'


def treeProducer(infiles, outdir, outfilename):
  # get the files
  print 'loading the file ...'
  files = glob.glob(infiles)
  events = Events(files)
  print 'nevts: ',events.size()
  print '... done!'

  # get the handles
  handles = OrderedDict()
  handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>'))

  # output file and tree gymnastics
  outfile = ROOT.TFile.Open('{}/{}'.format(outdir, outfilename), 'recreate')
  ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
  tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

  # get reconstruction weights
  bins_pt  = [(0,7), (7,10), (10,15),(15,30), (30,1000)] #rows
  bins_lxy = [(0,10),(10,30),(30,50),(50,100),(100,150),(150,300),(300,500),(500,1e10)] #columns, everything in mm
  reco_weights = pd.read_csv('reco_weights_updated.csv', sep=',', comment='#')

  # file ct weight
  #file_ct = ROOT.TFile.Open('weight_ct_m3_ctau1.root', 'READ')
  #hist_weight_ct = file_ct.Get('hist_weight')

  for i, event in enumerate(events):
    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())

    if i%1000==0:
        percentage = float(i)/events.size()*100.
        print '\t===> processing %d / %d event \t completed %.1f%s' %(i, events.size(), percentage, '%')

    #print '\n Event {a}'.format(a=i)
    #if float(i)>5000.: continue
        
    # get the heavy neutrino
    the_hns = [ip for ip in event.genP if abs(ip.pdgId())==9900015 and ip.isLastCopy()] 
    if len(the_hns):
       event.the_hn = the_hns[0] # one per event
       #print '0. found hnl of pdgId {a}'.format(a=event.the_hn.pdgId())
    else:
      print 'no hnl found'

    # find the B mother of the hnl 
    event.the_hn.mothers = [event.the_hn.mother(jj) for jj in range(event.the_hn.numberOfMothers())]
    the_b_mothers = sorted([ii for ii in event.the_hn.mothers if (abs(ii.pdgId())==521 or abs(ii.pdgId())==511 or abs(ii.pdgId())==531 or abs(ii.pdgId())==541)], key = lambda x : x.pt(), reverse=True)
    if len(the_b_mothers):
      event.the_b_mother = the_b_mothers[0]
      #if abs(event.the_b_mother.pdgId())==531: print '1. found a B of pdgId {a}'.format(a=event.the_b_mother.pdgId())
    else:
      event.the_b_mother = None
      print 'no mother found'

    # get the other daughters of the B meson
    event.the_b_mother.daughters = [event.the_b_mother.daughter(jj) for jj in range(event.the_b_mother.numberOfDaughters())]
    if len(event.the_b_mother.daughters)<2: print len(event.the_b_mother.daughters)

    ## first, the B muon
    the_pls = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt(), reverse=True)
    if len(the_pls):
      event.the_pl = the_pls[0]
      #print '3. found the B muon of pdgId {a}'.format(a=event.the_pl.pdgId())
      if event.the_pl.pt() > 6.8 and abs(event.the_pl.eta()) < 1.55: event.the_pl.satisfies_BParkHLT_cond = 1
      else: event.the_pl.satisfies_BParkHLT_cond = 0
    else:
      event.the_pl = None
      print 'no B muon found'

    #TODO check rate of different B decay channels: make it inclusive and per B species, with and without photon
    # then check for the accompagnying meson
    the_Xs = None
    event.the_X = False
    if len(event.the_b_mother.daughters)>2:
      the_Xs = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId()) not in [11, 13, 15, 9900015]], key = lambda x : x.pt(), reverse=True)
      if len(the_Xs):
        event.the_X = the_Xs[0]
        #if len(the_Xs) != 1 and abs(the_Xs[0].pdgId()) == 22:
        #  print '\n'
        #  for ix, X in enumerate(the_Xs):
        #    print the_Xs[ix].pdgId()
        #if abs(event.the_b_mother.pdgId()==531): print 'found accompagnying meson of pdgId {a}'.format(a=event.the_X.pdgId())
      else:
        event.the_X = None
        print 'no accompagnying meson found'

    # HNL daughters
    event.the_hn.initialdaus = [event.the_hn.daughter(jj) for jj in range(event.the_hn.numberOfDaughters())]

    ## the hnl lepton
    #TODO make sure that we only consider muon channel
    the_lep_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId()) in [11, 13, 15]], key = lambda x : x.pt(), reverse=True)
    if len(the_lep_daughters):
      event.the_hn.lep = the_lep_daughters[0]
      if event.the_hn.lep.pt() > 6.8 and abs(event.the_hn.lep.eta()) < 1.55: event.the_hn.lep.satisfies_BParkHLT_cond = 1
      else: event.the_hn.lep.satisfies_BParkHLT_cond = 0
       #print '6. found hnl daughter lepton of pdgId {a}'.format(a=event.the_hn.lep.pdgId())
    else: 
      event.the_hn.lep = None
      print 'hnl muon not found'
   
    # # the pion
    the_pi_daughters = sorted([ii for ii in event.the_hn.initialdaus if abs(ii.pdgId())==211], key = lambda x : x.pt(), reverse=True)
    if len(the_pi_daughters):
       event.the_hn.pi = the_pi_daughters[0]
       #print '7. found a daughter pion of pdgId {a}'.format(a=event.the_hn.pi.pdgId())
    else:
       event.the_hn.pi = None
       print 'hnl pion not found'

    # event topology
    if event.the_pl.satisfies_BParkHLT_cond + event.the_hn.lep.satisfies_BParkHLT_cond > 0:
      event.hnl_tag_side = 1
    else:
      event.hnl_tag_side = 0

    event.number_triggering_muon_hnlside = event.the_pl.satisfies_BParkHLT_cond + event.the_hn.lep.satisfies_BParkHLT_cond

    event.mu_fromHNL_only_satisfies_BParkHLT_cond_tagside = -99
    if event.hnl_tag_side: # either one or two triggering muons
      if event.the_pl.satisfies_BParkHLT_cond: 
        event.mu_fromHNL_only_satisfies_BParkHLT_cond_tagside = 0
      else:
        event.mu_fromHNL_only_satisfies_BParkHLT_cond_tagside = 1

    event.onetriggering_mu_fromB_tagside = -99
    event.onetriggering_mu_fromHNL_tagside = -99
    if event.number_triggering_muon_hnlside == 1 and event.the_pl.satisfies_BParkHLT_cond:
      event.onetriggering_mu_fromB_tagside = 1
      event.onetriggering_mu_fromHNL_tagside = 0
    elif event.number_triggering_muon_hnlside == 1 and event.the_hn.lep.satisfies_BParkHLT_cond:
      event.onetriggering_mu_fromB_tagside = 0
      event.onetriggering_mu_fromHNL_tagside = 1

    # checking that there is at least one triggering muon per event
    event.atleast_1triggeringmuon_otherside = -1
    if event.number_triggering_muon_hnlside == 0:
      the_trigger_muon = [ip for ip in event.genP if abs(ip.pdgId())==13 and ip.pt()>6.8 and abs(ip.eta())<1.55 and ip.isLastCopy()] 
      if len(the_trigger_muon): event.atleast_1triggeringmuon_otherside = 1
      else: event.atleast_1triggeringmuon_otherside = 0
      
    # invariant masses   
    # # to get the invariant mass of the HNL daughters
    if len(the_pi_daughters) and len(the_lep_daughters):   
      event.the_hnldaughters = event.the_hn.lep.p4() + event.the_hn.pi.p4()
    
    # # to get the invariant mass of the B daughters
    if the_Xs != None:
      event.the_bdaughters = event.the_hn.p4() + event.the_pl.p4() + event.the_X.p4()
    else:
      event.the_bdaughters = event.the_hn.p4() + event.the_pl.p4()
      
    # identify the primary vertex
    # for that, needs the B muon
    if len(the_pls):
      event.the_hn.the_pv = event.the_pl.vertex()
    
    # identify the secondary vertex
    if len(the_lep_daughters):
      event.the_hn.the_sv = event.the_hn.lep.vertex()
    
    # 2D transverse and 3D displacement, Pythagoras
    if len(the_pls) and len(the_lep_daughters):
      #print 'pl vx: {a} l1 vx: {b} hnl vx: {c}'.format(a=event.the_pl.vx(), b=event.the_hn.lep.vx(), c=event.the_hn.vx()) 
      event.Lxy  = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                           (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2) * 10 # in mm

      event.Lxyz = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
                           (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2 + \
                           (event.the_hn.the_pv.z() - event.the_hn.the_sv.z())**2) * 10 # in mm

      #event.Lxyz = np.sqrt((event.the_hn.vx() - event.the_hn.lep.vx())**2 + \
      #                  (event.the_hn.vy() - event.the_hn.lep.vy())**2 + \
      #                  (event.the_hn.vz() - event.the_hn.lep.vz())**2)

    # per event ct, as derived from the flight distance and Lorentz boost  
    event.the_hn.beta  = event.the_hn.p4().Beta()
    event.the_hn.gamma = event.the_hn.p4().Gamma()
    
    # we get the lifetime from the kinematics 
    if len(the_pls) and len(the_lep_daughters):
      event.the_hn.ct_reco = event.Lxyz / (event.the_hn.beta * event.the_hn.gamma)
      #print 'hnl ct reco: {a}'.format(a=event.the_hn.ct_reco)

    # temporary condition
    #if event.the_hn.ct_reco > 150: continue 

    # pointing angle    
    #if len(the_pls) and len(the_lep_daughters):
    #  hn_pt_vect = ROOT.math.XYZVector(event.the_hnldaughters.px(),
    #                                   event.the_hnldaughters.py(),
    #                                   0.)

    #  Lxy_vect = ROOT.GlobalPoint(-1*((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())), 
    #                              -1*((event.the_hn.the_pv.y() - event.the_hn.the_sv.y())),
    #                               0)

    #  vperptau = ROOT.math.XYZVector(Lxy_vect.x(), Lxy_vect.y(), 0.)
   
    #  event.cos_pointing = vperptau.Dot(hn_pt_vect)#/(vperptau.R()*hn_pt_vect.R())


    event.Lxyz_pl = np.sqrt((event.the_b_mother.vx() - event.the_pl.vx())**2 + \
                        (event.the_b_mother.vy() - event.the_pl.vy())**2 + \
                        (event.the_b_mother.vz() - event.the_pl.vz())**2) * 10 # in mm
    
    # get the lifetime of the B
    event.the_b_mother.beta  = event.the_b_mother.p4().Beta()
    event.the_b_mother.gamma = event.the_b_mother.p4().Gamma()
   
    event.the_b_mother.ct_reco = event.Lxyz_pl / (event.the_b_mother.beta * event.the_b_mother.gamma)

    # get reconstruction efficiency weight
    event.weight_reco = 1.
    for i,bin_pt in enumerate(bins_pt):
      for j,bin_lxy in enumerate(bins_lxy):
        if event.Lxy > bin_lxy[0] and event.Lxy < bin_lxy[1] and event.the_hn.pt() > bin_pt[0] and event.the_hn.pt() < bin_pt[1]:
          event.weight_reco = reco_weights.iloc[i][j] 

    #event.weight_ct = hist_weight_ct.GetBinContent(hist_weight_ct.GetXaxis().FindBin(event.the_hn.ct_reco/10.)) 
   
    # reset before filling
    for k, v in tofill.iteritems(): tofill[k] = -99. # initialise before filling

    #if event.the_pl.pt() > 9 and abs(event.the_pl.eta()) < 1.5 and event.the_hn.lep.pt() > 3 and abs(event.the_hn.lep.eta()) < 2.5: # and event.Lxy < 1:

    tofill['run'        ] = event.eventAuxiliary().run()
    tofill['lumi'       ] = event.eventAuxiliary().luminosityBlock()
    tofill['event'      ] = event.eventAuxiliary().event()
     
    if event.the_b_mother:
        tofill['b_pt'   ] = event.the_b_mother.pt()     
        tofill['b_eta'  ] = event.the_b_mother.eta()    
        tofill['b_phi'  ] = event.the_b_mother.phi()    
        tofill['b_mass' ] = event.the_b_mother.mass()   
        tofill['b_q'    ] = event.the_b_mother.charge()   
        tofill['b_pdgid'] = event.the_b_mother.pdgId()   
        tofill['b_ct_reco'] = event.the_b_mother.ct_reco   
     
    tofill['hnl_pt'     ] = event.the_hn.pt()     
    tofill['hnl_eta'    ] = event.the_hn.eta()    
    tofill['hnl_phi'    ] = event.the_hn.phi()    
    tofill['hnl_mass'   ] = event.the_hn.mass()   
    tofill['hnl_q'      ] = event.the_hn.lep.charge() + event.the_hn.pi.charge()
    #tofill['hnl_ct_lhe' ] = event.hnl_ct_lhe 
    if len(the_lep_daughters):
      tofill['hnl_ct_reco'] = event.the_hn.ct_reco #* 10. # convert cm to mm 
    tofill['hnl_beta'   ] = event.the_hn.beta  
    tofill['hnl_gamma'  ] = event.the_hn.gamma  
    tofill['hnl_pdgid'  ] = event.the_hn.pdgId()  

    if event.the_X:
        tofill['X_pt'   ] = event.the_X.pt()     
        tofill['X_eta'  ] = event.the_X.eta()    
        tofill['X_phi'  ] = event.the_X.phi()    
        tofill['X_mass' ] = event.the_X.mass()   
        tofill['X_q'    ] = event.the_X.charge()   
        tofill['X_pdgid'] = event.the_X.pdgId()   

    if event.the_pl:
       tofill['mu_fromB_pt'      ] = event.the_pl.pt()
       tofill['mu_fromB_eta'     ] = event.the_pl.eta()
       tofill['mu_fromB_phi'     ] = event.the_pl.phi()
       tofill['mu_fromB_mass'    ] = event.the_pl.mass()
       tofill['mu_fromB_q'       ] = event.the_pl.charge()
       tofill['mu_fromB_pdgid'   ] = event.the_pl.pdgId()
       tofill['mu_fromB_satisfies_BParkHLT_cond']   = event.the_pl.satisfies_BParkHLT_cond

    if event.the_hn.lep:
       tofill['mu_fromHNL_pt'      ] = event.the_hn.lep.pt()
       tofill['mu_fromHNL_eta'     ] = event.the_hn.lep.eta()
       tofill['mu_fromHNL_phi'     ] = event.the_hn.lep.phi()
       tofill['mu_fromHNL_mass'    ] = event.the_hn.lep.mass()
       tofill['mu_fromHNL_q'       ] = event.the_hn.lep.charge()
       tofill['mu_fromHNL_pdgid'   ] = event.the_hn.lep.pdgId()
       tofill['mu_fromHNL_satisfies_BParkHLT_cond'] = event.the_hn.lep.satisfies_BParkHLT_cond

    if event.the_hn.pi:   
       tofill['pi_fromHNL_pt'      ] = event.the_hn.pi.pt()
       tofill['pi_fromHNL_eta'     ] = event.the_hn.pi.eta()
       tofill['pi_fromHNL_phi'     ] = event.the_hn.pi.phi()
       tofill['pi_fromHNL_mass'    ] = event.the_hn.pi.mass()
       tofill['pi_fromHNL_q'       ] = event.the_hn.pi.charge()
       tofill['pi_fromHNL_pdgid'   ] = event.the_hn.pi.pdgId()

    # event topology
    tofill['hnl_tag_side'] = event.hnl_tag_side
    tofill['number_triggering_muon_hnlside'] = event.number_triggering_muon_hnlside
    tofill['atleast_1triggeringmuon_otherside'] = event.atleast_1triggeringmuon_otherside
    tofill['mu_fromHNL_only_satisfies_BParkHLT_cond_tagside'] = event.mu_fromHNL_only_satisfies_BParkHLT_cond_tagside
    tofill['onetriggering_mu_fromB_tagside'] = event.onetriggering_mu_fromB_tagside
    tofill['onetriggering_mu_fromHNL_tagside'] = event.onetriggering_mu_fromHNL_tagside

    # invariant mass
    tofill['mu_pi_invmass' ] = event.the_hnldaughters.mass()
    tofill['mu_HNL_X_invmass'] = event.the_bdaughters.mass()
    
    # hnl charge
    tofill['hnl_q'] = event.the_hn.lep.charge() + event.the_hn.pi.charge()   
    
    # displacement
    if len(the_pls) and len(the_lep_daughters):
      tofill['Lxy'] = event.Lxy
      tofill['Lxyz'] = event.Lxyz
      #tofill['Lxy_cos'] = event.cos_pointing
    tofill['Lxyz_l0'] = event.Lxyz_pl
      
    # reconstruction weight
    tofill['weight_reco'] = event.weight_reco

    # ct weight
    #tofill['weight_ct'] = event.weight_ct
    
    ntuple.Fill(array('f',tofill.values()))

  outfile.cd()
  ntuple.Write()
  outfile.Close()


if __name__ == "__main__":
  #version_label = 'V34_newfilter_genstudy_v3'
  #version_label = 'V34_newfilter_genstudy_Bc_v1'
  #version_label = 'V38_request_Bc'
  version_label = 'test_displacementFilter1e9mm'
  user = 'anlyon'

  indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}'.format(user, version_label)
  # get all the subdirectories (signal points)
  pointdirs = [f for f in glob.glob('{}/*'.format(indirectory))]

  for pointdir in pointdirs:
    print '----------------------'
    print ' Analysing {}'.format(pointdir)
    print '----------------------'

    #if pointdir != '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V38_request/mass4.5_ctau1.0': continue

    pointname = pointdir[pointdir.rfind('/')+1:].replace('.', 'p')

    infiles = '{}/step1*root'.format(pointdir)
    outdir = './outputfiles/{}'.format(version_label)
    outfilename = 'genTree_{}.root'.format(pointname)

    if not path.exists(outdir):
      os.system('mkdir -p {}'.format(outdir))

    treeProducer(infiles=infiles, outdir=outdir, outfilename=outfilename)

