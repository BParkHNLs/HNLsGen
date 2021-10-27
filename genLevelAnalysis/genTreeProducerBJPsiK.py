# http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf

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

    # the leading muon 1, daughter of JPsi
    'l1_pt',
    'l1_eta',
    'l1_phi',
    'l1_mass',
    'l1_q',
    'l1_pdgid',

    # the leading muon 2, daughter of JPsi
    'l2_pt',
    'l2_eta',
    'l2_phi',
    'l2_mass',
    'l2_q',
    'l2_pdgid',

    # the Kaon
    'k_pt',
    'k_eta',
    'k_phi',
    'k_mass',
    'k_q',
    'k_pdgid',
    
    # invariant masses
    'll_invmass',
    'llk_invmass',

]

###################
# Functions and classes
##################

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
   
    if event.the_b_mother==None: 
      print('B-mother in this event could not be found, skipping')
      continue

    b_p4 = event.the_b_mother.p4()
    event.the_b_mother.y = 0.5 * ROOT.TMath.Log( (b_p4.E()+b_p4.Pz()) / (b_p4.E()-b_p4.Pz()) ) if (b_p4.E()-b_p4.Pz())!=0 else -99
 
    # get the other daughters of the B meson
    event.the_b_mother.daughters = [event.the_b_mother.daughter(jj) for jj in range(event.the_b_mother.numberOfDaughters())]

    if opt.doDebug: 
      stampParticle(event.the_b_mother)
      stampDaugthers(event.the_b_mother)

    # the Kaon
    the_ks = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId()) == 321], key = lambda x : x.pt(), reverse=True)
    if len(the_ks) !=0:
      event.the_k = the_ks[0]

    # the JPsi
    the_JPsis = sorted([ii for ii in event.the_b_mother.daughters if abs(ii.pdgId()) == 443], key = lambda x : x.pt(), reverse=True) 
    if len(the_JPsis) !=0:
      event.the_JPsi = the_JPsis[0]

    # the JPsi daughter leptons
    the_JPsi_leptons = sorted([event.the_JPsi.daughter(jj) for jj in range(event.the_JPsi.numberOfDaughters()) \
                               if abs(event.the_JPsi.daughter(jj).pdgId()) in [11, 13, 15]], key = lambda x: x.pt(), reverse=True)
    
    if len(the_JPsi_leptons) !=0:
      event.the_l1 = the_JPsi_leptons[0]
    if len(the_JPsi_leptons) > 0:
      event.the_l2 = the_JPsi_leptons[1]
    
    ### Filling 
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
    
    if event.the_l1:
       tofill['l1_pt'      ] = event.the_l1.pt()
       tofill['l1_eta'     ] = event.the_l1.eta()
       tofill['l1_phi'     ] = event.the_l1.phi()
       tofill['l1_mass'    ] = event.the_l1.mass()
       tofill['l1_q'       ] = event.the_l1.charge()
       tofill['l1_pdgid'   ] = event.the_l1.pdgId()

    if event.the_l2:
       tofill['l2_pt'      ] = event.the_l2.pt()
       tofill['l2_eta'     ] = event.the_l2.eta()
       tofill['l2_phi'     ] = event.the_l2.phi()
       tofill['l2_mass'    ] = event.the_l2.mass()
       tofill['l2_q'       ] = event.the_l2.charge()
       tofill['l2_pdgid'   ] = event.the_l2.pdgId()
 
    if event.the_k:   
       tofill['k_pt'      ] = event.the_k.pt()
       tofill['k_eta'     ] = event.the_k.eta()
       tofill['k_phi'     ] = event.the_k.phi()
       tofill['k_mass'    ] = event.the_k.mass()
       tofill['k_q'       ] = event.the_k.charge()
       tofill['k_pdgid'   ] = event.the_k.pdgId()
 
    if event.the_l1 and event.the_l2: 
      tofill['ll_invmass' ] = (event.the_l1.p4() + event.the_l2.p4()).mass()
      if event.the_k:
        tofill['llk_invmass'] = (event.the_l1.p4() + event.the_l2.p4() + event.the_k.p4()).mass()  
    
    
    ##print(tofill.values())        
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

