import os
import sys
import numpy as np
from ROOT import TTree, TFile, TH1F, TEfficiency, TGraph, TCanvas, gROOT, TAxis, TMath, TLegend, gStyle, gPad, TLine, TGraphAsymmErrors, TGraphErrors, TChain
import ROOT
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *
from glob import glob
from array import array
import ratioplot as RP
# couplings to be tested, for which the reweight is run

#from genTreeProducer import new_vvs
#small_new_vvs = new_vvs[::4]  # only select one every three elements
#small_new_vvs.reverse()
#new_vvs.reverse()
#small_new_vvs = new_vvs[0:10]

from python.common import *
from python.decays import Decays

#from python.old_common import BR__B_l_N, BR__B_D0_l_N, BR__B_l_nu, BR__B_D0_l_nu, # BR_HNLmupion
global graph_saver
graph_saver=[]


def getOverflowedHisto(h):
  htemp = h.Clone()
  nbins = htemp.GetNbinsX()
  last_plus_overflow = htemp.GetBinContent(nbins) + htemp.GetBinContent(nbins+1)
  last_plus_overflow_error = TMath.Sqrt( htemp.GetBinError(nbins)*htemp.GetBinError(nbins)  + htemp.GetBinError(nbins+1)*htemp.GetBinError(nbins+1))
  htemp.SetBinContent(nbins,last_plus_overflow )
  htemp.SetBinError(nbins,last_plus_overflow_error)
  return htemp


class SampleQuantity(object):
  def __init__(self, name='vv', axis='x', title='|V|^{2}', log=True, Range=(0.,1), forceRange=False, err=False):
    self.name = name
    self.axis = axis
    self.title = title
    self.log = log
    self.Range = Range
    self.forceRange = forceRange
    self.err = err


class PlotOpt(object):
  def __init__(self, what, binning, xtitle, ytitle, logX, logY):
    self.what = what
    self.binning = binning
    self.xtitle = xtitle
    self.ytitle = ytitle
    self.logX = logX
    self.logY = logY

class Sample(object):
  '''
  Handles information of a single sample, i.e. mass,ctau , w/ or w/o reweighting
  '''
  def __init__(self,mass=-99, ctau=-99, vv=-99, infileNames=None, is_ctau_rw=False, is_reco_rw=False, orig_vv=None, label='V00', leglabel='', myfilterEff=None):
    self.mass = mass
    self.ctau = ctau
    self.vv = vv
    if vv==-99: 
      self.vv=getVV(mass=self.mass, ctau=self.ctau)
    if ctau==-99:
      self.ctau=getCtau(mass=self.mass, vv=self.vv)
    self.treeName='tree'
    self.infileNames=infileNames
    #if not os.path.isfile(self.infileName): raise RuntimeError('file %s not found' % self.infileName)
    self.is_ctau_rw=is_ctau_rw
    self.is_reco_rw=is_reco_rw
    if self.is_ctau_rw:
      self.orig_vv = orig_vv
      self.orig_ctau = getCtau(mass=self.mass,vv=orig_vv)
    else:
      self.orig_vv = self.vv
      self.orig_ctau = self.ctau
    self.label=label
    self.leglabel=leglabel
    self.ngenevts = -999  #int(label.split('_n')[1])
    self.njobs =    -999  #int(label.split('_njt')[1])
    #self.isMajorana
    self.name='bhnl_mass{m}_ctau{ctau}'.format(m=self.mass, ctau=self.ctau)
    #self.legname='{:3}: m={:.1f}GeV |V|^{{2}}={:.1e} c#tau={:.1f}mm {} {}'.format('rw' if self.is_ctau_rw else 'gen',self.mass,self.vv,self.ctau,'- orig |V|^{{2}}={:.1e}'.format(self.orig_vv) if self.is_ctau_rw else '', self.leglabel)
    #self.legname='{:3}: m={:.1f}GeV |V|^{{2}}={:.1e} c#tau={:.1f}mm'.format('rw' if self.is_ctau_rw else 'gen',self.mass,self.vv,self.ctau)
    #self.legname='{:3}: m={:.1f}GeV |V|^{{2}}={:.1e} {}'.format('rw' if self.is_ctau_rw else 'gen',self.mass,self.vv, self.leglabel)
    self.legname='{:3}: m={:.1f}GeV c#tau={:.1e}m {}'.format('rw' if self.is_ctau_rw else 'gen',self.mass,self.ctau/1000., self.leglabel)
    if self.is_ctau_rw:
      self.evt_w = '(weight_{vv})'.format(vv=str(self.vv).replace('-', 'm')) # reweight the tree (created with orig vv) to the vv of this sample
    elif self.is_reco_rw:
      #self.evt_w = '(weight_reco)*(0.4)'
      self.evt_w = '(weight_reco)' # removed the trigger efficiency because it is already taken into account
    else:
      self.evt_w = '(1)'
   
    self.myfilterEff = myfilterEff
 
    self.accB = {}
    self.accB_errup = {}
    self.accB_errdn = {}
    self.acc = -99
    self.acc_errup = -99
    self.acc_errdn = -99
    self.effFilter = -99
    self.effFilter_errup = -99
    self.effFilter_errdn = -99
    
    self.Bspecies = ['B', 'B0', 'Bs', 'Bc']
    self.BpdgIds  = [521, 511,  531,  541 ]
   
    self.histoDefs = {
    # b mother
    'b_pt'          : PlotOpt('b_pt', '(60,0,60)', 'B meson p_{T} [GeV]', 'a.u.', False, True),
    'b_signid'      : PlotOpt('b_pdgid/abs(b_pdgid)', '(2,-1.1,1.1)', 'B meson pdg Id sign', 'a.u.', False, False),   
    'b_pdgid'       : PlotOpt('abs(b_pdgid)', '(50,500,550)', 'B meson pdg Id', 'a.u.', False, False),   
    'b_eta'         : PlotOpt('b_eta', '(60,-6,6)', 'B meson #eta', 'a.u.', False, True),
    'b_ct'          : PlotOpt('b_ct_reco', '(50,0,100)', 'B meson ct [mm]', 'a.u.', False, True),
    'b_charge'      : PlotOpt('b_q', '(3,-1.5,1.5)', 'B meson charge', 'a.u.', False, False),   
    #'b_ct_large'    : PlotOpt('b_ct_reco', '(100,0,1000)', 'B meson ct [mm]', 'a.u.', False, True),     
    # daughters of the B
    ## the HNL
    'hnl_pzOverP'   : PlotOpt('abs(TMath::TanH(hnl_eta))', '(50,0,1)', 'HNL |p_{z}|/|p| [GeV]', 'a.u.', False, False),   
    'hnl_ptOverP'   : PlotOpt('1/abs(TMath::CosH(hnl_eta))', '(50,0,1)', 'HNL p_{T}/|p| [GeV]', 'a.u.', False, False),   
    'hnl_pz'        : PlotOpt('abs(hnl_pt * TMath::SinH(hnl_eta))', '(60,0,30)', 'HNL |pz| [GeV]', 'a.u.', False, True),   
    'hnl_pt'        : PlotOpt('hnl_pt', '(60,0,30)', 'HNL p_{T} [GeV]', 'a.u.', False, True),   
    'hnl_signid'    : PlotOpt('hnl_pdgid/abs(hnl_pdgid)', '(2,-1.1,1.1)', 'HNL pdg id sign', 'a.u.', False, False),   
    'hnl_eta'       : PlotOpt('hnl_eta', '(60,-6,6)', 'HNL #eta', 'a.u.', False, True),      
    'hnl_betagamma' : PlotOpt('hnl_beta*hnl_gamma', '(80,0,20)', 'HNL #beta#gamma', 'a.u.', False, True),
    'hnl_ct'        : PlotOpt('hnl_ct_reco', '(50,0,1000)', 'HNL ct [mm]', 'a.u.', False, True),
    'hnl_ct_large'  : PlotOpt('hnl_ct_reco', '(100,0,10000)', 'HNL ct [mm]', 'a.u.', False, True),    
    'hnl_Lz'        : PlotOpt('Lz', '(50,0,1000)', 'L_{z} [mm]', 'a.u.', False, True),
    'hnl_Lz_large'  : PlotOpt('Lz', '(100,0,10000)', 'L_{z} [mm]', 'a.u.', False, True),
    'hnl_Lxy'       : PlotOpt('Lxy', '(50,0,1000)', 'L_{xy} [mm]', 'a.u.', False, True),
    'hnl_Lxy_large' : PlotOpt('Lxy', '(100,0,10000)', 'L_{xy} [mm]', 'a.u.', False, True),    
    'hnl_Lxyz'       : PlotOpt('Lxyz', '(50,0,1000)', 'L_{xyz} [mm]', 'a.u.', False, True),
    'hnl_Lxyz_large' : PlotOpt('Lxyz', '(100,0,10000)', 'L_{xyz} [mm]', 'a.u.', False, True),    
    'hnl_charge'    : PlotOpt('hnl_q', '(3,-1.5,1.5)', 'HNL charge', 'a.u.', False, False),   

    ### the D meson
#    'd_pt'          : PlotOpt('d_pt', '(60,0,30)', 'D meson p_{T} [GeV]', 'a.u.', False, True),   
#    'd_eta'         : PlotOpt('d_eta', '(60,-6,6)', 'D meson #eta', 'a.u.', False, True),      

    ### the trigger lepton
    'mutrig_pt'     : PlotOpt('l0_pt', '(60,0,30)', '#mu^{trig} p_{T} [GeV]', 'a.u.', False, True),   
    'mutrig_eta'    : PlotOpt('l0_eta', '(60,-6,6)', '#mu^{trig} #eta', 'a.u.', False, True),      
    'mutrig_charge' : PlotOpt('l0_q', '(3,-1.5,1.5)', '#mu^{trig} charge', 'a.u.', False, False),   

    #### daughters of the D meson
    ###### the pion
#    'piD_pt'        : PlotOpt('pi_pt', '(60,0,30)', '#pi (from D) p_{T} [GeV]', 'a.u.', False, True),  
#    'piD_eta'       : PlotOpt('pi_eta', '(60,-6,6)', '#pi (from D) #eta', 'a.u.', False, True),      
    
    ###### the kaon
#    'k_pt'          : PlotOpt('k_pt', '(60,0,30)', 'K (from D) p_{T} [GeV]', 'a.u.', False, True),  
#    'k_eta'         : PlotOpt('k_eta', '(60,-6,6)', 'K (from D) #eta', 'a.u.', False, True),      
    
    #### daughters of the HNL
    ###### the lepton
    'mu_pt'         : PlotOpt('(l1_pt)*(abs(l1_pdgid)==13)', '(60,0,30)', '#mu (from HNL) p_{T} [GeV]', 'a.u.', False, True),  
    'mu_eta'        : PlotOpt('(l1_eta)*(abs(l1_pdgid)==13)', '(60,-6,6)', '#mu (from HNL) #eta', 'a.u.', False, True),      
    'mu_charge'     : PlotOpt('(l1_q)*(abs(l1_pdgid)==13)', '(3,-1.5,1.5)', '#mu (from HNL) charge', 'a.u.', False, False),   
    
    'el_pt'         : PlotOpt('(l1_pt)*(abs(l1_pdgid)==11)', '(60,0,30)', 'e (from HNL) p_{T} [GeV]', 'a.u.', False, True),  
    'el_eta'        : PlotOpt('(l1_eta)*(abs(l1_pdgid)==11)', '(60,-6,6)', 'e (from HNL) #eta', 'a.u.', False, True),      
    'el_charge'     : PlotOpt('(l1_q)*(abs(l1_pdgid)==11)', '(3,-1.5,1.5)', 'e (from HNL) charge', 'a.u.', False, False),   

    'l1_pdgid'      : PlotOpt('l1_pdgid','(31, -15.5, 15.5)', 'lep (from HNL) pdg Id', 'a.u.', False, False),

    ###### the pion
    'pi_pt'         : PlotOpt('pi1_pt', '(60,0,30)', '#pi (from HNL) p_{T} [GeV]', 'a.u.', False, True),   
    'pi_eta'        : PlotOpt('pi1_eta', '(60,-6,6)', '#pi (from HNL) #eta', 'a.u.', False, True),      
    'pi_charge'     : PlotOpt('pi1_q', '(3,-1.5,1.5)', '#pi (from HNL) charge', 'a.u.', False, False),   
  
    # charge of the two leptons
    'l1mutrig_charge': PlotOpt('l0_q*l1_q', '(3,-1.5,1.5)', '#mu^{trig} charge x #mu (from HNL) charge', 'a.u.', False, False),
    'l1pi_charge'    : PlotOpt('l1_q*pi1_q','(3,-1.5,1.5)', '#pi charge x #mu  charge', 'a.u.', False, False),

    # invariant masses
    'hnl_invmass'   : PlotOpt('lep_pi_invmass', '(50,0,5)', 'HNL invariant mass, m(#mu,#pi) [GeV]', 'a.u.', False, False),     
#    'd_invmass'     : PlotOpt('k_pi_invmass', '(50,0,5)', 'D meson invariant mass, m(K,#pi) [GeV]', 'a.u.', False, False),      
    #'b_invmass'     : PlotOpt('b_invmass', '(50,2,7)', 'B invariant mass [GeV]', 'a.u.', False, False),     
    'mulep_invmass'    : PlotOpt('lep_mu_invmass', '(70,0,7)', 'm(#ell,#mu^{trig}) [GeV]', 'a.u.', False, False),     
    'bpartial_invmass'     : PlotOpt('bpartial_invmass', '(70,0,7)', 'm(HNL,#mu^{trig}) [GeV]', 'a.u.', False, False),     
    'bpartial_invmass_log'     : PlotOpt('bpartial_invmass', '(70,0,7)', 'm(HNL,#mu^{trig}) [GeV]', 'a.u.', False, True),     
    #'Lxy_cos', # cosine of the pointing angle in the transverse plane
    #'Lxyz_b', #3D displacement of the B wrt to primary vertex
    #'Lxyz_l0' #3D displacement of the prompt lepton wrt to B vertex

    # deltaRs 
    'dR_l1pi'      : PlotOpt('sqrt(TVector2::Phi_mpi_pi(l1_phi-pi1_phi)*TVector2::Phi_mpi_pi(l1_phi-pi1_phi) + (l1_eta-pi1_eta)*(l1_eta-pi1_eta))', \
                              '(50,0.,5.)', '#DeltaR(#mu,#pi)', 'a.u.', False, False),
    'dR_mutrighnl' : PlotOpt('sqrt(TVector2::Phi_mpi_pi(hnl_phi-l0_phi)*TVector2::Phi_mpi_pi(hnl_phi-l0_phi) + (hnl_eta-l0_eta)*(hnl_eta-l0_eta))', \
                              '(50,0.,5.)', '#DeltaR(#mu^{trig},HNL)', 'a.u.', False, False),
    }
    self.histos = {}
   
  def stamp(self):
    '''
    This should print the basic information of the sample
    '''
    if debug: print('Info in Sample.stamp()')
    print('mass={m}GeV, ctau={ctau}mm VV={vv}, is_ctau_rw={is_ctau_rw}, orig_VV={ovv}, orig_ctau={octau}mm, acc={acc}, au={au}, ad={ad}, num={num}, den={den}, evt_w={ew} effFilter={ef}'.format( \
            m=self.mass,ctau=self.ctau,vv=self.vv,is_ctau_rw=self.is_ctau_rw,ovv=self.orig_vv,acc=self.acc,au=self.acc_errup,ad=self.acc_errdn,num=self.num,den=self.den,ew=self.evt_w,ef=self.effFilter,octau=self.orig_ctau))

    #print('TGraphErrors') 
    #print('mass={m}GeV, ctau={ctau}mm VV={vv}, is_ctau_rw={is_ctau_rw}, orig_VV={ovv}, orig_ctau={octau}mm, acc={acc}, au={au}, ad={ad}, num={num}, den={den}, evt_w={ew} effFilter={ef}'.format( \
    #        m=self.mass,ctau=self.ctau,vv=self.vv,is_ctau_rw=self.is_ctau_rw,ovv=self.orig_vv,acc=self.acc_tg,au=self.acc_errup_tg,ad=self.acc_errdn_tg,num=self.num,den=self.den,ew=self.evt_w,ef=self.effFilter,octau=self.orig_ctau))
 
  def fillHistos(self,sel='(1)',sellabel='noSel'):
    '''
    This is to fill the histograms
    '''
    if debug: print('Info in Sample.fillHistos()')
    chain = TChain(self.treeName)
    for fn in self.infileNames:
      chain.Add(fn)

    for name,spec in self.histoDefs.items():
      chain.Draw('%s>>%s%s' % (spec.what, name, spec.binning), sel+'*'+self.evt_w, 'goff') 
      h = chain.GetHistogram().Clone('%s_%s_%s' % (self.name, name, sellabel))
      h.SetDirectory(0)
      h.SetTitle(h.GetName())
      h.GetXaxis().SetTitle(spec.xtitle)
      h.GetXaxis().SetNdivisions(505)
      h.GetYaxis().SetTitle(spec.ytitle)
      h.GetYaxis().SetNdivisions(505)
      h.SetDirectory(0)
      self.histos[name] = h

    return True

  def fillHistoStats(self):
    '''
    To determine the stats
    '''
    if debug: print('Info in Sample.fillHistoStats()')
    chain = TChain(self.treeName)
    for fn in self.infileNames:
      chain.Add(fn)

    bins = [0.,10,50,100,300,500,1000] 
    hns = ROOT.TH1F('hns', 'w/o selection, w/o reco weight', len(bins)-1, array('d',bins))
    h   = ROOT.TH1F('h',   'w/  selection, w/o reco weight', len(bins)-1, array('d',bins))
    hw  = ROOT.TH1F('hw',  'w/  selection, w/  reco weight', len(bins)-1, array('d',bins))

    sel = '(l0_pt>7 && abs(l0_eta)<1.5 && pi1_pt > 0.7 && abs(pi1_eta) < 2. && l1_pt > 1.5 && abs(l1_eta) < 2.)'#
          # && sqrt(TVector2::Phi_mpi_pi(hnl_phi-l0_phi)*TVector2::Phi_mpi_pi(hnl_phi-l0_phi) + (hnl_eta-l0_eta)*(hnl_eta-l0_eta)) < 0.5'

    sel_and_weight = '({sel})*(weight_reco)*(0.4*0.8)'.format(sel=sel) # 0.8 pre-selection efficiency, 0.4 selection efficiency
    weight = '(weight_reco)*(0.4*0.8)'
    #sel_and_weight = '({sel})*(weight_reco)'.format(sel=sel) 
    #weight         = '(weight_reco)'                         
    #sel_and_weight = '({sel})*(weight_reco)'.format(sel=sel) 
    #weight         = '(weight_reco)'                         

    chain.Draw('Lxy>>hns', '(1)',          'goff')
    chain.Draw('Lxy>>h',   sel,            'goff')
    #chain.Draw('Lxy>>h',   weight,            'goff')
    chain.Draw('Lxy>>hw',  sel_and_weight, 'goff')

    hns.GetXaxis().SetTitle('L_{xy} (mm)')
    hns.GetYaxis().SetTitle('Entries')
    
    hns.SetLineColor(ROOT.kBlack)
    h.SetLineColor(ROOT.kBlue)
    hw.SetLineColor(ROOT.kRed)

    c = TCanvas(self.name, self.name, 800, 600)
    hns.Draw('histE')
    h.Draw('histEsame')
    hw.Draw('histEsame')
    c.BuildLegend()
    #c.SaveAs('./plots/' + self.label + suffix + '/stats_' + c.GetName()  + '.png')
    #c.SaveAs('./plots/' + self.label + suffix + '/stats_' + c.GetName()  + '.pdf')
    #c.SaveAs('./plots/' + self.label + suffix + '/stats_' + c.GetName()  + '.root')

    ## Print the values bin by bin,
    histos = [hns,h,hw]
    lines = []

    histos = [hw]
    for i,h in enumerate(histos):
      line = []
      tot_entries = h.Integral() + h.GetBinContent(h.GetNbinsX()+1)
      line.append('{:15.3f}'.format(tot_entries))
      for binx in range(1,h.GetNbinsX()+2):
        #line.append('{:15.1f}'.format(h.GetBinContent(binx)))
        line.append('{:15.3f}'.format(h.GetBinContent(binx)/float(tot_entries)))
      lines.append(' '.join(line))

    #print('{:15s} {:15s} {:15s} {:15s} {:15s} {:15s} {:15s}'.format('Integral', '[0,1]cm', '[1,5]cm', '[5,10]cm','[10,30]cm', '[30,100]cm', '[100,inf]cm')   )
    #for line in lines:
    #  print(line)


    tot_hns = hns.Integral() + hns.GetBinContent(hns.GetNbinsX()+1)
    #tot_hw  = hw.Integral()  + hw.Integral() + hw.GetBinContent(hw.GetNbinsX()+1)   ### ERROR!!!
    tot_hw  =  hw.Integral() + hw.GetBinContent(hw.GetNbinsX()+1)                    ### CORRECTED VERSION
    tot_eff = float(tot_hw) / float(tot_hns) 
    #print('{:15.1f} {:15.1f} {:15.1f}%'.format(tot_hns,tot_hw,tot_eff*100))
    print('{:15.4f}'.format(tot_eff))

    return True

  def saveHistos(self, norm):
    '''
    This is to save the histograms of the sample in a plot
    '''
    if debug: print('Info in Sample.saveHistos()')
    c = TCanvas(self.name, self.name, 800, 600)
    c.DivideSquare(len(self.histos.keys()))

    for i, what in enumerate(self.histos.keys()):
      h = self.histos[what]
      pad = c.cd(i + 1)
      pad.SetLogx(self.histoDefs[what].logX)
      pad.SetLogy(self.histoDefs[what].logY)
      if norm and h.Integral() != 0:
        #hdrawn = h.DrawNormalized('LPE')
        h.Scale(1/h.Integral(-1,h.GetNbinsX()+1))
        hdrawn = h.Draw('PLE')
      else:
        hdrawn = h.Draw('PLE')

      norm_suffix='_norm' if norm else ''

    c.SaveAs('./plots/' + self.label + suffix + '/' + c.GetName() + norm_suffix + '.png')
    c.SaveAs('./plots/' + self.label + suffix + '/' + c.GetName() + norm_suffix + '.pdf')
    c.SaveAs('./plots/' + self.label + suffix + '/' + c.GetName() + norm_suffix + '.C')

  def fillAcceptance(self):
    '''
    This is to calculate the acceptance of the sample
    '''
    if debug: print('Info in Sample.fillAcceptance()')
    #f = TFile.Open(self.infileName)
    #t = f.Get(self.treeName)
    #if not t:
    #  raise RuntimeError( 'ERROR: no tree in file %s' % self.infileName)
    chain = TChain(self.treeName)
    for fn in self.infileNames:
      chain.Add(fn)

    self.effnum  = ROOT.TH1F('effnum', 'effnum', 1, 0, 13000) 
    self.effden  = ROOT.TH1F('effden', 'effden', 1, 0, 13000)
    self.effnumB = dict([(b, ROOT.TH1F('effnum_%s' % b, 'effnum_%s' % b, 1, 0, 13000)) for b in self.Bspecies])
    self.effdenB = dict([(b, ROOT.TH1F('effden_%s' % b, 'effden_%s' % b, 1, 0, 13000)) for b in self.Bspecies])

    # preselectoin cuts
#    cutsnum = '(l0_pt>7 && abs(l0_eta)<1.5 && pi1_pt > 0.7 && abs(pi1_eta) < 2. && l1_pt > 1.5 && abs(l1_eta) < 2.)'#
#    cutsden = '(1)'
    ##### cutsden = '(l0_pt>7 && abs(l0_eta)<1.5)'

    
    ### cuts used for the july study
    cutsnum = '(l0_pt>{mp} && abs(l0_eta)<1.5'.format(mp=muTrigPt)
    if doSkipDispl:
      cutsnum += ''
    else:
      cutsnum += '&& Lxy < 500'
    if not doDisplZ:
      pass
    else:
      cutsnum += '&& Lz < 20'
    if doSkipHNLptEta:
      cutsnum += ''
    else:
      cutsnum += '&& l1_pt>3 && abs(l1_eta)<2.5 && pi1_pt>0.8 && abs(pi1_eta)<2.5'
    
    cutsnum += ')'
    cutsden = '(l0_pt>{mp} && abs(l0_eta)<1.5)'.format(mp=muTrigPt)

    
    ###### fill the total acceptance
    chain.Draw('hnl_pt>>effnum', cutsnum+'*'+self.evt_w, 'goff')
    chain.Draw('hnl_pt>>effden', cutsden, 'goff')
 
    if TEfficiency.CheckConsistency(self.effnum,self.effden): 
      peff = TEfficiency(self.effnum,self.effden)
     
      # check usage of TGraphAsymmetricErrors
      self.acc = peff.GetEfficiency(1)
      self.acc_errup = peff.GetEfficiencyErrorUp(1)
      self.acc_errdn = peff.GetEfficiencyErrorLow(1)

    # for debugging purposes
    self.num = self.effnum.GetEntries()
    if self.num==0: print('**** 0 entries for mass={}'.format(self.mass))
    self.den = self.effden.GetEntries()

    ###### fill the partial acceptances 
    for ib,b in enumerate(self.Bspecies):
      bsel   = '(abs(b_pdgid)=={bid})'.format(bid=self.BpdgIds[ib])
      selnum = '(' + cutsnum + '&&' + bsel + ')' 
      selden = '(' + cutsden + '&&' + bsel + ')'
      chain.Draw('hnl_pt>>effnum_{b}'.format(b=b), selnum + '*' + self.evt_w, 'goff')
      chain.Draw('hnl_pt>>effden_{b}'.format(b=b), selden, 'goff')

#      tgra = TGraphAsymmErrors()
#      tgra.BayesDivide(self.effnumB[b], self.effdenB[b])
#      print type(tgra.GetY()), self.effdenB[b]
#      if tgra.GetY():
#        self.accB[b] = tgra.GetY()[0]
#        self.accB_errup[b] = tgra.GetErrorYhigh(0)
#        self.accB_errdn[b] = tgra.GetErrorYlow(0)
#      else:
#        self.accB[b] = 0
#        self.accB_errup[b] = 0
#        self.accB_errdn[b] = 0

      if TEfficiency.CheckConsistency(self.effnumB[b],self.effdenB[b]):
         peff = TEfficiency(self.effnumB[b],self.effdenB[b])
        
         self.accB[b]       = peff.GetEfficiency(1)
         #self.accB_errup[b] = peff.GetEfficiencyErrorUp(1)
         #self.accB_errdn[b] = peff.GetEfficiencyErrorLow(1)
         self.accB_errup[b] = 0
         self.accB_errdn[b] = 0

      else:
         self.accB[b] = 0
         self.accB_errup[b] = 0
         self.accB_errdn[b] = 0

  def fillExpNevts(self):
    if debug: print('Info in Sample.fillExpNevts()')
    
    #N_B = 3.84E9 * 2. # 3.84E9, inclusive, = n(B+/-) * 2  
    # factor ~2 is for considering B0s as well, which are ~50% of the B-parking dataset
    
    N_B = 9.2001E9  # from slide 5 here http://cds.cern.ch/record/2704495/files/DP2019_043.pdf
        
    dec    = Decays(mass=p.mass, mixing_angle_square=1)
    dec_SM = Decays(mass=0., mixing_angle_square=1)

    B_w =        (dec.B_to_uHNL.BR+       \
                  dec.B_to_D0uHNL.BR+     \
                  dec.B_to_D0staruHNL.BR+ \
                  dec.B_to_pi0uHNL.BR+    \
                  dec.B_to_rho0uHNL.BR)   \
               / (dec_SM.B_to_uHNL.BR+       \
                  dec_SM.B_to_D0uHNL.BR+     \
                  dec_SM.B_to_D0staruHNL.BR+ \
                  dec_SM.B_to_pi0uHNL.BR+    \
                  dec_SM.B_to_rho0uHNL.BR)   

    B0_w =       (dec.B0_to_DuHNL.BR+     \
                  dec.B0_to_DstaruHNL.BR+ \
                  dec.B0_to_piuHNL.BR+    \
                  dec.B0_to_rhouHNL.BR)   \
               / (dec_SM.B0_to_DuHNL.BR    + \
                  dec_SM.B0_to_DstaruHNL.BR+ \
                  dec_SM.B0_to_piuHNL.BR   + \
                  dec_SM.B0_to_rhouHNL.BR)  \
   
    Bs_w =       (dec.Bs_to_DsuHNL.BR+ \
                  dec.Bs_to_DsstaruHNL.BR+ \
                  dec.Bs_to_KuHNL.BR+ \
                  dec.Bs_to_KstaruHNL.BR)  \
               / (dec_SM.Bs_to_DsuHNL.BR+ \
                  dec_SM.Bs_to_DsstaruHNL.BR+ \
                  dec_SM.Bs_to_KuHNL.BR+ \
                  dec_SM.Bs_to_KstaruHNL.BR) 
    
    Bc_w =     (dec.Bc_to_uHNL.BR) / (dec_SM.Bc_to_uHNL.BR)

    self.Bweights = [B_w, B0_w, Bs_w, Bc_w]
    self.Bfracs =   [0.4, 0.4, 0.1, 0.001]

    self.expNevts = 0
    self.expNevts_errup = 0
    self.expNevts_errdn = 0
    for ib,b in enumerate(self.Bspecies):
      self.expNevts        += N_B * self.Bfracs[ib] * self.Bweights[ib] * BR_HNLmupion(mass=self.mass) * self.vv * self.accB[b]
      self.expNevts_errup  += N_B * self.Bfracs[ib] * self.Bweights[ib] * BR_HNLmupion(mass=self.mass) * self.vv * self.accB_errup[b]
      self.expNevts_errdn  += N_B * self.Bfracs[ib] * self.Bweights[ib] * BR_HNLmupion(mass=self.mass) * self.vv * self.accB_errdn[b]
       
    #N_HNL_VV1 = N_B * float(BR_B_to_XuN) / float(BR_B_to_Xun)
    #self.expNevts       = N_HNL_VV1 * 2 * BR_HNLmupion(mass=self.mass) * self.vv * self.acc               # the factor 2 for BR_HNLmupion was needed in old_common.py
    #self.expNevts_errup = N_HNL_VV1 * 2 * BR_HNLmupion(mass=self.mass) * self.vv * self.acc_errup         # because of an inconsistency between gamma_part and gamma_tot
    #self.expNevts_errdn = N_HNL_VV1 * 2 * BR_HNLmupion(mass=self.mass) * self.vv * self.acc_errdn

  def fillExpNevts2(self):
    if debug: print('Info in Sample.fillExpNevts2()')
     
    #N_B = 9.2001E9  # from slide 5 here http://cds.cern.ch/record/2704495/files/DP2019_043.pdf
    # done with x-sec x lumi x filter_eff x acc x eff
    Lumi = 41.6 * 1E09 # 1/fb = 1E09 * 1/microbarn
    sigma_tot = 15.3 * 2 # microbarn, factor 2 for B+/B-  
    if self.myfilterEff:
      norm = Lumi * sigma_tot * self.myfilterEff
    else:
      print('WARNING: Filter eff not available, set to 1 for normalisation')
      norm = Lumi * sigma_tot

    # in microbarn, differentially in pT, not used
    sigma={}
    sigma[(10,13)] = 3.00      * 3    # eta < 1.45
    sigma[(13,17)] = 0.88      * 4    # ""
    sigma[(17,24)] = 0.39      * 6    # eta < 2.1
    sigma[(24,30)] = 0.12      * 6    # ""
    sigma[(30,40)] = 4.10E-02  * 10
    sigma[(40,50)] = 1.20E-02  * 10
    sigma[(50,60)] = 4.00E-03  * 10
    sigma[(60,70)] = 1.90E-03  * 10
    sigma[(70,100)]= 6.70E-04  * 30
   
    dec    = Decays(mass=p.mass, mixing_angle_square=1)
    dec_SM = Decays(mass=0., mixing_angle_square=1)

    B_w =        (dec.B_to_uHNL.BR+       \
                  dec.B_to_D0uHNL.BR+     \
                  dec.B_to_D0staruHNL.BR+ \
                  dec.B_to_pi0uHNL.BR+    \
                  dec.B_to_rho0uHNL.BR)  

    B0_w =       (dec.B0_to_DuHNL.BR+     \
                  dec.B0_to_DstaruHNL.BR+ \
                  dec.B0_to_piuHNL.BR+    \
                  dec.B0_to_rhouHNL.BR)   
   
    Bs_w =       (dec.Bs_to_DsuHNL.BR+ \
                  dec.Bs_to_DsstaruHNL.BR+ \
                  dec.Bs_to_KuHNL.BR+ \
                  dec.Bs_to_KstaruHNL.BR)  
    
    Bc_w =     (dec.Bc_to_uHNL.BR) 

    self.Bweights = [B_w, B0_w, Bs_w, Bc_w]
    self.Bfracs =   [0.4, 0.4, 0.1, 0.001]

    self.expNevts = 0
    self.expNevts_errup = 0
    self.expNevts_errdn = 0
    for ib,b in enumerate(self.Bspecies):
      #print norm, self.Bfracs[ib], self.Bweights[ib], BR_HNLmupion(mass=self.mass), self.vv, self.accB[b]
      self.expNevts        += norm * self.Bfracs[ib] / 0.4 * self.Bweights[ib] * BR_HNLmupion(mass=self.mass) * self.vv * self.accB[b]
      self.expNevts_errup  += norm * self.Bfracs[ib] / 0.4 * self.Bweights[ib] * BR_HNLmupion(mass=self.mass) * self.vv * self.accB_errup[b]
      self.expNevts_errdn  += norm * self.Bfracs[ib] / 0.4 * self.Bweights[ib] * BR_HNLmupion(mass=self.mass) * self.vv * self.accB_errdn[b]
       
  def fillFilterEff(self,dostamp=True):
    '''
    To retrieve and save the filter efficiency - from the minigentree # OUTDATED
    '''
    ##if debug: print('Info in Sample.fillFilterEff()')

    efffnum = ROOT.TH1F('efffnum', 'efffnum', 1, 0, 13000)
    efffden = ROOT.TH1F('efffden', 'efffden', 1, 0, 13000)

    chain = TChain(self.treeName)
    for fn in self.infileNames:
      chain.Add(fn)
    
    efffnum.SetBinContent(1, chain.GetEntries())
    efffnum.SetBinError(1, ROOT.TMath.Sqrt(chain.GetEntries()))
    
    # denominator = number of events that were run in the first place # access storage element... 
    #self.ngenevts_succ_afterfilter=0
    self.njobs_succ=0
    path = '/pnfs/psi.ch/cms/trivcat/store/user/{u}/BHNLsGen/{pl}/mass{m}_ctau{ctau}/step1*root'.format(u=os.environ['USER'],pl=self.label,m=self.mass,ctau=self.orig_ctau)
    for fname in glob(path):
        #if debug: print 'fname=',fname
        f = TFile.Open(fname)
        if f.GetListOfKeys().Contains('Events'):
          self.njobs_succ += 1

    #njobs_succ = len(glob(path))
    self.ngenevts_succ = float(self.ngenevts) / float(self.njobs) * float(self.njobs_succ)
    
    efffden.SetBinContent(1,self.ngenevts_succ)
    efffden.SetBinError(1,ROOT.TMath.Sqrt(self.ngenevts_succ))

    if TEfficiency.CheckConsistency(efffnum,efffden): 
      geneff = TEfficiency(efffnum,efffden) 
      self.effFilter = geneff.GetEfficiency(1)
      self.effFilter_errup = geneff.GetEfficiencyErrorUp(1)
      self.effFilter_errdn = geneff.GetEfficiencyErrorLow(1)

    # stamp basic info about filter eff
    if dostamp: print('mass={m}GeV, VV={vv:.1e}, nTot={nt}, nSucc={ns}, nSuccAfterFilter={nsf}, effFilter={ef:.3f}%, errup={eu:.3f}%, errdn={ed:.3f}% '.format( \
                m=self.mass,vv=self.vv,nt=self.ngenevts,ns=self.ngenevts_succ,nsf=efffnum.GetBinContent(1),ef=self.effFilter*100,eu=self.effFilter_errup*100,ed=self.effFilter_errdn*100))

class SampleList(object):
  '''
  Handles the plotting of several samples, both the histograms, and the general quantities (graphs)
  '''
  def __init__(self, name, samples, label=None):
    self.name = name
    self.samples = samples  # a list of objects of type Sample
    self.label = samples[0].label if label==None else label
    self.colors = [  ROOT.kOrange+1, ROOT.kRed, ROOT.kMagenta+2, ROOT.kViolet+8, ROOT.kAzure-8, ROOT.kAzure+6 ,
                      ROOT.kGreen+1, ROOT.kSpring+4, ROOT.kYellow -5, ROOT.kYellow -3, ROOT.kYellow, ROOT.kOrange,
                    ROOT.kBlack,  ROOT.kOrange+1, ROOT.kRed, ROOT.kMagenta+2, ROOT.kViolet+8, ROOT.kAzure-8, ROOT.kAzure+6 ,
                      ROOT.kGreen+1, ROOT.kSpring+4, ROOT.kYellow -5, ROOT.kYellow -3, ROOT.kYellow, ROOT.kOrange,
                  ]
    #self.colors = [ ROOT.kOrange+1, ROOT.kRed, ROOT.kMagenta+2, ROOT.kViolet+8, ROOT.kAzure-8, ROOT.kAzure+6 ,
    #                  ROOT.kGreen+1, ROOT.kSpring+4, ROOT.kYellow -5, ROOT.kYellow -3, ROOT.kYellow, ROOT.kOrange
    #              ]
    self.styles = [  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, ]
    self.graphs={}

    # x,y quantities for graphs
    self.quantities={
      #'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=True, Range=(0.,1), forceRange=True, err=True),
      'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=True, Range=(0.00001,2), forceRange=True, err=True),
      #'acc' : SampleQuantity(name='acc', axis='y', title='Acceptance', log=False, Range=(0.,1), forceRange=True, err=True),
      'expNevts': SampleQuantity(name='expNevts', axis='y', title='N_{#nu} x BR(HN#rightarrow#mu#pi) x Acc x V^{2}', log=True, err=True),
      'num' : SampleQuantity(name='num', axis='y', title='Numerator', log=True, Range=(1, 3000), forceRange=True, err=False),
      'den' : SampleQuantity(name='den', axis='y', title='Denominator', log=True, Range=(1, 3000), forceRange=True, err=False),
      #'filterEff' : 
      # 'xsec'
      'vv'  : SampleQuantity(name='vv', axis='x', title='|V|^{2}', log=True),
      'ctau': SampleQuantity(name='ctau', axis='x', title='c#tau [mm]', log=False),
      'mass': SampleQuantity(name='mass', axis='x', title='mass [GeV]', log=False),
    } 

    if len(self.samples) > len(self.colors):
      raise RuntimeError('need more colors (%d vs %d)' % (len(self.colors), len(self.samples)))
    self.checkConsistency()
  def add(self, sample):
    self.samples.insert(sample)
    if len(self.samples) > len(self.colors):
      raise RuntimeError('need more colors (%d vs %d)' % (len(self.colors), len(self.samples)))
    self.checkConsistency()
  def remove(self, sample):
    self.samples.remove(sample)
  def checkConsistency(self):
    '''
    Makes sure histograms of samples have same binning and histograms
    '''
    for sample1 in self.samples:
      for sample2 in self.samples:
        if len(sample1.histoDefs) != len(sample2.histoDefs):
          raise RuntimeError('number of PlotOpts differs for %s and %s' % (sample1.name, sample2.name))
        if len(sample1.histos) != len(sample2.histos):
          raise RuntimeError('number of histos differs for %s and %s' % (sample1.name, sample2.name))

  def plotHistos(self, norm=False, sameCanvas=False):
    ''' 
    Superimpose histograms for all samples, on the same canvas or on separate canvases
    '''
    if sameCanvas:
      c = TCanvas(self.name, self.name, 6400, 4800)
      tosave = { c.GetName() : c }
    else:
      c = None
      tosave = {}

    legends = []

    hMax = {}
    for j, what in enumerate(self.samples[0].histos.keys()):
      hMax[what] = []
      leg=defaultLegend(x1=0.15,y1=0.7,x2=0.95,y2=0.92,mult=1.2)
      legends.append(leg)
    # do the actual plotting
    for i, sample in enumerate(self.samples):
      if i == 0:
        opt = 'histE'
        if sameCanvas: c.DivideSquare(len(sample.histos.keys()))
      else: opt = 'histEsame'
      for j, what in enumerate(sample.histos.keys()):
        h = sample.histos[what]
        graph_saver.append(h)

        # kolmogorov test
        #if i == 0:
        #  h2 = self.samples[1].histos[what]
        #  h.KolmogorovTest(h2, 'OD')
        #

        if sameCanvas: pad = c.cd(j + 1)
        else:
          cname = '%s_%s' % (self.name, what)
          if cname not in tosave.keys():
            tosave[cname] = TCanvas(cname, cname, 700,600)
          pad = tosave[cname].cd()

        if i == 0:
          pad.SetLogx(sample.histoDefs[what].logX)
          pad.SetLogy(sample.histoDefs[what].logY)
        if norm and h.Integral()!=0:
          # new way with "Scale"
          h.Scale(1/h.Integral(-1,h.GetNbinsX()+1)) 
#          # usual way
#          hdrawn = h.DrawNormalized(opt)
#          hMax[what].append(hdrawn)
#          hdrawn.SetLineColor(self.colors[i])
#          hdrawn.SetMarkerColor(self.colors[i])
#          hdrawn.SetMarkerStyle(self.styles[i])
#          hdrawn.SetMarkerSize(3)
#          hdrawn.SetLineWidth(2)
#          legends[j].AddEntry(hdrawn,sample.legname, 'LP')
#          #legends[j].AddEntry(hdrawn,sample.legname + ' Mean={:.2f}'.format(hdrawn.GetMean()), 'LP')
#          #legends[j].AddEntry(hdrawn,sample.legname + ' Acc={:.2f}'.format(sample.acc), 'LP')

#        else:
        h.Draw(opt)
        hMax[what].append(h)
        h.SetLineColor(self.colors[i])
        #hMax[what][0].GetXaxis().SetTitle(h.GetXaxis().GetTitle())
        #hMax[what][0].GetYaxis().SetTitle(h.GetYaxis().GetTitle())
        h.SetMarkerColor(self.colors[i])
        h.SetMarkerStyle(self.styles[i])
        h.SetMarkerSize(4)
        h.SetLineWidth(2)
        #legends[j].AddEntry(h,sample.legname, 'LP')
        legends[j].AddEntry(h,sample.legname + ' Mean={:.2f}'.format(h.GetMean()), 'LP')
        #legends[j].AddEntry(h,sample.legname + ' Acc={:.2f}'.format(sample.acc), 'LP')

        if i == len(self.samples) - 1:
          legends[j].Draw('same')

    # compute the max of norm histo and set the range accordingly
    for j, what in enumerate(self.samples[0].histos.keys()):
      if len(hMax[what])>0:
        max = 0
        for ih in hMax[what]:
          if(max < ih.GetMaximum()):
            max=ih.GetMaximum()
        if (sample.histoDefs[what].logY==True):
          #hMax[name][0].GetYaxis().SetRangeUser(0.01, max * 7)
          hMax[what][0].SetMaximum(max*30)
          #print hMax[what][0].GetMaximum()
        else:
          #hMax[name][0].GetYaxis().SetRangeUser(0.01, max * 1.5)
          hMax[what][0].SetMaximum(max * 1.7)

    norm_suffix='_norm' if norm else ''

    for cname, canv in tosave.items():
      canv.SaveAs('./plots/' + self.label + suffix + '/' + cname + norm_suffix + '.png')
      canv.SaveAs('./plots/' + self.label + suffix + '/' + cname + norm_suffix + '.pdf')
      canv.SaveAs('./plots/' + self.label + suffix + '/' + cname + norm_suffix + '.C')
      canv.SaveAs('./plots/' + self.label + suffix + '/' + cname + norm_suffix + '.root')
  
  def plotRatios(self):
    ''' 
    On the same canva, plot ratio for the first two distributions
    ''' 
    snum = self.samples[0]
    sden = self.samples[1]
    sden2 = self.samples[2]

    for j, what in enumerate(snum.histos.keys()):     
      hNum = snum.histos[what]
      hDen = sden.histos[what]      
      graph_saver.append(hNum)
      graph_saver.append(hDen)
      if sden2:
        hDen2 = sden2.histos[what]
        graph_saver.append(hDen2)

      if norm and hNum.Integral()!=0 and hDen.Integral()!=0:
        hNum.Scale(1/hNum.Integral(-1,hNum.GetNbinsX()+1)) 
        hDen.Scale(1/hDen.Integral(-1,hDen.GetNbinsX()+1)) 
      if norm and sden2 and hDen2.Integral()!=0:
        hDen2.Scale(1/hDen2.Integral(-1,hDen2.GetNbinsX()+1))    
  
      norm_suffix='_norm' if norm else ''
      outDir = './plots/' + self.label + suffix
      cname = '%s_%s' % (self.name, what)
      plotName =  cname + norm_suffix + '_ratio'

      nameNum = snum.legname + ' Mean={:.2f}'.format(hNum.GetMean())
      nameDen = sden.legname + ' Mean={:.2f}'.format(hDen.GetMean())
      nameDen2 = sden2.legname + ' Mean={:.2f}'.format(hDen2.GetMean()) if sden2 else ""

      hDen2_to_pass = hDen2 if sden2 else ""
      RP.makeRatioPlot(hNum, hDen, hDen2=hDen2_to_pass, nameNum=nameNum, nameDen=nameDen, nameDen2=nameDen2, \
                    xtitle=snum.histoDefs[what].xtitle,ytitle=snum.histoDefs[what].ytitle, \
                    ratiotitle="Ratio", norm=False, log=snum.histoDefs[what].logY,plotName=plotName,outDir=outDir,ratioyrange=(0,2.)) 
      # handle normalisation before RP.makeRatioPlot()

  def plotGraph(self, x='vv', y='acc'): 
    '''
    Plot a graph with specified quantities on x and y axes , and saves the graph
    '''

    if (x not in self.quantities.keys()) or (y not in self.quantities.keys()):
      raise RuntimeError('selected quantities not available, available quantities are: \n{}'.format(self.quantities.keys()))

    xq = self.quantities[x]
    yq = self.quantities[y]

    #graph = TGraphAsymmErrors()
    #graph = TGraph()
    graph = TGraphErrors()
    for i,s in enumerate(self.samples):
      graph.SetPoint(i,getattr(s, xq.name), getattr(s, yq.name) )
      #if xq.err: 
      #  graph.SetPointEXhigh(i, getattr(s, xq.name+'_errup'))   # errup errdn
      #  graph.SetPointEXlow (i, getattr(s, xq.name+'_errdn'))
      if yq.err: 
        graph.SetPointError(i, 0, getattr(s, yq.name+'_errup'))
      #  graph.SetPointEYhigh(i, getattr(s, yq.name+'_errup'))  
      #  graph.SetPointEYlow (i, getattr(s, yq.name+'_errdn'))

    c = TCanvas()
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(22)
    graph.SetTitle(';{x};{y}'.format(y=yq.title,x=xq.title))
    graph.Draw('APLE')

    if yq.forceRange:
      graph.SetMinimum(yq.Range[0])
      graph.SetMaximum(yq.Range[1])

    gPad.Modified()
    gPad.Update()
    if yq.name=='expNevts':
      line = TLine(gPad.GetUxmin(),3,gPad.GetUxmax(),3)
      line.SetLineColor(ROOT.kBlue)
      line.Draw('same')
      #graph.SetMinimum(0.01)
      #graph.SetMaximum(1E06)

    if xq.log: c.SetLogx()
    if yq.log: c.SetLogy()
    c.SetGridx()
    c.SetGridy()
    c.SaveAs('./plots/{}{}/{}_{}VS{}.pdf'.format(self.label,suffix,self.name,yq.name,xq.name))
    c.SaveAs('./plots/{}{}/{}_{}VS{}.C'.format(self.label,suffix,self.name,yq.name,xq.name))
    c.SaveAs('./plots/{}{}/{}_{}VS{}.png'.format(self.label,suffix,self.name,yq.name,xq.name))

    self.graphs['{}VS{}'.format(yq.name,xq.name)] = graph
    # add the graph container for memory?
    graph_saver.append(graph)

def doAnalysis(path,pl,points,name,leglabel='',path2=None,pl2=None,points2=None,leglabel2='',path3=None,pl3=None,points3=None,leglabel3=''):
  '''
  Perform plotting of samples lists
  Optionally allows to create sample list out of 2 different productions
  TODO: quantities to plot in graph as an option...
  '''
  print('\n*************************************************')
  print('  => Going to do plotting for production={}, name={}'.format(pl,name))
  samples = []
  for p in points:
    fns = glob(path.format(m=p.mass,ctau=p.orig_ctau))
    s = Sample(mass=p.mass, ctau=p.ctau, vv=p.vv, infileNames=fns, is_ctau_rw=p.is_ctau_rw, is_reco_rw=p.is_reco_rw, orig_vv=p.orig_vv, 
               label=pl, leglabel=leglabel,myfilterEff=p.myfilterEff)
    if doRwAnalysis:
      s.fillHistos(sel='(Lxy<1000)', sellabel='')
    elif doGridAnalysis:
      s.fillHistoStats()
      s.fillHistos()
    elif doFixedMassAnalysis:
      pass 
      # do not fill histograms, only acceptance is needed
    else:
      s.fillHistos() # no selection

    if not doGridAnalysis:
      s.fillAcceptance()
      #s.fillExpNevts2()
      s.fillExpNevts()
    #s.fillFilterEff()
    #s.stamp()
    samples.append(s)

  if path2 and pl2 and points2:
    print('\n*************************************************')
    print('  => Going to do plotting for production={}, name={}'.format(pl2,name))
    for p in points2:
      fns = glob(path2.format(m=p.mass,ctau=p.orig_ctau))
      s = Sample(mass=p.mass, ctau=p.ctau, vv=p.vv, infileNames=fns, is_ctau_rw=p.is_ctau_rw, orig_vv=p.orig_vv, label=pl2, leglabel=leglabel2) 
      s.fillHistos()
      s.fillAcceptance()
      s.fillExpNevts2()
      #s.fillFilterEff()
      s.stamp()
      samples.append(s)

  if path3 and pl3 and points3:
    print('\n*************************************************')
    print('  => Going to do plotting for production={}, name={}'.format(pl3,name))
    for p in points3:
      fns = glob(path3.format(m=p.mass,ctau=p.orig_ctau))
      s = Sample(mass=p.mass, ctau=p.ctau, vv=p.vv, infileNames=fns, is_ctau_rw=p.is_ctau_rw, orig_vv=p.orig_vv, label=pl3, leglabel=leglabel3) 
      s.fillHistos()
      s.fillAcceptance()
      s.fillExpNevts2()
      #s.fillFilterEff()
      s.stamp()
      samples.append(s)


  label = pl if not pl2 else ('{}_VS_{}'.format(pl,pl2) if not pl3 else '{}_VS_{}_VS_{}'.format(pl,pl2,pl3))

  slist = SampleList(name=name, samples=samples, label=label)
  if doCompareAnalysis or doRwAnalysis:
    slist.plotRatios()
    #slist.plotHistos(norm=norm, sameCanvas=False)
  else:
    slist.plotHistos(norm=norm, sameCanvas=False)
  #if 'fixedM' in name or 'closure' in name: 
  if 'fixedM' in name:
    slist.plotGraph(x='vv',y='acc')
    slist.plotGraph(x='vv',y='expNevts')
  elif 'fixedVV' in name:
    slist.plotGraph(x='mass',y='acc')
    slist.plotGraph(x='mass',y='num')
    slist.plotGraph(x='mass',y='den')
    ###slist.plotGraph(x='mass',y='expNevts')

  return slist


def doGraphComparison(list1,list2,what):

  c = TCanvas()
  list1.graphs[what].Draw('APL')
  list2.graphs[what].Draw('same')
  c.SaveAs(what+'.pdf')
  c.SaveAs(what+'.png')

def doLimitGraph(slists,what, label):

  c = TCanvas()
  graph = TGraph()
  for ilist,slist in enumerate(slists):
    # get the inverted graph
    inverted_graph = TGraph()
    xs = slist.graphs[what].GetX()
    ys = slist.graphs[what].GetY()
    N =  slist.graphs[what].GetN()
    xs_arr = np.ndarray(N,'d',xs) 
    ys_arr = np.ndarray(N,'d',ys) 
    inverted_graph = TGraph(N,ys,xs)
    # use Eval
    limit = inverted_graph.Eval(3.) # linear interpolation between closest points
    graph.SetPoint(ilist,slist.samples[0].mass,limit)
  graph.SetLineWidth(2)
  #if doInclusive: 
  graph.SetMarkerStyle(22)
  #else: graph.SetMarkerStyle(24)
  graph.SetTitle(';mass [GeV]; exp. limit on V^{2}')
  graph.SetMaximum(1E-02)
  graph.SetMinimum(1E-08)
  graph.Draw('APL')
  graph.GetXaxis().SetLimits(0.1, 10)
  graph.Draw('APL')
  c.Update()

  c.SetLogx() 
  c.SetLogy()
  #c.SetGridx()
  c.SetGridy()
  c.SaveAs('./plots/{}{}/limit_vv.pdf'.format(label,suffix))

def checkFiles(path,points):
  '''
  Check existence of files for analysis, return existing points
  '''
  for p in points:
    expr = path.format(m=p.mass,ctau=p.orig_ctau) # check the file where the tree is stored
    print('**********',expr)
    nfiles = 0
    for fn in glob(expr):
      print('**********',glob(expr))
      if os.path.isfile(fn):
        nfiles+=1
    if nfiles==0:
      p.missing=True
    else: 
      p.missing=False
 
  existing_points = [p for p in points if not p.missing]
  
  if len(existing_points)>0:
    return existing_points 
  else:
    raise RuntimeError('no files for analysis')
    return []


  print('  => Going to do analysis on following points:')
  print('  ')

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='', add_help=True)
  parser.add_argument('--pl', type=str, dest='pl', help='production label of input', default='V_blabla')  
  parser.add_argument('--pl2', type=str, dest='pl2', help='production label of second input', default=None)  
  parser.add_argument('--pl3', type=str, dest='pl3', help='production label of second input', default=None)  
  return parser.parse_args()

if __name__ == "__main__":

  #### ROOT Options
  gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()
  ROOT.TH1.StatOverflows(ROOT.kTRUE) # consider overflows for mean and rms calculation
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gStyle.SetTitleXOffset(1.1);
  gStyle.SetTitleYOffset(1.45);

  ## globals 
  global debug
  #global doInclusive
  global doSkipDispl
  global doSkipDisplZ
  global suffix
  global muTrigPt
  global norm

  #### options
  debug = False
  norm = True
  #doInclusive = True # 
  doSkipDispl = False #
  doDisplZ = False #
  doSkipHNLptEta = False
  doCompareAnalysis = False
  doTestAnalysis = False
  doFixedMassAnalysis = False
  doOld = False
  doRwAnalysis = False
  doFixedVVAnalysis = False
  doBenchmarkAnalysis = True
  doGridAnalysis = False
  doGridBc = False
  muTrigPt = 9 # 0 1 2 5 7 9
  ####

  suffix='' 
  if doSkipDispl: suffix += '_skipDispl'
  if doDisplZ: suffix += '_wLz'
  if doSkipHNLptEta: suffix += '_skipHNLptEta'
  #suffix += '_muTrigPt{mp}'.format(mp=muTrigPt)
  opt = getOptions()
  
  complabel = '' if not opt.pl2 else ('_VS_{}'.format(opt.pl2) if not opt.pl3 else '_VS_{}_VS_{}'.format(opt.pl2,opt.pl3))
  os.system('mkdir -p ./plots/{}{}{}'.format(opt.pl,complabel,suffix))

  path = './outputfiles/' + opt.pl + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  #path = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGenDump/' + opt.pl + "/mass{m}_ctau{ctau}"  + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  if opt.pl2: 
    path2 = './outputfiles/' + opt.pl2 + '/mass{m}_ctau{ctau}_miniGenTree*.root'
    #path2 = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGenDump/' + opt.pl2 + "/mass{m}_ctau{ctau}"  + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  if opt.pl3: 
    path3 = './outputfiles/' + opt.pl3 + '/mass{m}_ctau{ctau}_miniGenTree*.root'
    #path3 = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGenDump/' + opt.pl3 + "/mass{m}_ctau{ctau}"  + '/mass{m}_ctau{ctau}_miniGenTree*.root'
  

  if doCompareAnalysis:
    #points = [Point(mass=1.5,ctau=None,vv=1e-03,is_ctau_rw=False)]
    #points2 = [Point(mass=1.5,ctau=None,vv=1e-03,is_ctau_rw=False)]
    #points3 = [Point(mass=1.5,ctau=None,vv=1e-03,is_ctau_rw=False)]  
    #points = [Point(mass=2.0,ctau=None,vv=1.5e-05,is_ctau_rw=False)]
    #points2 = [Point(mass=2.0,ctau=None,vv=1.5e-05,is_ctau_rw=False)]
    #points  = [Point(mass=3.0,ctau=184.256851021,vv=None,is_ctau_rw=False)]
    points  = [Point(mass=3.0,ctau=184.0,vv=None,is_ctau_rw=False)]
    points2 = [Point(mass=3.0,ctau=184.0,vv=None,is_ctau_rw=False)]
    points3 = [Point(mass=3.0,ctau=184.0,vv=None,is_ctau_rw=False)]
    #points = [ Point(mass=2.0, ctau=100.0,  vv=None,is_reco_rw=False,),]
    #points2 =[ Point(mass=2.0, ctau=100.0,  vv=None,is_reco_rw=False,),]

    for p in points:
      p.stamp()
    existing_points=checkFiles(path=path,points=points)
    for p in points2:
      p.stamp()
    existing_points2=checkFiles(path=path2,points=points2)
    if not opt.pl3: 
      doAnalysis(path=path,pl=opt.pl,points=existing_points,name='comp_BvsBc',path2=path2,pl2=opt.pl2,points2=existing_points2,leglabel='B^{#pm},B^{0},B_{s}', leglabel2='B_{c}')
    else:
      for p in points3:
        p.stamp()
      existing_points3=checkFiles(path=path3,points=points3)
      doAnalysis(path=path,pl=opt.pl,points=existing_points,name='comp_valid',path2=path2,pl2=opt.pl2,points2=existing_points2,path3=path3,pl3=opt.pl3,points3=existing_points3, leglabel='prod:sara', leglabel2='prod:camilla', leglabel3='prod:mg')

  if doTestAnalysis:

    gStyle.SetOptStat("mreuo")
    #points = [Point(mass=1.0,ctau=None,vv=5e-04,is_ctau_rw=False)]
    #points = [Point(mass=1.5,ctau=None,vv=1e-03,is_ctau_rw=False)]
    #points = [Point(mass=3.0,ctau=None,vv=5e-05,is_ctau_rw=False)]
    #points = [Point(mass=2.0,ctau=None,vv=1.5e-05,is_ctau_rw=False)]
    #points = [Point(mass=2.0,ctau=None,vv=5.0e-05,is_ctau_rw=False)]
    points = [ Point(mass=3.0,ctau=None,vv=1.0e-05)] # same settings as V11_3_points.py

    for p in points:
      p.stamp()
    existing_points=checkFiles(path=path,points=points)
    doAnalysis(path=path,pl=opt.pl,points=existing_points,name='testpoint_norw')

  if doFixedMassAnalysis and doOld:
   
    slists_fixedMass = []

    ################
    points = [
      Point(mass=0.5,ctau=None,vv=1e-03,is_ctau_rw=True,orig_vv=8.6e-02),
      Point(mass=0.5,ctau=None,vv=5e-04,is_ctau_rw=True,orig_vv=8.6e-02),
      Point(mass=0.5,ctau=None,vv=1e-04,is_ctau_rw=True,orig_vv=8.6e-02),
      Point(mass=0.5,ctau=None,vv=5e-05,is_ctau_rw=True,orig_vv=8.6e-02),
      Point(mass=0.5,ctau=None,vv=1e-05,is_ctau_rw=True,orig_vv=8.6e-02),
      Point(mass=0.5,ctau=None,vv=5e-06,is_ctau_rw=True,orig_vv=8.6e-02),
      Point(mass=0.5,ctau=None,vv=1e-06,is_ctau_rw=True,orig_vv=8.6e-02),
      Point(mass=0.5,ctau=None,vv=5e-07,is_ctau_rw=True,orig_vv=8.6e-02),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass0.5_norw'))
    ################
    points = [
      Point(mass=1.0,ctau=None,vv=1e-03,is_ctau_rw=True,orig_vv=2.7e-03),
      Point(mass=1.0,ctau=None,vv=5e-04,is_ctau_rw=True,orig_vv=2.7e-03),
      Point(mass=1.0,ctau=None,vv=1e-04,is_ctau_rw=True,orig_vv=2.7e-03),
      Point(mass=1.0,ctau=None,vv=5e-05,is_ctau_rw=True,orig_vv=2.7e-03),
      Point(mass=1.0,ctau=None,vv=1e-05,is_ctau_rw=True,orig_vv=2.7e-03),
      Point(mass=1.0,ctau=None,vv=5e-06,is_ctau_rw=True,orig_vv=2.7e-03),
      Point(mass=1.0,ctau=None,vv=1e-06,is_ctau_rw=True,orig_vv=2.7e-03),
      Point(mass=1.0,ctau=None,vv=5e-07,is_ctau_rw=True,orig_vv=2.7e-03),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.0_norw'))
    ################
    points = [
      Point(mass=1.5,ctau=None,vv=1e-03,is_ctau_rw=True,orig_vv=3.5e-04),
      Point(mass=1.5,ctau=None,vv=5e-04,is_ctau_rw=True,orig_vv=3.5e-04),
      Point(mass=1.5,ctau=None,vv=1e-04,is_ctau_rw=True,orig_vv=3.5e-04),
      Point(mass=1.5,ctau=None,vv=5e-05,is_ctau_rw=True,orig_vv=3.5e-04),
      Point(mass=1.5,ctau=None,vv=1e-05,is_ctau_rw=True,orig_vv=3.5e-04),
      Point(mass=1.5,ctau=None,vv=5e-06,is_ctau_rw=True,orig_vv=3.5e-04),
      Point(mass=1.5,ctau=None,vv=1e-06,is_ctau_rw=True,orig_vv=3.5e-04),
      Point(mass=1.5,ctau=None,vv=5e-07,is_ctau_rw=True,orig_vv=3.5e-04),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.5_norw'))
    ################
    points = [
      Point(mass=2.0,ctau=None,vv=5e-04,is_ctau_rw=True,orig_vv=8.4e-05),
      Point(mass=2.0,ctau=None,vv=1e-04,is_ctau_rw=True,orig_vv=8.4e-05),
      Point(mass=2.0,ctau=None,vv=5e-05,is_ctau_rw=True,orig_vv=8.4e-05),
      Point(mass=2.0,ctau=None,vv=1e-05,is_ctau_rw=True,orig_vv=8.4e-05),
      Point(mass=2.0,ctau=None,vv=5e-06,is_ctau_rw=True,orig_vv=8.4e-05),
      Point(mass=2.0,ctau=None,vv=1e-06,is_ctau_rw=True,orig_vv=8.4e-05),
      Point(mass=2.0,ctau=None,vv=5e-07,is_ctau_rw=True,orig_vv=8.4e-05),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass2.0_norw'))
    ################
    points = [
      Point(mass=3.0,ctau=None,vv=1e-03,is_ctau_rw=True,orig_vv=1.1e-05),
      Point(mass=3.0,ctau=None,vv=5e-04,is_ctau_rw=True,orig_vv=1.1e-05),
      Point(mass=3.0,ctau=None,vv=1e-04,is_ctau_rw=True,orig_vv=1.1e-05),
      Point(mass=3.0,ctau=None,vv=5e-05,is_ctau_rw=True,orig_vv=1.1e-05),
      Point(mass=3.0,ctau=None,vv=1e-05,is_ctau_rw=True,orig_vv=1.1e-05),
      Point(mass=3.0,ctau=None,vv=5e-06,is_ctau_rw=True,orig_vv=1.1e-05),
      Point(mass=3.0,ctau=None,vv=1e-06,is_ctau_rw=True,orig_vv=1.1e-05),
      Point(mass=3.0,ctau=None,vv=5e-07,is_ctau_rw=True,orig_vv=1.1e-05),
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slist_norw = doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass3.0_norw')
    slists_fixedMass.append(slist_norw)
    ################
    points = [
      Point(mass=4.5,ctau=None,vv=1e-03,is_ctau_rw=True,orig_vv=2.0e-04),
      Point(mass=4.5,ctau=None,vv=5e-04,is_ctau_rw=True,orig_vv=2.0e-04),
      Point(mass=4.5,ctau=None,vv=1e-04,is_ctau_rw=True,orig_vv=2.0e-04),
      Point(mass=4.5,ctau=None,vv=5e-05,is_ctau_rw=True,orig_vv=2.0e-04),
      Point(mass=4.5,ctau=None,vv=1e-05,is_ctau_rw=True,orig_vv=2.0e-04),
      Point(mass=4.5,ctau=None,vv=5e-06,is_ctau_rw=True,orig_vv=2.0e-04),
      Point(mass=4.5,ctau=None,vv=1e-06,is_ctau_rw=True,orig_vv=2.0e-04),
      Point(mass=4.5,ctau=None,vv=5e-07,is_ctau_rw=True,orig_vv=2.0e-04),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass4.5_norw'))

    doLimitGraph(slists_fixedMass,'expNevtsVSvv', label=opt.pl)

  if doFixedMassAnalysis and not doOld:
   
    slists_fixedMass = []

    ################
    points = [
        Point(mass=1.0, ctau=  10000.0,  vv=None,is_reco_rw=True, myfilterEff=2.24E-04), # filter effs are useless...
        Point(mass=1.0, ctau=   1000.0,  vv=None,is_reco_rw=True, myfilterEff=1.79E-03),
        Point(mass=1.0, ctau=    100.0,  vv=None,is_reco_rw=True, myfilterEff=5.72E-03),
        Point(mass=1.0, ctau=     10.0,  vv=None,is_reco_rw=True, myfilterEff=6.57E-03),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.0_norw'))
    ################
    points = [
        Point(mass=1.5, ctau=  10000.0,  vv=None,is_reco_rw=True, myfilterEff=1.07E-04),
        Point(mass=1.5, ctau=   1000.0,  vv=None,is_reco_rw=True, myfilterEff=8.83E-04),
        Point(mass=1.5, ctau=    100.0,  vv=None,is_reco_rw=True, myfilterEff=3.09E-03),
        Point(mass=1.5, ctau=     10.0,  vv=None,is_reco_rw=True, myfilterEff=3.56E-03),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.5_norw'))
    ################
    points = [
        Point(mass=2.0, ctau=  10000.0,  vv=None,is_reco_rw=True, myfilterEff=4.40E-05),
        Point(mass=2.0, ctau=   1000.0,  vv=None,is_reco_rw=True, myfilterEff=3.66E-04),
        Point(mass=2.0, ctau=    100.0,  vv=None,is_reco_rw=True, myfilterEff=1.34E-03),
        Point(mass=2.0, ctau=     10.0,  vv=None,is_reco_rw=True, myfilterEff=1.57E-03),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass2.0_norw'))
    ################
    points = [
        Point(mass=3.0, ctau=   1000.0,  vv=None,is_reco_rw=True, myfilterEff=1.93E-03),
        Point(mass=3.0, ctau=    100.0,  vv=None,is_reco_rw=True, myfilterEff=5.15E-03),
        Point(mass=3.0, ctau=     10.0,  vv=None,is_reco_rw=True, myfilterEff=5.46E-03),
        Point(mass=3.0, ctau=      1.0,  vv=None,is_reco_rw=True, myfilterEff=5.47E-03),
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slist_norw = doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass3.0_norw')
    slists_fixedMass.append(slist_norw)
    ################
    points = [
        Point(mass=4.5, ctau=    100.0,  vv=None,is_reco_rw=True, myfilterEff=3.53E-04),
        Point(mass=4.5, ctau=     10.0,  vv=None,is_reco_rw=True, myfilterEff=4.59E-04),
        Point(mass=4.5, ctau=      1.0,  vv=None,is_reco_rw=True, myfilterEff=4.59E-04),
        Point(mass=4.5, ctau=      0.1,  vv=None,is_reco_rw=True, myfilterEff=4.58E-04  ),
      
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slists_fixedMass.append(doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass4.5_norw'))

    doLimitGraph(slists_fixedMass,'expNevtsVSvv', label=opt.pl)

  ################
  if doRwAnalysis:

    points = [
      #Point(mass=2.0,ctau=None,vv=1.5e-05,is_ctau_rw=False),
      Point(mass=3.0,ctau=None,vv=1e-06,is_ctau_rw=True,orig_vv=2.0e-06),
      Point(mass=3.0,ctau=None,vv=3e-06,is_ctau_rw=True,orig_vv=2.0e-06),
      Point(mass=3.0,ctau=None,vv=2e-05,is_ctau_rw=True,orig_vv=2.0e-06),
      Point(mass=3.0,ctau=None,vv=7e-05,is_ctau_rw=True,orig_vv=2.0e-06),
      Point(mass=3.0,ctau=None,vv=3e-04,is_ctau_rw=True,orig_vv=2.0e-06),
      #Point(mass=2.0,ctau=None,vv=0.0002,is_ctau_rw=True,orig_vv=1.5e-05),
    ]
    for p in points:
      p.stamp()
    existing_points=checkFiles(path=path,points=points)
    slist_rw = doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedMass1.0')
  
    ## compare the graphs w/ and w/o reweighting
    #doGraphComparison(slist_norw,slist_rw,what='accVSvv')
    #doGraphComparison(slist_norw,slist_rw,what='expNevtsVSvv')

    ###############
    points = [
      #Point(mass=2.0,ctau=None,vv=3e-06,is_ctau_rw=True,orig_vv=1e-04),
      #Point(mass=2.0,ctau=None,vv=3e-06,is_ctau_rw=False),
      Point(mass=2.0,ctau=None,vv=3.0e-06,is_ctau_rw=True,orig_vv=5.0e-05),
      Point(mass=2.0,ctau=None,vv=3.0e-06,is_ctau_rw=False),
    ]
    existing_points = checkFiles(path,points)
    #doAnalysis(path=path,pl=opt.pl,points=existing_points,name='closureRw_Mass2.0_VV3em06')


  if doFixedVVAnalysis:

    ################
    points = [
      # normal B species
      Point(mass=0.5,ctau=None,vv=1e-05,orig_vv=8.6e-02,is_ctau_rw=True),
      Point(mass=1.0,ctau=None,vv=1e-05,orig_vv=2.7e-03,is_ctau_rw=True),
      Point(mass=1.5,ctau=None,vv=1e-05,orig_vv=3.5e-04,is_ctau_rw=True),
      Point(mass=2.0,ctau=None,vv=1e-05,orig_vv=8.4e-05,is_ctau_rw=True),
      Point(mass=3.0,ctau=None,vv=1e-05,orig_vv=1.1e-05,is_ctau_rw=True),
      Point(mass=4.5,ctau=None,vv=1e-05,orig_vv=2.0e-04,is_ctau_rw=True),
      
      # Bc points
      #Point(mass=2.0,ctau=None,vv=8.4e-05,is_ctau_rw=False),
      #Point(mass=3.0,ctau=None,vv=1.1e-05,is_ctau_rw=False),
      ##Point(mass=3.0,ctau=None,vv=4.0e-04,is_ctau_rw=False),
      #Point(mass=4.5,ctau=None,vv=2.0e-04,is_ctau_rw=False),
      ##Point(mass=4.5,ctau=None,vv=6.0e-03,is_ctau_rw=False),
      #Point(mass=5.5,ctau=None,vv=1.0e-03,is_ctau_rw=False),
      #Point(mass=5.5,ctau=None,vv=3.0e-02,is_ctau_rw=False),

      #Point(mass=1.0,ctau=None,vv=1e-05,is_ctau_rw=False),
      #Point(mass=2.0,ctau=None,vv=1e-05,is_ctau_rw=False),
      #Point(mass=3.0,ctau=None,vv=1e-05,is_ctau_rw=False),
      #Point(mass=3.0,ctau=None,vv=1e-05,orig_vv=5e-05,is_ctau_rw=True),
    ]
    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    #doAnalysis(path=path,pl=opt.pl,points=existing_points,name='grid')
    doAnalysis(path=path,pl=opt.pl,points=existing_points,name='fixedVV1em05')
  

  if doBenchmarkAnalysis:
    ################
    points = [
    
##      Point(mass=0.5,ctau=None,vv=2e-05,orig_vv=8.6e-02,is_ctau_rw=True),
##      Point(mass=1.0,ctau=None,vv=8e-06,orig_vv=2.7e-03,is_ctau_rw=True),
#      Point(mass=1.5,ctau=None,vv=2e-06,orig_vv=3.5e-04,is_ctau_rw=True),
##      Point(mass=2.0,ctau=None,vv=2e-06,orig_vv=8.4e-05,is_ctau_rw=True),
##      Point(mass=3.0,ctau=None,vv=1e-05,orig_vv=1.1e-05,is_ctau_rw=True),
##      Point(mass=4.5,ctau=None,vv=1e-04,orig_vv=2.0e-04,is_ctau_rw=True),
#      Point(mass=0.5,ctau=None,vv=2e-05,is_ctau_rw=False),
#      Point(mass=1.0,ctau=None,vv=8e-06,is_ctau_rw=False),
#      Point(mass=1.5,ctau=None,vv=2e-06,is_ctau_rw=False),
#      Point(mass=2.0,ctau=None,vv=2e-06,is_ctau_rw=False),
#      Point(mass=3.0,ctau=None,vv=1e-05,is_ctau_rw=False),
#      Point(mass=4.5,ctau=None,vv=1e-04,is_ctau_rw=False),
      Point(mass=1.0,ctau=10000.0,vv=None,is_ctau_rw=False),
      Point(mass=3.0,ctau=100.0,vv=None,is_ctau_rw=False),
      Point(mass=4.5,ctau=1.0,vv=None,is_ctau_rw=False),

    ]

    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    doAnalysis(path=path,pl=opt.pl,points=existing_points,name='benchmark')
  
  if doGridAnalysis:
    ################

    if not doGridBc:
      points = [
    
##        Point(mass=0.5, ctau= 100000.0,  vv=None,is_reco_rw=True, myfilterEff=3.51E-05),
##        Point(mass=0.5, ctau=  10000.0,  vv=None,is_reco_rw=True, myfilterEff=3.34E-04),
##        Point(mass=0.5, ctau=   1000.0,  vv=None,is_reco_rw=True, myfilterEff=2.50E-03),
##        Point(mass=1.0, ctau= 100000.0,  vv=None,is_reco_rw=True, myfilterEff=2.34E-05),
        Point(mass=1.0, ctau=  10000.0,  vv=None,is_reco_rw=False, myfilterEff=2.24E-04),   # FILTER EFF are OBSOLETE
        Point(mass=1.0, ctau=   1000.0,  vv=None,is_reco_rw=False, myfilterEff=1.79E-03),
        Point(mass=1.0, ctau=    100.0,  vv=None,is_reco_rw=False, myfilterEff=5.72E-03),
        Point(mass=1.0, ctau=     10.0,  vv=None,is_reco_rw=False, myfilterEff=6.57E-03),
        Point(mass=1.5, ctau=  10000.0,  vv=None,is_reco_rw=False, myfilterEff=1.07E-04),
        Point(mass=1.5, ctau=   1000.0,  vv=None,is_reco_rw=False, myfilterEff=8.83E-04),
        Point(mass=1.5, ctau=    100.0,  vv=None,is_reco_rw=False, myfilterEff=3.09E-03),
        Point(mass=1.5, ctau=     10.0,  vv=None,is_reco_rw=False, myfilterEff=3.56E-03),
        Point(mass=2.0, ctau=  10000.0,  vv=None,is_reco_rw=False, myfilterEff=4.40E-05),
        Point(mass=2.0, ctau=   1000.0,  vv=None,is_reco_rw=False, myfilterEff=3.66E-04),
        Point(mass=2.0, ctau=    100.0,  vv=None,is_reco_rw=False, myfilterEff=1.34E-03),
        Point(mass=2.0, ctau=     10.0,  vv=None,is_reco_rw=False, myfilterEff=1.57E-03),
        Point(mass=3.0, ctau=   1000.0,  vv=None,is_reco_rw=False, myfilterEff=1.93E-03),
        Point(mass=3.0, ctau=    100.0,  vv=None,is_reco_rw=False, myfilterEff=5.15E-03),
        Point(mass=3.0, ctau=     10.0,  vv=None,is_reco_rw=False, myfilterEff=5.46E-03),
        Point(mass=3.0, ctau=      1.0,  vv=None,is_reco_rw=False, myfilterEff=5.47E-03),
        Point(mass=4.5, ctau=    100.0,  vv=None,is_reco_rw=False, myfilterEff=3.53E-04),
        Point(mass=4.5, ctau=     10.0,  vv=None,is_reco_rw=False, myfilterEff=4.59E-04),
        Point(mass=4.5, ctau=      1.0,  vv=None,is_reco_rw=False, myfilterEff=4.59E-04),
        Point(mass=4.5, ctau=      0.1,  vv=None,is_reco_rw=False, myfilterEff=4.58E-04  ),
      ]
    else:
       points = [

         Point(mass=2.0, ctau=   10000.0,  vv=None,is_reco_rw=False, myfilterEff=1.61E-02   ),  # FILTER EFF are OBSOLETE
         Point(mass=2.0, ctau=    1000.0,  vv=None,is_reco_rw=False, myfilterEff=2.21E-01   ),
         Point(mass=2.0, ctau=     100.0,  vv=None,is_reco_rw=False, myfilterEff=2.21E-01   ),
         Point(mass=2.0, ctau=      10.0,  vv=None,is_reco_rw=False, myfilterEff=2.34E-01   ),
         Point(mass=3.0, ctau=    1000.0,  vv=None,is_reco_rw=False, myfilterEff=7.99E-02   ),
         Point(mass=3.0, ctau=     100.0,  vv=None,is_reco_rw=False, myfilterEff=1.68E-01   ),
         Point(mass=3.0, ctau=      10.0,  vv=None,is_reco_rw=False, myfilterEff=1.69E-01   ),
         Point(mass=3.0, ctau=       1.0,  vv=None,is_reco_rw=False, myfilterEff=1.69E-01   ),
         Point(mass=4.5, ctau=     100.0,  vv=None,is_reco_rw=False, myfilterEff=3.63E-02   ),
         Point(mass=4.5, ctau=      10.0,  vv=None,is_reco_rw=False, myfilterEff=3.77E-02   ),
         Point(mass=4.5, ctau=       1.0,  vv=None,is_reco_rw=False, myfilterEff=3.81E-02   ),
         Point(mass=4.5, ctau=       0.1,  vv=None,is_reco_rw=False, myfilterEff=3.79E-02   ),
         Point(mass=5.5, ctau=      10.0,  vv=None,is_reco_rw=False, myfilterEff=1.28E-03   ),
         Point(mass=5.5, ctau=       1.0,  vv=None,is_reco_rw=False, myfilterEff=1.28E-03   ),
         Point(mass=5.5, ctau=       0.1,  vv=None,is_reco_rw=False, myfilterEff=1.28E-03   ),
         Point(mass=5.5, ctau=      0.01,  vv=None,is_reco_rw=False, myfilterEff=1.28E-03   ),
      ]                               

    for p in points:
     p.stamp()
    existing_points=checkFiles(path=path,points=points)
    doAnalysis(path=path,pl=opt.pl,points=existing_points,name='grid')


