import ROOT
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import sys
sys.path.append('../python')
from decays import Decays
from common import getVV, BR_HNLmupion


class Quantity(object):
  def __init__(self, name='', title='', label='', nbins=0, bin_min=0., bin_max=0.):
    self.name = name
    self.title = title
    self.label = label
    self.nbins = nbins
    self.bin_min = bin_min
    self.bin_max = bin_max


class SignalPoint(object):
  def __init__(self, mass, ctau):
    self.mass = mass
    self.ctau = ctau


class Sample(object):
  def __init__(self, filename, mass, ctau, filter_efficiency):
    self.filename = filename
    self.mass = mass
    self.ctau = ctau
    self.filter_efficiency = filter_efficiency


def plotTagRate(version_label):
  print '----------------------'
  print 'Will plot tag rate for ntuples {}'.format(version_label)
  print '----------------------'

  hnl_tag_side = Quantity('hnl_tag_side', 'HNL on the tag side', 'hnl_tag_side', 2, 0, 2)

  filenames = [f for f in glob.glob('./outputfiles/{}/*'.format(version_label))]

  # get the masses
  masses = []
  for filename in filenames:
    mass = filename[filename.rfind('mass')+4:filename.rfind('_ctau')]
    if mass not in masses:
      masses.append(mass)

  # get the signal points
  signal_points = []
  for mass in masses:
    for filename in filenames:
      if 'mass'+mass not in filename: continue
      ctau = filename[filename.rfind('ctau')+4:filename.rfind('.root')]
      signal_points.append(SignalPoint(mass=mass, ctau=ctau))
      
  # get the tag and probe rates
  tag_rates = {}
  for mass in masses:
    for signal_point in signal_points:
      if signal_point.mass != mass: continue
      ctau = signal_point.ctau

      for filename in filenames:
        if mass not in filename or ctau not in filename: continue

        f = ROOT.TFile.Open(filename, 'READ')
        tree = f.Get('tree')

        hist_name = 'hist_tag_side_{}_{}'.format(mass, ctau)
        hist = ROOT.TH1D(hist_name, hist_name, hnl_tag_side.nbins, hnl_tag_side.bin_min, hnl_tag_side.bin_max)
        tree.Project(hist_name, 'hnl_tag_side', '')
        hist.SetDirectory(0)

        tag_rates['mass{}_ctau{}'.format(mass, ctau)] = [hist.GetMean(), 1-hist.GetMean()]

  # plot the rates
  the_df = pd.DataFrame(tag_rates, index=['tag_side', 'probe_side'])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='Pastel2', figsize=(7, 10))#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #ax.set_title('Mass 4.5 GeV')
  for ictau, tag in enumerate(tag_rates):
    widths_tag = the_df['tag_side']
    widths_probe = the_df['probe_side']
    xpos_tag = 0.2 # widths_tag[ictau]/2.
    xpos_probe = 1 - widths_probe[ictau]/2.
    ypos = ictau
    ax.text(xpos_tag, ypos, round(widths_tag[ictau], 2), ha='center', va='center',color='black')
    #ax.text(xpos_probe, ypos, round(widths_probe[ictau], 2), ha='center', va='center',color='black')
  ax.figure.savefig('plots/tag_side_{}.png'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('plots/tag_side_{}.pdf'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> plots/tag_side_{}.png created'.format(version_label)


def plotMuFromHNLTriggeringRate(version_label):
  print '----------------------'
  print 'Will plot triggering rate for HNL mu for ntuples {}'.format(version_label)
  print '----------------------'

  triggering_rate_mu_fromHNL = Quantity('mu_fromHNL_only_satisfies_BParkHLT_cond_tagside', 'Mu from HNL triggers', 'triggering_rate_mu_fromHNL', 2, 0, 2)

  filenames = [f for f in glob.glob('./outputfiles/{}/*'.format(version_label))]

  # get the masses
  masses = []
  for filename in filenames:
    mass = filename[filename.rfind('mass')+4:filename.rfind('_ctau')]
    if mass not in masses:
      masses.append(mass)

  # get the signal points
  signal_points = []
  for mass in masses:
    for filename in filenames:
      if 'mass'+mass not in filename: continue
      ctau = filename[filename.rfind('ctau')+4:filename.rfind('.root')]
      signal_points.append(SignalPoint(mass=mass, ctau=ctau))
      
  # get the tag and probe rates
  triggering_rates = {}
  for mass in masses:
    for signal_point in signal_points:
      if signal_point.mass != mass: continue
      ctau = signal_point.ctau

      for filename in filenames:
        if mass not in filename or ctau not in filename: continue

        f = ROOT.TFile.Open(filename, 'READ')
        tree = f.Get('tree')

        hist_name = 'hist_tag_side_{}_{}'.format(mass, ctau)
        hist = ROOT.TH1D(hist_name, hist_name, triggering_rate_mu_fromHNL.nbins, triggering_rate_mu_fromHNL.bin_min, triggering_rate_mu_fromHNL.bin_max)
        tree.Project(hist_name, triggering_rate_mu_fromHNL.name, '')
        hist.SetDirectory(0)

        triggering_rates['mass{}_ctau{}'.format(mass, ctau)] = [hist.GetMean(), 1-hist.GetMean()]

  # plot the rates
  the_df = pd.DataFrame(triggering_rates, index=[r'$\mu$'+ ' from HNL', r'$\mu$'+ ' from B'])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='Pastel1', figsize=(7, 10))#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  ax.set_title('Triggering rates on hnl side')
  for ictau, tag in enumerate(triggering_rates):
    widths_hnl = the_df[r'$\mu$'+ ' from HNL']
    widths_b = the_df[r'$\mu$'+ ' from B']
    xpos_hnl = 0.1 # widths_hnl[ictau]/2.
    xpos_b = 1 - widths_b[ictau]/2.
    ypos = ictau
    ax.text(xpos_hnl, ypos, round(widths_hnl[ictau], 2), ha='center', va='center',color='black')
    #ax.text(xpos_b, ypos, round(widths_probe[ictau], 2), ha='center', va='center',color='black')
  ax.figure.savefig('plots/triggering_rate_mu_fromHNL_{}.png'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('plots/triggering_rate_mu_fromHNL_{}.pdf'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> plots/triggering_rate_mu_fromHNL_{}.png created'.format(version_label)


def plotNumberTriggeringMuons(version_label):
  print '----------------------'
  print 'Will plot the number of triggering muons in the hnl side for ntuples {}'.format(version_label)
  print '----------------------'

  number_triggering_muons = Quantity('number_triggering_muon_hnlside', 'Number of triggering muons', 'number_triggering_muons', 30, 0, 3)

  filenames = [f for f in glob.glob('./outputfiles/{}/*'.format(version_label))]

  # get the masses
  masses = []
  for filename in filenames:
    mass = filename[filename.rfind('mass')+4:filename.rfind('_ctau')]
    if mass not in masses:
      masses.append(mass)

  # get the signal points
  signal_points = []
  for mass in masses:
    for filename in filenames:
      if 'mass'+mass not in filename: continue
      ctau = filename[filename.rfind('ctau')+4:filename.rfind('.root')]
      signal_points.append(SignalPoint(mass=mass, ctau=ctau))
      
  # get the tag and probe rates
  triggering_rates = {}
  for mass in masses:
    for signal_point in signal_points:
      if signal_point.mass != mass: continue
      ctau = signal_point.ctau

      for filename in filenames:
        if mass not in filename or ctau not in filename: continue

        f = ROOT.TFile.Open(filename, 'READ')
        tree = f.Get('tree')

        hist_name_nb0 = 'hist_tag_side_{}_{}_nb0'.format(mass, ctau)
        hist_nb0 = ROOT.TH1D(hist_name_nb0, hist_name_nb0, number_triggering_muons.nbins, number_triggering_muons.bin_min, number_triggering_muons.bin_max)
        tree.Project(hist_name_nb0, number_triggering_muons.name, 'number_triggering_muon_hnlside<0.5')
        hist_nb0.SetDirectory(0)

        hist_name_nb1 = 'hist_tag_side_{}_{}_nb1'.format(mass, ctau)
        hist_nb1 = ROOT.TH1D(hist_name_nb1, hist_name_nb1, number_triggering_muons.nbins, number_triggering_muons.bin_min, number_triggering_muons.bin_max)
        tree.Project(hist_name_nb1, number_triggering_muons.name, 'number_triggering_muon_hnlside>0.5 && number_triggering_muon_hnlside<1.5')
        hist_nb1.SetDirectory(0)

        hist_name_nb2 = 'hist_tag_side_{}_{}_nb2'.format(mass, ctau)
        hist_nb2 = ROOT.TH1D(hist_name_nb2, hist_name_nb2, number_triggering_muons.nbins, number_triggering_muons.bin_min, number_triggering_muons.bin_max)
        tree.Project(hist_name_nb2, number_triggering_muons.name, 'number_triggering_muon_hnlside>1.5')
        hist_nb2.SetDirectory(0)

        entries_tot = hist_nb0.GetEntries() + hist_nb1.GetEntries() + hist_nb2.GetEntries()
        triggering_rates['mass{}_ctau{}'.format(mass, ctau)] = [hist_nb0.GetEntries()/entries_tot, hist_nb1.GetEntries()/entries_tot, hist_nb2.GetEntries()/entries_tot]

  # plot the rates
  the_df = pd.DataFrame(triggering_rates, index=['0 muon triggering', '1 muon triggering', '2 muons triggering'])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='Pastel1', figsize=(7, 10))#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  ax.set_title('Number of triggering muons on hnl side')
  ax.figure.savefig('plots/number_triggering_muons_{}.png'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('plots/number_triggering_muons_{}.pdf'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> plots/number_triggering_muons_{}.png created'.format(version_label)


def plotBspecies(version_label):
  print '----------------------'
  print 'Will plot B species for ntuples {}'.format(version_label)
  print '----------------------'

  b_pdgid = Quantity('fabs(b_pdgid)', 'B pdgId', 'b_pdgid', 600, 0, 600)

  filenames = [f for f in glob.glob('./outputfiles/{}/*'.format(version_label))]

  # get the masses
  masses = []
  for filename in filenames:
    mass = filename[filename.rfind('mass')+4:filename.rfind('_ctau')]
    if mass not in masses:
      masses.append(mass)

  # get the signal points
  signal_points = []
  for mass in masses:
    for filename in filenames:
      if 'mass'+mass not in filename: continue
      ctau = filename[filename.rfind('ctau')+4:filename.rfind('.root')]
      signal_points.append(SignalPoint(mass=mass, ctau=ctau))
      
  # get the tag and probe rates
  b_species = {}
  for mass in masses:
    for signal_point in signal_points:
      if signal_point.mass != mass: continue
      ctau = signal_point.ctau

      for filename in filenames:
        if mass not in filename or ctau not in filename: continue

        f = ROOT.TFile.Open(filename, 'READ')
        tree = f.Get('tree')

        hist_name = 'hist_b_pdgid_{}_{}'.format(mass, ctau)
        hist = ROOT.TH1D(hist_name, hist_name, b_pdgid.nbins, b_pdgid.bin_min, b_pdgid.bin_max)
        tree.Project(hist_name, 'fabs(b_pdgid)', '')
        #tree.Project(hist_name, 'fabs(b_pdgid)', 'hnl_tag_side==1')
        hist.SetDirectory(0)

        b_isB = hist.GetBinContent(hist.GetXaxis().FindBin(521)) / hist.GetEntries()
        b_isB0 = hist.GetBinContent(hist.GetXaxis().FindBin(511)) / hist.GetEntries()
        b_isBs = hist.GetBinContent(hist.GetXaxis().FindBin(531)) / hist.GetEntries()
        b_species['mass{}_ctau{}'.format(mass, ctau)] = [b_isB, b_isB0, b_isBs]

  # plot the rates
  the_df = pd.DataFrame(b_species, index=['B', 'B0', 'Bs'])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='Set3', figsize=(7, 10))#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  #ax.set_title('Mass 4.5 GeV')
  for ictau, tag in enumerate(b_species):
    widths_B = the_df['B']
    widths_B0 = the_df['B0']
    widths_Bs = the_df['Bs']
    #xpos_tag = 0.2 # widths_tag[ictau]/2.
    #xpos_probe = 1 - widths_probe[ictau]/2.
    xpos_B = 0.25
    xpos_B0 = 0.75
    xpos_Bs = 0.96
    ypos = ictau
    ax.text(xpos_B, ypos, round(widths_B[ictau], 2), ha='center', va='center',color='black')
    ax.text(xpos_B0, ypos, round(widths_B0[ictau], 2), ha='center', va='center',color='black')
    ax.text(xpos_Bs, ypos, round(widths_Bs[ictau], 2), ha='center', va='center',color='black')
    #ax.text(xpos_tag, ypos, round(widths_tag[ictau], 2), ha='center', va='center',color='black')
    #ax.text(xpos_probe, ypos, round(widths_probe[ictau], 2), ha='center', va='center',color='black')
  ax.figure.savefig('plots/b_species_{}.png'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('plots/b_species_{}.pdf'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> plots/b_species_{}.png created'.format(version_label)


def plotAccompagnyingMeson(version_label):
  print '----------------------'
  print 'Will plot X pdgId for ntuples {}'.format(version_label)
  print '----------------------'

  X_pdgid = Quantity('fabs(X_pdgid)', 'X pdgid', 'X_pdgid', 600, 0, 600)

  filenames = [f for f in glob.glob('./outputfiles/{}/*'.format(version_label))]

  # get the masses
  masses = []
  for filename in filenames:
    mass = filename[filename.rfind('mass')+4:filename.rfind('_ctau')]
    if mass not in masses:
      masses.append(mass)

  # get the signal points
  signal_points = []
  for mass in masses:
    for filename in filenames:
      if 'mass'+mass not in filename: continue
      ctau = filename[filename.rfind('ctau')+4:filename.rfind('.root')]
      signal_points.append(SignalPoint(mass=mass, ctau=ctau))
      
  # get the tag and probe rates
  X_pdgids = {}
  for mass in masses:
    for signal_point in signal_points:
      if signal_point.mass != mass: continue
      ctau = signal_point.ctau

      for filename in filenames:
        if mass not in filename or ctau not in filename: continue

        f = ROOT.TFile.Open(filename, 'READ')
        tree = f.Get('tree')

        hist_name = 'hist_X_pdgid_{}_{}'.format(mass, ctau)
        hist = ROOT.TH1D(hist_name, hist_name, X_pdgid.nbins, X_pdgid.bin_min, X_pdgid.bin_max)
        #tree.Project(hist_name, 'fabs(X_pdgid)', '')
        tree.Project(hist_name, 'fabs(X_pdgid)', 'hnl_tag_side==1')
        hist.SetDirectory(0)

        X_isD0 = hist.GetBinContent(hist.GetXaxis().FindBin(421)) / hist.GetEntries()
        X_isD0star = hist.GetBinContent(hist.GetXaxis().FindBin(423)) / hist.GetEntries()
        X_isD = hist.GetBinContent(hist.GetXaxis().FindBin(411)) / hist.GetEntries()
        X_isDstar = hist.GetBinContent(hist.GetXaxis().FindBin(413)) / hist.GetEntries()
        X_isDs = hist.GetBinContent(hist.GetXaxis().FindBin(431)) / hist.GetEntries()
        X_isDsstar = hist.GetBinContent(hist.GetXaxis().FindBin(433)) / hist.GetEntries()
        X_isK = hist.GetBinContent(hist.GetXaxis().FindBin(321)) / hist.GetEntries()
        X_isKstar = hist.GetBinContent(hist.GetXaxis().FindBin(323)) / hist.GetEntries()
        X_ispi = hist.GetBinContent(hist.GetXaxis().FindBin(211)) / hist.GetEntries()
        X_ispi0 = hist.GetBinContent(hist.GetXaxis().FindBin(111)) / hist.GetEntries()
        X_isrho = hist.GetBinContent(hist.GetXaxis().FindBin(213)) / hist.GetEntries()
        X_isrho0 = hist.GetBinContent(hist.GetXaxis().FindBin(113)) / hist.GetEntries()
        X_isgamma = hist.GetBinContent(hist.GetXaxis().FindBin(22)) / hist.GetEntries()
        X_isNothing = hist.GetBinContent(hist.GetXaxis().FindBin(99)) / hist.GetEntries()

        X_pdgids['mass{}_ctau{}'.format(mass, ctau)] = [X_isD0, X_isD0star, X_isD, X_isDstar, X_isDs, X_isDsstar, X_isK, X_isKstar, X_ispi, X_ispi0, X_isrho, X_isrho0, X_isgamma, X_isNothing]
        #print X_isD0 + X_isD0star + X_isD + X_isDstar + X_isDs + X_isDsstar + X_isK + X_isKstar + X_ispi + X_ispi0 + X_isrho + X_isrho0 + X_isgamma + X_isNothing 

  # plot the rates
  the_df = pd.DataFrame(X_pdgids, index=['D0', 'D0star', 'D', 'Dstar', 'Ds', 'Dsstar', 'K', 'Kstar', 'pi', 'pi0', 'rho', 'rho0', 'gamma', '-'])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='PiYG', figsize=(7, 10))#RdBu#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  ax.set_title('Inclusive B-decay (tag side)')
  #for ictau, tag in enumerate(b_species):
  #  widths_B = the_df['B']
  #  widths_B0 = the_df['B0']
  #  widths_Bs = the_df['Bs']
  #  #xpos_tag = 0.2 # widths_tag[ictau]/2.
  #  #xpos_probe = 1 - widths_probe[ictau]/2.
  #  xpos_B = 0.25
  #  xpos_B0 = 0.75
  #  xpos_Bs = 0.96
  #  ypos = ictau
  #  ax.text(xpos_B, ypos, round(widths_B[ictau], 2), ha='center', va='center',color='black')
  #  ax.text(xpos_B0, ypos, round(widths_B0[ictau], 2), ha='center', va='center',color='black')
  #  ax.text(xpos_Bs, ypos, round(widths_Bs[ictau], 2), ha='center', va='center',color='black')
  #  #ax.text(xpos_tag, ypos, round(widths_tag[ictau], 2), ha='center', va='center',color='black')
  #  #ax.text(xpos_probe, ypos, round(widths_probe[ictau], 2), ha='center', va='center',color='black')
  plt.xlim([0, 1])
  ax.figure.savefig('plots/accompagnying_meson_pdgId_{}.png'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('plots/accompagnying_meson_pdgId_{}.pdf'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> plots/accompagnying_meson_pdgId_{}.png created'.format(version_label)


def computeSignalYields(sample, do_old, do_bc):
  '''
  yields = sigma * lumi * filter_eff * acc_eff
  acc_eff comprises reco_eff, tag_eff
  '''

  #filename = './outputfiles/V34_newfilter_genstudy_v3/genTree_mass1p0_ctau100p0.root'
  #filename = './outputfiles/V35_oldfilter_genstudy_v1/genTree_mass1p0_ctau100p0.root'
  #filename = './outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau10p0.root'
  mass = sample.mass
  ctau = sample.ctau
  filter_efficiency = sample.filter_efficiency #.47e-03 #1.50e-03 #5.72E-03 #7.9e-3

  sigma_B = 472.8 * 1e9
  if not do_bc: BR_prod = Decays(mass=mass, mixing_angle_square=1).BR_tot_mu
  else: BR_prod = Decays(mass=mass, mixing_angle_square=1).BR_tot_mu_Bc
  BR_decay = BR_HNLmupion(mass=mass)
  v_square = getVV(mass=mass, ctau=ctau)
  print 'v_square {}'.format(v_square)
  sigma_HNL = sigma_B / 0.4 * BR_prod * BR_decay * v_square

  lumi = 41.6

  preselection_efficiency = 0.6

  # compute acceptance
  f = ROOT.TFile.Open(sample.filename, 'READ')
  tree = f.Get('tree')

  hist_num = ROOT.TH1D('hist_num', 'hist_num', 100, 0, 10)
  hist_deno = ROOT.TH1D('hist_deno', 'hist_deno', 100, 0, 10)
  #hist_num.SetDirectory(0)
  #hist_deno.SetDirectory(0)
        
  #tree.Project('hist_num', 'hnl_mass', '((mu_fromB_satisfies_BParkHLT_cond || mu_fromHNL_satisfies_BParkHLT_cond) && mu_fromHNL_pt>1.5 && abs(mu_fromHNL_eta)<2.5 && pi_fromHNL_pt>0.8 && abs(pi_fromHNL_eta)<2.5) * (weight_reco)')
  #tree.Project('hist_num', 'hnl_mass', '(mu_fromB_pt>9 && abs(mu_fromB_eta)<1.5 && mu_fromHNL_pt>1.5 && abs(mu_fromHNL_eta)<2.5 && pi_fromHNL_pt>0.8 && abs(pi_fromHNL_eta)<2.5) * (weight_reco) * (0.6)')
  #tree.Project('hist_deno', 'hnl_mass', 'mu_fromB_pt>9 && abs(mu_fromB_eta)<1.5')
  if do_old:
    tree.Project('hist_num', 'hnl_mass', '(mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5 && mu_fromHNL_pt>1.5 && abs(mu_fromHNL_eta)<2.5 && pi_fromHNL_pt>0.8 && abs(pi_fromHNL_eta)<2.5) * (weight_reco) * ({})'.format(preselection_efficiency))
    tree.Project('hist_deno', 'hnl_mass', 'mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5')
  else:
    tree.Project('hist_num', 'hnl_mass', '(((mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5) || (mu_fromHNL_pt>7 && abs(mu_fromHNL_eta)<1.5)) && mu_fromHNL_pt>1.5 && abs(mu_fromHNL_eta)<2.5 && pi_fromHNL_pt>0.8 && abs(pi_fromHNL_eta)<2.5) * (weight_reco) * ({})'.format(preselection_efficiency))
    tree.Project('hist_deno', 'hnl_mass', '((mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5) || (mu_fromHNL_pt>7 && abs(mu_fromHNL_eta)<1.5))')
  #tree.Project('hist_num', 'hnl_mass', '(mu_fromB_satisfies_BParkHLT_cond==1 && mu_fromB_pt>9 && abs(mu_fromB_eta)<1.5 && mu_fromHNL_pt>1.5 && abs(mu_fromHNL_eta)<2.5 && pi_fromHNL_pt>0.8 && abs(pi_fromHNL_eta)<2.5) * (weight_reco) * ({})'.format(preselection_efficiency))
  #tree.Project('hist_deno', 'hnl_mass', 'hnl_tag_side==1 && mu_fromB_satisfies_BParkHLT_cond==1 && mu_fromB_pt>9 && abs(mu_fromB_eta)<1.5')
  #tree.Project('hist_num', 'hnl_mass', '(l0_pt>7 && abs(l0_eta)<1.5 && l1_pt>1.5 && abs(l1_eta)<2.5 && pi1_pt>0.8 && abs(pi1_eta)<2.5) * (weight_reco)')
  #tree.Project('hist_num', 'hnl_mass', '(l0_pt>7 && abs(l0_eta)<1.5 && l1_pt>1.5 && abs(l1_eta)<2.5 && pi1_pt>0.8 && abs(pi1_eta)<2.5)')
  #tree.Project('hist_deno', 'hnl_mass', 'l0_pt>7 && abs(l0_eta)<1.5')

  acceptance = hist_num.Integral() / hist_deno.Integral()

  yields = sigma_HNL * lumi * filter_efficiency * acceptance
  print '{:.2e}'.format(yields)

  return v_square, yields


def plotSignalYields():

  ROOT.gROOT.SetBatch(True)

  samples_m1_new = [
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass1p0_ctau1000p0.root', 1.0, 1000.0, 1.61e-3),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass1p0_ctau100p0.root', 1.0, 100.0, 6.20e-3),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass1p0_ctau10p0.root', 1.0, 10.0, 7.8e-3),
  ]

  samples_m3_new = [
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 4.96e-3),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 1.4e-2),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 1.48e-2),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 1.48e-2),
  ]

  samples_m4p5_new = [
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass4p5_ctau100p0.root', 4.5, 100.0, 2.32e-2),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass4p5_ctau10p0.root', 4.5, 10.0, 2.38e-2),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass4p5_ctau1p0.root', 4.5, 1.0, 2.38e-2),
    Sample('./outputfiles/V34_newfilter_genstudy_v3/genTree_mass4p5_ctau0p1.root', 4.5, 0.1, 2.38e-2),
  ]

  samples_m1_old = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau1000p0.root', 1.0, 1000.0, 1.59e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau100p0.root', 1.0, 100.0, 5.17e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau10p0.root', 1.0, 10.0, 5.85e-3),
  ]

  samples_m3_old = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 2.4e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 5.58e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 5.87e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 5.87e-3),
  ]

  samples_m4p5_old = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau100p0.root', 4.5, 100.0, 6.1e-4),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau10p0.root', 4.5, 10.0, 8.0e-4),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau1p0.root', 4.5, 1.0, 8.35e-4),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau0p1.root', 4.5, 0.1, 8.62e-4),
  ]

  canv = ROOT.TCanvas('canv', 'canv', 900, 800)
  canv.SetLogx()
  canv.SetLogy()
  canv.SetGrid()

  graph_m1_new = ROOT.TGraph()
  graph_m3_new = ROOT.TGraph()
  graph_m4p5_new = ROOT.TGraph()
  graph_m1_old = ROOT.TGraph()
  graph_m3_old = ROOT.TGraph()
  graph_m4p5_old = ROOT.TGraph()

  coupling_m1_new = []
  yields_m1_new = []
  for sample in samples_m1_new:
    coupling, signal_yield = computeSignalYields(sample, do_old=False)
    coupling_m1_new.append(coupling)
    yields_m1_new.append(signal_yield)

  for pt in range(0, len(coupling_m1_new)):
    point = graph_m1_new.GetN()
    graph_m1_new.SetPoint(point, coupling_m1_new[pt], yields_m1_new[pt])

  graph_m1_new.SetMarkerStyle(24)
  graph_m1_new.SetMarkerSize(2)
  graph_m1_new.SetMarkerColor(ROOT.kOrange+0)
  graph_m1_new.SetLineStyle(9)
  graph_m1_new.SetLineWidth(2)
  graph_m1_new.SetLineColor(ROOT.kOrange+0)

  coupling_m3_new = []
  yields_m3_new = []
  for sample in samples_m3_new:
    coupling, signal_yield = computeSignalYields(sample, do_old=False)
    coupling_m3_new.append(coupling)
    yields_m3_new.append(signal_yield)

  for pt in range(0, len(coupling_m3_new)):
    point = graph_m3_new.GetN()
    graph_m3_new.SetPoint(point, coupling_m3_new[pt], yields_m3_new[pt])

  graph_m3_new.SetMarkerStyle(24)
  graph_m3_new.SetMarkerSize(2)
  graph_m3_new.SetMarkerColor(ROOT.kRed+1)
  graph_m3_new.SetLineStyle(9)
  graph_m3_new.SetLineWidth(2)
  graph_m3_new.SetLineColor(ROOT.kRed+1)

  coupling_m4p5_new = []
  yields_m4p5_new = []
  for sample in samples_m4p5_new:
    coupling, signal_yield = computeSignalYields(sample, do_old=False)
    coupling_m4p5_new.append(coupling)
    yields_m4p5_new.append(signal_yield)

  for pt in range(0, len(coupling_m4p5_new)):
    point = graph_m4p5_new.GetN()
    graph_m4p5_new.SetPoint(point, coupling_m4p5_new[pt], yields_m4p5_new[pt])

  graph_m4p5_new.SetMarkerStyle(24)
  graph_m4p5_new.SetMarkerSize(2)
  graph_m4p5_new.SetMarkerColor(ROOT.kRed+4)
  graph_m4p5_new.SetLineStyle(9)
  graph_m4p5_new.SetLineWidth(2)
  graph_m4p5_new.SetLineColor(ROOT.kRed+4)

  coupling_m1_old = []
  yields_m1_old = []
  for sample in samples_m1_old:
    coupling, signal_yield = computeSignalYields(sample, do_old=True)
    coupling_m1_old.append(coupling)
    yields_m1_old.append(signal_yield)

  for pt in range(0, len(coupling_m1_old)):
    point = graph_m1_old.GetN()
    graph_m1_old.SetPoint(point, coupling_m1_old[pt], yields_m1_old[pt])

  graph_m1_old.SetMarkerStyle(20)
  graph_m1_old.SetMarkerSize(2)
  graph_m1_old.SetMarkerColor(ROOT.kOrange+0)
  graph_m1_old.SetLineStyle(9)
  graph_m1_old.SetLineWidth(2)
  graph_m1_old.SetLineColor(ROOT.kOrange+0)

  coupling_m3_old = []
  yields_m3_old = []
  for sample in samples_m3_old:
    coupling, signal_yield = computeSignalYields(sample, do_old=True)
    coupling_m3_old.append(coupling)
    yields_m3_old.append(signal_yield)

  for pt in range(0, len(coupling_m3_old)):
    point = graph_m3_old.GetN()
    graph_m3_old.SetPoint(point, coupling_m3_old[pt], yields_m3_old[pt])

  graph_m3_old.SetMarkerStyle(20)
  graph_m3_old.SetMarkerSize(2)
  graph_m3_old.SetMarkerColor(ROOT.kRed+1)
  graph_m3_old.SetLineStyle(9)
  graph_m3_old.SetLineWidth(2)
  graph_m3_old.SetLineColor(ROOT.kRed+1)

  coupling_m4p5_old = []
  yields_m4p5_old = []
  for sample in samples_m4p5_old:
    coupling, signal_yield = computeSignalYields(sample, do_old=True)
    coupling_m4p5_old.append(coupling)
    yields_m4p5_old.append(signal_yield)

  for pt in range(0, len(coupling_m4p5_old)):
    point = graph_m4p5_old.GetN()
    graph_m4p5_old.SetPoint(point, coupling_m4p5_old[pt], yields_m4p5_old[pt])

  graph_m4p5_old.SetMarkerStyle(20)
  graph_m4p5_old.SetMarkerSize(2)
  graph_m4p5_old.SetMarkerColor(ROOT.kRed+4)
  graph_m4p5_old.SetLineStyle(9)
  graph_m4p5_old.SetLineWidth(2)
  graph_m4p5_old.SetLineColor(ROOT.kRed+4)

  graph_dummy = ROOT.TGraph()
  graph_dummy.SetPoint(0, 1e-6, 1e-7)
  graph_dummy.SetPoint(1, 1, 1e10)
  graph_dummy.SetMarkerStyle(0)
  graph_dummy.SetMarkerSize(0)
  graph_dummy.SetMarkerColor(0)
  graph_dummy.GetXaxis().SetTitle('|V^{2}|')
  graph_dummy.GetXaxis().SetLabelSize(0.037)
  graph_dummy.GetXaxis().SetTitleSize(0.042)
  graph_dummy.GetXaxis().SetTitleOffset(1.1)
  graph_dummy.GetYaxis().SetTitle('Yields')
  graph_dummy.GetYaxis().SetLabelSize(0.037)
  graph_dummy.GetYaxis().SetTitleSize(0.042)
  graph_dummy.GetYaxis().SetTitleOffset(1.1)

  graph_dummy.Draw('AP')

  graph_m1_new.Draw('PL same')
  graph_m3_new.Draw('PL same')
  graph_m4p5_new.Draw('PL same')

  graph_m1_old.Draw('PL same')
  graph_m3_old.Draw('PL same')
  graph_m4p5_old.Draw('PL same')

  legend = ROOT.TLegend(0.55, 0.2, 0.8, 0.45)
  legend.SetTextSize(0.03)
  legend.SetLineColor(0)
  legend.SetFillColor(0)
  legend.SetBorderSize(0)
  legend.AddEntry(graph_m1_new, 'm=1GeV, new filter')
  legend.AddEntry(graph_m1_old, 'm=1GeV, old filter')
  legend.AddEntry(graph_m3_new, 'm=3GeV, new filter')
  legend.AddEntry(graph_m3_old, 'm=3GeV, old filter')
  legend.AddEntry(graph_m4p5_new, 'm=4.5GeV, new filter')
  legend.AddEntry(graph_m4p5_old, 'm=4.5GeV, old filter')
  legend.Draw()

  canv.SaveAs('./plots/signal_yields_V34_vs_V33.png')


def plotSignalYieldsBc():

  ROOT.gROOT.SetBatch(True)

  samples_m5p5_new = [
    #Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass5p5_ctau100p0.root', 5.5, 100.0, 1.41e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass5p5_ctau10p0.root', 5.5, 10.0, 1.41e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass5p5_ctau1p0.root', 5.5, 1.0, 1.41e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass5p5_ctau0p1.root', 5.5, 0.1, 1.41e-1),
  ]

  samples_m3_new = [
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 8.64e-2),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 1.88e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 1.93e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 1.93e-1),
  ]

  samples_m4p5_new = [
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass4p5_ctau100p0.root', 4.5, 100.0, 1.30e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass4p5_ctau10p0.root', 4.5, 10.0, 1.32e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass4p5_ctau1p0.root', 4.5, 1.0, 1.31e-1),
    Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass4p5_ctau0p1.root', 4.5, 0.1, 1.31e-1),
  ]

  samples_m5p5_old = [
    #Sample('./outputfiles/V34_newfilter_genstudy_Bc_v1/genTree_mass5p5_ctau100p0.root', 5.5, 100.0, 1.22e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass5p5_ctau10p0.root', 5.5, 10.0, 1.22e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass5p5_ctau1p0.root', 5.5, 1.0, 1.23e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass5p5_ctau0p1.root', 5.5, 0.1, 1.27e-3),
  ]

  samples_m3_old = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 8.52e-2),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 1.51e-1),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 1.54e-1),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 1.54e-1),
  ]

  samples_m4p5_old = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass4p5_ctau100p0.root', 4.5, 100.0, 3.40e-2),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass4p5_ctau10p0.root', 4.5, 10.0, 3.49e-2),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass4p5_ctau1p0.root', 4.5, 1.0, 3.47e-2),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc/genTree_mass4p5_ctau0p1.root', 4.5, 0.1, 3.47e-2),
  ]

  canv = ROOT.TCanvas('canv', 'canv', 900, 800)
  canv.SetLogx()
  canv.SetLogy()
  canv.SetGrid()

  graph_m5p5_new = ROOT.TGraph()
  graph_m3_new = ROOT.TGraph()
  graph_m4p5_new = ROOT.TGraph()
  graph_m5p5_old = ROOT.TGraph()
  graph_m3_old = ROOT.TGraph()
  graph_m4p5_old = ROOT.TGraph()

  coupling_m5p5_new = []
  yields_m5p5_new = []
  for sample in samples_m5p5_new:
    coupling, signal_yield = computeSignalYields(sample, do_old=False, do_bc=True)
    coupling_m5p5_new.append(coupling)
    yields_m5p5_new.append(signal_yield)

  for pt in range(0, len(coupling_m5p5_new)):
    point = graph_m5p5_new.GetN()
    graph_m5p5_new.SetPoint(point, coupling_m5p5_new[pt], yields_m5p5_new[pt])

  graph_m5p5_new.SetMarkerStyle(24)
  graph_m5p5_new.SetMarkerSize(2)
  graph_m5p5_new.SetMarkerColor(ROOT.kBlue+0)
  graph_m5p5_new.SetLineStyle(9)
  graph_m5p5_new.SetLineWidth(2)
  graph_m5p5_new.SetLineColor(ROOT.kBlue+0)

  coupling_m3_new = []
  yields_m3_new = []
  for sample in samples_m3_new:
    coupling, signal_yield = computeSignalYields(sample, do_old=False, do_bc=True)
    coupling_m3_new.append(coupling)
    yields_m3_new.append(signal_yield)

  for pt in range(0, len(coupling_m3_new)):
    point = graph_m3_new.GetN()
    graph_m3_new.SetPoint(point, coupling_m3_new[pt], yields_m3_new[pt])

  graph_m3_new.SetMarkerStyle(24)
  graph_m3_new.SetMarkerSize(2)
  graph_m3_new.SetMarkerColor(ROOT.kRed+1)
  graph_m3_new.SetLineStyle(9)
  graph_m3_new.SetLineWidth(2)
  graph_m3_new.SetLineColor(ROOT.kRed+1)

  coupling_m4p5_new = []
  yields_m4p5_new = []
  for sample in samples_m4p5_new:
    coupling, signal_yield = computeSignalYields(sample, do_old=False, do_bc=True)
    coupling_m4p5_new.append(coupling)
    yields_m4p5_new.append(signal_yield)

  for pt in range(0, len(coupling_m4p5_new)):
    point = graph_m4p5_new.GetN()
    graph_m4p5_new.SetPoint(point, coupling_m4p5_new[pt], yields_m4p5_new[pt])

  graph_m4p5_new.SetMarkerStyle(24)
  graph_m4p5_new.SetMarkerSize(2)
  graph_m4p5_new.SetMarkerColor(ROOT.kRed+4)
  graph_m4p5_new.SetLineStyle(9)
  graph_m4p5_new.SetLineWidth(2)
  graph_m4p5_new.SetLineColor(ROOT.kRed+4)

  coupling_m5p5_old = []
  yields_m5p5_old = []
  for sample in samples_m5p5_old:
    coupling, signal_yield = computeSignalYields(sample, do_old=True, do_bc=True)
    coupling_m5p5_old.append(coupling)
    yields_m5p5_old.append(signal_yield)

  for pt in range(0, len(coupling_m5p5_old)):
    point = graph_m5p5_old.GetN()
    graph_m5p5_old.SetPoint(point, coupling_m5p5_old[pt], yields_m5p5_old[pt])

  graph_m5p5_old.SetMarkerStyle(20)
  graph_m5p5_old.SetMarkerSize(2)
  graph_m5p5_old.SetMarkerColor(ROOT.kBlue)
  graph_m5p5_old.SetLineStyle(9)
  graph_m5p5_old.SetLineWidth(2)
  graph_m5p5_old.SetLineColor(ROOT.kBlue)

  coupling_m3_old = []
  yields_m3_old = []
  for sample in samples_m3_old:
    coupling, signal_yield = computeSignalYields(sample, do_old=True, do_bc=True)
    coupling_m3_old.append(coupling)
    yields_m3_old.append(signal_yield)

  for pt in range(0, len(coupling_m3_old)):
    point = graph_m3_old.GetN()
    graph_m3_old.SetPoint(point, coupling_m3_old[pt], yields_m3_old[pt])

  graph_m3_old.SetMarkerStyle(20)
  graph_m3_old.SetMarkerSize(2)
  graph_m3_old.SetMarkerColor(ROOT.kRed+1)
  graph_m3_old.SetLineStyle(9)
  graph_m3_old.SetLineWidth(2)
  graph_m3_old.SetLineColor(ROOT.kRed+1)

  coupling_m4p5_old = []
  yields_m4p5_old = []
  for sample in samples_m4p5_old:
    coupling, signal_yield = computeSignalYields(sample, do_old=True, do_bc=True)
    coupling_m4p5_old.append(coupling)
    yields_m4p5_old.append(signal_yield)

  for pt in range(0, len(coupling_m4p5_old)):
    point = graph_m4p5_old.GetN()
    graph_m4p5_old.SetPoint(point, coupling_m4p5_old[pt], yields_m4p5_old[pt])

  graph_m4p5_old.SetMarkerStyle(20)
  graph_m4p5_old.SetMarkerSize(2)
  graph_m4p5_old.SetMarkerColor(ROOT.kRed+4)
  graph_m4p5_old.SetLineStyle(9)
  graph_m4p5_old.SetLineWidth(2)
  graph_m4p5_old.SetLineColor(ROOT.kRed+4)

  graph_dummy = ROOT.TGraph()
  graph_dummy.SetPoint(0, 1e-6, 1e-7)
  graph_dummy.SetPoint(1, 1, 1e10)
  graph_dummy.SetMarkerStyle(0)
  graph_dummy.SetMarkerSize(0)
  graph_dummy.SetMarkerColor(0)
  graph_dummy.GetXaxis().SetTitle('|V^{2}|')
  graph_dummy.GetXaxis().SetLabelSize(0.037)
  graph_dummy.GetXaxis().SetTitleSize(0.042)
  graph_dummy.GetXaxis().SetTitleOffset(1.1)
  graph_dummy.GetYaxis().SetTitle('Yields')
  graph_dummy.GetYaxis().SetLabelSize(0.037)
  graph_dummy.GetYaxis().SetTitleSize(0.042)
  graph_dummy.GetYaxis().SetTitleOffset(1.1)

  graph_dummy.Draw('AP')

  graph_m5p5_new.Draw('PL same')
  graph_m3_new.Draw('PL same')
  graph_m4p5_new.Draw('PL same')

  graph_m5p5_old.Draw('PL same')
  graph_m3_old.Draw('PL same')
  graph_m4p5_old.Draw('PL same')

  legend = ROOT.TLegend(0.55, 0.2, 0.8, 0.45)
  legend.SetTextSize(0.03)
  legend.SetLineColor(0)
  legend.SetFillColor(0)
  legend.SetBorderSize(0)
  legend.AddEntry(graph_m3_new, 'm=3GeV, new filter')
  legend.AddEntry(graph_m3_old, 'm=3GeV, old filter')
  legend.AddEntry(graph_m4p5_new, 'm=4.5GeV, new filter')
  legend.AddEntry(graph_m4p5_old, 'm=4.5GeV, old filter')
  legend.AddEntry(graph_m5p5_new, 'm=5.5GeV, new filter')
  legend.AddEntry(graph_m5p5_old, 'm=5.5GeV, old filter')
  legend.Draw()

  canv.SaveAs('./plots/signal_yields_V34_vs_V33_Bc.png')



if __name__ == "__main__":

  version_label = 'V38_request_Bc'
  #version_label = 'V34_newfilter_genstudy_Bc_v1'
  #version_label = 'test_modfilter_v3_n10000000_njt500'

  plotTagRate(version_label=version_label)
  plotMuFromHNLTriggeringRate(version_label=version_label)
  plotNumberTriggeringMuons(version_label=version_label)
  plotBspecies(version_label=version_label)
  plotAccompagnyingMeson(version_label=version_label)

  #plotSignalYields()
  #plotSignalYieldsBc()


