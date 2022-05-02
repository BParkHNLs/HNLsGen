import os
import sys
import ROOT
sys.path.append('../python')
from decays import Decays
from common import getVV, BR_HNLmupion

'''
  This script performs the lifetime reweighting studies at gen-level
'''

class Sample(object):
  def __init__(self, filename, mass, ctau, filter_efficiency):
    self.filename = filename
    self.mass = mass
    self.ctau = ctau
    self.filter_efficiency = filter_efficiency

def computeSignalYields(sample, sample_list='', do_old='', do_bc=''):
  '''
  yields = sigma * lumi * filter_eff * acc_eff
  acc_eff comprises reco_eff, tag_eff
  '''

  if '*' in sample.filename: do_reweighting = True
  else: do_reweighting = False

  samples_for_reweighting = []

  if do_reweighting:
    do_exclusive_reweighting_fromlargerctau = False
    do_exclusive_reweighting_fromsmallerctau = False
    do_inclusive_reweighting = True

    if do_exclusive_reweighting_fromlargerctau:
      # exclusive reweighting from sample having larger ctau
      for the_signal_sample in sample_list:
        if str(sample.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if '*' in the_signal_sample.filename: continue
        if the_signal_sample.ctau >= sample.ctau:
          samples_for_reweighting.append(the_signal_sample)
        elif the_signal_sample.ctau < sample.ctau and len(samples_for_reweighting) == 0: 
          # for mass 1, we reweight to even larger ctaus than the ones generated
          samples_for_reweighting.append(the_signal_sample)

      if len(samples_for_reweighting) > 1:
        # only keep last element of the list
        samples_for_reweighting_tmp = []
        samples_for_reweighting_tmp.append(samples_for_reweighting.pop())
        samples_for_reweighting = samples_for_reweighting_tmp

    if do_exclusive_reweighting_fromsmallerctau:
      # exclusive reweighting from sample having smaller ctau
      for the_signal_sample in sample_list:
        if str(sample.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if '*' in the_signal_sample.filename: continue
        if the_signal_sample.ctau <= sample.ctau:
          samples_for_reweighting.append(the_signal_sample)
        elif the_signal_sample.ctau < sample.ctau and len(samples_for_reweighting) == 0: 
          # for mass 1, we reweight to even larger ctaus than the ones generated
          samples_for_reweighting.append(the_signal_sample)

      if len(samples_for_reweighting) > 1:
        # only keep last element of the list
        samples_for_reweighting_tmp = []
        samples_for_reweighting_tmp.append(samples_for_reweighting[0])
        samples_for_reweighting = samples_for_reweighting_tmp

    if do_inclusive_reweighting:
      for the_signal_sample in sample_list:
        if str(sample.mass).replace('.', 'p') not in the_signal_sample.filename: continue
        if '*' in the_signal_sample.filename: continue
        samples_for_reweighting.append(the_signal_sample)

  print '\n {}'.format(sample.ctau)
  for the_signal_sample in samples_for_reweighting:
    print the_signal_sample.filename

  mass = sample.mass
  ctau = sample.ctau

  if not do_reweighting:
    filter_efficiency = sample.filter_efficiency #TODO
  else:
    filter_efficiency = 0.
    n_events_tot = 0
    for the_signal_sample in samples_for_reweighting:
      the_file = ROOT.TFile.Open(the_signal_sample.filename, 'READ')
      the_tree = the_file.Get('tree')
      n_events = the_tree.GetEntries()
      n_events_tot += n_events
      filter_efficiency += n_events * the_signal_sample.filter_efficiency
    filter_efficiency = filter_efficiency / n_events_tot

  sigma_B = 472.8 * 1e9
  if not do_bc: BR_prod = Decays(mass=mass, mixing_angle_square=1).BR_tot_mu
  else: BR_prod = Decays(mass=mass, mixing_angle_square=1).BR_tot_mu_Bc
  BR_decay = BR_HNLmupion(mass=mass)
  v_square = getVV(mass=mass, ctau=ctau)
  sigma_HNL = sigma_B / 0.4 * BR_prod * BR_decay * v_square

  lumi = 41.6

  preselection_efficiency = 0.6

  # compute acceptance
  if not do_reweighting:
    f = ROOT.TFile.Open(sample.filename, 'READ')
    tree = f.Get('tree')
  else:
    tree = ROOT.TChain('tree')
    for the_sample in samples_for_reweighting:
      tree.Add(the_sample.filename)

  # get lifetime weight
  if not do_reweighting:
    weight_ctau = '1.'
  else:
    deno_weight = ''
    for ifile, the_signal_sample in enumerate(samples_for_reweighting):
      the_file = ROOT.TFile.Open(the_signal_sample.filename, 'READ')
      the_tree = the_file.Get('tree')
      n_events = the_tree.GetEntries()
      if ifile == 0:
        deno_weight += ' {n0} / {ctau0} * exp(-hnl_ct_reco / {ctau0})'.format(
              n0 = n_events,
              ctau0 = the_signal_sample.ctau,
              )
      else:
        deno_weight += ' + {n0} / {ctau0} * exp(-hnl_ct_reco / {ctau0})'.format(
              n0 = n_events,
              ctau0 = the_signal_sample.ctau,
              )
    weight_ctau = '({ntot} / {ctau1} * exp(-hnl_ct_reco / {ctau1}) * (1. / ({deno_weight})))'.format(
        ntot = n_events_tot,
        ctau1 = sample.ctau,
        deno_weight = deno_weight,
        )

  if do_old:
    acceptance_cuts_num = 'mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5 && mu_fromHNL_pt>1.5 && abs(mu_fromHNL_eta)<2.5 && pi_fromHNL_pt>0.8 && abs(pi_fromHNL_eta)<2.5 && Lxyz<1200'
    acceptance_cuts_deno = 'mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5'
  else:
    acceptance_cuts_num = '((mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5) || (mu_fromHNL_pt>7 && abs(mu_fromHNL_eta)<1.5)) && mu_fromHNL_pt>1.5 && abs(mu_fromHNL_eta)<2.5 && pi_fromHNL_pt>0.8 && abs(pi_fromHNL_eta)<2.5 && Lxyz<1200'
    acceptance_cuts_deno = '((mu_fromB_pt>7 && abs(mu_fromB_eta)<1.5) || (mu_fromHNL_pt>7 && abs(mu_fromHNL_eta)<1.5))'

  weight_tot_num = '({}) * (weight_reco) * ({}) * ({})'.format(acceptance_cuts_num, preselection_efficiency, weight_ctau)
  weight_tot_deno = '({})'.format(acceptance_cuts_deno)

  hist_num_name = 'hist_num_{}_{}'.format(sample.mass, sample.ctau)
  hist_deno_name = 'hist_deno_{}_{}'.format(sample.mass, sample.ctau)
  hist_num = ROOT.TH1D(hist_num_name, hist_num_name, 100, 0, 10)
  hist_deno = ROOT.TH1D(hist_deno_name, hist_deno_name, 100, 0, 10)
        
  tree.Project(hist_num_name, 'hnl_mass', weight_tot_num)
  tree.Project(hist_deno_name, 'hnl_mass', weight_tot_deno)

  acceptance = hist_num.Integral() / hist_deno.Integral()

  yields = sigma_HNL * lumi * filter_efficiency * acceptance

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

  samples_m1_old_large = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau1000p0.root', 1.0, 1000.0, 1.59e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau100p0.root', 1.0, 100.0, 5.17e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau10p0.root', 1.0, 10.0, 5.85e-3),
    #Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 10000.0, -99.),
    #Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 7000.0, -99.),
    #Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 5000.0, -99.),
    #Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 4000.0, -99.),
    #Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 3000.0, -99.),
    #Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 2000.0, -99.),
    #Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 1500.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 700.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 500.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 400.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 300.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 200.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 150.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 70.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 50.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 40.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 30.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 20.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass1p0_ctau*.root', 1.0, 15.0, -99.),
  ]

  samples_m3_old_large = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 2.4e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 5.58e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 5.87e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 5.87e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 700.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 500.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 400.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 300.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 200.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 150.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 70.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 50.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 40.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 30.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 20.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 15.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 7.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 5.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 4.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 3.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 2.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass3p0_ctau*.root', 3.0, 1.5, -99.),
  ]

  samples_m4p5_old_large = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau100p0.root', 4.5, 100.0, 6.1e-4),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau10p0.root', 4.5, 10.0, 8.0e-4),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau1p0.root', 4.5, 1.0, 8.35e-4),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau0p1.root', 4.5, 0.1, 8.62e-4),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 70.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 50.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 40.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 30.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 20.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 15.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 7.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 5.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 4.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 3.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 2.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 1.5, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 0.7, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 0.5, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 0.4, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 0.3, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 0.2, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/genTree_mass4p5_ctau*.root', 4.5, 0.15, -99.),
  ]

  samples_m3_old_ctcut_large = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 2.4e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 5.58e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 5.87e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 5.87e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 700.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 500.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 400.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 300.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 200.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 150.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 70.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 50.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 40.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 30.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 20.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 15.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 7.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 5.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 4.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 3.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 2.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_ctsmaller150mm/genTree_mass3p0_ctau*.root', 3.0, 1.5, -99.),
  ]

  samples_m3_old_ctweightedtonano_large = [
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 2.4e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 5.58e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 5.87e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 5.87e-3),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 700.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 500.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 400.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 300.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 200.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 150.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 70.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 50.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 40.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 30.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 20.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 15.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 7.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 5.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 4.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 3.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 2.0, -99.),
    Sample('./outputfiles/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_v2/genTree_mass3p0_ctau*.root', 3.0, 1.5, -99.),
  ]

  samples_m1_V38_request_large = [
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau1000p0.root', 1.0, 1000.0, 1.59e-3),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau100p0.root', 1.0, 100.0, 6.24e-3),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau10p0.root', 1.0, 10.0, 7.75e-3),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 700.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 500.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 400.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 300.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 200.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 150.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 70.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 50.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 40.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 30.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 20.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass1p0_ctau*.root', 1.0, 15.0, -99.),
  ]

  samples_m3_V38_request_large = [
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 4.97e-3),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 1.40e-2),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 1.48e-2),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 1.49e-2),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 700.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 500.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 400.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 300.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 200.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 150.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 70.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 50.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 40.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 30.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 20.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 15.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 7.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 5.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 4.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 3.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 2.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass3p0_ctau*.root', 3.0, 1.5, -99.),
  ]

  samples_m4p5_V38_request_large = [
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau100p0.root', 4.5, 100.0, 2.30e-2),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau10p0.root', 4.5, 10.0, 2.38e-2),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau1p0.root', 4.5, 1.0, 2.38e-2),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau0p1.root', 4.5, 0.1, 2.38e-2),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 70.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 50.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 40.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 30.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 20.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 15.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 7.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 5.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 4.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 3.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 2.0, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 1.5, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 0.7, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 0.5, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 0.4, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 0.3, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 0.2, -99.),
    Sample('./outputfiles/V38_request/genTree_mass4p5_ctau*.root', 4.5, 0.15, -99.),
  ]

  samples_m1_noDisplacementFilter_large = [
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau1000p0.root', 1.0, 1000.0, 7.77e-03),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau100p0.root', 1.0, 100.0, 7.77e-03),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau10p0.root', 1.0, 10.0, 7.81e-03),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 700.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 500.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 400.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 300.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 200.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 150.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 70.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 50.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 40.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 30.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 20.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass1p0_ctau*.root', 1.0, 15.0, -99.),
  ]

  samples_m3_noDisplacementFilter_large = [
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau1000p0.root', 3.0, 1000.0, 1.49e-02),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau100p0.root', 3.0, 100.0, 1.48e-02),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau10p0.root', 3.0, 10.0, 1.48e-02),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau1p0.root', 3.0, 1.0, 1.48e-02),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 700.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 500.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 400.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 300.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 200.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 150.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 70.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 50.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 40.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 30.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 20.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 15.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 7.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 5.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 4.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 3.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 2.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass3p0_ctau*.root', 3.0, 1.5, -99.),
  ]
  
  samples_m4p5_noDisplacementFilter_large = [
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau100p0.root', 4.5, 100.0, 2.39e-2),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau10p0.root', 4.5, 10.0, 2.38e-2),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau1p0.root', 4.5, 1.0, 2.38e-2),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau0p1.root', 4.5, 0.1, 2.38e-2),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 70.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 50.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 40.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 30.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 20.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 15.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 7.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 5.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 4.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 3.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 2.0, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 1.5, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 0.7, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 0.5, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 0.4, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 0.3, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 0.2, -99.),
    Sample('./outputfiles/test_displacementFilter1e9mm/genTree_mass4p5_ctau*.root', 4.5, 0.15, -99.),
  ]

  samples_m1 = samples_m1_old_large
  samples_m3 = samples_m3_old_ctweightedtonano_large
  samples_m4p5 = samples_m4p5_old_large

  # old means samples produced with filter strategy a la July 21
  do_old = True

  canv = ROOT.TCanvas('canv', 'canv', 1200, 800)
  canv.SetLogx()
  canv.SetLogy()
  canv.SetGrid()

  graph_m1 = ROOT.TGraph()
  graph_m3 = ROOT.TGraph()
  graph_m4p5 = ROOT.TGraph()

  coupling_m1 = []
  yields_m1 = []
  for sample in samples_m1:
    coupling, signal_yield = computeSignalYields(sample, sample_list=samples_m1, do_old=do_old, do_bc=False)
    coupling_m1.append(coupling)
    yields_m1.append(signal_yield)

  for pt in range(0, len(coupling_m1)):
    point = graph_m1.GetN()
    graph_m1.SetPoint(point, coupling_m1[pt], yields_m1[pt])

  graph_m1.SetMarkerStyle(20)
  graph_m1.SetMarkerSize(2)
  graph_m1.SetMarkerColor(ROOT.kOrange+0)
  graph_m1.SetLineStyle(9)
  graph_m1.SetLineWidth(0)
  graph_m1.SetLineColor(ROOT.kOrange+0)

  coupling_m3 = []
  yields_m3 = []
  for sample in samples_m3:
    coupling, signal_yield = computeSignalYields(sample, sample_list=samples_m3, do_old=do_old, do_bc=False)
    coupling_m3.append(coupling)
    yields_m3.append(signal_yield)

  for pt in range(0, len(coupling_m3)):
    point = graph_m3.GetN()
    graph_m3.SetPoint(point, coupling_m3[pt], yields_m3[pt])

  graph_m3.SetMarkerStyle(20)
  graph_m3.SetMarkerSize(2)
  graph_m3.SetMarkerColor(ROOT.kRed+1)
  graph_m3.SetLineStyle(9)
  graph_m3.SetLineWidth(0)
  graph_m3.SetLineColor(ROOT.kRed+1)

  coupling_m4p5 = []
  yields_m4p5 = []
  for sample in samples_m4p5:
    coupling, signal_yield = computeSignalYields(sample, sample_list=samples_m4p5, do_old=do_old, do_bc=False)
    coupling_m4p5.append(coupling)
    yields_m4p5.append(signal_yield)

  for pt in range(0, len(coupling_m4p5)):
    point = graph_m4p5.GetN()
    graph_m4p5.SetPoint(point, coupling_m4p5[pt], yields_m4p5[pt])

  graph_m4p5.SetMarkerStyle(20)
  graph_m4p5.SetMarkerSize(2)
  graph_m4p5.SetMarkerColor(ROOT.kRed+4)
  graph_m4p5.SetLineStyle(9)
  graph_m4p5.SetLineWidth(0)
  graph_m4p5.SetLineColor(ROOT.kRed+4)

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

  #graph_m1.Draw('PL same')
  graph_m3.Draw('PL same')
  #graph_m4p5.Draw('PL same')

  legend = ROOT.TLegend(0.55, 0.2, 0.8, 0.45)
  legend.SetTextSize(0.03)
  legend.SetLineColor(0)
  legend.SetFillColor(0)
  legend.SetBorderSize(0)
  #legend.AddEntry(graph_m1, 'm=1GeV')
  legend.AddEntry(graph_m3, 'm=3GeV')
  #legend.AddEntry(graph_m4p5, 'm=4.5GeV')
  legend.Draw()

  canv.SaveAs('./plots/signal_yields_lifetime_reweighting.png')


if __name__ == "__main__":

  plotSignalYields()

