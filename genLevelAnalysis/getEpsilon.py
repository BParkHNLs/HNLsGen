import ROOT

b_fiducial_cuts = 'b_pt > 10 && b_pt < 17 && abs(b_y) < 1.45'
gen_filter_cuts = 'l1_pt > 6.8 && abs(l1_eta) < 1.55'


### Epsilon 2 
#### from "unfiltered" sample (i.e. no muon filter, but B filter with fiducial)

pl = 'outputfiles/V15_control_BfilterNoMufilter'
fname = 'mass999_ctau999_miniGenTree.root'
f = ROOT.TFile.Open(pl+'/'+fname,'read')
t = f.Get('tree')

num_cuts = b_fiducial_cuts + ' && ' + gen_filter_cuts
den_cuts = b_fiducial_cuts

hnum = ROOT.TH1F('hnum', 'hnum', 1, 0., 20.)
hden = ROOT.TH1F('hden', 'hden', 1, 0., 20.)

t.Draw('b_pt>>hnum', num_cuts, 'goff')
t.Draw('b_pt>>hden', den_cuts, 'goff')


peff = ROOT.TEfficiency(hnum,hden)

print('===> Resulting efficiency Epsilon2')
print('eff = {:.4f} + {:.4f} - {:.4f}'.format(peff.GetEfficiency(1), peff.GetEfficiencyErrorUp(1), peff.GetEfficiencyErrorLow(1)))


### Epsilon 1, 
#### denominator from filtered sample (i.e. with muon filter)
#### numerator from filtered sample, but at RECO level
pl = 'outputfiles/V15_control'
fname = 'mass999_ctau999_miniGenTree.root'
f = ROOT.TFile.Open(pl+'/'+fname,'read')
freco = ROOT.TFile.Open('/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_Oct20_Oct20.root','read')
t = f.Get('tree')
treco = freco.Get('control_tree')

reco_cuts = 'b_cos2d > 0.995 && dimu_mass > (3.097-0.05) && dimu_mass < (3.097+0.05) && b_mass > 5  && sv_lxy > 0.1  && sv_prob > 0. && k_pt > 1.5  && abs(k_eta) < 2.  && l1_pt > 9.  && l2_pt > 4.  && abs(l2_eta) < 2  && hlt_mu9_ip6==1'

num_cuts = b_fiducial_cuts + ' && ' + reco_cuts
den_cuts = b_fiducial_cuts + ' && ' + gen_filter_cuts

hnum = ROOT.TH1F('hnum', 'hnum', 1, 0., 20.)
hden = ROOT.TH1F('hden', 'hden', 1, 0., 20.)

treco.Draw('b_pt>>hnum', num_cuts, 'godff')
t.Draw('b_pt>>hden', den_cuts, 'goff')

peff = ROOT.TEfficiency(hnum,hden)

print('===> Resulting efficiency Epsilon1')
print('eff = {:.4f} + {:.4f} - {:.4f}'.format(peff.GetEfficiency(1), peff.GetEfficiencyErrorUp(1), peff.GetEfficiencyErrorLow(1)))



