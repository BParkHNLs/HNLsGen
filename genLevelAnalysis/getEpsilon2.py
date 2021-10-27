import ROOT

pl = 'outputfiles/V15_control_BfilterNoMufilter'
fname = 'mass999_ctau999_miniGenTree.root'
f = ROOT.TFile.Open(pl+'/'+fname,'read')
t = f.Get('tree')

b_fiducial_cuts = 'b_pt > 10 && b_pt < 17 && abs(b_y) < 1.45'
gen_filter_cuts = 'l1_pt > 6.8 && abs(l1_eta) < 1.55'

num_cuts = b_fiducial_cuts + ' && ' + gen_filter_cuts
den_cuts = b_fiducial_cuts

hnum = ROOT.TH1F('hnum', 'hnum', 1, 0., 20.)
hden = ROOT.TH1F('hden', 'hden', 1, 0., 20.)

t.Draw('b_pt>>hnum', num_cuts, 'goff')
t.Draw('b_pt>>hden', den_cuts, 'goff')


peff = ROOT.TEfficiency(hnum,hden)

print('===> Resulting efficiency')
print('eff = {:.4f} + {:.4f} - {:.4f}'.format(peff.GetEfficiency(1), peff.GetEfficiencyErrorUp(1), peff.GetEfficiencyErrorLow(1)))
