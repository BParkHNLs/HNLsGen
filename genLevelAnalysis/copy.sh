dirname="blabla"

# assumes eos mounted on t3home
#dirname="testManyChan_wMuonFilter_n2000000_njt100_incl_muTrigPt9"
#dirname="testManyChan_wMuonFilter_n2000000_njt100_incl_muTrigPt9_noLxyCut"
#dirname="testManyChan_wMuonFilter_n2000000_njt100_incl_muTrigPt9_Lxy1000"
#dirname="testManyChan_wMuonFilter_n2000000_njt100_incl_muTrigPt9_Lxy500"
#dirname="testManyChan_wMuonFilter_n2000000_njt100_incl_muTrigPt9_Lxy300"
dirname="compare3_pt5_eta1p6"
dirname="testManyChan_wMuonFilter_n20000000_njt200_VS_testBc_n250000_njt100_incl_muTrigPt9"
dirname="testManyChan_wMuonFilter_n20000000_njt200_incl_muTrigPt9"
dirname="testManyChan_wMuonFilter_n20000000_njt200_incl_muTrigPt9_Lxy1000_clipw100"
dirname="V11_inclB_n6700000_njt200_VS_V11_inclB_n25000000_njt200_VS_V11_inclB_n4200000_njt200_incl_muTrigPt9"
dirname="testManyChan_wMuonFilter_n20000000_njt200_VS_testManyChan_wMuonFilter_SoftQCD_n17000000_njt250_incl_muTrigPt9"
dirname="testManyChan_wMuonFilter_n20000000_njt200_incl_muTrigPt9_Lxy1000"
dirname="testNoDispl_n5000_njt20_VS_testDispl_n10000_njt20_incl_muTrigPt9"
dirname="testDirac_n5000_njt20_VS_testMajorana_n5000_njt20_incl_muTrigPt9"
dirname="testDirac_n5000_njt20_VS_testDirac_n5000_njt20_incl_muTrigPt9"
dirname="testDirac_n5000_njt20_VS_testMajorana_n5000_njt20_incl_muTrigPt9"
dirname="pilot_masses_noMuon_noDispl_V13_incl_muTrigPt9"
dirname="pilotV15_incl_muTrigPt9"
dirname="pilotV15_Bc_incl_muTrigPt9"
dirname="pilotV15_incl_muTrigPt9"
dirname="V20_emu_incl_muTrigPt9"
dirname="pilotV15_incl_muTrigPt9"
dirname="V15_full_incl_muTrigPt9"
dirname="V23_scaleToFilter_1_VS_V23_scaleToFilter_2_VS_V23_scaleToFilter_4_incl_muTrigPt9"
dirname="V23_scaleToFilter_1_VS_V23_scaleToFilter_5_VS_V23_scaleToFilter_10_incl_muTrigPt9"
dirname="V23_scaleToFilter_1_VS_V23_scaleToFilter_10_incl_muTrigPt9"
dirname="V23_scaleToFilter_1_VS_V23_scaleToFilter_5_incl_muTrigPt9"
dirname="V15_full_VS_V23_scaleToFilter_5_incl_muTrigPt9"
dirname="V30_nofilter_incl_muTrigPt9"
dirname="V31_stats_incl"
dirname="V31_stats_Bxsec_recoW"
dirname="V31_stats_NBs_recoW"
#dirname="V31_stats_NBs_NOrecoW"
dirname="V30_nofilter_benchmark_incl"
#dirname="V15_full_VS_V23_scaleToFilter_10_incl_muTrigPt9"
#dirname="V12_multiB_n4200000_njt200_incl_muTrigPt9"
dirname="pilot_V32_stats_Lxy1300_tkPt400MeV_lepPt500MeV_incl"
dirname="V32_stats_Lxy1300_incl"
rm -rf /eos/home-m/mratti/www/BHNL/genAnalysis/$dirname
cp -r plots/$dirname /eos/home-m/mratti/www/BHNL/genAnalysis/$dirname
cp HTACCESS /eos/home-m/mratti/www/BHNL/genAnalysis/$dirname/.htaccess
