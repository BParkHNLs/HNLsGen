dirname="blabla"

# assumes eos mounted on t3home
dirname="V30_nofilter_benchmark"
dirname="BvsBc_comparison"
dirname="benchmark_V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV"
#dirname="benchmark_V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV_Bc"
#dirname="V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV"
rm -rf /eos/home-m/mratti/www/BHNL/genAnalysisNov2021/$dirname
cp -r plots/$dirname /eos/home-m/mratti/www/BHNL/genAnalysisNov2021/$dirname
cp HTACCESS /eos/home-m/mratti/www/BHNL/genAnalysisNov2021/$dirname/.htaccess

