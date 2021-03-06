## Instructions to set up environment for HNL generation

First Installation
```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init

git cms-merge-topic mgratti:BHNL

# *commented lines are not needed*
# git cms-addpkg GeneratorInterface/EvtGenInterface
# git cms-addpkg GeneratorInterface/ExternalDecays
# git cms-addpkg GeneratorInterface/Pythia8Interface

git clone git@github.com:BParkHNLs/HNLsGen.git

export CMSSW_SEARCH_PATH=$CMSSW_BASE/src/HNLsGen/evtGenData/:$CMSSW_SEARCH_PATH  # needed to use local evt_xx.pdl file
scram b

cd HNLsGen
git checkout -b mybranch

export PYTHONPATH=$PYTHONPATH:$PWD 

```

After first installation:
```
cd CMSSW_10_2_15/src/HNLsGen
cmsenv
export PYTHONPATH=$PYTHONPATH:$PWD 
```

## NOT NEEDED ANYMORE - Instructions to set up a different version of Pythia within CMSSW

This needs to be started from clean CMSSW directory, before cmsenv

* On my laptop, cloned sonia's version of pythia
* compiled there (make)
* copied to t3
* then followed https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pythia8Interface#How_to_setup_the_SCRAM_tool_with 
* however does not compile due to imcompatibility between versions 8.230 which is what is used by CMSSSW_10_2_0 and 8.240_sonia 

Therefore adopt different strategy:
* start from 8230 version (as obtained from pythia8 website)
* compile on laptop => but then removed in /lib the file
* copy it somewhere on t3
* copy there from sonia's pythia the files she added / changed
* the follow again above mentioend twiki page...


## Drivers 
```
BPH_start_cfg.py                  => mod tau->3mu  with Fall18, starting point
BPH_mod_cfg.py                    => FROZEN DRIVER, for sharing (version shared with Rome)
step1.py                          => HNL from B species other than Bc
step1_Bc.py                       => HNL from Bcs 
step1_control.py                  => control channel generation                
```

## Production
### vacuum ==> GEN-SIM ==> miniAOD
Submission to slurm:
```
cd slurm
```
Create a VXX_points.py files starting from points.py, then submit
```
python prodHelper.py --help
```
The ```genHelper.py``` is now deprecated.

### Notes for Bc production
We currently shower Bc+ events, starting from LHE files produced by someone else.

The LHE => ROOT conversion is done via ```cmsDrivers/BHNL_Bc_LHEtoRoot_TEMPLATE.py``` and submitted via ```slurm/submitter_lhegen.py```

For showering, instead of using the ```GeneratorFilter```, we use the ```HadronizerFilter```, details in ```cmsDrivers/step1_Bc.py```

## Analyze
To visualize the decay chain in a tree (printout to screen), using ```vector<reco::genParticles>```
```
cd genLevelAnalysis
cmsRun test_ParticleTreeDrawer.py maxEvents=1 inputFiles=file:/work/mratti/GEN_HNL/CMSSW_10_2_15/src/HNLsGen/genSimFiles/BPH-test_HardQCDon.root
```

To get the gen-level ntuples
```
cd genLevelAnalysis
python genTreeProducer.py --pl V05_muFromB_pt5_eta1p6_n600000_njt30 --points V05_points.py
```

To analyze the gen-level ntuples, edit possible options and then run:
```
python genAnalysis.py --pl V05_muFromB_pt5_eta1p6_n600000_njt30
```


***OBOLETE***
Proto-analyzer of ```edm::HepMCProduct```
```
cd genLevelAnalysis
cmsRun test_EvtGenTestAnalyzer.py
```

To analyze particles' distributions
```
cd genLevelAnalysis
python hnl_signal_reweighter.py
```

