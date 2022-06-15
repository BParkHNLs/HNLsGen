## Instructions to set up environment for HNL generation

First Installation
```
cmsrel CMSSW_10_2_28_patch1
cd CMSSW_10_2_28_patch1/src
cmsenv
git cms-init

git clone git@github.com:BParkHNLs/HNLsGen.git

export CMSSW_SEARCH_PATH=$CMSSW_BASE/src/HNLsGen/evtGenData/:$CMSSW_SEARCH_PATH  # needed to use local evt_xx.pdl file
scram b

cd HNLsGen
git checkout -b mybranch

export PYTHONPATH=$PYTHONPATH:$PWD 

```

After first installation:
```
cd CMSSW_10_2_28_patch1/src/HNLsGen
source setup.sh
```

## The [./python](./python) directory
Besides the classes for HNL production and decay, it includes several plotting scripts, denoted with prefix “study”. Each script has a short comment at the beginning for documentation purposes.

## The [./cmsDrivers](./cmsDrivers) directory
It contains the cmsdrivers used for steering the production

[step1.py](./cmsDrivers/step1.py)                          => HNL from B species other than Bc

[step1_Bc.py](./cmsDrivers/step1_Bc.py)                       => HNL from Bcs 

[step1_control.py](./cmsDrivers/step1_control.py)                  => control channel generation                

## The [./evtGenData](./evtGenData) directory
It contains files needed by evtgen for the definition of the HNL and its properties, as well as for the decay chain including the HNL.

## Production, the [./slurm](./slurm) directory
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
We currently shower Bc+ events, starting from existing LHE files.

The LHE => ROOT conversion is done via ```cmsDrivers/BHNL_Bc_LHEtoRoot_TEMPLATE.py``` and submitted via ```slurm/submitter_lhegen.py```

For showering, instead of using the ```GeneratorFilter```, we use the ```HadronizerFilter```, details in ```cmsDrivers/step1_Bc.py```

## Analyze, the [./genLevelAnalysis](./genLevelAnalysis) directory
To visualize the decay chain in a tree (printout to screen), using ```vector<reco::genParticles>```
```
cd genLevelAnalysis
cmsRun test_ParticleTreeDrawer.py maxEvents=1 inputFiles=file:/work/mratti/GEN_HNL/CMSSW_10_2_28_patch1/src/HNLsGen/genSimFiles/BPH-test_HardQCDon.root
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

