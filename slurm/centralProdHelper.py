'''
Production helper for generating private samples for B-initiated HNLs
a version of prodHelper.py adapted to run on fragments rather than cmsdrivers
'''

import sys
import os
import subprocess

from python.common import Point
from python.decays import Decays

class Job(object):
  def __init__(self,opt):

    self.opt = opt
    for k,v in sorted(vars(opt).items()):
      setattr(self,k,v)

    ps = __import__(self.pointFile.split('.py')[0])
    self.points = ps.points
    self.njobs = self.njobs if self.domultijob else 1
    self.nevtsjob = self.nevts if not self.domultijob else self.nevts/self.njobs
    self.prodLabel = '{v}_n{n}_njt{nj}'.format(v=self.ver,n=self.nevts,nj=self.njobs)
    self.nthr = 8 if self.domultithread else 1

    self.user = os.environ["USER"]
#    if self.dobc:
#      self.jop1_in = 'step1_Bc.py'
#    elif self.docontrol:
#      self.jop1_in = 'step1_control.py'
#    else:
#      self.jop1_in = 'step1.py' 
#    self.jop1 = 'step1.py'
#    self.jop2 = 'step2.py'
#    self.jop3 = 'step3.py'
#    self.jop4 = 'step4.py'
   
    # run checks
    if not os.path.isfile(self.pointFile): raise RuntimeError('Provided file for points to scan does not exist, {}'.format(self.pointFile))
    if self.domultijob and self.njobs <= 1: raise RuntimeError('when running multiple jobs, the number of parallel jobs should be larger than 1')
    if self.domultijob and self.nevts % self.njobs != 0: raise RuntimeError('cannot split events in njobs evenly, please change njobs / nevts')
    if self.domultijob and self.domultithread: raise RuntimeError('either multijob or multithread, choose, otherwise seed for generation will repeat')
    if self.dobc and self.nevtsjob > 1000000: raise RuntimeError('Not enough events in the Bc LHE->ROOT files, please reduce number of events per job')
    if self.dobc and self.njobs > 107: raise RuntimeError('Currently we access only 107 M Bc events, either find more Bc events or reduce the total number of jobs')
    # TODO: raise a warning if nevtsjob * filter_eff > npremixfiles * 1200
   
    self.override = (False not in [p.cfg is not None for p in self.points]) # override only if every point has the cfg set
    if self.override:
      print('===> Will override several job configurations for all points:')
      print('       njobs, nevtsjob, time, prodLabel')
      self.prodLabel = '{v}'.format(v=self.ver)


  def makeProdDir(self):
    if not os.path.isdir(self.prodLabel):
      os.system('mkdir -p ./{}'.format(self.prodLabel))
    # otherwise will overwrite 
    if not os.path.isdir(self.prodLabel+'/logs'):
      os.system('mkdir {}/logs'.format(self.prodLabel)) 
    print('===> Created directory for submission {}\n'.format(self.prodLabel))

    print('===> Points to be run')
    for p in self.points:
      p.stamp_simpli()
      if p.cfg is not None: p.cfg.stamp()
    print('')


  def makeTimeStamp(self):
    timestamp=[
        'RUNTIME_step1=$((DATE_END_step1-DATE_START_step1))',
        'echo "Wallclock running time: $RUNTIME_step1 s"'
    ]
    return '\n'.join(timestamp)


  def makeTemplates(self):
    for p in self.points:
      nevtsjob_toset = self.nevtsjob if not self.override else p.cfg.nevtsjob 
      template = [
        '#!/bin/bash',
        '',
        '#SBATCH -J prod_m{m}_ctau{ctau}',
        '#SBATCH -o logs/prod_mass{m}_ctau{ctau}_%a.log', 
        '#SBATCH -e logs/prod_mass{m}_ctau{ctau}_%a.log',
        '#SBATCH -p standard',
        '#SBATCH -t {hh}:00:00',
        '#SBATCH --mem {mem}',
        '#SBATCH --array={arr}',
        '#SBATCH --ntasks=1',
        '#SBATCH --account=t3',
        '',
        'DIRNAME="{pl}"/mass{m}_ctau{ctau}/',
        'STARTDIR=$CMSSW_BASE/src/HNLsGen/', # Will take the cmssw version used at submissio time
        'TOPWORKDIR="/scratch/{user}/"',
        'JOBDIR="gen_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}"', # MIND THE PARENTHESIS
        'WORKDIR=$TOPWORKDIR/$JOBDIR',
        'INSEPREFIX="root://t3dcachedb.psi.ch:1094/"',
        'OUTSEPREFIX="root://t3dcachedb.psi.ch:1094/"',
        'SERESULTDIR="/pnfs/psi.ch/cms/trivcat/store/user/{user}/BHNLsGen/"$DIRNAME'
        '',
        'source $VO_CMS_SW_DIR/cmsset_default.sh',
        'shopt -s expand_aliases',
        'echo ""',
        'echo "Going to set up cms environment"',
        'cd $STARTDIR',
        'cmsenv',
        'echo ""',
        '',
        'echo "Going to create work dir"',
        'mkdir -p $WORKDIR',
        'echo "workdir: "',
        'echo $WORKDIR',
        'echo ""',
        '',
        'echo "Going to create the output dir"',
        'echo "May give an error if the directory already exists, which can be ignored"', # once python bindings to interact with SE are available, should be easier...
        'xrdfs $T3CACHE mkdir -p $SERESULTDIR',
        'echo ""',
        '',
        'echo "Going to copy fragment to config dir"',
        'mkdir -p $CMSSW_BASE/Configuration/GenProduction/python/.'
        'cp $STARTDIR/slurm/{lbldir}/{frag} $CMSSW_BASE/Configuration/GenProduction/python/.',
        'echo ""',
        '',
        'cd $WORKDIR',
        'pwd',
        'echo "Going to run step1"',
        'DATE_START_step1=`date +%s`',
        #'cmsRun {jop1} maxEvents={nevtsjob} nThr={nthr} mass={m} ctau={ctau} outputFile=BPH-step1.root seedOffset=$SLURM_ARRAY_TASK_ID doSkipMuonFilter={dsmf} doDisplFilter={ddf} doMajorana={dmj} doElectron={de} scaleToFilter={stf} minTrackPt={mtpt} minLeptonPt={mlpt} maxDisplacement={mtdd}',
        'cmsDriver.py Configuration/GenProduction/python/{frag} --fileout file:BPH-step1.root --mc --eventcontent FEVTDEBUG --datatier GEN-SIM --conditions 102X_upgrade2018_realistic_v11 --beamspot Realistic25ns13TeVEarly2018Collision --step GEN,SIM --geometry DB:Extended --era Run2_2018 --python_filename step1.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n {nevtsjob} --mc --customise_commands process.RandomNumberGeneratorService.eventSeedOffset="cms.untracked.uint32($SLURM_ARRAY_TASK_ID)" {bcadd}',
        'cmsRun step1.py',
        'DATE_END_step1=`date +%s`',
        'if [ $? -eq 0 ]; then echo "Successfully run step 1"; else exit $?; fi',
        'echo "Finished running step1"',
        'echo "Content of current directory"',
        'ls -al',
        'echo ""',
        '',
        'echo "Going to copy output to result directory"',
        'xrdcp -f $WORKDIR/BPH-step1.root $OUTSEPREFIX/$SERESULTDIR/step1_nj$SLURM_ARRAY_TASK_ID".root"',
        'if [ $? -eq 0 ]; then echo "Successfully copied step1 file"; else exit $?; fi',
        '',
        'echo ""',
        'echo "Cleaning up $WORKDIR"',
        'rm -rf $WORKDIR',
        '{timestamp}',
        'cd $STARTDIR',
      ]
      template = '\n'.join(template)
      template = template.format(
          m=p.mass,
          ctau=p.ctau,
          hh=self.time if not self.override else p.cfg.timejob,
          mem=self.mem,
          lbldir=self.prodLabel,
          arr='1-{}'.format(self.njobs if not self.override else p.cfg.njobs),
          pl=self.prodLabel,
          user=self.user,
          frag=p.fragname,
          nevtsjob=nevtsjob_toset,
          bcadd='--filein root://eosuser.cern.ch//eos/user/m/mratti/BHNL_Bc_LHEtoRoot_step0_nj$SLURM_ARRAY_TASK_ID.root',
          timestamp=self.makeTimeStamp()
          )
      launcherFile = '{pl}/slurm_mass{m}_ctau{ctau}_prod.sh'.format(pl=self.prodLabel,m=p.mass,ctau=p.ctau)
      with open(launcherFile, 'w') as f:
        f.write(template)
    
    print('===> Created templates for batch submission\n')
    
  def writeCfg(self):
    with open('{}/cfg.txt'.format(self.prodLabel), 'w') as f:
      f.write('Run prodHelper.py with following options\n')
      for k,v in sorted(vars(self.opt).items()):
        f.write('{:15s}: {:10s}\n'.format(str(k),str(v)))
    #os.system('cp ../cmsDrivers/{jop_in} ./{pl}/{jop}'.format(jop=self.jop1,jop_in=self.jop1_in,pl=self.prodLabel))
    #if not self.dogenonly:
    #  os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop2,pl=self.prodLabel))
    #  os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop3,pl=self.prodLabel))
    #  os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop4,pl=self.prodLabel))
      

  def submit(self):
    with open('{}/jobs.txt'.format(self.prodLabel), 'w') as f:
      os.chdir(self.prodLabel)
      for p in self.points:
        command = 'sbatch slurm_mass{m}_ctau{ctau}_prod.sh'.format(m=p.mass,ctau=p.ctau)
        out = subprocess.check_output(command, shell=True)
        jobn = out.split('\n')[0].split(' ')[-1]      
        f.write('mass{m}_ctau{ctau} {j}\n'.format(m=p.mass,ctau=p.ctau,j=jobn))
      os.chdir('../')
    print('')
    print('===> Submitted {n} job arrays for {pl}\n'.format(n=len(self.points),pl=self.prodLabel))
  

  def writeFragments(self):
    for i,p in enumerate(self.points):
      fragname = 'BToNMuX_NToEMuPi_SoftQCD_b_mN{:.1f}_ctau{:.1f}mm_TuneCP5_13TeV_pythia8-evtgen_cfi.py'.format(p.mass,p.ctau)
      p.fragname = fragname
      with open('{}/{}'.format(self.prodLabel,p.fragname), 'w') as f:
        tobewritten = '''
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Production Info
configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('B -> mu N X, with long-lived N, m={MASS:.1f}GeV, ctau={CTAU:.1f}mm'),
    name = cms.untracked.string('B -> mu N X, with long-lived N, m={MASS:.1f}GeV, ctau={CTAU:.1f}mm'),
    version = cms.untracked.string('$1.0$')
)

BFilter = cms.EDFilter("MCMultiParticleFilter",
   NumRequired = cms.int32(1),
   AcceptMore = cms.bool(True),
   ParticleID = cms.vint32(521,511,531),
   PtMin = cms.vdouble(0.,0.,0.),
   EtaMax = cms.vdouble(10.,10.,10.),
   Status = cms.vint32(0,0,0), 
)

BToNMuXFilter = cms.EDFilter("PythiaFilterMotherSister", 
    MaxEta = cms.untracked.double(1.55),
    MinEta = cms.untracked.double(-1.55),
    MinPt = cms.untracked.double(6.8), 
    ParticleID = cms.untracked.int32(13),
    MotherIDs = cms.untracked.vint32(521, 511, 531), 
    SisterID = cms.untracked.int32(9900015), 
    MaxSisterDisplacement = cms.untracked.double(1300.), # max Lxy displacement to generate in mm, -1 for no max
    NephewIDs = cms.untracked.vint32(11,13,211), # ids of the nephews you want to check the pt of
    MinNephewPts = cms.untracked.vdouble(0.4,0.4,0.5), 
)

generator = cms.EDFilter("Pythia8GeneratorFilter",
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            convertPythiaCodes = cms.untracked.bool(False),
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            
            list_forced_decays = cms.vstring(       
                'myB+', 
                'myB-',
                'myB0',
                'myB0bar',
                'myB0s',
                'myB0sbar',
            ),
            
            operates_on_particles = cms.vint32(521, -521, 511, -511, 531, -531), 
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass{MASS:.1f}_ctau{CTAU:.1f}_maj.pdl'),
            user_decay_file = cms.vstring('GeneratorInterface/EvtGenInterface/data/HNLdecay_mass{MASS:.1f}_maj_emu.DEC'),
        ),
        parameterSets = cms.vstring('EvtGen130'),
    ),

    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring('SoftQCD:nonDiffractive = on',
                                        'PTFilter:filter = on',      
                                        'PTFilter:quarkToFilter = 5', 
                                        'PTFilter:scaleToFilter = 5.0',
                                        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'processParameters',
                                    ),
    ), 

    comEnergy = cms.double(13000.0),
    filterEfficiency = cms.untracked.double(-1),      # this will not be used by Pythia, only saved in GenInfo 
    maxEventsToPrint = cms.untracked.int32(0),        
    pythiaHepMCVerbosity = cms.untracked.bool(False), 
    pythiaPylistVerbosity = cms.untracked.int32(0),
)

ProductionFilterSequence = cms.Sequence(generator+BFilter+BToNMuXFilter)

'''.format(MASS=p.mass,CTAU=p.ctau)
        f.write(tobewritten)
    print('')
    print('===> Wrote the fragments \n')

  def writeFragmentsBc(self):
    for p in self.points:
      if p.ctau!=0.01:
       fragname = 'BcToNMuX_NToEMuPi_SoftQCD_b_mN{:.1f}_ctau{:.1f}mm_TuneCP5_13TeV_pythia8-evtgen_cfi.py'.format(p.mass,p.ctau)
      else:
       fragname = 'BcToNMuX_NToEMuPi_SoftQCD_b_mN{:.1f}_ctau{:.2f}mm_TuneCP5_13TeV_pythia8-evtgen_cfi.py'.format(p.mass,p.ctau)
      p.fragname = fragname
      with open('{}/{}'.format(self.prodLabel,p.fragname), 'w') as f:
        tobewritten = '''
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Production Info
configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Bc -> mu N X, with long-lived N, m={MASS:.1f}GeV, ctau={CTAU:.1f}mm'),
    name = cms.untracked.string('Bc -> mu N X, with long-lived N, m={MASS:.1f}GeV, ctau={CTAU:.1f}mm'),
    version = cms.untracked.string('$1.0$')
)

BFilter = cms.EDFilter("MCMultiParticleFilter",
   NumRequired = cms.int32(1),
   AcceptMore = cms.bool(True),
   ParticleID = cms.vint32(541),
   PtMin = cms.vdouble(0.),
   EtaMax = cms.vdouble(10.),
   Status = cms.vint32(0), 
)

BToNMuXFilter = cms.EDFilter("PythiaFilterMotherSister", 
    MaxEta = cms.untracked.double(1.55),
    MinEta = cms.untracked.double(-1.55),
    MinPt = cms.untracked.double(6.8), 
    ParticleID = cms.untracked.int32(13),
    MotherIDs = cms.untracked.vint32(541), 
    SisterID = cms.untracked.int32(9900015), 
    MaxSisterDisplacement = cms.untracked.double(1300.), # max Lxy displacement to generate in mm, -1 for no max
    NephewIDs = cms.untracked.vint32(11,13,211), # ids of the nephews you want to check the pt of
    MinNephewPts = cms.untracked.vdouble(0.4,0.4,0.5), 
)

generator = cms.EDFilter("Pythia8HadronizerFilter",
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            convertPythiaCodes = cms.untracked.bool(False),
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            
            list_forced_decays = cms.vstring(       
                'myBc+', 
                'myBc-',
            ),
            
            operates_on_particles = cms.vint32(541, -541), 
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass{MASS:.1f}_ctau{CTAU:.1f}_maj.pdl'),
            user_decay_file = cms.vstring('GeneratorInterface/EvtGenInterface/data/HNLdecay_mass{MASS:.1f}_maj_emu_Bc.DEC'),
        ),
        parameterSets = cms.vstring('EvtGen130'),
    ),

    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring('541:m0 = 6.275',
                                        '541:tau0 = 0.153',
                                        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'processParameters',
                                    ),
    ),

    comEnergy = cms.double(13000.0),
    filterEfficiency = cms.untracked.double(-1),      # this will not be used by Pythia, only saved in GenInfo 
    maxEventsToPrint = cms.untracked.int32(0),        
    pythiaHepMCVerbosity = cms.untracked.bool(False), 
    pythiaPylistVerbosity = cms.untracked.int32(0),
)


ProductionFilterSequence = cms.Sequence(generator+BFilter+BToNMuXFilter)

'''.format(MASS=p.mass,CTAU=p.ctau)
        f.write(tobewritten)
    print('')
    print('===> Wrote the fragments \n')



def getOptions():

  # convention: no capital letters

  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for B-initiated HNL signals', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  parser.add_argument('-n','--nevts', type=int, dest='nevts', help='total number of events to be generated', default=10)
  parser.add_argument('--time', type=str, dest='time', help='allowed time for each job', default='08')
  parser.add_argument('--mem', type=str, dest='mem', help='allowed memory for each job in [MB]', default='4500')
  parser.add_argument('--njobs', type=int, dest='njobs', help='number of parallel jobs to submit', default=10)
  parser.add_argument('--points', type=str, dest='pointFile', help='name of file contaning information on scan to be run', default='points.py')
  parser.add_argument('--domultithread', dest='domultithread', help='run multithreaded', action='store_true', default=False)
  parser.add_argument('--domultijob', dest='domultijob', help='run several separate jobs', action='store_true', default=False)
  parser.add_argument('--dosubmit', dest='dosubmit', help='submit to slurm', action='store_true', default=False)
  #parser.add_argument('--dogenonly', dest='dogenonly', help='produce sample until gen', action='store_true', default=False)
  #parser.add_argument('--doskipmuonfilter', dest='doskipmuonfilter', help='skip the muon filter', action='store_true', default=False)
  #parser.add_argument('--dodisplfilter', dest='dodisplfilter', help='add a filter on the HNL displacement, Lxyz<1.5m', action='store_true', default=False)
  parser.add_argument('--dobc', dest='dobc', help='do the Bc generation instead of other B species', action='store_true', default=False)
  #parser.add_argument('--docontrol', dest='docontrol', help='do the generation for the control channel B->JpsiK', action='store_true', default=False)
  #parser.add_argument('--domajorana', dest='domajorana', help='consider the HNL as a Majorana particle instead of Dirac', action='store_true', default=False)
  #parser.add_argument('--doelectron', dest='doelectron', help='do electron decay in addition to muon', action='store_true', default=False)
  #parser.add_argument('--pythiascale', type=float, dest='pythiascale', help='a parameter in Pythia to scale the pt of the quark (?)', default=5.0)
  #parser.add_argument('--maxdisplacement', type=float, dest='maxdisplacement', help='maximum 2D displacement in mm', default=1300)
  #parser.add_argument('--mintrackpt', type=float, dest='mintrackpt', help='minimum track pt', default=0.0)
  #parser.add_argument('--minleptonpt', type=float, dest='minleptonpt', help='minimum lepton pt', default=0.0)
  #parser.add_argument('--dofragments', dest='dofragments', help='write the relevant fragments', action='store_true', default=False)


  return parser.parse_args()

if __name__ == "__main__":

  opt = getOptions()

  job = Job(opt)

  job.makeProdDir()

  if opt.dobc:
    job.writeFragmentsBc()
  else:
    job.writeFragments()

  #if not opt.docontrol:
  #  job.makeEvtGenData()

  #if opt.dobc:
  #  job.makeEvtGenDecayBc()
  #elif opt.docontrol:
  #  job.makeEvtGenDecayControl()
  #else:
  #  job.makeEvtGenDecay()

  job.makeTemplates()

  job.writeCfg()   

  #if opt.dofragments:
  #  if opt.dobc:
  #    job.writeFragmentsBc()
  #  else:
  #    job.writeFragments()
   

  if opt.dosubmit:
    job.submit()

