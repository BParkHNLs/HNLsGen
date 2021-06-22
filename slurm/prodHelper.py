'''
Production helper for generating private samples for B-initiated HNLs
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
    if self.dobc:
      self.jop1_in = 'step1_Bc.py'
    elif self.docontrol:
      self.jop1_in = 'step1_control.py'
    else:
      self.jop1_in = 'step1.py' 
    self.jop1 = 'step1.py'
    self.jop2 = 'step2.py'
    self.jop3 = 'step3.py'
    self.jop4 = 'step4.py'
   
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


  def makeEvtGenData(self):
    for p in self.points:      
      if self.domajorana:
        hnl_lines = 'add  p Particle  hnl                          9900015  {:.7e}  0.0000000e+00  0.0000000e+00     0     1  {:.7e}    9900015\n'.format(p.mass,p.ctau,p.mass,p.ctau)
      else:
        hnl_lines = 'add  p Particle  hnl                          9900015  {:.7e}  0.0000000e+00  0.0000000e+00     0     1  {:.7e}    9900015\nadd  p Particle  anti_hnl                    -9900015  {:.7e}  0.0000000e+00  0.0000000e+00     0     1  {:.7e}   -9900015\n'.format(p.mass,p.ctau,p.mass,p.ctau)

      with open('../evtGenData/evt_2014.pdl', 'r') as fin:
        contents = fin.readlines()
        contents.insert(4, hnl_lines)
        contents = ''.join(contents)
      with open('../evtGenData/evt_2014_mass{m}_ctau{ctau}_{dm}.pdl'.format(m=p.mass, ctau=p.ctau, dm='maj' if self.domajorana else 'dirac'), 'w') as fout:
        fout.write(contents)
    print('===> Created evtGen particle property files\n')

  def makeEvtGenDecayControl(self):
    for p in self.points: # it should be only one...
      decay_table = [
       'Alias myB+ B+',
       'Alias myB- B-',
       #'Alias myB0 B0',
       #'Alias myB0bar anti-B0',
       'Alias myJpsi J/psi',
       '',
       'ChargeConj myB+ myB-',
       #'ChargeConj myB0bar myB0',
       '',
       'Decay myB+',
       '1.0     myJpsi  K+  SVS;', 
       'Enddecay',
       'CDecay myB-',
       '',
       #'Decay myB0bar',
       #'1.0     myJpsi  anti-K*0  SVV_HELAMP PKHminus PKphHminus PKHzero PKphHzero PKHplus PKphHplus;', 
       #'Enddecay',
       #'CDecay myB0bar',
       #'',
       'Decay myJpsi',
       '1.0    mu-    mu+   PHOTOS  VLL;',
       'Enddecay',
       '',
       'End',      
       '',
      ]
      
      decay_table = '\n'.join(decay_table)

      with open('../evtGenData/ControlChannel.DEC', 'w') as fout:
        fout.write(decay_table)
    print('===> Created evtGen decay files for Control Channel \n')

  def makeEvtGenDecayBc(self):
    for p in self.points:
      decay_table = [
       'Alias myBc+ B_c+',
       'Alias myBc- B_c-',
       '',
       'ChargeConj myBc+ myBc-',
       '{cconj}',
       '',
       'Decay myBc+',
       '{Bc_br0:.10f}               mu+    hnl    PHSP;',
       '',
       'Enddecay',
       'CDecay myBc-',
       '',
       'Decay hnl',
       '0.5     mu-    pi+    PHSP;',
       '{maj_decay}',
       '{el_decay}',
       '{el_maj_decay}',
       'Enddecay',
       '{cdec}',
       '',
       'End',      
       '',      
      ]
      
      decay_table = '\n'.join(decay_table)
      dec = Decays(mass=p.mass, mixing_angle_square=1)

      decay_table = decay_table.format(
                         Bc_br0=dec.Bc_to_uHNL.BR*1000,
                         #pfx='anti_' if not self.domajorana else '',

                         cconj = 'ChargeConj hnl anti_hnl' if not self.domajorana else '',
                         cdec = 'CDecay anti_hnl' if not self.domajorana else '',
                         maj_decay    =  '0.5     mu+    pi-    PHSP;' if self.domajorana else '',
                         el_decay     =  '0.5     e-     pi+    PHSP;' if self.doelectron else '',
                         el_maj_decay =  '0.5     e+     pi-    PHSP;' if self.doelectron and self.domajorana else '',

                         )

      with open('../evtGenData/HNLdecay_mass{m}_{dm}_{de}_Bc.DEC'.format(m=p.mass, 
                                                                         dm='maj' if self.domajorana else 'dirac', 
                                                                         de='emu' if self.doelectron else 'mu'), 'w') as fout:
        fout.write(decay_table)
    print('===> Created evtGen decay files for Bc \n')

  def makeEvtGenDecay(self):

    for p in self.points:
      decay_table = [
       'Alias myB+ B+',
       'Alias myB- B-',
       'Alias myB0 B0',
       'Alias myB0bar anti-B0',
       'Alias myB0s B_s0',
       'Alias myB0sbar anti-B_s0',
       '',
       'ChargeConj myB+ myB-',
       'ChargeConj myB0 myB0bar',
       'ChargeConj myB0s myB0sbar', 
       '{cconj}',
       '',
       'Decay myB+',
       '{Bp_br0:.10f}               mu+    hnl    PHSP;',
       '{Bp_br1:.10f}    anti-D0    mu+    hnl    PHSP;',
       '{Bp_br2:.10f}    anti-D*0   mu+    hnl    PHSP;',
       '{Bp_br3:.10f}    pi0        mu+    hnl    PHSP;',
       '{Bp_br4:.10f}    rho0       mu+    hnl    PHSP;',
       'Enddecay',
       'CDecay myB-',
       '',
       'Decay myB0',
       '{B0_br1:.10f}    D-    mu+    hnl    PHSP;',
       '{B0_br2:.10f}    D*-   mu+    hnl    PHSP;',
       '{B0_br3:.10f}    pi-   mu+    hnl    PHSP;',
       '{B0_br4:.10f}   rho-   mu+    hnl    PHSP;',
       'Enddecay',
       'CDecay myB0bar',
       '',
       'Decay myB0s',
       '{B0s_br1:.10f}    D_s-    mu+    hnl    PHSP;',
       '{B0s_br2:.10f}    D_s*-   mu+    hnl    PHSP;',
       '{B0s_br3:.10f}    K-      mu+    hnl    PHSP;',
       '{B0s_br4:.10f}    K*-     mu+    hnl    PHSP;',
       'Enddecay',
       'CDecay myB0sbar',
       '',
       'Decay hnl',
       '0.5     mu-    pi+    PHSP;',
       '{maj_decay}',
       '{el_decay}',
       '{el_maj_decay}',
       'Enddecay',
       '{cdec}',
       '',
       'End',      
       '',
      ]

      decay_table = '\n'.join(decay_table)
      dec = Decays(mass=p.mass, mixing_angle_square=1)

      decay_table = decay_table.format(
                         Bp_br0=dec.B_to_uHNL.BR*1000,
                         Bp_br1=dec.B_to_D0uHNL.BR*1000,
                         Bp_br2=dec.B_to_D0staruHNL.BR*1000,
                         Bp_br3=dec.B_to_pi0uHNL.BR*1000,
                         Bp_br4=dec.B_to_rho0uHNL.BR*1000,

                         B0_br1=dec.B0_to_DuHNL.BR*1000,
                         B0_br2=dec.B0_to_DstaruHNL.BR*1000,
                         B0_br3=dec.B0_to_piuHNL.BR*1000,
                         B0_br4=dec.B0_to_rhouHNL.BR*1000,

                         B0s_br1=dec.Bs_to_DsuHNL.BR*1000,
                         B0s_br2=dec.Bs_to_DsstaruHNL.BR*1000,
                         B0s_br3=dec.Bs_to_KuHNL.BR*1000,
                         B0s_br4=dec.Bs_to_KstaruHNL.BR*1000,

                         cconj = 'ChargeConj hnl anti_hnl' if not self.domajorana else '',
                         cdec = 'CDecay anti_hnl' if not self.domajorana else '',
                         maj_decay    =  '0.5     mu+    pi-    PHSP;' if self.domajorana else '',
                         el_decay     =  '0.5     e-     pi+    PHSP;' if self.doelectron else '',
                         el_maj_decay =  '0.5     e+     pi-    PHSP;' if self.doelectron and self.domajorana else '',
                         )

      with open('../evtGenData/HNLdecay_mass{m}_{dm}_{de}.DEC'.format(m=p.mass, 
                                                                      dm='maj' if self.domajorana else 'dirac', 
                                                                      de='emu' if self.doelectron else 'mu'), 'w') as fout:
        fout.write(decay_table)
    print('===> Created evtGen decay files\n')


  def appendTemplate(self, jopa, jopb, nthr, nevtsjob):
    if self.dogenonly:
      addlines = ''
    else:
      labela = jopa[0:len(jopa)-3]
      labelb = jopb[0:len(jopb)-3]
      if jopa == 'step2.py':
        command = 'cmsRun {jopa} randomizePremix=1 nPremixFiles={nprx} nThr={nthr} inputFile=BPH-{lblb}.root outputFile=BPH-{lbla}.root seedOffset=$SLURM_ARRAY_TASK_ID'.format(jopa=jopa, nprx=self.npremixfiles, nthr=nthr, lblb=labelb, nevtsjob=nevtsjob, lbla=labela)
      else:
        command = 'cmsRun {jopa} nThr={nthr} inputFile=BPH-{lblb}.root outputFile=BPH-{lbla}.root seedOffset=$SLURM_ARRAY_TASK_ID'.format(jopa=jopa, nthr=1, lblb=labelb, lbla=labela)

      addlines = [
        '',
        '### {lbla} ###',
        'echo "Going to copy cmsdriver to work dir"',
        'cp $STARTDIR/slurm/{lbldir}/{jopa} $WORKDIR/. ',
        'echo "Copying pileup profile file to work dir"',
        'cp $STARTDIR/data/pileup_2018.root $WORKDIR/.',
        'echo "Going to run {lbla}"',
        'DATE_START_{lbla}=`date +%s`',
        '{command}',
        'DATE_END_{lbla}=`date +%s`',
        'if [ $? -eq 0 ]; then echo "Successfully run step"; else exit $?; fi',
        'echo "Finished running {lbla}"',
        'echo "Content of current directory"',
        'ls -al',
        'echo ""',
        'RUNTIME_{lbla}=$((DATE_END_{lbla}-DATE_START_{lbla}))',
        'echo "Intermediate wallclock running time {lbla}: $RUNTIME_{lbla} s"',
        '',
      ]
      if jopa == 'step4.py':
        addlines += [
        '',
        'echo "Going to copy output to result directory"',
        'xrdcp -f $WORKDIR/BPH-{lbla}.root $OUTSEPREFIX/$SERESULTDIR/{lbla}_nj$SLURM_ARRAY_TASK_ID".root"',
        'if [ $? -eq 0 ]; then echo "Successfully copied step4 file"; else exit $?; fi',
        '',
        ]

      addlines = '\n'.join(addlines)
      addlines = addlines.format(jopa=jopa, nthr=nthr, lblb=labelb, nevtsjob=nevtsjob, lbldir=self.prodLabel, lbla=labela, command=command)

    return addlines


  def makeTimeStamp(self):
    if self.dogenonly:
      timestamp=[
        'RUNTIME_step1=$((DATE_END_step1-DATE_START_step1))',
        'echo "Wallclock running time: $RUNTIME_step1 s"'
      ]
    else:
      timestamp=[
        'RUNTIME_step1=$((DATE_END_step1-DATE_START_step1))',
        'RUNTIME_step2=$((DATE_END_step2-DATE_START_step2))',
        'RUNTIME_step3=$((DATE_END_step3-DATE_START_step3))',
        'RUNTIME_step4=$((DATE_END_step4-DATE_START_step4))',
        'RUNTIME_tot=$((DATE_END_step4-DATE_START_step1))',
        'echo "Wallclock running time step1: $RUNTIME_step1 s"',
        'echo "                       step2: $RUNTIME_step2 s"',
        'echo "                       step3: $RUNTIME_step3 s"',
        'echo "                       step4: $RUNTIME_step4 s"',
        'echo "                         tot: $RUNTIME_tot s"'
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
        'echo "Going to copy cmsdriver to work dir"',
        'cp $STARTDIR/slurm/{lbldir}/{jop1} $WORKDIR/. ',
        'echo ""',
        '',
        #'echo "Going to copy the evtGen particle data file"',  # currently not needed, cannot use local copy of file
        #'mkdir -p $WORKDIR/HNLsGen/evtGenData/',
        #'cp ../../evtGenData/evt_2014_mass{m}_ctau{ctau}.pdl $WORKDIR/HNLsGen/evtGenData/.',
        #'echo ""',
        #'',
        'cd $WORKDIR',
        'pwd',
        'echo "Going to run step1"',
        'DATE_START_step1=`date +%s`',
        'cmsRun {jop1} maxEvents={nevtsjob} nThr={nthr} mass={m} ctau={ctau} outputFile=BPH-step1.root seedOffset=$SLURM_ARRAY_TASK_ID doSkipMuonFilter={dsmf} doDisplFilter={ddf} doMajorana={dmj} doElectron={de} scaleToFilter={stf} minTrackPt={mtpt} minLeptonPt={mlpt} maxDisplacement={mtdd}',
        'DATE_END_step1=`date +%s`',
        'if [ $? -eq 0 ]; then echo "Successfully run step 1"; else exit $?; fi',
        'echo "Finished running step1"',
        'mv $WORKDIR/BPH-step1_numEvent{nevtsjob}.root $WORKDIR/BPH-step1.root',
        'echo "Content of current directory"',
        'ls -al',
        'echo ""',
        '',
        'echo "Going to copy output to result directory"',
        'xrdcp -f $WORKDIR/BPH-step1.root $OUTSEPREFIX/$SERESULTDIR/step1_nj$SLURM_ARRAY_TASK_ID".root"',
        'if [ $? -eq 0 ]; then echo "Successfully copied step1 file"; else exit $?; fi',
        '',
        '{addstep2}',
        '{addstep3}',
        '{addstep4}',
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
          jop1=self.jop1,
          dsmf=self.doskipmuonfilter,
          ddf=self.dodisplfilter,
          dmj=self.domajorana,
          de=self.doelectron,
          stf=self.pythiascale,
          mtpt=self.mintrackpt,
          mlpt=self.minleptonpt,
          mtdd=self.maxdisplacement,
          nevtsjob=nevtsjob_toset,
          nthr=self.nthr,
          jop2=self.jop2,
          jop3=self.jop3,
          jop4=self.jop4,
          addstep2=self.appendTemplate(self.jop2,self.jop1,self.nthr,nevtsjob_toset),
          addstep3=self.appendTemplate(self.jop3,self.jop2,self.nthr,nevtsjob_toset),
          addstep4=self.appendTemplate(self.jop4,self.jop3,self.nthr,nevtsjob_toset),

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
    os.system('cp ../cmsDrivers/{jop_in} ./{pl}/{jop}'.format(jop=self.jop1,jop_in=self.jop1_in,pl=self.prodLabel))
    if not self.dogenonly:
      os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop2,pl=self.prodLabel))
      os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop3,pl=self.prodLabel))
      os.system('cp ../cmsDrivers/{jop} ./{pl}/.'.format(jop=self.jop4,pl=self.prodLabel))
      

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
    for p in self.points:
      fragname = 'BToNMuX_NToEMuPi_SoftQCD_b_mN{:.1f}_ctau{:.1f}mm_TuneCP5_13TeV_pythia8-evtgen_cfi.py'.format(p.mass,p.ctau)
      with open('{}/{}'.format(self.prodLabel,fragname), 'w') as f:
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
            particle_property_file = cms.FileInPath('McRequest/evtGenData/evt_2014_mass{MASS:.1f}_ctau{CTAU:.1f}_maj.pdl'),
            user_decay_file = cms.vstring('McRequest/evtGenData/HNLdecay_mass{MASS:.1f}_maj_emu.DEC'),
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
      fragname = 'BcToNMuX_NToEMuPi_SoftQCD_b_mN{:.1f}_ctau{:.1f}mm_TuneCP5_13TeV_pythia8-evtgen_cfi.py'.format(p.mass,p.ctau)
      with open('{}/{}'.format(self.prodLabel,fragname), 'w') as f:
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
            particle_property_file = cms.FileInPath('McRequest/evtGenData/evt_2014_mass{MASS:.1f}_ctau{CTAU:.1f}_maj.pdl'),
            user_decay_file = cms.vstring('McRequest/evtGenData/HNLdecay_mass{MASS:.1f}_maj_emu_Bc.DEC'),
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
  parser.add_argument('--npremixfiles', type=str, dest='npremixfiles', help='number of premixing files to be randomly chosen', default=20)
  parser.add_argument('--domultithread', dest='domultithread', help='run multithreaded', action='store_true', default=False)
  parser.add_argument('--domultijob', dest='domultijob', help='run several separate jobs', action='store_true', default=False)
  parser.add_argument('--dosubmit', dest='dosubmit', help='submit to slurm', action='store_true', default=False)
  parser.add_argument('--dogenonly', dest='dogenonly', help='produce sample until gen', action='store_true', default=False)
  parser.add_argument('--doskipmuonfilter', dest='doskipmuonfilter', help='skip the muon filter', action='store_true', default=False)
  parser.add_argument('--dodisplfilter', dest='dodisplfilter', help='add a filter on the HNL displacement, Lxyz<1.5m', action='store_true', default=False)
  parser.add_argument('--dobc', dest='dobc', help='do the Bc generation instead of other B species', action='store_true', default=False)
  parser.add_argument('--docontrol', dest='docontrol', help='do the generation for the control channel B->JpsiK', action='store_true', default=False)
  parser.add_argument('--domajorana', dest='domajorana', help='consider the HNL as a Majorana particle instead of Dirac', action='store_true', default=False)
  parser.add_argument('--doelectron', dest='doelectron', help='do electron decay in addition to muon', action='store_true', default=False)
  parser.add_argument('--pythiascale', type=float, dest='pythiascale', help='a parameter in Pythia to scale the pt of the quark (?)', default=5.0)
  parser.add_argument('--maxdisplacement', type=float, dest='maxdisplacement', help='maximum 2D displacement in mm', default=1300)
  parser.add_argument('--mintrackpt', type=float, dest='mintrackpt', help='minimum track pt', default=0.0)
  parser.add_argument('--minleptonpt', type=float, dest='minleptonpt', help='minimum lepton pt', default=0.0)
  parser.add_argument('--dofragments', dest='dofragments', help='write the relevant fragments', action='store_true', default=False)


  return parser.parse_args()

if __name__ == "__main__":

  opt = getOptions()

  job = Job(opt)

  job.makeProdDir()

  if not opt.docontrol:
    job.makeEvtGenData()

  if opt.dobc:
    job.makeEvtGenDecayBc()
  elif opt.docontrol:
    job.makeEvtGenDecayControl()
  else:
    job.makeEvtGenDecay()

  job.makeTemplates()

  job.writeCfg()   

  if opt.dofragments:
    if opt.dobc:
      job.writeFragmentsBc()
    else:
      job.writeFragments()
   

  if opt.dosubmit:
    job.submit()

