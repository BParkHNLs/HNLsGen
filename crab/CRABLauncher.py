import os
import sys
from os import path
import random

#sys.path.append('../python')
from python.common import Point


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='CRAB Launcher', add_help=True)

  parser.add_argument('--pl'      , type=str, dest='pl'         , help='production label'                 ,                      default='V00_v00')
  parser.add_argument('--points'  , type=str, dest='points_file', help='point file'                       ,                      default='points.py')
  parser.add_argument('--nevents' , type=str, dest='nevents'    , help='requested events (analysis level)',                      default='20000')
  parser.add_argument('--dosubmit',           dest='dosubmit'   , help='submit to slurm'                  , action='store_true', default=False)
  parser.add_argument('--dobc'    ,           dest='dobc'       , help='Bc grid'                          , action='store_true', default=False)

  return parser.parse_args()


class CRABLauncher(object):
  def __init__(self, options):
    self.opt = options
    for k,v in sorted(vars(opt).items()):
      setattr(self,k,v)

    ps = __import__(self.points_file.split('.py')[0])
    self.points = ps.points

    # some fixed parameters
    self.nevents_perminiaod = 1500
    self.eff_nanoaod = 0.1


  def createOuputDir(self):
    for point in self.points:
      dirname = '{}/mass{mass:.2f}_ctau{ctau:.1f}'.format(self.pl, mass=point.mass, ctau=point.ctau)
      if not path.exists(dirname):
        os.system('mkdir -p {}'.format(dirname))


  def getRandomLHEfiles(self, prefix=''):
    random_lhe_files = []
    # one lhe file contains 200k events, we input 5 of them for each point
    for i in range(0, 5):
      lhe_file_list = open('../data/lhe_files/full_list_lhe_files.txt')
      lhe_files = lhe_file_list.readlines()
      lhe_file_idx = random.randint(1,17127) # there are 17127 lhe files in total
      for ifile, lhe_file in enumerate(lhe_files):
        if ifile != lhe_file_idx: continue
        random_lhe_file = lhe_file[:lhe_file.rfind('.xz')+3]
      random_lhe_files.append(random_lhe_file)
    for ifile, random_lhe_file in enumerate(random_lhe_files):
      if ifile == 0:
        random_lhe_files_str = prefix + random_lhe_file
      else:
        random_lhe_files_str += ',' + prefix + random_lhe_file
    return random_lhe_files_str


  def createCRABConfig(self):
    for point in self.points:
      self.nevents_togenerate = float(self.nevents) / (float(point.myfilterEff) * float(self.eff_nanoaod))
      self.nevents_perjob = int(self.nevents_perminiaod / float(point.myfilterEff))
      if self.nevents_perjob > 75000:
        raise RuntimeError('WARNING - the number of events per job ({}) is larger than 75k. Number of events per LHE file is 200k. Set n_miniaod to a lower number')
      if self.nevents_togenerate > self.nevents_perjob:
        self.njobs = int(self.nevents_togenerate / self.nevents_perjob)
      else:
        self.nevents_perjob = int(self.nevents_togenerate)
        self.njobs = 1

      config = [
        'from CRABClient.UserUtilities import config',
        'config = config()',
        '',
        'import datetime, time',
        'ts = time.time()',
        'st = datetime.datetime.fromtimestamp(ts).strftime("%Y%m%d_%H%M%S")',
        'config.General.requestName = "{pl}_mass{MASS}_ctau{CTAU}_"+ st',
        'config.General.transferOutputs = True',
        'config.General.transferLogs = True',
        'config.General.workArea = "crab_workdir"',
        '',
        'config.JobType.pluginName = "PrivateMC"',
        '{addgnrt}',
        'config.JobType.psetName = "step1.py"',
        'config.JobType.inputFiles = ["../../data/FrameworkJobReport.xml", "../../pdl_files/evt_BHNL_mass{mass:.2f}_ctau{ctau:.1f}_maj.pdl", "step1.py", "../../cmsDrivers/step2.py", "../../data/pileup_2018.root", "../../cmsDrivers/step3.py", "../../cmsDrivers/step4.py"]',
        'config.JobType.outputFiles = ["step4.root"]',
        'config.JobType.scriptExe = "submitter.sh"',
        'config.JobType.disableAutomaticOutputCollection = True',
        'ncores = 1',
        'mempecore = 3000',
        'config.JobType.maxMemoryMB  = ncores*mempecore',                                                                                                                                  
        'config.JobType.maxJobRuntimeMin = {time}',
        'config.JobType.allowUndistributedCMSSW = True',
        #'config.JobType.priority = 100',
        '',
        'config.Data.outputPrimaryDataset = "mass{MASS}_ctau{CTAU}"',
        'config.Data.outLFNDirBase = "/store/user/anlyon/BHNLsGen/{pl}"',
        'NUMEVENTSPERJOB={nevtsjob}',
        'NUMJOBS={njobs}',
        'config.Data.totalUnits = NUMEVENTSPERJOB*NUMJOBS',
        'config.Data.publication = False',
        'config.Data.unitsPerJob = NUMEVENTSPERJOB',
        'config.Data.splitting = "EventBased"',
        'config.Data.ignoreLocality = False',
        '',
        'config.Site.storageSite = "T3_CH_PSI"',
        'config.Site.whitelist = ["T2_CH_CSCS"]',
        #'config.Site.whitelist = ["T1_US_FNAL", "T2_US_Caltech", "T2_CH_CSCS", "T2_US_MIT"]',
        ]

      config = '\n'.join(config)
      config = config.format(
          MASS = str(point.mass).replace('.', 'p'),
          CTAU = str(point.ctau).replace('.', 'p'),
          mass = point.mass,
          ctau = point.ctau,
          addgnrt = 'config.JobType.generator = "lhe"' if self.dobc else '',
          pl = self.pl,
          time = 1800 if float(point.mass) < 3 else 3000,
          nevtsjob = self.nevents_perjob,
          njobs = self.njobs,
          )

      config_filename = '{pl}/mass{mass:.2f}_ctau{ctau:.1f}/crab_config.py'.format(pl=self.pl, mass=point.mass, ctau=point.ctau)
      config_file = open(config_filename, 'w+')
      config_file.write(config)
      config_file.close()


  def createSubmitter(self):
    for point in self.points:
      submitter = [
        '#!/bin/bash',
        'echo "will copy pdl file to /srv/CMSSW_10_2_28_patch1/src"',
        'cp ./CMSSW_10_2_28_patch1/evt_BHNL_mass{mass:.2f}_ctau{ctau:.1f}_maj.pdl /srv/CMSSW_10_2_28_patch1/src',
        'echo "end copying"',
        'echo " "',
        'echo "content of /srv/CMSSW_10_2_28_patch1/src"',
        'ls -al /srv/CMSSW_10_2_28_patch1/src',
        'echo " "',
        'echo "will run step1"',
        'cmsRun -j step1.log step1.py',
        'echo "end run step1"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "will run step2"',
        'cmsRun -j step2.log step2.py',
        'echo "end run step2"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "will run step3"',
        'cmsRun -j step3.log step3.py',
        'echo "end run step3"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "will run step4"',
        'cmsRun -e -j FrameworkJobReport.xml step4.py',
        'echo "end run step4"',
        'echo "Done"',
        ]

      submitter = '\n'.join(submitter)
      submitter = submitter.format(
          mass = point.mass,
          ctau = point.ctau,
          )

      submitter_filename = '{pl}/mass{mass:.2f}_ctau{ctau:.1f}/submitter.sh'.format(pl=self.pl, mass=point.mass, ctau=point.ctau)
      submitter_file = open(submitter_filename, 'w+')
      submitter_file.write(submitter)
      submitter_file.close()


  def createDriver(self):
    for ipoint, point in enumerate(self.points):
      if not self.dobc:
        fragment_name = 'BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL{mass:.2f}_ctau{ctau:.1f}mm_TuneCP5_13TeV_pythia8-evtgen_cfi.py'.format(
            mass = point.mass,
            ctau = point.ctau,
            )
      else:
        fragment_name = 'BcToNMuX_NToEMuPi_SoftQCD_b_mN{mass:.2f}_ctau{ctau:.1f}mm_TuneCP5_13TeV_pythia8-evtgen_cfi.py'.format(
            mass = point.mass,
            ctau = point.ctau,
            )

      if not self.dobc:
        command = 'cmsDriver.py Configuration/GenProduction/python/fragment.py --fileout file:step1.root --mc --eventcontent FEVTDEBUG --datatier GEN-SIM --conditions 102X_upgrade2018_realistic_v11 --beamspot Realistic25ns13TeVEarly2018Collision --step GEN,SIM --geometry DB:Extended --era Run2_2018 --python_filename step1.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n {nevts} --mc --customise_commands "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate()"'

        command = command.format(
              nevts = self.nevents_perjob,
              )
        command_replace = ''
      else:
        command = 'cmsDriver.py Configuration/GenProduction/python/fragment.py --fileout file:step1.root --mc --eventcontent FEVTDEBUG --datatier GEN-SIM --conditions 102X_upgrade2018_realistic_v11 --beamspot Realistic25ns13TeVEarly2018Collision --step GEN,SIM --geometry DB:Extended --era Run2_2018 --python_filename step1.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n {nevts} --mc --customise_commands "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate()" --filein {lhefile}'
        
        command = command.format(
              nevts = self.nevents_perjob,
              lhefile = self.getRandomLHEfiles(prefix='gsiftp://storage01.lcg.cscs.ch/'),
              )
        command_replace = 'sed -i "s/PoolSource/LHESource/g" step1.py'

      if not path.exists('{pl}/mass{mass:.2f}_ctau{ctau:.1f}/step1.py'.format(pl=self.pl, mass=point.mass, ctau=point.ctau)):
        submitter_tmp = [
          '#!/bin/bash',
          'STARTDIR=$PWD',
          'PRODDIR=$PWD/{pl}',
          'WORKDIR=$CMSSW_BASE/src',
          'cd $WORKDIR',
          'mkdir -p Configuration/GenProduction/python/',
          'cp $STARTDIR/fragments/{frgmt} Configuration/GenProduction/python/fragment.py',
          'scram b -j 8',
          '{cmd}',
          '{cmd_rpl}',
          'cp step1.py $PRODDIR/mass{mass:.2f}_ctau{ctau:.1f}/step1.py',
          'cd $STARTDIR',
          ]

        submitter_tmp = '\n'.join(submitter_tmp)
        submitter_tmp = submitter_tmp.format(
            pl = self.pl,
            frgmt = fragment_name,
            cmd = command,
            cmd_rpl = command_replace,
            mass = point.mass,
            ctau = point.ctau,
            )

        submitter_tmp_file = open('submitter_tmp.sh', 'w+')
        submitter_tmp_file.write(submitter_tmp)
        submitter_tmp_file.close()
        
        os.system('sh submitter_tmp.sh > log.txt')

        percentage = float(ipoint+1) / float(len(self.points)) * 100
        print '   ---> {}% completed'.format(round(percentage, 1))

        os.system('rm submitter_tmp.sh') 
        os.system('rm log.txt') 


  def submit(self):
    for ipoint, point in enumerate(self.points):
      dirname = '{}/mass{mass:.2f}_ctau{ctau:.1f}'.format(self.pl, mass=point.mass, ctau=point.ctau)
      os.chdir(dirname)
      command_submit = 'crab submit -c crab_config.py > job_submit.log'
      os.system(command_submit)
      percentage = float(ipoint+1) / float(len(self.points)) * 100
      print '   ---> {}% completed'.format(round(percentage, 1))
      os.chdir('../../')


  def process(self):
    print '------------------------------------------'
    print '              CRAB Launcher               '
    print '------------------------------------------'
    print ' '
    print 'Will submit production with:'
    print '           - points:     {}'.format(self.points_file) 
    print '           - prod label: {}'.format(self.pl) 
    print '------------------------------------------'
    print ' ' 

    print ' -> Creating output directories'
    self.createOuputDir()

    print '\n -> Creating CRAB configuration files'
    self.createCRABConfig()
      
    print '\n -> Creating submitter files'
    self.createSubmitter()

    print '\n -> Creating cmsDrivers'
    self.createDriver()

    if self.dosubmit:
      print '\n --> Submitting...'
      self.submit()

    print '\nDone'




if __name__ == "__main__":
  
  opt = getOptions()

  launcher = CRABLauncher(options=opt)
  launcher.process()



