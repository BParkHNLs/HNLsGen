'''
Job option for the B-initiated HNL generation - specialised for the Bc case
'''

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# define the defaults here, changed from command line
options.maxEvents = -1 # -1 means all events, maxEvents considers the total over files considered
options.outputFile = 'BPH-test.root'
# add costum parameters
options.register ('severityLevel',
                  'ERROR', # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.string,          # string, int, or float
                  'severity level for log messages, DEBUG, INFO, WARNING, ERROR')
options.register('nThr',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Number of threads')
options.register('seedOffset',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Seed offset')
options.register('mass',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 'mass of the HNL')
options.register('ctau',
                 100,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 'ctau of the HNL [mm]')
options.register('doSkipMuonFilter',    
                  False, 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.bool, 
                  'Skip the muon filter' )
options.register('doDisplFilter',    
                  False, 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.bool, 
                  'In muon filter, add a cut on the HNL displacement' )
options.register('doMajorana',    
                  False, 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.bool, 
                  'HNL is majorana particle, otherwise dirac' )
options.register('doElectron',    
                  False, 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.bool, 
                  'HNL can decay to epi too, otherwise mupi only' )
#options.register ("doDirac",
#                  1, # default value
#                  VarParsing.multiplicity.singleton, # singleton or list
#                  VarParsing.varType.int,          # string, int, or float
#                  "do Dirac HNL? otherwise Majorana")
options.parseArguments()
print options

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2018Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
### LHE->root files, each file has 1 M events, _0.root does not work, start from 1
LHErootfile = 'root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/BHNL_Bc_LHEGEN_v0/BHNL_Bc_LHEtoRoot_step0_nj{ijob}.root'.format(ijob=options.seedOffset)
process.source = cms.Source("PoolSource",
  dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
  fileNames = cms.untracked.vstring(LHErootfile),
  #skipEvents=cms.untracked.uint32(2), ## not needed
  inputCommands = cms.untracked.vstring(
      'keep *', 
      'drop LHEXMLStringProduct_*_*_*'
  ),
  secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('test'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('file:{}'.format(options.outputFile)),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Message logger
#process.MessageLogger = cms.Service("MessageLogger",
#    cout = cms.untracked.PSet(
#         threshold  = cms.untracked.string(options.severityLevel) 
#    ),
#    destinations = cms.untracked.vstring('cout')
#)

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v11', '')

process.MuFilter = cms.EDFilter("MCParticlePairFilter",
    MaxEta = cms.untracked.vdouble(2.45, 2.45),
    MinEta = cms.untracked.vdouble(-2.45, -2.45),
    MinPt = cms.untracked.vdouble(2.7, 2.7),
    ParticleID1 = cms.untracked.vint32(13),
    ParticleID2 = cms.untracked.vint32(13)
)


### Operates on all particles in the HepMC::GenEvent
### accpects events if:
###  - there is at least one particle with specified pdgID in the entire HepMC::GenEvent
###  - any status (but can be specified)
process.BpFilter = cms.EDFilter("PythiaFilter",
    ParticleID = cms.untracked.int32(541) # Bc+ Bc- filter 
)

process.BFilter = cms.EDFilter("MCMultiParticleFilter",
   NumRequired = cms.int32(1),
   AcceptMore = cms.bool(True),
   #ParticleID = cms.vint32(521,511,531), # abs not needed
   ParticleID = cms.vint32(541), # abs already taken into account, Bc+, Bc-
   PtMin = cms.vdouble(0.),
   EtaMax = cms.vdouble(10.),
   Status = cms.vint32(0), 
)


if options.doDisplFilter:
  maxDispl = cms.untracked.double(1500) 
else:
  maxDispl = cms.untracked.double(-1)

process.SingleMuFilter = cms.EDFilter("PythiaFilterMotherSister", 
    #MaxEta = cms.untracked.double(6),
    #MinEta = cms.untracked.double(-6),
    #MinPt = cms.untracked.double(0.0), # <=== keep it a bit lower than the pt cut at reco level... 
    MaxEta = cms.untracked.double(1.55),
    MinEta = cms.untracked.double(-1.55),
    MinPt = cms.untracked.double(6.8), # <=== keep it a bit lower than the pt cut at reco level... #### FIXME should be raised to 6.5 - 7
    ParticleID = cms.untracked.int32(13), # abs value is taken
    #Status = cms.untracked.int32(1),
    MotherIDs = cms.untracked.vint32(541), # require muon to come from Bc+/Bc- decay
    SisterID = cms.untracked.int32(9900015), # require HNL sister
    MaxSisterDisplacement = maxDispl, # max Lxyz displacement to generate in mm, -1 for no max
)

process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            ### for info, see https://twiki.cern.ch/twiki/bin/view/CMS/EvtGenInterface
            convertPythiaCodes = cms.untracked.bool(False),
            
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            
            ### the list of particles that are aliased and forced to be decayed by EvtGen
            list_forced_decays = cms.vstring(       
                'myBc+', 
                'myBc-',
            ),
            
            ### the list of particles that remain undecayed by Pythia for EvtGen to operate on. 
            ### If the vector has a size 0 or size of 1 with a value of 0, the default list is used. 
            ### These are are hard-coded in: GeneratorInterface/EvtGenInterface/plugins/EvtGen/EvtGenInterface.cc., in the function SetDefault_m_PDGs().            
            operates_on_particles = cms.vint32(541, -541),  # 541 is Bc+

            ### The file with properties of all particles
            particle_property_file = cms.FileInPath('HNLsGen/evtGenData/evt_2014_mass{m}_ctau{ctau}_{dm}.pdl'.format(\
                                                       m=options.mass,ctau=options.ctau,dm='maj' if options.doMajorana else 'dirac')), 
            #https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEdmFileInPath 

            ### The decay file 
            user_decay_file = cms.vstring('HNLsGen/evtGenData/HNLdecay_mass{m}_{dm}_{de}_Bc.DEC'.format(\
                                             m=options.mass, 
                                             dm='maj' if options.doMajorana else 'dirac',
                                             de='emu' if options.doElectron else 'mu')),

        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CUEP8M1Settings',  # pythia8CP5Settings ?
                                       # pythia8PSweightsSettings ?
            'processParameters'
        ),
        processParameters = cms.vstring(
            '541:m0 = 6.275',
            '541:tau0 = 0.153'
        ),
        pythia8CUEP8M1Settings = cms.vstring( # these probably remain the same
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024',    # default is 2.28000 
            'MultipartonInteractions:ecmPow=0.25208',   # default is 0.21500
            'MultipartonInteractions:expPow=1.6'        # default is 1.85000
        ),
        pythia8CommonSettings = cms.vstring( 
            'Tune:preferLHAPDF = 2',                    # default is 1 
            'Main:timesAllowErrors = 10000',            # default is 10
            'Check:epTolErr = 0.01',                    # default is 1.0000e-04
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on',                         # default is 100
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on',            # default is false
            'ParticleDecays:tau0Max = 10',   
            'ParticleDecays:allowPhotonRadiation = on'  # default is false
        ) 
        # do we want pythia8PSweightsSettings ?
    ),
    comEnergy = cms.double(13000.0),
    #filterEfficiency = cms.untracked.double(0.0013),  # this will not be used by Pythia, only saved in GenInfo
    maxEventsToPrint = cms.untracked.int32(0),        # max events to print the complete event list information
    pythiaHepMCVerbosity = cms.untracked.bool(False), # to display HepMC information: vertices and particles (not interesting)
    pythiaPylistVerbosity = cms.untracked.int32(0)    # 1 for "normal" verbosity, 11 to display all Pythia Settings
)


if options.doSkipMuonFilter:
  process.ProductionFilterSequence = cms.Sequence(process.generator+process.BFilter) 
else:
  process.ProductionFilterSequence = cms.Sequence(process.generator+process.BFilter+process.SingleMuFilter)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(options.nThr)
process.options.numberOfStreams=cms.untracked.uint32(0)

# set a different offset seed, if you run multiple jobs 
process.RandomNumberGeneratorService.eventSeedOffset=cms.untracked.uint32(options.seedOffset)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions

# Customisation from command line
process.MessageLogger.cerr.FwkReport.reportEvery=100

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
