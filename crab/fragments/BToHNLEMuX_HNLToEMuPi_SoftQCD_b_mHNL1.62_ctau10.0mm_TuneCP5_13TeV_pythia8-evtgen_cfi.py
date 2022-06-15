
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Production Info
configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('B -> mu N X, with long-lived N, m=1.62GeV, ctau=10.0mm'),
    name = cms.untracked.string('B -> mu N X, with long-lived N, m=1.62GeV, ctau=10.0mm'),
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

DoubleLeptonFilter = cms.EDFilter("MCParticlePairFilter",
    MaxEta = cms.untracked.vdouble(1.55, 2.45),
    MinEta = cms.untracked.vdouble(-1.55, -2.45),
    MinPt = cms.untracked.vdouble(6.8, 1.0),
    ParticleID1 = cms.untracked.vint32(-13, 13),
    ParticleID2 = cms.untracked.vint32(-11, 11,-13,13),
    MaxInvMass = cms.untracked.double(10.),
    Status = cms.untracked.vint32(1, 1),
)

TriggerMuonFilter = cms.EDFilter("PythiaFilterMultiMother", 
    MaxEta = cms.untracked.double(1.55),
    MinEta = cms.untracked.double(-1.55),
    MinPt = cms.untracked.double(6.8), 
    ParticleID = cms.untracked.int32(13),
    MotherIDs = cms.untracked.vint32(521, 511, 531, 9900015),
)

HNLPionFilter = cms.EDFilter("PythiaFilterMultiMother", 
    MaxEta = cms.untracked.double(10.),
    MinEta = cms.untracked.double(-10.),
    MinPt = cms.untracked.double(0.5), 
    ParticleID = cms.untracked.int32(211),
    MotherIDs = cms.untracked.vint32(9900015),
)

HNLDisplacementFilter = cms.EDFilter("MCDisplacementFilter",
    ParticleIDs = cms.vint32(9900015),
    LengMin = cms.double(0), 
    LengMax = cms.double(1300), 
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
                'myHNL_mu',
                'myHNL_e',
            ),
            
            operates_on_particles = cms.vint32(521, -521, 511, -511, 531, -531, 9900015), 
            #particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.62_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('McRequest/evtGenData/evt_BHNL_mass1.62_ctau10.0_maj.pdl'),
            particle_property_file = cms.FileInPath('evt_BHNL_mass1.62_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('/GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.62_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('evt_BHNL_mass1.62_ctau10.0_maj.pdl'),
            user_decay_embedded = cms.vstring(
              
               'Alias myB+ B+',
               'Alias myB- B-',
               'Alias myB0 B0',
               'Alias myB0bar anti-B0',
               'Alias myB0s B_s0',
               'Alias myB0sbar anti-B_s0',
               'Alias myHNL_mu hnl',
               'Alias myHNL_e hnl',
               'ChargeConj myB+ myB-',
               'ChargeConj myB0 myB0bar',
               'ChargeConj myB0s myB0sbar', 
               'Decay myB+',
               '0.0901837932               mu+    myHNL_mu    PHSP;',
               '9.2211008009    anti-D0    mu+    myHNL_mu    PHSP;',
               '19.4631687928    anti-D*0   mu+    myHNL_mu    PHSP;',
               '0.0492541084    pi0        mu+    myHNL_mu    PHSP;',
               '0.1118948257    rho0       mu+    myHNL_mu    PHSP;',
               '0.0901837932               mu+    myHNL_e     PHSP;',
               '9.2211008009    anti-D0    mu+    myHNL_e     PHSP;',
               '19.4631687928    anti-D*0   mu+    myHNL_e     PHSP;',
               '0.0492541084    pi0        mu+    myHNL_e     PHSP;',
               '0.1118948257    rho0       mu+    myHNL_e     PHSP;',
               '0.0897314412                e+     myHNL_mu    PHSP;',
               '9.3051403943     anti-D0    e+     myHNL_mu    PHSP;',
               '19.6726940566     anti-D*0   e+     myHNL_mu    PHSP;',
               '0.0493028512     pi0        e+     myHNL_mu    PHSP;',
               '0.1122271629     rho0       e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB-',
               'Decay myB0',
               '8.4876582131    D-     mu+    myHNL_mu    PHSP;',
               '18.2791164523    D*-    mu+    myHNL_mu    PHSP;',
               '0.0929862529    pi-    mu+    myHNL_mu    PHSP;',
               '0.2076372345    rho-   mu+    myHNL_mu    PHSP;',
               '8.4876582131    D-     mu+    myHNL_e     PHSP;',
               '18.2791164523    D*-    mu+    myHNL_e     PHSP;',
               '0.0929862529    pi-    mu+    myHNL_e     PHSP;',
               '0.2076372345    rho-   mu+    myHNL_e     PHSP;',
               '8.5654355341     D-     e+     myHNL_mu    PHSP;',
               '18.4744995476     D*-    e+     myHNL_mu    PHSP;',
               '0.0930808668     pi-    e+     myHNL_mu    PHSP;',
               '0.2082536416     rho-   e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB0bar',
               'Decay myB0s',
               '8.6092989744    D_s-    mu+    myHNL_mu    PHSP;',
               '18.1943834768    D_s*-   mu+    myHNL_mu    PHSP;',
               '0.1210807074    K-      mu+    myHNL_mu    PHSP;',
               '0.2211221675    K*-     mu+    myHNL_mu    PHSP;',
               '8.6092989744    D_s-    mu+    myHNL_e    PHSP;',
               '18.1943834768    D_s*-   mu+    myHNL_e    PHSP;',
               '0.1210807074    K-      mu+    myHNL_e    PHSP;',
               '0.2211221675    K*-     mu+    myHNL_e    PHSP;',
               '8.6900292952    D_s-     e+     myHNL_mu    PHSP;',
               '18.3945204390    D_s*-    e+     myHNL_mu    PHSP;',
               '0.1213007993    K-       e+     myHNL_mu    PHSP;',
               '0.2218080584    K*-      e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB0sbar',
               'Decay myHNL_mu',
               '0.5     mu-    pi+    PHSP;',
               '0.5     mu+    pi-    PHSP;',
               'Enddecay',
               'Decay myHNL_e',
               '0.5     e-    pi+    PHSP;',
               '0.5     e+    pi-    PHSP;',
               'Enddecay',
               'End',      

            )
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
    maxEventsToPrint = cms.untracked.int32(0),        
    pythiaHepMCVerbosity = cms.untracked.bool(False), 
    pythiaPylistVerbosity = cms.untracked.int32(0),
)

ProductionFilterSequence = cms.Sequence(generator+BFilter+DoubleLeptonFilter+TriggerMuonFilter+HNLPionFilter+HNLDisplacementFilter)

