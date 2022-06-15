
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Production Info
configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('B -> mu N X, with long-lived N, m=1.59GeV, ctau=10.0mm'),
    name = cms.untracked.string('B -> mu N X, with long-lived N, m=1.59GeV, ctau=10.0mm'),
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
            #particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.59_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('McRequest/evtGenData/evt_BHNL_mass1.59_ctau10.0_maj.pdl'),
            particle_property_file = cms.FileInPath('evt_BHNL_mass1.59_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('/GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.59_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('evt_BHNL_mass1.59_ctau10.0_maj.pdl'),
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
               '0.0875527399               mu+    myHNL_mu    PHSP;',
               '9.6099116884    anti-D0    mu+    myHNL_mu    PHSP;',
               '20.3872929095    anti-D*0   mu+    myHNL_mu    PHSP;',
               '0.0497044520    pi0        mu+    myHNL_mu    PHSP;',
               '0.1140293176    rho0       mu+    myHNL_mu    PHSP;',
               '0.0875527399               mu+    myHNL_e     PHSP;',
               '9.6099116884    anti-D0    mu+    myHNL_e     PHSP;',
               '20.3872929095    anti-D*0   mu+    myHNL_e     PHSP;',
               '0.0497044520    pi0        mu+    myHNL_e     PHSP;',
               '0.1140293176    rho0       mu+    myHNL_e     PHSP;',
               '0.0870994934                e+     myHNL_mu    PHSP;',
               '9.6951380612     anti-D0    e+     myHNL_mu    PHSP;',
               '20.6000209850     anti-D*0   e+     myHNL_mu    PHSP;',
               '0.0497526896     pi0        e+     myHNL_mu    PHSP;',
               '0.1143621453     rho0       e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB-',
               'Decay myB0',
               '8.8468934079    D-     mu+    myHNL_mu    PHSP;',
               '19.1427463827    D*-    mu+    myHNL_mu    PHSP;',
               '0.0938491901    pi-    mu+    myHNL_mu    PHSP;',
               '0.2115967121    rho-   mu+    myHNL_mu    PHSP;',
               '8.8468934079    D-     mu+    myHNL_e     PHSP;',
               '19.1427463827    D*-    mu+    myHNL_e     PHSP;',
               '0.0938491901    pi-    mu+    myHNL_e     PHSP;',
               '0.2115967121    rho-   mu+    myHNL_e     PHSP;',
               '8.9257747581     D-     e+     myHNL_mu    PHSP;',
               '19.3410981367     D*-    e+     myHNL_mu    PHSP;',
               '0.0939428832     pi-    e+     myHNL_mu    PHSP;',
               '0.2122140269     rho-   e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB0bar',
               'Decay myB0s',
               '8.9791002787    D_s-    mu+    myHNL_mu    PHSP;',
               '19.0715872030    D_s*-   mu+    myHNL_mu    PHSP;',
               '0.1227009261    K-      mu+    myHNL_mu    PHSP;',
               '0.2253821402    K*-     mu+    myHNL_mu    PHSP;',
               '8.9791002787    D_s-    mu+    myHNL_e    PHSP;',
               '19.0715872030    D_s*-   mu+    myHNL_e    PHSP;',
               '0.1227009261    K-      mu+    myHNL_e    PHSP;',
               '0.2253821402    K*-     mu+    myHNL_e    PHSP;',
               '9.0610153199    D_s-     e+     myHNL_mu    PHSP;',
               '19.2748183361    D_s*-    e+     myHNL_mu    PHSP;',
               '0.1229205165    K-       e+     myHNL_mu    PHSP;',
               '0.2260710881    K*-      e+     myHNL_mu    PHSP;',
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

