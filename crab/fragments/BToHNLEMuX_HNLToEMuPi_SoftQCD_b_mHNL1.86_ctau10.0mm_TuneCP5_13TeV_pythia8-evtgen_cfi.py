
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Production Info
configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('B -> mu N X, with long-lived N, m=1.86GeV, ctau=10.0mm'),
    name = cms.untracked.string('B -> mu N X, with long-lived N, m=1.86GeV, ctau=10.0mm'),
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
            #particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.86_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('McRequest/evtGenData/evt_BHNL_mass1.86_ctau10.0_maj.pdl'),
            particle_property_file = cms.FileInPath('evt_BHNL_mass1.86_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('/GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.86_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('evt_BHNL_mass1.86_ctau10.0_maj.pdl'),
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
               '0.1110343958               mu+    myHNL_mu    PHSP;',
               '6.3451291268    anti-D0    mu+    myHNL_mu    PHSP;',
               '12.7484334041    anti-D*0   mu+    myHNL_mu    PHSP;',
               '0.0453147989    pi0        mu+    myHNL_mu    PHSP;',
               '0.0945145916    rho0       mu+    myHNL_mu    PHSP;',
               '0.1110343958               mu+    myHNL_e     PHSP;',
               '6.3451291268    anti-D0    mu+    myHNL_e     PHSP;',
               '12.7484334041    anti-D*0   mu+    myHNL_e     PHSP;',
               '0.0453147989    pi0        mu+    myHNL_e     PHSP;',
               '0.0945145916    rho0       mu+    myHNL_e     PHSP;',
               '0.1105912010                e+     myHNL_mu    PHSP;',
               '6.4190512895     anti-D0    e+     myHNL_mu    PHSP;',
               '12.9299812470     anti-D*0   e+     myHNL_mu    PHSP;',
               '0.0453677968     pi0        e+     myHNL_mu    PHSP;',
               '0.0948415503     rho0       e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB-',
               'Decay myB0',
               '5.8320239744    D-     mu+    myHNL_mu    PHSP;',
               '11.9988075734    D*-    mu+    myHNL_mu    PHSP;',
               '0.0854530865    pi-    mu+    myHNL_mu    PHSP;',
               '0.1753959468    rho-   mu+    myHNL_mu    PHSP;',
               '5.8320239744    D-     mu+    myHNL_e     PHSP;',
               '11.9988075734    D*-    mu+    myHNL_e     PHSP;',
               '0.0854530865    pi-    mu+    myHNL_e     PHSP;',
               '0.1753959468    rho-   mu+    myHNL_e     PHSP;',
               '5.9003913753     D-     e+     myHNL_mu    PHSP;',
               '12.1682532484     D*-    e+     myHNL_mu    PHSP;',
               '0.0855554768     pi-    e+     myHNL_mu    PHSP;',
               '0.1760023963     rho-   e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB0bar',
               'Decay myB0s',
               '5.8833691007    D_s-    mu+    myHNL_mu    PHSP;',
               '11.8359071901    D_s*-   mu+    myHNL_mu    PHSP;',
               '0.1074779124    K-      mu+    myHNL_mu    PHSP;',
               '0.1866963910    K*-     mu+    myHNL_mu    PHSP;',
               '5.8833691007    D_s-    mu+    myHNL_e    PHSP;',
               '11.8359071901    D_s*-   mu+    myHNL_e    PHSP;',
               '0.1074779124    K-      mu+    myHNL_e    PHSP;',
               '0.1866963910    K*-     mu+    myHNL_e    PHSP;',
               '5.9540352881    D_s-     e+     myHNL_mu    PHSP;',
               '12.0089530870    D_s*-    e+     myHNL_mu    PHSP;',
               '0.1077017334    K-       e+     myHNL_mu    PHSP;',
               '0.1873563067    K*-      e+     myHNL_mu    PHSP;',
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

