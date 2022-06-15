
import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

# Production Info
configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('B -> mu N X, with long-lived N, m=1.18GeV, ctau=10.0mm'),
    name = cms.untracked.string('B -> mu N X, with long-lived N, m=1.18GeV, ctau=10.0mm'),
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
            #particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.18_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('McRequest/evtGenData/evt_BHNL_mass1.18_ctau10.0_maj.pdl'),
            particle_property_file = cms.FileInPath('evt_BHNL_mass1.18_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('/GeneratorInterface/EvtGenInterface/data/evt_BHNL_mass1.18_ctau10.0_maj.pdl'),
            #particle_property_file = cms.FileInPath('evt_BHNL_mass1.18_ctau10.0_maj.pdl'),
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
               '0.0528288230               mu+    myHNL_mu    PHSP;',
               '15.3059351103    anti-D0    mu+    myHNL_mu    PHSP;',
               '34.2046330433    anti-D*0   mu+    myHNL_mu    PHSP;',
               '0.0548625789    pi0        mu+    myHNL_mu    PHSP;',
               '0.1413200919    rho0       mu+    myHNL_mu    PHSP;',
               '0.0528288230               mu+    myHNL_e     PHSP;',
               '15.3059351103    anti-D0    mu+    myHNL_e     PHSP;',
               '34.2046330433    anti-D*0   mu+    myHNL_e     PHSP;',
               '0.0548625789    pi0        mu+    myHNL_e     PHSP;',
               '0.1413200919    rho0       mu+    myHNL_e     PHSP;',
               '0.0523675452                e+     myHNL_mu    PHSP;',
               '15.4056192464     anti-D0    e+     myHNL_mu    PHSP;',
               '34.4542433396     anti-D*0   e+     myHNL_mu    PHSP;',
               '0.0549046977     pi0        e+     myHNL_mu    PHSP;',
               '0.1416566523     rho0       e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB-',
               'Decay myB0',
               '14.1138647364    D-     mu+    myHNL_mu    PHSP;',
               '32.0419090204    D*-    mu+    myHNL_mu    PHSP;',
               '0.1037695316    pi-    mu+    myHNL_mu    PHSP;',
               '0.2622188560    rho-   mu+    myHNL_mu    PHSP;',
               '14.1138647364    D-     mu+    myHNL_e     PHSP;',
               '32.0419090204    D*-    mu+    myHNL_e     PHSP;',
               '0.1037695316    pi-    mu+    myHNL_e     PHSP;',
               '0.2622188560    rho-   mu+    myHNL_e     PHSP;',
               '14.2062015119     D-     e+     myHNL_mu    PHSP;',
               '32.2744399303     D*-    e+     myHNL_mu    PHSP;',
               '0.1038521727     pi-    e+     myHNL_mu    PHSP;',
               '0.2628430703     rho-   e+     myHNL_mu    PHSP;',
               'Enddecay',
               'CDecay myB0bar',
               'Decay myB0s',
               '14.4242314332    D_s-    mu+    myHNL_mu    PHSP;',
               '32.2238264934    D_s*-   mu+    myHNL_mu    PHSP;',
               '0.1425216353    K-      mu+    myHNL_mu    PHSP;',
               '0.2807170168    K*-     mu+    myHNL_mu    PHSP;',
               '14.4242314332    D_s-    mu+    myHNL_e    PHSP;',
               '32.2238264934    D_s*-   mu+    myHNL_e    PHSP;',
               '0.1425216353    K-      mu+    myHNL_e    PHSP;',
               '0.2807170168    K*-     mu+    myHNL_e    PHSP;',
               '14.5206790573    D_s-     e+     myHNL_mu    PHSP;',
               '32.4624474685    D_s*-    e+     myHNL_mu    PHSP;',
               '0.1427342475    K-       e+     myHNL_mu    PHSP;',
               '0.2814434058    K*-      e+     myHNL_mu    PHSP;',
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

