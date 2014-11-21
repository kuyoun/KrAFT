import FWCore.ParameterSet.Config as cms

partons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop *",
        "drop pt <= 0",
#        "keep status = 3", # For the pythia
        "keep abs(pdgId) < 10",
        "keep abs(pdgId) == 24 || abs(pdgId) == 23",	
    ),
)

