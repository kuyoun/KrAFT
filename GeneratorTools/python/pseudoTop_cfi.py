import FWCore.ParameterSet.Config as cms

pseudoTop = cms.EDProducer("PseudoTopObjectProducer",
    src = cms.InputTag("genParticles"),
    leptonMinPt = cms.double(5.0),
    jetMinPt = cms.double(20),
    leptonConeSize = cms.double(0.1),
    jetConeSize = cms.double(0.5),
)
