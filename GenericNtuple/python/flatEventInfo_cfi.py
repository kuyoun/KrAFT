import FWCore.ParameterSet.Config as cms

flatEventInfo = cms.EDProducer("FlatEventInfoProducer",
    genInfo = cms.InputTag("generator"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    hltPaths = cms.vstring(),
)

