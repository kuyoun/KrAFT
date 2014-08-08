import FWCore.ParameterSet.Config as cms

flatEventInfo = cms.EDProducer("FlatEventInfoProducer",
    genInfo = cms.InputTag("generator"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    hltProcessName = cms.string("HLT"),
    HLT = cms.PSet(
        DoubleMu = cms.vstring(),
        DoubleElectron = cms.vstring(),
        MuEG = cms.vstring(),
    ),
)

