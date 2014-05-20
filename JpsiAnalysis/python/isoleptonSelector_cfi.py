import FWCore.ParameterSet.Config as cms

goodMuonsForJpsi = cms.EDFilter("KIsoMuonSelector",
    src = cms.InputTag("goodMuons"),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
)

goodElectronsForJpsi = cms.EDFilter("KIsoElectronSelector",
    src = cms.InputTag("goodElectrons"),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
)

