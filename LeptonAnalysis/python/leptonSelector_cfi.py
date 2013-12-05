import FWCore.ParameterSet.Config as cms

goodMuons = cms.EDFilter("KCMuonSelector",
    rho = cms.InputTag("kt6PFJets", "rho"),
    src = cms.InputTag("patMuonsWithTrigger"),
    precut = cms.string("pt > 20 && abs(eta) < 2.5"),
    cut = cms.string(""),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)

goodElectrons = cms.EDFilter("KCElectronSelector",
    rho = cms.InputTag("kt6PFJets", "rho"),
    src = cms.InputTag("patElectronsWithTrigger"),
    precut = cms.string("pt > 20 && abs(eta) < 2.5"),
    cut = cms.string(""),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)

