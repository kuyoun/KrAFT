import FWCore.ParameterSet.Config as cms

goodMuons = cms.EDFilter("KMuonSelector",
    rho = cms.InputTag("kt6PFJets", "rho"),
    src = cms.InputTag("patMuonsPFlow"),
    precut = cms.string("pt > 20 && abs(eta) < 2.5"),
    cut = cms.string(""),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)

goodElectrons = cms.EDFilter("KElectronSelector",
    rho = cms.InputTag("kt6PFJets", "rho"),
    src = cms.InputTag("patElectronsPFlow"),
    precut = cms.string("pt > 20 && abs(eta) < 2.5"),
    cut = cms.string(""),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)

