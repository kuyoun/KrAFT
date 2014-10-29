import FWCore.ParameterSet.Config as cms

goodMuons = cms.EDFilter("KMuonSelector",
    rho = cms.InputTag("fixedGridRhoAll"),
    src = cms.InputTag("patMuonsPFlow"),

    precut = cms.string("pt > 20 && abs(eta) < 2.5"),

    cut = cms.string(""),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    maxDz = cms.double(1e9), ## 0.5 for tight muons
    maxDxy = cms.double(1e9), ## 0.2 for tight muons

    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)

goodElectrons = cms.EDFilter("KElectronSelector",
    rho = cms.InputTag("fixedGridRhoAll"),
    src = cms.InputTag("patElectronsPFlow"),

    precut = cms.string("pt > 20 && abs(eta) < 2.5"),

    cut = cms.string(""),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    maxDz = cms.double(1e9), ## for cut based WP: 0.2/0.2/0.1/0.1 for veto/loose/medium/tight electrons
    maxDxy = cms.double(1e9), ## for cut based WP: 0.04/0.02/0.02/0.02 for veto/loose/medium/tight electrons

    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    coneSize = cms.double(0.3),
)

