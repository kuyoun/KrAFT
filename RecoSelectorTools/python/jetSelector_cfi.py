import FWCore.ParameterSet.Config as cms

cleanJets = cms.EDFilter("TopCleanJetSelector",
    debug = cms.untracked.bool(False),
    src = cms.InputTag("cmgPFJetSelCHS"),
    cut = cms.string(
        " abs(eta) < 2.5 && pt > 30 && nConstituents > 1"
        " && component(5).fraction < 0.99 && component(4).fraction < 0.99"
        " && (abs(eta) >= 2.4 || component(2).fraction < 0.99 )"
        " && (abs(eta) >= 2.4 || component(1).fraction > 0 )"
        " && (abs(eta) >= 2.4 || component(1).number > 0 ) "
    ),
    overlapCands = cms.VInputTag(
        cms.InputTag("cmgMuonSel"),
        cms.InputTag("cmgElectronSel"),
    ),
    overlapDeltaR = cms.double(0.5),
    #cleanMethod = cms.untracked.string("subtract"),
    cleanMethod = cms.untracked.string("subtractAndRestore"),
    #cleanMethod = cms.untracked.string(""),
    minPt = cms.untracked.double(30),
    maxEta = cms.untracked.double(2.5),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
)

