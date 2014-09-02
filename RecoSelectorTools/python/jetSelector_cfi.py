import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

jetUncertainties = cms.EDProducer("KJetMetUncProducer",
    jet = cms.InputTag("patJetsPFlow"),
    met = cms.InputTag("patMETsPFlow"),

    jecFile = cms.string("KrAFT/RecoSelectorTools/data/JEC/Summer13_V4/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt"),
    jecLevels = cms.vstring("Total", "Absolute", 
        "RelativeJEREC1", "RelativeJEREC2",
        "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativeFSR",
        "PileUpDataMC", "PileUpPtBB", "PileUpPtEC", "PileUpPtHF",
        "SubTotalPileUp", "SubTotalRelative",
    ),
)

goodJets = cms.EDFilter("KCleanJetSelector",
    jet = cms.InputTag("patJetsPFlow"),
    jes = cms.VInputTag(),

    selection = cms.PSet(
        jetId = pfJetIDSelector,
        minPt = cms.double(30),
        maxEta = cms.double(2.5),
    ),

    cleaning = cms.PSet(
        doClean = cms.bool(False),
        overlapDeltaR = cms.double(0.5),
        overlapCands = cms.VInputTag(
            cms.InputTag("goodMuons"),
            cms.InputTag("goodElectrons"),
        ),
    ),

    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
)

