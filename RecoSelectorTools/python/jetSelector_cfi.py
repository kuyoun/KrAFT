import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

goodJets = cms.EDFilter("KJetSelector",
    debug = cms.untracked.bool(False),

    isMC = cms.bool(False),

    jet = cms.InputTag("patJetsPFlow"),
    met = cms.InputTag("patMETsPFlow"),

    selection = cms.PSet(
        jetId = pfJetIDSelector,
        minPt = cms.double(30),
        maxEta = cms.double(2.5),
    ),

    cleaning = cms.PSet(
        overlapDeltaR = cms.double(0.5),
        overlapCands = cms.VInputTag(
            cms.InputTag("goodMuons"),
            cms.InputTag("goodElectrons"),
        ),
        #cleanMethod = cms.string("subtract"),
        cleanMethod = cms.string("subtractAndRestore"),
        #cleanMethod = cms.string("cleanAll"),
    ),
    jecFileRD = cms.string("KCMSAnalyses/RecoSelectorTools/data/JEC/Summer13_V4/Summer13_V4_DATA_Uncertainty_AK5PFchs.txt"),
    jecFileMC = cms.string("KCMSAnalyses/RecoSelectorTools/data/JEC/Summer13_V4/Summer13_V4_MC_Uncertainty_AK5PFchs.txt"),

    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
)

