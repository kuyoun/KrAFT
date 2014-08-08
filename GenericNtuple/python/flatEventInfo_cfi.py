import FWCore.ParameterSet.Config as cms

flatEventInfo = cms.EDProducer("FlatEventInfoProducer",
    genInfo = cms.InputTag("generator"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    hltProcessName = cms.string("HLT"),
    HLT = cms.PSet(
        DoubleMu = cms.vstring(
            "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"
        ),
        DoubleElectron = cms.vstring(
            "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
        ),
        MuEG = cms.vstring(
            "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
            "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
        ),
        MuJet = cms.vstring(
            "HLT_IsoMu17_eta2p1_*Central*", "HLT_Mu17_eta2p1_*Central*",
        ),
        ElJet = cms.vstring(
            "HLT_Ele25_CaloIdVL_*", "HLT_Ele25_CaloIdVT_*",
        ),
    ),
)

