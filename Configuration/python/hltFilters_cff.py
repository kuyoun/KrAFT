import FWCore.ParameterSet.Config as cms

#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#
#hltHighLevel.throw = False

hltElEl = [
    "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
]
hltMuMu = [
    "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"
]
hltMuEl = [
    "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
    "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
]
hltMuJets = [
    "HLT_IsoMu17_eta2p1_*Central*", "HLT_Mu17_eta2p1_*Central*",
]
hltElJets = [
    "HLT_Ele25_CaloIdVL_*", "HLT_Ele25_CaloIdVT_*", 
]

