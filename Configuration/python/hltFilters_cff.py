import FWCore.ParameterSet.Config as cms

from HLTrigger.HLTfilters.hltHighLevel_cfi import *

hltHighLevel.throw = False

hltElEl = hltHighLevel.clone(HLTPaths = cms.vstring(
    "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
))
hltMuMu = hltHighLevel.clone(HLTPaths = cms.vstring(
    "HLT_Mu17_Mu8_v*", "HLT_Mu17_TkMu8_v*"
))
hltMuEl = hltHighLevel.clone(HLTPaths = cms.vstring(
    "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
    "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",
))
hltMuJets = hltHighLevel.clone(HLTPaths = cms.vstring(
    "HLT_IsoMu17_eta2p1_*Central*", "HLT_Mu17_eta2p1_*Central*",
))
hltElJets = hltHighLevel.clone(HLTPaths = cms.vstring(
    "HLT_Ele25_CaloIdVL_*", "HLT_Ele25_CaloIdVT_*", 
))

