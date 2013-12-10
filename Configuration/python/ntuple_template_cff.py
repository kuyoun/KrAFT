import FWCore.ParameterSet.Config as cms

from KCMSAnalyses.GeneratorTools.pileupWeight_cff import *
from KCMSAnalyses.RecoSelectorTools.leptonSelector_cfi import *
from KCMSAnalyses.RecoSelectorTools.jetSelector_cfi import *

TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

from KCMSAnalyses.FlatTree.flatTreeMaker_cfi import *
event.eventCounters = ["nEventsTotal", "nEventsClean", "nEventsPAT",]

MuMu = event.clone()
ElEl = event.clone()
MuEl = event.clone()
MuMu.muon.minNumber = 2
ElEl.electron.minNumber = 2
MuEl.muon.minNumber = 1
MuEl.electron.minNumber = 1

MuJets = event.clone()
MuJets.muon.minNumber = 1
#MuJets.muon.maxNumber = 1
#MuJets.electron.maxNumber = 0

ElJets = event.clone()
ElJets.electron.minNumber = 1
#ElJets.electron.maxNumber = 1
#ElJets.muon.maxNumber = 0

nEventsNtupleElEl = cms.EDProducer("EventCountProducer")
nEventsNtupleMuMu = cms.EDProducer("EventCountProducer")
nEventsNtupleMuEl = cms.EDProducer("EventCountProducer")
nEventsNtupleMuJets = cms.EDProducer("EventCountProducer")
nEventsNtupleElJets = cms.EDProducer("EventCountProducer")

ntupleSequenceElEl = cms.Sequence(
    pileupWeight
  + nEventsNtupleElEl
  + goodMuons + goodElectrons * goodJets
  * ElEl
)

ntupleSequenceMuMu = cms.Sequence(
    pileupWeight
  + nEventsNtupleMuMu
  + goodMuons + goodElectrons * goodJets
  * MuMu
)

ntupleSequenceMuEl = cms.Sequence(
    pileupWeight
  + nEventsNtupleMuEl
  + goodMuons + goodElectrons * goodJets
  * MuEl
)

ntupleSequenceMuJets = cms.Sequence(
    pileupWeight
  + nEventsNtupleMuJets
  + goodMuons + goodElectrons * goodJets
  * MuJets
)

ntupleSequenceElJets = cms.Sequence(
    pileupWeight
  + nEventsNtupleElJets
  + goodMuons + goodElectrons * goodJets
  * ElJets
)

