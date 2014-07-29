import FWCore.ParameterSet.Config as cms

from KrAFT.GeneratorTools.pileupWeight_cff import *
from KrAFT.GeneratorTools.pdfWeight_cff import *
from KrAFT.RecoSelectorTools.leptonSelector_cfi import *
from KrAFT.RecoSelectorTools.jetSelector_cfi import *
from KrAFT.RecoSelectorTools.jpsiSelector_cfi import *

TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

from KrAFT.GenericNtuple.genericNtupleMaker_cfi import *
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
MuJets.jetMET.minNumber = 3

ElJets = event.clone()
ElJets.electron.minNumber = 1
ElJets.jetMET.minNumber = 3

nEventsNtupleElEl = cms.EDProducer("EventCountProducer")
nEventsNtupleMuMu = cms.EDProducer("EventCountProducer")
nEventsNtupleMuEl = cms.EDProducer("EventCountProducer")
nEventsNtupleMuJets = cms.EDProducer("EventCountProducer")
nEventsNtupleElJets = cms.EDProducer("EventCountProducer")

ntupleSequenceElEl = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleElEl
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu+jpsiToElEl
  * ElEl
)

ntupleSequenceMuMu = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleMuMu
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu+jpsiToElEl
  * MuMu
)

ntupleSequenceMuEl = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleMuEl
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu+jpsiToElEl
  * MuEl
)

ntupleSequenceMuJets = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleMuJets
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu+jpsiToElEl
  * MuJets
)

ntupleSequenceElJets = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleElJets
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu+jpsiToElEl
  * ElJets
)

