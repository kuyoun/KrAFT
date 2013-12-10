import FWCore.ParameterSet.Config as cms

from KCMSAnalyses.GeneratorTools.pileupWeight_cff import *
from KCMSAnalyses.RecoSelectorTools.leptonSelector_cfi import *
from KCMSAnalyses.RecoSelectorTools.jetSelector_cfi import *

TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

from KCMSAnalyses.FlatTree.flatTreeMaker_cfi import *
event.eventCounters = ["nEventsTotal", "nEventsClean", "nEventsPAT", "nEventsNtuple"]

MuMu = event.clone()
ElEl = event.clone()
MuEl = event.clone()
MuMu.muon.minNumber = 2
ElEl.electron.minNumber = 2
MuEl.muon.minNumber = 1
MuEl.electron.minNumber = 1

MuJet = event.clone()
MuJet.muon.minNumber = 1
MuJet.muon.maxNumber = 1
MuJet.electron.maxNumber = 0

ElJet = event.clone()
ElJet.muon.maxNumber = 0
ElJet.electron.minNumber = 1
ElJet.electron.maxNumber = 1

nEventsNtuple = cms.EDProducer("EventCountProducer")

ntupleSequenceDilepton = cms.Sequence(
    pileupWeight
  + nEventsNtuple
  + goodMuons + goodElectrons * goodJets
  * MuMu + ElEl + MuEl
)

ntupleSequenceMuJet = cms.Sequence(
    pileupWeight
  + nEventsNtuple
  + goodMuons + goodElectrons * goodJets
  * MuJet
)

ntupleSequenceElJet = cms.Sequence(
    pileupWeight
  + nEventsNtuple
  + goodMuons + goodElectrons * goodJets
  * ElJet
)

