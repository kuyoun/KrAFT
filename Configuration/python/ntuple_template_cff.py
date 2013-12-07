import FWCore.ParameterSet.Config as cms

from KCMSAnalyses.GeneratorTools.pileupWeight_cff import *
from KCMSAnalyses.RecoSelectorTools.leptonSelector_cfi import *
from KCMSAnalyses.RecoSelectorTools.jetSelector_cfi import *

goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    src = cms.InputTag('offlinePrimaryVertices'),
    filterParams =  cms.PSet(
        minNdof = cms.double(4.),
        maxZ    = cms.double(24.),
        maxRho  = cms.double(2.)
    ),
    filter = cms.bool(True),
)

TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

from KCMSAnalyses.FlatTree.flatTreeMaker_cfi import *

MuMu = event.clone()
ElEl = event.clone()
MuEl = event.clone()
MuMu.muon.minNumber = 2
ElEl.electron.minNumber = 2
MuEl.muon.minNumber = 1
MuEl.electron.minNumber = 1

ntupleSequence = cms.Sequence(
    goodOfflinePrimaryVertices * pileupWeight
  + goodMuons + goodElectrons
  * goodJets
  * MuMu + ElEl + MuEl
)
