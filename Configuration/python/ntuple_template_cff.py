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

event = cms.EDAnalyzer("KFlatTreeMaker",
    isMC = cms.bool(False),

    gen = cms.InputTag("genParticles"),
    genJet = cms.InputTag("ak5GenJets"),
    recoToGenJetMap = cms.InputTag("recoToGenJetMap"),
    genJetToPartonsMap = cms.InputTag("genJetToPartonsMap"),

    weight = cms.string("pileupWeight"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    eventCounters = cms.vstring(
    ),

    electron = cms.PSet(
        src = cms.InputTag("goodElectrons"),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    muon = cms.PSet(
        src = cms.InputTag("goodMuons"),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    jet = cms.PSet(
        src = cms.string("goodJets"),
        leptonDeltaR = cms.double(0.5),
        bTagType = cms.string("combinedSecondaryVertexBJetTags"),
    ),
    met = cms.PSet(
        src = cms.string("goodJets"),
    ),
)

ntupleSequence = cms.Sequence(
    goodOfflinePrimaryVertices * pileupWeight
  + goodMuons + goodElectrons 
  * goodJets
  * event
)
