import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/relval/CMSSW_5_3_11/RelValTTbar/GEN-SIM-RECO/START53_LV3_Alca7TeV_14Jun2013-v1/00000/905DC97D-F3D4-E211-9B9E-003048F1C9CA.root',
    '/store/relval/CMSSW_5_3_11/RelValTTbar/GEN-SIM-RECO/START53_LV3_Alca7TeV_14Jun2013-v1/00000/A4086D2F-F3D4-E211-A474-003048F1C5FA.root',
    '/store/relval/CMSSW_5_3_11/RelValTTbar/GEN-SIM-RECO/START53_LV3_Alca7TeV_14Jun2013-v1/00000/CA8F250D-F1D4-E211-895D-001E67398B2E.root',
]

process.load("KrAFT.GeneratorTools.genJetAssociation_cff")
process.load("KrAFT.GeneratorTools.pileupWeight_cff")
process.load("KrAFT.GeneratorTools.pdfWeight_cff")
process.load("KrAFT.GeneratorTools.pseudoTop_cfi")

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_Ana",
    )
)

#process.outPath = cms.EndPath(process.out)

process.wlnu = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("status == 3 && (abs(pdgId) == 11 || abs(pdgId) == 13) && abs(mother.pdgId) == 24"),
)

process.dileptonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("wlnu"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
)

process.genAna = cms.EDAnalyzer("PseudoTopAnalyzer",
    genParticles = cms.InputTag("genParticles"),
    genJets = cms.InputTag("ak5GenJets"),
    pseudoTop = cms.InputTag("pseudoTop"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

process.p = cms.Path(
    process.wlnu * process.dileptonFilter
  + process.pseudoTop * process.genAna
)

