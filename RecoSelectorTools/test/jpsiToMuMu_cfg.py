import FWCore.ParameterSet.Config as cms
runOnMC = True

process = cms.Process("JPSI")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['startup']
#else: process.GlobalTag.globaltag = autoCond['com10']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file://jpsi.root"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.genJpsi = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("abs(pdgId) == 443 && numberOfDaughters == 2 && abs(daughter(0).pdgId) == 13 && abs(daughter(1).pdgId) == 13"),
    minNumber = cms.uint32(1),
    filter = cms.bool(True),
)

process.filterJpsi = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("genJpsi"),
    minNumber = cms.uint32(1),
)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("jpsi.root"),
#    outputCommands = cms.untracked.vstring(
#        "keep *",
#        "drop *_*_*_FILTER",
#    ),
#    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring('p') )
#)
#process.outPath = cms.EndPath(process.out)


process.load("KrAFT.RecoSelectorTools.jpsiToMuMu_cfi")
process.p = cms.Path(
    process.genJpsi * process.filterJpsi
  + process.jpsiToMuMu
)

#process.maxEvents.input = 100
