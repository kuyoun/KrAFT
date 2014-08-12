import FWCore.ParameterSet.Config as cms

runOnMC = False

process = cms.Process("KrAFT")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
from Configuration.AlCa.autoCond import autoCond
if runOnMC: process.GlobalTag.globaltag = autoCond['startup']
else: process.GlobalTag.globaltag = autoCond['com10']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012D/SingleMu/AOD/22Jan2013-v1/10000/FE618190-93A7-E211-B50C-90E6BA19A215.root'
    ),
)

process.out = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_TriggerResults_*_HLT',
        'keep *_flat*_*_*',
        'keep *_TriggerResults_*_%s' % process.process,
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("CANDSEL"),
    ),
)

from KrAFT.Configuration.customise_cff import *
customisePAT(process, runOnMC=runOnMC, outputModules=[])

process.load("KrAFT.Configuration.flatEDM_RD_cff")
process.load("KrAFT.Configuration.commonFilters_RD_cff")

process.CANDSEL = cms.Path(
    process.preFilterSequence
  + process.patPF2PATSequencePFlow
  + process.analysisObjectSequence
)

process.output = cms.EndPath(process.out)

