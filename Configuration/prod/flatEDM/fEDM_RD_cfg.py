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

process.load("TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi")
process.goodOfflinePrimaryVertices.filter = True

process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
if runOnMC: process.eventCleaning += process.eventCleaningMC
else: process.eventCleaning += process.eventCleaningData

process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsClean = cms.EDProducer("EventCountProducer")
process.nEventsPAT   = cms.EDProducer("EventCountProducer")

## Default path
process.load("KrAFT.GeneratorTools.pileupWeight_cff")
process.load("KrAFT.GeneratorTools.pdfWeight_cff")
process.load("KrAFT.GeneratorTools.pseudoTop_cfi")
process.load("KrAFT.RecoSelectorTools.leptonSelector_cfi")
process.load("KrAFT.RecoSelectorTools.jetSelector_cfi")
process.load("KrAFT.RecoSelectorTools.jpsiSelector_cfi")
process.load("KrAFT.GenericNtuple.flatEventInfo_cfi")
process.load("KrAFT.GenericNtuple.flatCands_cfi")
process.goodJets.isMC = runOnMC
delattr(process.flatJets.variables, 'res')
delattr(process.flatJets.variables, 'resUp')
delattr(process.flatJets.variables, 'resDn')
delattr(process, 'flatMETsRes')
delattr(process, 'flatMETsResUp')
delattr(process, 'flatMETsResDn')

process.analysisObjectSequence = cms.Sequence(
    process.goodMuons + process.goodElectrons * process.goodJets
  * process.jpsiToMuMu# + process.jpsiToElEl

  + process.flatEventInfo
  * process.flatMuons + process.flatElectrons + process.flatJets
  + process.flatMETs + process.flatMETsUp + process.flatMETsDn
  + process.flatJpsiMuMu# + process.flatJpsiElEl
)

process.CANDSEL = cms.Path(
    process.nEventsTotal
  + process.goodOfflinePrimaryVertices + process.eventCleaning + process.nEventsClean
  + process.patPF2PATSequencePFlow + process.nEventsPAT
  + process.analysisObjectSequence
)

process.output = cms.EndPath(process.out)

