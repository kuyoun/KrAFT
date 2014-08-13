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
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
from Configuration.AlCa.autoCond import autoCond
if runOnMC: process.GlobalTag.globaltag = autoCond['startup']
else: process.GlobalTag.globaltag = autoCond['com10']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
				'/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10000/0014C5C0-1A8F-E211-AD90-0026189438A9.root',
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
        'keep edmMergeableCounter_*_*_*',
        'keep *_partons_*_*',
        #'keep *_pseudoTop_*_*',
        'keep *_pileupWeight_*_*',
        'keep *_pdfWeight_*_*',
        'keep *_flat*_*_*',
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

process.partons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop *",
        "drop pt <= 0",
        "keep status = 3", # For the pythia
				"+keep pdgId = 443", 
				"keep+ pdgId = 443" 
    ),
)

process.analysisObjectSequence = cms.Sequence(
    process.pileupWeight + process.pdfWeight
  + process.goodMuons + process.goodElectrons * process.goodJets
  * process.jpsiToMuMu + process.jpsiToElEl

  + process.flatEventInfo
  * process.flatMuons + process.flatElectrons + process.flatJets
  + process.flatMETs + process.flatMETsUp + process.flatMETsDn
  + process.flatMETsRes + process.flatMETsResUp + process.flatMETsResDn
  + process.flatJpsiMuMu + process.flatJpsiElEl
)

if runOnMC is True : 

	process.pGen = cms.Path(
  	  process.pseudoTop
	  + process.partons
	  * process.flatPseudoTopLepton + process.flatPseudoTopNu + process.flatPseudoTopJet
	)

process.p = cms.Path(
    process.nEventsTotal
  + process.goodOfflinePrimaryVertices + process.eventCleaning + process.nEventsClean
  + process.patPF2PATSequencePFlow + process.nEventsPAT
  + process.analysisObjectSequence
)

process.output = cms.EndPath(process.out)

