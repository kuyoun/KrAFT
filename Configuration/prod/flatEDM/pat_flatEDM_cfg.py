import FWCore.ParameterSet.Config as cms

runOnMC = True

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

#process.Tracer = cms.Service("Tracer")
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_7_0_6/RelValTTbarLepton_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/4EA445A2-3BFA-E311-B066-0026189438BA.root',
        '/store/relval/CMSSW_7_0_6/RelValTTbarLepton_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/F0AADC82-39FA-E311-8050-002354EF3BD0.root',
        '/store/relval/CMSSW_7_0_6/RelValTTbarLepton_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/FC321902-35FA-E311-B1FE-002618943957.root',
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
#process.load("KrAFT.Configuration.commonFilters_cff")
process.load("KrAFT.GeneratorTools.pileupWeight_cff")
process.load("KrAFT.GeneratorTools.pdfWeight_cff")
process.load("KrAFT.GeneratorTools.pseudoTop_cfi")
process.load("KrAFT.RecoSelectorTools.leptonSelector_cfi")
process.load("KrAFT.RecoSelectorTools.jetSelector_cfi")
process.load("KrAFT.RecoSelectorTools.jpsiSelector_cfi")
process.load("KrAFT.GenericNtuple.flatEventInfo_cfi")
process.load("KrAFT.GenericNtuple.flatCands_cfi")

process.goodMuons.rho = "fixedGridRhoFastjetAll"
process.goodElectrons.rho = "fixedGridRhoFastjetAll"
process.goodJets.isMC = runOnMC

process.partons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop *",
        "drop pt <= 0",
        "keep status = 3", # For the pythia
    ),
)

"""
process.analysisObjectSequence = cms.Sequence(
    process.pileupWeight + process.pdfWeight
  + process.goodMuons + process.goodElectrons * process.goodJets
  * process.jpsiToMuMu# + process.jpsiToElEl

  + process.flatEventInfo
  * process.flatMuons + process.flatElectrons + process.flatJets
  + process.flatMETs + process.flatMETsUp + process.flatMETsDn
  + process.flatMETsRes + process.flatMETsResUp + process.flatMETsResDn
  + process.flatJpsiMuMu# + process.flatJpsiElEl
)

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
"""

process.output = cms.EndPath(process.out)
