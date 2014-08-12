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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
from Configuration.AlCa.autoCond import autoCond
if runOnMC: process.GlobalTag.globaltag = autoCond['startup']
else: process.GlobalTag.globaltag = autoCond['com10']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_5_3_12_patch2/RelValProdTTbar/GEN-SIM-RECO/START53_LV2-v1/00000/5E865D62-AA2B-E311-AA04-002618943962.root',
        '/store/relval/CMSSW_5_3_12_patch2/RelValProdTTbar/GEN-SIM-RECO/START53_LV2-v1/00000/92EB24DF-C72B-E311-8AA2-00261894390E.root',
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
        'keep *_pseudoTop_*_*',
        'keep *_pileupWeight_*_*',
        'keep *_pdfWeight_*_*',
        'keep *_flat*_*_*',
        'keep *_TriggerResults_*_%s' % process.process,
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("GEN", "CANDSEL"),
    ),
)

from KrAFT.Configuration.customise_cff import *
customisePAT(process, runOnMC=runOnMC, outputModules=[])

process.load("KrAFT.Configuration.flatEDM_MC_cff")
process.load("KrAFT.Configuration.commonFilters_MC_cff")

process.partons = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop *",
        "drop pt <= 0",
        "keep status = 3", # For the pythia
    ),
)

process.GEN = cms.Path(
    process.pseudoTop
  + process.partons
  * process.flatPseudoTopLepton + process.flatPseudoTopNu + process.flatPseudoTopJet
)

process.CANDSEL = cms.Path(
    process.preFilterSequence
  + process.patPF2PATSequencePFlow
  + process.analysisObjectSequence
)

process.output = cms.EndPath(process.out)

