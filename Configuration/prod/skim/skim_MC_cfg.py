import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
#    '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10000/FEEEC639-4A98-E211-BE1C-002618943919.root'
     'file:mc_aod.root'
]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("skim.root"),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_selectedMuons_*_SKIM",
        "drop *_selectedElectrons_*_SKIM",
        "drop *_z*Cands_*_SKIM",
        "drop *_selectedJets_*_SKIM",
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pElEl', 'pMuMu', 'pMuEl', 'pSingleMu', 'pSingleElectron') 
    ),
)
process.outPath = cms.EndPath(process.out)

process.load("KrAFT.Configuration.commonFilters_cff")

process.pElEl = cms.Path(
    process.commonSequenceForMC
  + process.filterDoubleElectronSequence
)

process.pMuMu = cms.Path(
    process.commonSequenceForMC
  + process.filterDoubleMuSequence
)

process.pMuEl = cms.Path(
    process.commonSequenceForMC
  + process.filterMuEGSequence
)

process.pSingleMu = cms.Path(
    process.commonSequenceForMC
  + process.filterSingleMuSequence
)

process.pSingleElectron = cms.Path(
    process.commonSequenceForMC
  + process.filterSingleElectronSequence
)
