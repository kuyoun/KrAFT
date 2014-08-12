import FWCore.ParameterSet.Config as cms

runOnMC = False

from KrAFT.Configuration.customise_cff import *
process = initialize(runOnMC)
customisePAT(process, runOnMC=runOnMC, outputModules=[])

process.source.fileNames = [
    '/store/data/Run2012D/SingleMu/AOD/22Jan2013-v1/10000/FE618190-93A7-E211-B50C-90E6BA19A215.root',
]

process.load("KrAFT.Configuration.flatEDM_RD_cff")
process.load("KrAFT.Configuration.commonFilters_RD_cff")

process.CANDSEL = cms.Path(
    process.preFilterSequence
  + process.patPF2PATSequencePFlow
  + process.analysisObjectSequence
)

process.output = cms.EndPath(process.out)

