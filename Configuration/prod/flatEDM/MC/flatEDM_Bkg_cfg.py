import FWCore.ParameterSet.Config as cms

runOnMC = True

from KrAFT.Configuration.customise_cff import *
process = initialize(runOnMC)
customisePAT(process, runOnMC=runOnMC, outputModules=[])

process.source.fileNames = [
    '/store/relval/CMSSW_5_3_12_patch2/RelValZMM/GEN-SIM-RECO/START53_LV2-v1/00000/28D552C4-A82B-E311-952C-002590596486.root',
    '/store/relval/CMSSW_5_3_12_patch2/RelValZMM/GEN-SIM-RECO/START53_LV2-v1/00000/EA8973F4-A92B-E311-A0A6-00261894396D.root',
]

process.load("KrAFT.Configuration.flatEDM_MC_cff")
process.load("KrAFT.Configuration.commonFilters_MC_cff")

process.CANDSEL = cms.Path(
    process.preFilterSequence
  + process.patPF2PATSequencePFlow
  + process.analysisObjectSequence
  + process.partons
)

process.output = cms.EndPath(process.out)

