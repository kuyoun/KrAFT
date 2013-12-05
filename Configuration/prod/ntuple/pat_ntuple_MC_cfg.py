import FWCore.ParameterSet.Config as cms
process = cms.Process("PAT")
runOnMC = True

from KCMSAnalyses.Configuration.customise_cff import customise
customise(process, runOnMC=runOnMC)

process.source.fileNames = [
    '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10000/000560C1-FD97-E211-9F33-00304867924E.root'
]

process.maxEvents.input = 100

