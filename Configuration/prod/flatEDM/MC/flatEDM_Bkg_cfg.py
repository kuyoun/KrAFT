import FWCore.ParameterSet.Config as cms

runOnMC = True

from KrAFT.Configuration.customise_cff import *
process = initialize(runOnMC)
process.options.allowUnscheduled = cms.untracked.bool(True)
customisePAT(process, runOnMC=runOnMC, outputModules=[])

process.source.fileNames = [
    '/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/BC91BA37-E2F2-E311-A317-0025905A612E.root',
]

process.load("KrAFT.Configuration.flatEDM_MC_cff")
process.load("KrAFT.Configuration.commonFilters_MC_cff")

process.CANDSEL = cms.Path(
    process.preFilterSequence
#  + process.patPF2PATSequencePFlow
  + process.analysisObjectSequence
  + process.partons
)

process.output = cms.EndPath(process.out)

