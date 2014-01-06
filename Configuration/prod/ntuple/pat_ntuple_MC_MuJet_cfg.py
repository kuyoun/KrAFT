import FWCore.ParameterSet.Config as cms
runOnMC = True

from KrAFT.Configuration.customise_cff import *
process = initialise(decayMode="MuJets", runOnMC=runOnMC)
addNtupleStep(process, runOnMC=runOnMC)

process.source.fileNames = [
    '/store/mc/Summer12_DR53X/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v1/20000/FED34040-AD84-E211-B269-782BCB536A50.root'
]

process.maxEvents.input = 100

