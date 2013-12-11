import FWCore.ParameterSet.Config as cms
runOnMC = False

from KCMSAnalyses.Configuration.customise_cff import *
process = initialise(decayMode="ElEl", runOnMC=runOnMC)
addNtupleStep(process, runOnMC=runOnMC)

process.source.fileNames = [
    '/store/user/jhgoh/MuEG/Run2012A-22Jan2013-v1-KCMSSkim20131027_1/4407ef23eed415918ae815f01ecb7627/skim_1_1_2W6.root',
]

