import FWCore.ParameterSet.Config as cms
runOnMC = False

from KrAFT.Configuration.customise_cff import *
process = initialise(decayMode="MuEl", runOnMC=runOnMC)
addNtupleStep(process, runOnMC=runOnMC)

import os
hostName = os.environ['HOSTNAME']
if 'cern.ch' in hostName:
    process.source.fileNames = [
        '/store/data/Run2012D/MuEG/AOD/22Jan2013-v1/20000/FE5E0B33-0D82-E211-AB78-00237DA15C7C.root',
    ]
elif 'uos.ac.kr' in hostName:
    process.source.fileNames = [
        '/store/user/jhgoh/MuEG/Run2012A-22Jan2013-v1-KCMSSkim20131027_1/4407ef23eed415918ae815f01ecb7627/skim_1_1_2W6.root',
    ]

process.maxEvents.input = 100

