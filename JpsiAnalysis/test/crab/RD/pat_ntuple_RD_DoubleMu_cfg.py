import FWCore.ParameterSet.Config as cms
runOnMC = False

from KrAFT.Configuration.customise_cff import *
process = initialise(decayMode="MuMu", runOnMC=runOnMC)
addNtupleStep(process, runOnMC=runOnMC)

import os
hostName = os.environ['HOSTNAME']
if 'cern.ch' in hostName:
    process.source.fileNames = [
        '/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/20000/FEF2A4B5-2082-E211-BC3E-E0CB4EA0A904.root',
    ]
elif 'uos.ac.kr' in hostName:
    process.source.fileNames = [
        '/store/user/jhgoh/MuEG/Run2012A-22Jan2013-v1-KCMSSkim20131027_1/4407ef23eed415918ae815f01ecb7627/skim_1_1_2W6.root',
    ]

process.maxEvents.input = -1

process.pMuMu.remove(process.MuMu)
index = process.pMuMu.index(process.jpsiToMuMu)
process.pMuMu.insert(index+1, process.jpsiToElEl)
process.pMuMu.insert(index+2, process.MuMu)
