import FWCore.ParameterSet.Config as cms
runOnMC = False

from KrAFT.Configuration.customise_cff import *
process = initialise(decayMode="ElEl", runOnMC=runOnMC)
addNtupleStep(process, runOnMC=runOnMC)
process.event.muon.src=cms.InputTag("goodMuonsForJpsi")
process.event.muon.src=cms.InputTag("goodElectronsForJpsi")

import os
hostName = os.environ['HOSTNAME']
if 'cern.ch' in hostName:
    process.source.fileNames = [
        '/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10000/FEED5E9F-6A8F-E211-91C4-00261894391C.root',
    ]
elif 'uos.ac.kr' in hostName:
    process.source.fileNames = [
        '/store/user/jhgoh/MuEG/Run2012A-22Jan2013-v1-KCMSSkim20131027_1/4407ef23eed415918ae815f01ecb7627/skim_1_1_2W6.root',
    ]

process.maxEvents.input = -1
process.pElEl.remove(process.ElEl)
"""
index = process.pElEl.index(process.jpsiToMuMu)
process.pElEl.insert(index+1, process.jpsiToElEl)
process.pElEl.insert(index+2, process.ElEl)
"""

