import FWCore.ParameterSet.Config as cms
runOnMC = True

from KrAFT.Configuration.customise_cff import *
process = initialise(decayMode="dilepton", runOnMC=runOnMC)
addNtupleStep(process, runOnMC=runOnMC)

import os
hostName = os.environ['HOSTNAME']
if 'cern.ch' in hostName:
    process.source.fileNames = [
        '/store/relval//CMSSW_5_3_12_patch2/RelValProdTTbar/GEN-SIM-RECO/START53_LV2-v1/00000/5E865D62-AA2B-E311-AA04-002618943962.root',
        '/store/relval//CMSSW_5_3_12_patch2/RelValProdTTbar/GEN-SIM-RECO/START53_LV2-v1/00000/92EB24DF-C72B-E311-8AA2-00261894390E.root',
    ]
elif 'uos.ac.kr' in hostName
    process.source.fileNames = [
        '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10000/000560C1-FD97-E211-9F33-00304867924E.root'
]

process.maxEvents.input = 100

