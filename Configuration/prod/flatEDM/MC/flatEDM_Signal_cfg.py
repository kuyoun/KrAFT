import FWCore.ParameterSet.Config as cms

runOnMC = True

from KrAFT.Configuration.customise_cff import *
process = initialize(runOnMC)
process.options.allowUnscheduled = cms.untracked.bool(True)
customisePAT(process, runOnMC=runOnMC, outputModules=[])

process.source.fileNames = [
	'/store/relval/CMSSW_7_2_0/RelValTTbar_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v2/00000/3EEC24CE-9D59-E411-9662-002618943875.root',
	'/store/relval/CMSSW_7_2_0/RelValTTbar_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v2/00000/BC30DFD3-9D59-E411-87FE-0025905A612E.root',
	'/store/relval/CMSSW_7_2_0/RelValTTbar_13/GEN-SIM-RECO/PHYS14_25_V1_Phys14-v2/00000/EC548A8B-9659-E411-B063-0025905A607E.root',    

#	'file:/afs/cern.ch/work/y/youn/work/BC91BA37-E2F2-E311-A317-0025905A612E.root',

#    '/store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/BC91BA37-E2F2-E311-A317-0025905A612E.root',
#	'/store/relval/CMSSW_7_2_0/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v2/00000/D61BBE8C-8359-E411-9FF7-00261894386A.root',
]

process.load("KrAFT.Configuration.flatEDM_MC_cff")
process.load("KrAFT.Configuration.commonFilters_MC_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.GEN = cms.Path(
    process.pseudoTop
  + process.partons
  * process.flatPseudoTopLepton + process.flatPseudoTopNu + process.flatPseudoTopJet + process.flatPseudoTopWdau
)

process.CANDSEL = cms.Path(
    process.preFilterSequence
  #  process.patPF2PATSequencePFlow
  + process.analysisObjectSequence
)

process.out.SelectEvents.SelectEvents.append("GEN")
process.output = cms.EndPath(process.out)

