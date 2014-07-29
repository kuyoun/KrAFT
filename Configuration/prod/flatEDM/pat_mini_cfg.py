from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
jetAlgo="AK5"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix, outputModules=[])

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning

# top projections in PF2PAT:
getattr(process,"pfNoPileUpJME"+postfix).enable = True
getattr(process,"pfNoMuonJME"+postfix).enable = True
getattr(process,"pfNoElectronJME"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True
#getattr(process,"pfNoTau"+postfix).enable = True
getattr(process,"pfNoMuonJME"+postfix).verbose = False
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = cms.bool(True)

# load KrAFT
process.load("KrAFT.Configuration.commonFilters_cff")
process.load("KrAFT.RecoSelectorTools.leptonSelector_cfi")
process.load("KrAFT.RecoSelectorTools.jetSelector_cfi")
process.load("KrAFT.GenericNtuple.flatCands_cfi")
process.goodMuons.rho = "fixedGridRhoAll"
process.goodElectrons.rho = "fixedGridRhoAll"

process.options.wantSummary = False
process.maxEvents.input = -1
process.source.fileNames = [
    '/store/relval/CMSSW_7_0_6_patch3/RelValTTbar_13/MINIAODSIM/PUpmx50ns_PLS170_V7AN1-v2/00000/329775B2-BB12-E411-A99D-0025905A609E.root',
    '/store/relval/CMSSW_7_0_6_patch3/RelValTTbar_13/MINIAODSIM/PUpmx50ns_PLS170_V7AN1-v2/00000/90974C6E-BB12-E411-A947-002618943834.root'
]
process.out.outputCommands = ['drop *', 'keep *_goodMuons_*_*',]
process.out.fileName = "out.root"
#process.outPath = cms.EndPath()

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string("ntuple.root"),
#)
#process.p = cms.Path(process.fEvent)
