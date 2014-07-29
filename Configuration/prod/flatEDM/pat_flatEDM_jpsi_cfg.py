from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
jetAlgo="AK5"
isMC = True
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=isMC, postfix=postfix, outputModules=[])

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
process.load("KrAFT.RecoSelectorTools.jpsiToMuMu_cfi")

process.load("KrAFT.GenericNtuple.flatCands_cfi")
process.goodMuons.rho = "fixedGridRhoFastjetAll"
process.goodElectrons.rho = "fixedGridRhoFastjetAll"
process.goodJets.isMC = cms.bool(isMC)

process.options.wantSummary = False
process.maxEvents.input = -1
process.source.fileNames = [
    '/store/relval/CMSSW_7_0_6_patch1/RelValJpsiMM_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/DC160579-8802-E411-AEB6-002618943978.root',
    '/store/relval/CMSSW_7_0_6_patch1/RelValJpsiMM_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/FA2F148F-8B02-E411-B8CB-0025905A6090.root',
]
process.out.outputCommands = ['drop *','keep *_flat*_*_*']
process.out.fileName = "out.root"
process.outPath = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)
process.p = cms.Path(process.fEvent)
