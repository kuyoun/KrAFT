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
process.goodMuons.rho = "ak5PFJets:rho"
process.goodElectrons.rho = "ak5PFJets:rho"

process.options.wantSummary = False
process.maxEvents.input = 100
process.source.fileNames = [
    '/store/relval/CMSSW_7_0_6/RelValTTbarLepton_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/4EA445A2-3BFA-E311-B066-0026189438BA.root',
    '/store/relval/CMSSW_7_0_6/RelValTTbarLepton_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/F0AADC82-39FA-E311-8050-002354EF3BD0.root',
    '/store/relval/CMSSW_7_0_6/RelValTTbarLepton_13/GEN-SIM-RECO/PLS170_V7AN1-v1/00000/FC321902-35FA-E311-B1FE-002618943957.root',

]
process.out.outputCommands = ['drop *', 'keep *_goodMuons_*_*',]
process.out.fileName = "out.root"
#process.outPath = cms.EndPath()

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string("ntuple.root"),
#)
#process.p = cms.Path(process.fEvent)
