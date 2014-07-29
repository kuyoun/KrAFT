#from PhysicsTools.PatAlgos.patTemplate_cfg import *
import FWCore.ParameterSet.Config as cms
process = cms.Process("KrAFT")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

#from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
jetAlgo="AK5"
isMC = True
"""
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=isMC, postfix=postfix, outputModules=[])

# top projections in PF2PAT:
getattr(process,"pfNoPileUpJME"+postfix).enable = True
getattr(process,"pfNoMuonJME"+postfix).enable = True
getattr(process,"pfNoElectronJME"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True
#getattr(process,"pfNoTau"+postfix).enable = True
getattr(process,"pfNoMuonJME"+postfix).verbose = False
getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = cms.bool(True)
"""
# load KrAFT
process.load("KrAFT.Configuration.commonFilters_cff")
process.load("KrAFT.RecoSelectorTools.leptonSelector_cfi")
process.load("KrAFT.RecoSelectorTools.jetSelector_cfi")
process.load("KrAFT.RecoSelectorTools.jpsiToMuMu_cfi")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'PLS170_V6AN1::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.load("KrAFT.GenericNtuple.flatCands_cfi")



process.goodMuons.src = "slimmedMuons"
process.goodMuons.rho = "fixedGridRhoAll"

process.goodElectrons.src = "slimmedElectrons"
process.goodElectrons.rho = "fixedGridRhoAll"

process.goodJets.jet = "slimmedJets"
process.goodJets.met = "slimmedMETs"
process.goodJets.isMC= isMC

process.jpsiToMuMu.muonSrc = "slimmedMuons"
process.jpsiToMuMu.electronSrc = "slimmedElectrons"

process.goodOfflinePrimaryVertices.src= cms.InputTag("offlineSlimmedPrimaryVertices")

process.options.wantSummary = False
process.source.fileNames = [
    '/store/relval/CMSSW_7_0_6_patch1/RelValJpsiMM_13/MINIAODSIM/PLS170_V6AN1-v1/00000/AE482444-9502-E411-B5E7-002618943810.root'
]
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               ## save only events passing the full path
                               #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               ## save PAT output; you need a '*' to unpack the list of commands
                               ## 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *')
                               )



process.out.outputCommands = ['drop *', 'keep *_flat*_*_*',]
process.out.fileName = "out.root"
process.outPath = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)
process.p = cms.Path(process.fEvent)
