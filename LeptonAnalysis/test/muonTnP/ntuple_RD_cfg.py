import FWCore.ParameterSet.Config as cms

runOnMC = False

tightMuonCut = "&&".join([
    "isGlobalMuon && isPFMuon",
    "globalTrack().normalizedChi2()<10.0",
    "globalTrack().hitPattern().numberOfValidMuonHits()>0",
    "numberOfMatchedStations>1",
    "dB < 0.2",
    "innerTrack().hitPattern().numberOfValidPixelHits() > 0",
    "track().hitPattern().trackerLayersWithMeasurement() > 5",
])
isoMuonCut = "(chargedHadronIso+max(0.,neutralHadronIso+photonIso-0.50*puChargedHadronIso)) < 0.12*pt"

process = cms.Process("PAT")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#process.GlobalTag.globaltag = "START53_V15::All"
process.GlobalTag.globaltag = "FT_P_V42_AN3::All"
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.suppressError = cms.untracked.vstring('patTriggerFull')
process.MessageLogger.suppressWarning = cms.untracked.vstring('patTriggerFull')
process.MessageLogger.suppressInfo = cms.untracked.vstring('patTriggerFull')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#        "/store/relval/CMSSW_5_3_6/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START53_V14-v1/0003/2668E6A9-E92C-E211-BA7D-003048D37666.root",
"/store/data/Run2012D/SingleMu/AOD/22Jan2013-v1/10000/FE618190-93A7-E211-B50C-90E6BA19A215.root"
    )
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.throw = False
process.hltHighLevel.HLTPaths = ["HLT_IsoMu24_v*", "HLT_IsoMu24_eta2p1_v*",]

outputModules = []

## Load PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Apply MVA
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

## Load trigger matching
## Trigger matching with PAT
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1

## Apply PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
if runOnMC: jecLevels = ['L1FastJet','L2Relative','L3Absolute']
else: jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

#usePFBRECO(process,runPFBRECO=True,
usePF2PAT(process, runPF2PAT=True,
          runOnMC=runOnMC, outputModules = outputModules, postfix="PFlow",
          jetAlgo="AK5", jetCorrections=("AK5PFchs", jecLevels),
          typeIMetCorrections=True)

# top projections in PF2PAT:
process.pfNoPileUpPFlow.enable = True
process.pfNoMuonPFlow.enable = True
process.pfNoElectronPFlow.enable = True
process.pfNoTauPFlow.enable = False
process.pfNoJetPFlow.enable = True

# verbose flags for the PF2PAT modules
process.pfNoMuonPFlow.verbose = False

"""
# Change DR cone size to 0.3
process.pfIsolatedMuonsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PFlow'))
process.pfIsolatedMuonsPFlow.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PFlow')
process.pfIsolatedMuonsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PFlow'),
                                                                     cms.InputTag('muPFIsoValueGamma03PFlow'),)
process.pfMuonsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PFlow') )
process.pfMuonsPFlow.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PFlow')
process.pfMuonsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PFlow'),
                                                             cms.InputTag('muPFIsoValueGamma03PFlow'),)
process.patMuonsPFlow.isolationValues.pfNeutralHadrons   = cms.InputTag('muPFIsoValueNeutral03PFlow')
process.patMuonsPFlow.isolationValues.pfChargedAll       = cms.InputTag('muPFIsoValueChargedAll03PFlow')
process.patMuonsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag('muPFIsoValuePU03PFlow')
process.patMuonsPFlow.isolationValues.pfPhotons          = cms.InputTag('muPFIsoValueGamma03PFlow')
process.patMuonsPFlow.isolationValues.pfChargedHadrons   = cms.InputTag('muPFIsoValueCharged03PFlow')
"""

process.pfIsolatedElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFlow'))
process.pfIsolatedElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
process.pfIsolatedElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPFlow'),
                                                                         cms.InputTag('elPFIsoValueGamma03PFIdPFlow'))
process.pfElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFlow'))
process.pfElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
process.pfElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPFlow'),
                                                                 cms.InputTag('elPFIsoValueGamma03PFIdPFlow'))
process.patElectronsPFlow.isolationValues.pfNeutralHadrons   = cms.InputTag('elPFIsoValueNeutral03PFIdPFlow')
process.patElectronsPFlow.isolationValues.pfChargedAll       = cms.InputTag('elPFIsoValueChargedAll03PFIdPFlow')
process.patElectronsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
process.patElectronsPFlow.isolationValues.pfPhotons          = cms.InputTag('elPFIsoValueGamma03PFIdPFlow')
process.patElectronsPFlow.isolationValues.pfChargedHadrons   = cms.InputTag('elPFIsoValueCharged03PFIdPFlow')

#from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
#from PhysicsTools.PatAlgos.tools.trigTools import *
#switchOnTriggerMatchEmbedding(process, outputModule="")

## Add common filters
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
process.goodOfflinePrimaryVertices.filter = True

process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
if runOnMC: process.eventCleaning += process.eventCleaningMC
else: process.eventCleaning += process.eventCleaningData

# event counters
process.nEventsTotal = cms.EDProducer("EventCountProducer")
process.nEventsClean = cms.EDProducer("EventCountProducer")
process.nEventsPAT   = cms.EDProducer("EventCountProducer")

process.commonFilterSequence = cms.Sequence(
    process.goodOfflinePrimaryVertices
  * process.eventCleaning
  + process.nEventsClean
)

process.patSequenceComplete = cms.Sequence(
#  + process.patDefaultSequence
#  + process.patPFBRECOSequencePFlow
    process.patPF2PATSequencePFlow
  + process.nEventsPAT
)

## Tag and probes
#process.load("KrAFT.RecoSelectorTools.leptonSelector_cfi")
#process.goodMuons.src = "patMuonsWithTrigger"
#process.goodMuons.cut = tightMuonCut
#process.goodMuons.coneSize = 0.4
#process.goodMuons.maxDz = 0.5
#process.goodMuons.maxDxy = 0.2

process.goodMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(tightMuonCut + " && pt > 20 && abs(eta) < 2.5"),
    filter = cms.bool(False),
)

process.isoMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("goodMuons"),
    cut = cms.string(isoMuonCut),
    #cut = cms.string("userIsolation('User2Iso') < 0.12"),
    filter = cms.bool(False),
)

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("isoMuons"),
    cut = cms.string(
        "pt > 25"
        "&& (!triggerObjectMatchesByPath('HLT_IsoMu24_v*').empty()"
        "  || !triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*').empty() )"
    ),
    filter = cms.bool(True),
)

process.tagMuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("tagMuons"),
    minNumber = cms.uint32(1),
)

process.goodTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string("pt>20 && abs(eta)<2.5"),
)

process.trackProbes = cms.EDProducer("ConcreteChargedCandidateProducer",
    src  = cms.InputTag("goodTracks"),
    particleType = cms.string("mu+"),
)

## Matchings
process.matchTightMuons = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("trackProbes"),
    matched = cms.InputTag("goodMuons"),
    algorithm = cms.string("byDirectComparison"),
    srcTrack = cms.string("tracker"),
    srcState = cms.string("atVertex"),
    matchedTrack = cms.string("tracker"),
    matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(0.01),
    maxDeltaLocalPos = cms.double(0.01),
    maxDeltaPtRel    = cms.double(0.01),
    sortBy           = cms.string("deltaR"),
)

process.matchIsoMuons = process.matchTightMuons.clone(
    matched = cms.InputTag("isoMuons")
)

## Passing probe candidates
process.passTightMuons = cms.EDProducer("MatchedCandidateSelector",
    src = cms.InputTag("trackProbes"),
    match = cms.InputTag("matchTightMuons"),
)

process.passIsoMuons = cms.EDProducer("MatchedCandidateSelector",
    src = cms.InputTag("trackProbes"),
    match = cms.InputTag("matchIsoMuons"),
)

process.muonSelectionSequence = cms.Sequence(
    process.patMuonsWithTriggerSequence
  * process.goodMuons * process.isoMuons 
  * process.tagMuons * process.tagMuonFilter
  + process.goodTracks * process.trackProbes
  + process.matchTightMuons * process.passTightMuons
  * process.matchIsoMuons * process.passIsoMuons
)

## Build Tag-Probe pair and tree
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ trackProbes@-"),
    cut   = cms.string("60 < mass < 120"),
)

process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        pt  = cms.string("pt"),
        abseta    = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        tightMuons   = cms.InputTag("passTightMuons"),
        isoMuons     = cms.InputTag("passIsoMuons"  ),
    ),
    isMC = cms.bool(False),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("tnpTree.root"),
)

process.load("KrAFT.RecoSelectorTools.jetSelector_cfi")

## Defile paths
process.pMu = cms.Path(
    process.nEventsTotal
  + process.commonFilterSequence
  + process.patSequenceComplete
  + process.hltHighLevel
  + process.muonSelectionSequence
#  * process.goodElectrons + process.goodJets
#  * process.tpPairs * process.muonEffs
)

