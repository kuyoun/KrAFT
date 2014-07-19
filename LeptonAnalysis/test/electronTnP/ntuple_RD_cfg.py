import FWCore.ParameterSet.Config as cms

runOnMC = False

process = cms.Process("PAT")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

if runOnMC:
    process.GlobalTag.globaltag = "START53_V27::All"
else:
    process.GlobalTag.globaltag = "FT53_V21A_AN6::All"
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.suppressError = cms.untracked.vstring('patTriggerFull')
process.MessageLogger.suppressWarning = cms.untracked.vstring('patTriggerFull')
process.MessageLogger.suppressInfo = cms.untracked.vstring('patTriggerFull')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#        "/store/relval/CMSSW_5_3_6/RelValZmumuJets_Pt_20_300/GEN-SIM-RECO/PU_START53_V14-v1/0003/2668E6A9-E92C-E211-BA7D-003048D37666.root",
        "/store/data/Run2012D/SingleElectron/AOD/22Jan2013-v1/10000/FEE2508B-1193-E211-9327-002590593920.root",
    )
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.throw = False
process.hltHighLevel.HLTPaths = ["HLT_Ele27_WP80_v*", ]

## Load PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Apply MVA
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

## Apply PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
if runOnMC: jecLevels = ['L1FastJet','L2Relative','L3Absolute']
else: jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

#usePFBRECO(process,runPFBRECO=True,
usePF2PAT(process, runPF2PAT=True,
          runOnMC=runOnMC, outputModules = [], postfix="PFlow",
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

## Load trigger matching
## Trigger matching with PAT
from PhysicsTools.PatAlgos.tools.trigTools import *
process.patElectronsWithTrigger = cms.EDProducer("PATTriggerMatchElectronEmbedder",
    src     = cms.InputTag("patElectronsPFlow"),
    matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT')),#, cms.InputTag('eleIdTriggerMatchHLT'))
)
process.eleTriggerMatchHLT = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    src = cms.InputTag( "patElectronsPFlow" ),
    matched = cms.InputTag( "patTrigger" ),
    maxDPtRel = cms.double(999),
    maxDeltaR = cms.double(0.3),
    resolveAmbiguities = cms.bool( True ),
    matchedCuts = cms.string(
        'filter("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter") || '
        'filter("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter") || '
        'filter("hltEle27WP80TrackIsoFilter")'
    ),
    resolveByMatchQuality = cms.bool( True )
)

#switchOnTrigger(process,sequence='patDefaultSequencePFlow',hltProcess = '*',outputModule="")
#from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = [ 'eleTriggerMatchHLT' ], 
                              sequence='patDefaultSequencePFlow', hltProcess = '*', outputModule="")

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

## Tag and probes
#process.load("KrAFT.RecoSelectorTools.leptonSelector_cfi")

process.patElectronsWithIso = cms.EDFilter("KElectronSelector",
    rho = cms.InputTag("kt6PFJets", "rho"),
    src = cms.InputTag("patElectronsPFlow"),
    precut = cms.string("pt >= 4 && abs(eta) < 2.5 && passConversionVeto && gsfTrack.isNonnull"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    cut = cms.string(""),
    coneSize = cms.double(0.3),
    maxDz = cms.double(0.5),
    maxDxy = cms.double(0.2),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
    filter = cms.bool(False),
)

process.patElectronsWithIsoForTag = process.patElectronsWithIso.clone(
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(
        "pt > 30 && ( !triggerObjectMatchesByPath('HLT_Ele27_WP80_v*').empty() )"
    ),
)

process.tagElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithIsoForTag"),
    cut = cms.string("userIsolation('User3Iso') < 0.15"),
    filter = cms.bool(True),
)

process.tagElectronFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("tagElectrons"),
    minNumber = cms.uint32(1),
)

## Passing probe candidates
process.leptonSelectionSequence = cms.Sequence(
    process.eleTriggerMatchHLT * process.patElectronsWithTrigger
  * process.patElectronsWithIsoForTag * process.tagElectrons * process.tagElectronFilter
  * process.patElectronsWithIso
)

## Build Tag-Probe pair and tree
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagElectrons@+ patElectronsWithIso@-"),
    cut   = cms.string("60 < mass < 120"),
)

process.tp = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairs"),
    jets = cms.InputTag("goodJets"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        scEta = cms.string("superCluster.eta"),
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        pt  = cms.string("pt"),
        abseta    = cms.string("abs(eta)"),
        mva = cms.string("electronID('mvaNonTrigV0')"),
        nHit = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits"),
        isoRho = cms.string("userIsolation('User3Iso')"),
        isoDbeta = cms.string("userIsolation('User2Iso')"),
    ),
    flags = cms.PSet(
        pf = cms.string("isPF"),
        baseCut = cms.string("isPF && dB < 0.02"),

        mva00 = cms.string("electronID('mvaNonTrigV0') > 0.0"),
        mva05 = cms.string("electronID('mvaNonTrigV0') > 0.5"),
        mva09 = cms.string("electronID('mvaNonTrigV0') > 0.9"),

        iso15 = cms.string("userIsolation('User3Iso') < 0.15"),
        qCheck = cms.string("isGsfCtfScPixChargeConsistent"),
    ),
    isMC = cms.bool(runOnMC),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("tnpTree.root"),
)

process.load("KrAFT.RecoSelectorTools.jetSelector_cfi")

## Defile paths
process.patSequenceComplete = cms.Sequence(
#  + process.patDefaultSequence
#  + process.patPFBRECOSequencePFlow
    process.patPF2PATSequencePFlow
  + process.nEventsPAT
)

process.p = cms.Path(
    process.nEventsTotal
  + process.commonFilterSequence
  + process.patSequenceComplete
  + process.hltHighLevel
  + process.leptonSelectionSequence
  + process.goodJets
  * process.tpPairs * process.tp
)

