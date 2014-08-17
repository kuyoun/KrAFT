import FWCore.ParameterSet.Config as cms

from TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi import *
goodOfflinePrimaryVertices.filter = True

from TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff import *
trackingFailureFilter.VertexSource = 'goodOfflinePrimaryVertices'
eventCleaning += eventCleaningMC

nEventsTotal = cms.EDProducer("EventCountProducer")
nEventsClean = cms.EDProducer("EventCountProducer")

preFilterSequence = cms.Sequence(
    nEventsTotal
  + goodOfflinePrimaryVertices + eventCleaning
  + nEventsClean
)

"""

selectedMuons = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("muons"),
    cut = cms.string(
        "abs(eta) < 2.6 && pt > 17"
        " && isPFMuon && (isGlobalMuon || isTrackerMuon)"),
    filter = cms.bool(True),
)

selectedElectrons = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string(
      "abs(eta) < 2.6 && pt > 17"
      " && gsfTrack.isNonnull"
      #" && passConversionVeto "
      " && gsfTrack.trackerExpectedHitsInner.numberOfHits<=0"),
    filter = cms.bool(True),
)

zMuMuCands = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedMuons@+ selectedMuons@-"),
    checkCharge = cms.bool(False),
    cut = cms.string("10 < mass"),
)

zElElCands = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedElectrons@+ selectedElectrons@-"),
    checkCharge = cms.bool(False),
    cut = cms.string("10 < mass"),
)

zMuElCands = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("selectedMuons@+ selectedElectrons@-"),
    checkCharge = cms.bool(False),
    cut = cms.string("10 < mass"),
)

nZMuMuCands = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zMuMuCands"),
    minNumber = cms.uint32(1),
)
nZElElCands = nZMuMuCands.clone(src = cms.InputTag("zElElCands"))
nZMuElCands = nZMuMuCands.clone(src = cms.InputTag("zMuElCands"))

selectedSingleMuon = selectedMuons.clone(
    cut = cms.string(
        "abs(eta) < 2.6 && pt > 24 && isPFMuon && isGlobalMuon"
        " && globalTrack.normalizedChi2 < 10"
        " && globalTrack.hitPattern.numberOfValidMuonHits > 0"
        " && numberOfMatchedStations > 1"
        " && innerTrack.hitPattern.numberOfValidPixelHits > 0"
        " && track.hitPattern.trackerLayersWithMeasurement > 5"
    ),
)

selectedSingleElectron = selectedElectrons.clone(
    cut = cms.string(
      "abs(eta) < 2.6 && pt > 27"
      " && gsfTrack.isNonnull"
      #" && passConversionVeto "
      " && gsfTrack.trackerExpectedHitsInner.numberOfHits<=0"
      " && !(1.4442 < abs(superCluster.eta) && abs(superCluster.eta) < 1.5660)"
    ),
)

nMuonFilterSingleLepton = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("selectedSingleMuon"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
)

nElectronFilterSingleLepton = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("selectedSingleElectron"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1),
)

selectedJets = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("ak5PFJets"),
    cut = cms.string(
        "abs(eta) < 2.6 && pt > 20"
        " && numberOfDaughters > 1"
        " && neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99"
        " && (abs(eta) >= 2.4 || chargedEmEnergyFraction < 0.99)"
        " && (abs(eta) >= 2.4 || chargedHadronEnergyFraction > 0.)"
        " && (abs(eta) >= 2.4 || chargedMultiplicity > 0)"),
)

nJetFilterSingleLepton = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("selectedJets"),
    minNumber = cms.uint32(3),
)

nEventsTotal = cms.EDProducer("EventCountProducer")

commonSequenceForData = cms.Sequence(
    goodOfflinePrimaryVertices
  + noscraping
)

commonSequenceForMC = cms.Sequence(
    nEventsTotal
  + goodOfflinePrimaryVertices
)

filterDoubleMuSequence = cms.Sequence(
    selectedMuons * zMuMuCands * nZMuMuCands
)

filterDoubleElectronSequence = cms.Sequence(
    selectedElectrons * zElElCands * nZElElCands
)

filterMuEGSequence = cms.Sequence(
    selectedMuons * selectedElectrons * zMuElCands * nZMuElCands
)

filterSingleMuSequence = cms.Sequence(
    selectedSingleMuon * nMuonFilterSingleLepton
  + selectedJets * nJetFilterSingleLepton
)

filterSingleElectronSequence = cms.Sequence(
    selectedSingleElectron * nElectronFilterSingleLepton
  + selectedJets * nJetFilterSingleLepton
)
"""
