import FWCore.ParameterSet.Config as cms

jpsiToMuMu = cms.EDFilter("KVertexToMuMuProducer",
    track = cms.PSet(
        src = cms.InputTag("generalTracks"),
        minPt = cms.double(4.0),
        maxEta = cms.double(2.5),
        chi2 = cms.double(5.),
        nHit = cms.int32(6),
        signif = cms.double(-5),
        DCA = cms.double(1.),
    ),
    vertex = cms.PSet(
        chi2 = cms.double(7.),
        minLxy = cms.double(-40),
        maxLxy = cms.double(40),
        signif = cms.double(-5.0),
    ),
    muon = cms.PSet(
        src = cms.InputTag("muons"),
        dPtRel = cms.double(0.01),
        dR     = cms.double(0.01),
    ),
    pdgId = cms.uint32(443),
    leg1Id = cms.uint32(13),
    leg2Id = cms.uint32(13),
    rawMassMin = cms.double(2.75),
    rawMassMax = cms.double(3.45),
    massMin = cms.double(2.80),
    massMax = cms.double(3.40),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(100),
)

