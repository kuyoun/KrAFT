import FWCore.ParameterSet.Config as cms

jpsiToElEl = cms.EDFilter("KJpsiElElProducer",
    src = cms.InputTag("patElectronsPFlow"),
    track = cms.PSet(
        minPt = cms.double(1.0),
        maxEta = cms.double(2.5),
        chi2 = cms.double(5.),
        nHit = cms.int32(6),
        signif = cms.double(-5),
        DCA = cms.double(1.),
    ),
    vertex = cms.PSet(
        chi2 = cms.double(7.),
        minLxy = cms.double(-100),
        maxLxy = cms.double(100),
        signif = cms.double(-5.0),
    ),
    pdgId = cms.uint32(443),
    leg1Id = cms.uint32(11),
    leg2Id = cms.uint32(11),
    rawMassMin = cms.double(2),
    rawMassMax = cms.double(4),
    massMin = cms.double(2.80),
    massMax = cms.double(3.40),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(100),
)

