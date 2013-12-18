import FWCore.ParameterSet.Config as cms

kshortVertex = cms.EDFilter("VertexCandProducer",
    track = cms.PSet(
        src = cms.InputTag("generalTracks"),
        minPt = cms.double(4.0),
        maxEta = cms.double(2.5),
        chi2 = cms.double(5.),
        nHit = cms.int32(6),
        signif = cms.double(0.5),
        DCA = cms.double(1.),
    ),
    vertex = cms.PSet(
        chi2 = cms.double(7.),
        minLxy = cms.double(0.0),
        maxLxy = cms.double(40),
        signif = cms.double(5.0),
    ),
    pdgId = cms.uint32(310),
    leg1Id = cms.uint32(211),
    leg2Id = cms.uint32(211),
    rawMassMin = cms.double(0.35),
    rawMassMax = cms.double(0.65),
    massMin = cms.double(0.40),
    massMax = cms.double(0.60),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(100),
)

lambdaVertex = kshortVertex.clone(
    pdgId = cms.uint32(3122),
    leg1Id = cms.uint32(2212),
    leg2Id = cms.uint32(211),
    rawMassMin = cms.double(1.04),
    rawMassMax = cms.double(1.24),
    massMin = cms.double(1.06),
    massMax = cms.double(1.22),
)

phiVertex = kshortVertex.clone(
    pdgId = cms.uint32(333),
    leg1Id = cms.uint32(321),
    leg2Id = cms.uint32(321),
    rawMassMin = cms.double(0.96),
    rawMassMax = cms.double(1.08),
    massMin = cms.double(0.98),
    massMax = cms.double(1.06),
)
phiVertex.track.signif  = -5
phiVertex.vertex.signif = -5
phiVertex.vertex.minLxy = -40
phiVertex.vertex.maxLxy =  40

jpsiVertex = kshortVertex.clone(
    pdgId = cms.uint32(443),
    leg1Id = cms.uint32(13),
    leg2Id = cms.uint32(13),
    rawMassMin = cms.double(2.75),
    rawMassMax = cms.double(3.45),
    massMin = cms.double(2.80),
    massMax = cms.double(3.40),
)
jpsiVertex.track.signif  = -5
jpsiVertex.vertex.signif = -5
jpsiVertex.vertex.minLxy = -40
jpsiVertex.vertex.maxLxy =  40
