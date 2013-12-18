import FWCore.ParameterSet.Config as cms

event = cms.EDAnalyzer("KFlatTreeMaker",
    isMC = cms.bool(False),

    genEventInfo = cms.InputTag("generator"),
    genParticle = cms.InputTag("genParticles"),
    genJet = cms.InputTag("ak5GenJets"),
    recoToGenJetMap = cms.InputTag("recoToGenJetMap"),
    genJetToPartonsMap = cms.InputTag("genJetToPartonsMap"),

    puWeight = cms.InputTag("pileupWeight"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    eventCounters = cms.vstring(),

    electron = cms.PSet(
        src = cms.InputTag("goodElectrons"),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    muon = cms.PSet(
        src = cms.InputTag("goodMuons"),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    jet = cms.PSet(
        src = cms.InputTag("goodJets"),
        leptonDeltaR = cms.double(0.5),
        bTagType = cms.string("combinedSecondaryVertexBJetTags"),
    ),
    met = cms.PSet(
        src = cms.InputTag("goodJets"),
    ),
    jpsi = cms.PSet(
        src = cms.InputTag("jpsiToMuMu"),
    ),
)
