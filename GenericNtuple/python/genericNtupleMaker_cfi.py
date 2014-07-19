import FWCore.ParameterSet.Config as cms

event = cms.EDAnalyzer("KGenericNtupleMaker",
    isMC = cms.bool(False),

    genEventInfo = cms.InputTag("generator"),
    genParticle = cms.InputTag("genParticles"),
    genJet = cms.InputTag("ak5GenJets"),
    recoToGenJetMap = cms.InputTag("recoToGenJetMap"),
    genJetToPartonsMap = cms.InputTag("genJetToPartonsMap"),

    pdfWeights = cms.InputTag("pdfWeight"),
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
    jetMET = cms.PSet(
        jet = cms.InputTag("goodJets"),
        met = cms.InputTag("patMETsPFlow"),
        unc = cms.InputTag("goodJets"),
        minNumber = cms.uint32(0),
        leptonDeltaR = cms.double(0.5),
        bTagType = cms.string("combinedSecondaryVertexBJetTags"),
    ),
    jpsi = cms.PSet(
        src = cms.InputTag("jpsiToMuMu"),
    ),
)
