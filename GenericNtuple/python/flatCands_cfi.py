import FWCore.ParameterSet.Config as cms

flatMuons = cms.EDProducer("FlatCandProducer",
    type = cms.string("muon"),
    src = cms.InputTag("goodMuons"),
    vmaps = cms.VInputTag(),
)

flatElectrons = cms.EDProducer("FlatCandProducer",
    type = cms.string("electron"),
    src = cms.InputTag("goodElectrons"),
    vmaps = cms.VInputTag(),
)

flatJets = cms.EDProducer("FlatCandProducer",
    type = cms.string("jet"),
    src = cms.InputTag("goodJets"),
    vmaps = cms.VInputTag(
        cms.InputTag("goodJets", "up"),
        cms.InputTag("goodJets", "dn"),
        cms.InputTag("goodJets", "res"),
        cms.InputTag("goodJets", "resUp"),
        cms.InputTag("goodJets", "resDn"),
    ),
)

fEvent = cms.EDAnalyzer("FlatCandToNtupleMaker",
    srcs = cms.VPSet(
        cms.PSet(
            src = cms.InputTag("flatMuons"),
            vmaps = cms.VInputTag(
                cms.InputTag("flatMuons", "isTight"),
                cms.InputTag("flatMuons", "isLoose"),
                cms.InputTag("flatMuons", "relIso"),
                cms.InputTag("flatMuons", "dxy"),
            ),
        ),
        cms.PSet(
            src = cms.InputTag("flatElectrons"),
            vmaps = cms.VInputTag(
                cms.InputTag("flatElectrons", "mva"),
                cms.InputTag("flatElectrons", "relIso"),
                cms.InputTag("flatElectrons", "scEta"),
                cms.InputTag("flatElectrons", "dxy"),
                cms.InputTag("flatElectrons", "chargeID"),
            ),
        ),
        cms.PSet(
            src = cms.InputTag("flatJets"),
            vmaps = cms.VInputTag(
                cms.InputTag("flatJets", "bTagCSV"),
                cms.InputTag("flatJets", "up"),
                cms.InputTag("flatJets", "dn"),
                cms.InputTag("flatJets", "res"),
                cms.InputTag("flatJets", "resUp"),
                cms.InputTag("flatJets", "resDn"),
            ),
        ),
    ),
)
