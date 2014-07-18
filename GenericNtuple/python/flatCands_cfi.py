import FWCore.ParameterSet.Config as cms

flatMuonCands = cms.EDProducer("FlatCandProducer",
    type = cms.string("muon"),
    src = cms.InputTag("goodMuons"),
    vmaps = cms.VInputTag(),
)

flatElectronCands = cms.EDProducer("FlatCandProducer",
    type = cms.string("electron"),
    src = cms.InputTag("goodElectrons"),
    vmaps = cms.VInputTag(),
)

flatJetCands = cms.EDProducer("FlatCandProducer",
    type = cms.string("jet"),
    src = cms.InputTag("patJetsPFlow"),
    vmaps = cms.VInputTag(
        cms.InputTag("jetUnc", "up"),
        cms.InputTag("jetUnc", "dn"),
        cms.InputTag("jetUnc", "res"),
        cms.InputTag("jetUnc", "resUp"),
        cms.InputTag("jetUnc", "resDn"),
    ),
)

flatCandNtuple = cms.EDAnalyzer("FlatCandToNtupleMaker",
    srcs = cms.VPSet(
        cms.PSet(
            src = cms.InputTag("flatMuonCands"),
            vmaps = cms.VInputTag(
                cms.InputTag("flatMuonCands", "isTight"),
                cms.InputTag("flatMuonCands", "isLoose"),
                cms.InputTag("flatMuonCands", "relIsoDbeta03"),
                cms.InputTag("flatMuonCands", "relIsoDbeta04"),
            ),
        ),
        cms.PSet(
            src = cms.InputTag("flatElectronCands"),
            vmaps = cms.VInputTag(
                cms.InputTag("flatElectronCands", "mva"),
                cms.InputTag("flatElectronCands", "relIsoDbeta03"),
                cms.InputTag("flatElectronCands", "relIsoDbeta04"),
                cms.InputTag("flatElectronCands", "relIsoRho03"),
                cms.InputTag("flatElectronCands", "relIsoRho04"),
            ),
        ),
        cms.PSet(
            src = cms.InputTag("flatJetCands"),
            vmaps = cms.VInputTag(
                cms.InputTag("flatJetCands", "bTagCSV"),
            ),
        ),
    ),
)
