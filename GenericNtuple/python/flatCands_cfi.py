import FWCore.ParameterSet.Config as cms

flatMuonCands = cms.EDProducer("FlatCandProducer",
    type = cms.string("muon"),
    src = cms.InputTag("selectedMuonsPFlow"),
)

flatElectronCands = cms.EDProducer("FlatCandProducer",
    type = cms.string("electron"),
    src = cms.InputTag("selectedElectronsPFlow"),
)

flatJetCands = cms.EDProducer("FlatCandProducer",
    type = cms.string("jet"),
    src = cms.InputTag("selectedJetsPFlow"),
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
