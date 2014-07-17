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

