import FWCore.ParameterSet.Config as cms

fEvent = cms.EDAnalyzer("FlatCandToNtupleMaker",
    weight = cms.PSet(
        puWeight   = cms.PSet(src = cms.InputTag("pileupWeight")),
        puWeightUp = cms.PSet(src = cms.InputTag("pileupWeight", "up")),
        puWeightDn = cms.PSet(src = cms.InputTag("pileupWeight", "dn")),
    ),
    vWeight = cms.PSet(
        pdfWeight = cms.PSet(src = cms.InputTag("pdfWeight")),
    ),
    cands = cms.PSet(
        muons = cms.PSet(
            src = cms.InputTag("flatMuons"),
            vmaps = cms.vstring("isTight", "isLoose", "relIso", "dxy", "dz", ),
        ),
        electrons = cms.PSet(
            src = cms.InputTag("flatElectrons"),
            vmaps = cms.vstring("mva", "scEta", "relIso", "dxy", "dz", "chargeIDFull",),
        ),
        jets = cms.PSet(
            src = cms.InputTag("flatJets"),
            vmaps = cms.vstring("bTagCSV", "up", "dn", "res", "resUp", "resDn"),
        ),
        jpsiMuMu = cms.PSet(
            src = cms.InputTag("flatJpsiMuMu"),
            vmaps = cms.vstring("lxy", "l3D", "jetDR", "vProb"),
        ),
        jpsiElEl = cms.PSet(
            src = cms.InputTag("flatJpsiElEl"),
            vmaps = cms.vstring("lxy", "l3D", "jetDR", "vProb"),
        ),
        pseudoTopLepton = cms.PSet(
            src = cms.InputTag("flatPseudoTopLepton"),
            vmaps = cms.vstring(),
        ),
        pseudoTopNu = cms.PSet(
            src = cms.InputTag("flatPseudoTopNu"),
            vmaps = cms.vstring(),
        ),
        pseudoTopJet = cms.PSet(
            src = cms.InputTag("flatPseudoTopJet"),
            vmaps = cms.vstring(),
        ),
    ),
)
