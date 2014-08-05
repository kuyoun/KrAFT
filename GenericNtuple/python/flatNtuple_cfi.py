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
            fillPdg = cms.untracked.bool(False),
            vmaps = cms.untracked.vstring("isTight", "isLoose", "relIso", "dxy", "dz", ),
        ),
        electrons = cms.PSet(
            src = cms.InputTag("flatElectrons"),
            fillPdg = cms.untracked.bool(False),
            vmaps = cms.untracked.vstring("mva", "scEta", "relIso", "dxy", "dz", "chargeIDFull",),
        ),
        jets = cms.PSet(
            src = cms.InputTag("flatJets"),
            fillPdg = cms.untracked.bool(False),
            vmaps = cms.untracked.vstring("bTagCSV", "up", "dn", "res", "resUp", "resDn"),
        ),
        mets = cms.PSet(
            src = cms.InputTag("flatMETs"),
            fillEta = cms.untracked.bool(False),
            fillM   = cms.untracked.bool(False),
            fillQ   = cms.untracked.bool(False),
            fillPdg = cms.untracked.bool(False),
        ),
        metsUp = cms.PSet(
            src = cms.InputTag("flatMETsUp"),
            fillEta = cms.untracked.bool(False),
            fillM   = cms.untracked.bool(False),
            fillQ   = cms.untracked.bool(False),
            fillPdg = cms.untracked.bool(False),
        ),
        metsDn = cms.PSet(
            src = cms.InputTag("flatMETsDn"),
            fillEta = cms.untracked.bool(False),
            fillM   = cms.untracked.bool(False),
            fillQ   = cms.untracked.bool(False),
            fillPdg = cms.untracked.bool(False),
        ),
        metsRes = cms.PSet(
            src = cms.InputTag("flatMETsRes"),
            fillEta = cms.untracked.bool(False),
            fillM   = cms.untracked.bool(False),
            fillQ   = cms.untracked.bool(False),
            fillPdg = cms.untracked.bool(False),
        ),
        metsResUp = cms.PSet(
            src = cms.InputTag("flatMETsResUp"),
            fillEta = cms.untracked.bool(False),
            fillM   = cms.untracked.bool(False),
            fillQ   = cms.untracked.bool(False),
            fillPdg = cms.untracked.bool(False),
        ),
        metsResDn = cms.PSet(
            src = cms.InputTag("flatMETsResDn"),
            fillEta = cms.untracked.bool(False),
            fillM   = cms.untracked.bool(False),
            fillQ   = cms.untracked.bool(False),
            fillPdg = cms.untracked.bool(False),
        ),
        jpsiMuMu = cms.PSet(
            src = cms.InputTag("flatJpsiMuMu"),
            fillPdg = cms.untracked.bool(False),
            vmaps = cms.untracked.vstring("lxy", "l3D", "jetDR", "vProb"),
        ),
        jpsiElEl = cms.PSet(
            src = cms.InputTag("flatJpsiElEl"),
            fillPdg = cms.untracked.bool(False),
            vmaps = cms.untracked.vstring("lxy", "l3D", "jetDR", "vProb"),
        ),
        pseudoTopLepton = cms.PSet(
            src = cms.InputTag("flatPseudoTopLepton"),
        ),
        pseudoTopNu = cms.PSet(
            src = cms.InputTag("flatPseudoTopNu"),
            fillQ = cms.untracked.bool(False),
        ),
        pseudoTopJet = cms.PSet(
            src = cms.InputTag("flatPseudoTopJet"),
        ),
    ),
)
