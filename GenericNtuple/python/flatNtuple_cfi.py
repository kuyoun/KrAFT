import FWCore.ParameterSet.Config as cms

p4Set = cms.untracked.PSet(
    pt = cms.string("pt"),
    eta = cms.string("eta"),
    phi = cms.string("phi"),
    m = cms.string("mass"),
    q = cms.string("charge"),
    pdgId = cms.string("pdgId"),
)

fEvent = cms.EDAnalyzer("FlatCandToNtupleMaker",
    failureMode = cms.untracked.string("keep"), # choose one among keep/skip/error
    eventCounters = cms.vstring("nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        nVertex = cms.PSet(src = cms.InputTag("flatEventInfo", "pvN")),
    ),
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
            exprs = p4Set.clone(),
            vmaps = cms.untracked.vstring("isTight", "isLoose", "relIso", "dxy", "dz", ),
        ),
        electrons = cms.PSet(
            src = cms.InputTag("flatElectrons"),
            exprs = p4Set.clone(),
            vmaps = cms.untracked.vstring("mva", "scEta", "relIso", "dxy", "dz", "chargeIDFull","conversionVeto","isPF"),
        ),
        jets = cms.PSet(
            src = cms.InputTag("flatJets"),
            exprs = p4Set.clone(),
            vmaps = cms.untracked.vstring("bTagCSV", "up", "dn", "res", "resUp", "resDn"),
        ),
        mets = cms.PSet(
            src = cms.InputTag("flatMETs"),
            exprs = cms.untracked.PSet(pt = cms.string("pt"), phi = cms.string("phi")),
        ),
        metsUp = cms.PSet(
            src = cms.InputTag("flatMETsUp"),
            exprs = cms.untracked.PSet(pt = cms.string("pt"), phi = cms.string("phi")),
        ),
        metsDn = cms.PSet(
            src = cms.InputTag("flatMETsDn"),
            exprs = cms.untracked.PSet(pt = cms.string("pt"), phi = cms.string("phi")),
        ),
        metsRes = cms.PSet(
            src = cms.InputTag("flatMETsRes"),
            exprs = cms.untracked.PSet(pt = cms.string("pt"), phi = cms.string("phi")),
        ),
        metsResUp = cms.PSet(
            src = cms.InputTag("flatMETsResUp"),
            exprs = cms.untracked.PSet(pt = cms.string("pt"), phi = cms.string("phi")),
        ),
        metsResDn = cms.PSet(
            src = cms.InputTag("flatMETsResDn"),
        ),
        jpsiMuMu = cms.PSet(
            src = cms.InputTag("flatJpsiMuMu"),
            exprs = p4Set.clone(),
            vmaps = cms.untracked.vstring("lxy", "l3D", "jetDR", "vProb"),
        ),
        jpsiElEl = cms.PSet(
            src = cms.InputTag("flatJpsiElEl"),
            exprs = p4Set.clone(),
            vmaps = cms.untracked.vstring("lxy", "l3D", "jetDR", "vProb"),
        ),
        pseudoTopLepton = cms.PSet(
            src = cms.InputTag("flatPseudoTopLepton"),
            exprs = p4Set.clone(),
        ),
        pseudoTopNu = cms.PSet(
            src = cms.InputTag("flatPseudoTopNu"),
            exprs = p4Set.clone(),
        ),
        pseudoTopJet = cms.PSet(
            src = cms.InputTag("flatPseudoTopJet"),
            exprs = p4Set.clone(),
        ),
    ),
)

delattr(fEvent.cands.muons.exprs, 'pdgId')
delattr(fEvent.cands.electrons.exprs, 'pdgId')
delattr(fEvent.cands.jpsiMuMu.exprs, 'pdgId')
delattr(fEvent.cands.jpsiMuMu.exprs, 'q')
delattr(fEvent.cands.jpsiElEl.exprs, 'pdgId')
delattr(fEvent.cands.jpsiElEl.exprs, 'q')
delattr(fEvent.cands.pseudoTopNu.exprs, 'm')
delattr(fEvent.cands.pseudoTopNu.exprs, 'q')
delattr(fEvent.cands.pseudoTopJet.exprs, 'q')
