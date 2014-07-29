import FWCore.ParameterSet.Config as cms

flatMuons = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("goodMuons"),
    variables = cms.PSet(
        #isTight = cms.string("userInt('isTight')"),
        #isLoose = cms.string("userInt('isLoose')"),
        relIso = cms.string("userIso(1)"),
        dxy = cms.string("dB"),
        dz = cms.string("userFloat('dz')"),
    ),
)

flatElectrons = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("goodElectrons"),
    variables = cms.PSet(
        #mva = cms.string("electronID('mvaTrigV0')"),  # FIXME for miniAOD 
        relIso = cms.string("userIso(2)"),
        scEta = cms.string("superCluster.eta"),
        dxy = cms.string("dB"),
        dz = cms.string("userFloat('dz')"),
        #chargeID = cms.string("isGsfCtfScPixChargeConsistent ? 3 : isGsfScPixChargeConsistent ? 2 : isGsfCtfChargeConsistent ? 1 : 0"),
    ),
)

flatJets = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("goodJets"),
    variables = cms.PSet(
        bTagCSV = cms.string("bDiscriminator('combinedSecondaryVertexBJetTags')"),
        up = cms.InputTag("goodJets", "up"),
        dn = cms.InputTag("goodJets", "dn"),
        res = cms.InputTag("goodJets", "res"),
        resUp = cms.InputTag("goodJets", "resUp"),
        resDn = cms.InputTag("goodJets", "resDn"),
    ),
)
flatJpsiMuMu = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("jpsiToMuMu"),
    variables = cms.PSet(
      lxy = cms.InputTag("jpsiToMuMu", "lxy"),
      l3D = cms.InputTag("jpsiToMuMu", "l3D"),
      jetDR = cms.InputTag("jpsiToMuMu", "jetDR"),
      vProb = cms.InputTag("jpsiToMuMu", "vProb"),
    ),
)


fEvent = cms.EDAnalyzer("FlatCandToNtupleMaker",
    cands = cms.PSet(
        muons = cms.PSet(
            src = cms.InputTag("flatMuons"),
            vmaps = cms.vstring(
                #"isTight", "isLoose",
                "relIso", "dxy", "dz",
            ),
        ),
        electrons = cms.PSet(
            src = cms.InputTag("flatElectrons"),
            vmaps = cms.vstring(
                #"mva",   # FIXME for miniAOD
                "scEta",
                "relIso", "dxy", "dz",
                #"chargeID",
            ),
        ),
        jets = cms.PSet(
            src = cms.InputTag("flatJets"),
            vmaps = cms.vstring(
                "bTagCSV", "up", "dn", "res", "resUp", "resDn",
            ),
        ),
        jpsis = cms.PSet(
            src = cms.InputTag("flatJpsiMuMu"),
            vmaps = cms.vstring(
                "lxy", "l3D", "jetDR", "vProb",
            ),
        ),
    ),
)
