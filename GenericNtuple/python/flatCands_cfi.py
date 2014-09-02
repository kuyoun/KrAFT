import FWCore.ParameterSet.Config as cms

flatDummy = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag(""),
    variables = cms.PSet(),
    selections = cms.PSet(),
)

flatMuons = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("goodMuons"),
    variables = cms.PSet(
        relIso = cms.string("userIso(1)"),
        dxy = cms.string("dB"),
        dz = cms.string("userFloat('dz')"),
    ),
    selections = cms.PSet(
        isTight = cms.string("isPFMuon && isGlobalMuon && globalTrack.normalizedChi2<10.0 && innerTrack.hitPattern.numberOfValidHits > 0 && track.hitPattern.trackerLayersWithMeasurement > 5 && numberOfMatchedStations > 1 && globalTrack.hitPattern.numberOfValidMuonHits > 0"),
        isLoose = cms.string("isPFMuon && (isTrackerMuon || isGlobalMuon)"),
    ),
)

flatElectrons = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("goodElectrons"),
    variables = cms.PSet(
        mva = cms.string("electronID('mvaTrigV0')"),
        relIso = cms.string("userIso(2)"),
        scEta = cms.string("superCluster.eta"),
        dxy = cms.string("dB"),
        dz = cms.string("userFloat('dz')"),
    ),
    selections = cms.PSet(
        chargeIDFull = cms.string("isGsfCtfScPixChargeConsistent"),
        #isGsfScPixChargeConsistent isGsfCtfChargeConsistent),
    ),
)

flatJets = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("goodJets"),
    variables = cms.PSet(
        bTagCSV = cms.string("bDiscriminator('combinedSecondaryVertexBJetTags')"),
        up = cms.InputTag("goodJets", "TotalUp"),
        dn = cms.InputTag("goodJets", "TotalDn"),
        res = cms.InputTag("goodJets", "res"),
        resUp = cms.InputTag("goodJets", "resUp"),
        resDn = cms.InputTag("goodJets", "resDn"),
    ),
    selections = cms.PSet(),
)

flatMETs    = flatDummy.clone(src = cms.InputTag("patMETsPFlow"))
flatMETsUp  = flatDummy.clone(src = cms.InputTag("jetUncertainties", "TotalUp"))
flatMETsDn  = flatDummy.clone(src = cms.InputTag("jetUncertainties", "TotalDn"))
flatMETsRes = flatDummy.clone(src = cms.InputTag("jetUncertainties", "res"))
flatMETsResUp = flatDummy.clone(src = cms.InputTag("jetUncertainties", "resDn"))
flatMETsResDn = flatDummy.clone(src = cms.InputTag("jetUncertainties", "resUp"))

flatJpsiMuMu = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("jpsiToMuMu"),
    variables = cms.PSet(
      lxy = cms.InputTag("jpsiToMuMu", "lxy"),
      l3D = cms.InputTag("jpsiToMuMu", "l3D"),
      jetDR = cms.InputTag("jpsiToMuMu", "jetDR"),
      vProb = cms.InputTag("jpsiToMuMu", "vProb"),
    ),
    selections = cms.PSet(),
)
flatJpsiElEl = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("jpsiToElEl"),
    variables = cms.PSet(
      lxy = cms.InputTag("jpsiToElEl", "lxy"),
      l3D = cms.InputTag("jpsiToElEl", "l3D"),
      jetDR = cms.InputTag("jpsiToElEl", "jetDR"),
      vProb = cms.InputTag("jpsiToElEl", "vProb"),
    ),
    selections = cms.PSet(),
)

flatPseudoTopLepton = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("pseudoTop", "leptons"),
    variables = cms.PSet(),
    selections = cms.PSet(),
)

flatPseudoTopNu = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("pseudoTop", "neutrinos"),
    variables = cms.PSet(),
    selections = cms.PSet(),
)

flatPseudoTopJet = cms.EDProducer("FlatCandProducer",
    src = cms.InputTag("pseudoTop", "jets"),
    variables = cms.PSet(),
    selections = cms.PSet(),
)

