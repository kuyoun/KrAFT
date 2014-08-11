import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:out.root',
    ),
)

process.load("KrAFT.GenericNtuple.flatNtuple_cfi")
delattr(process.fEvent.cands, 'jpsiElEl')

process.passFEDM = cms.EDFilter("HLTHighLevel",
    eventSetupPathsKey = cms.string(''),
    TriggerResultsTag = cms.InputTag("TriggerResults","","KrAFT"),
    HLTPaths = cms.vstring(
        "CANDSEL",
    ),
    throw = cms.bool(False),
    andOr = cms.bool(True)
)

#process.p = cms.Path(process.passFEDM+process.fEvent)
process.p = cms.Path(process.fEvent)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

