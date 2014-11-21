import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:MC/out.root',
        'file:MC/crab_0_141118_160859/res/out_1_1_tO4.root',  
        'file:MC/crab_0_141118_160859/res/out_3_1_uER.root',
        'file:MC/crab_0_141118_160859/res/out_2_1_6Y8.root',
    ),
)

process.load("KrAFT.GenericNtuple.flatNtuple_cfi")

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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

