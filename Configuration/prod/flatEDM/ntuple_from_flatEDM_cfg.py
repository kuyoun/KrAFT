import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:out.root',
    ),
)

process.load("KrAFT.GenericNtuple.flatNtuple_cfi")
delattr(process.fEvent.cands, 'jpsiElEl')

process.p = cms.Path(process.fEvent)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

