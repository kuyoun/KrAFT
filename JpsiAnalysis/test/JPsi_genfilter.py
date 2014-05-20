#made by j.h.Goh.
import FWCore.ParameterSet.Config as cms
import sys, os
hostName = os.environ["HOSTNAME"]
process = cms.Process("FILTER")

process.genJpsi = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("abs(pdgId) == 443 && numberOfDaughters == 2 && abs(daughter(0).pdgId) == 13 && abs(daughter(1).pdgId) == 13"),
    minNumber = cms.uint32(1),
    filter = cms.bool(True),
)

process.filterJpsi = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("genJpsi"),
    minNumber = cms.uint32(1),
)

process.p = cms.Path(
    process.genJpsi * process.filterJpsi
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("jpsi.root"),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_*_*_FILTER",
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring('p') )
)
process.outPath = cms.EndPath(process.out)
process.source.fileNames = [
        'file:a.root',
    ]

