from KrAFT.Configuration.customise_cff import *
runOnMC = True
process = initialise(decayMode="dilepton", runOnMC=runOnMC)
process.load("KrAFT.GenericNtuple.flatCands_cfi")
addNtupleStep(process, runOnMC=runOnMC)
process.options.wantSummary = False
process.maxEvents.input = 100
process.source.fileNames = [
        '/store/relval//CMSSW_5_3_12_patch2/RelValProdTTbar/GEN-SIM-RECO/START53_LV2-v1/00000/5E865D62-AA2B-E311-AA04-002618943962.root',
        '/store/relval//CMSSW_5_3_12_patch2/RelValProdTTbar/GEN-SIM-RECO/START53_LV2-v1/00000/92EB24DF-C72B-E311-8AA2-00261894390E.root',
]
process.out = cms.OutputModule("PoolOutputModule",
                fileName = cms.untracked.string("out.root"),
                outputCommands = cms.untracked.vstring(['drop *', 'keep *_flat*_*_*'])
  )

process.maxEvents.input = 100

process.pMuMu += cms.Sequence(process.flatMuons+process.flatJpsiMuMu+process.flatJpsiElEl)
process.output = cms.EndPath(process.out)

