import FWCore.ParameterSet.Config as cms

def customise(process, runOnMC):
    process.load("Configuration.StandardSequences.Services_cff")
    process.load("Configuration.Geometry.GeometryDB_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    #process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = 10000

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    from Configuration.AlCa.autoCond import autoCond
    if runOnMC:
        process.GlobalTag.globaltag = autoCond['startup']
    else:
        process.GlobalTag.globaltag = autoCond['com10']
    process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

    from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("out.root"),
        outputCommands = cms.untracked.vstring(
            'drop *',
            'keep recoPFCandidates_particleFlow_*_*',
            *patEventContentNoCleaning
        )
    )
    #process.outPath = cms.EndPath(process.out)

    process.load("PhysicsTools.PatAlgos.patSequences_cff")
    ## Apply MVA
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT

    postfix = "PFlow"
    #usePFBRECO(process,runPFBRECO=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix)
    usePF2PAT(process, runOnMC=runOnMC, postfix=postfix)

    # top projections in PF2PAT:
    getattr(process,"pfNoPileUp"+postfix).enable = True
    getattr(process,"pfNoMuon"+postfix).enable = True
    getattr(process,"pfNoElectron"+postfix).enable = True
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoJet"+postfix).enable = True

    # verbose flags for the PF2PAT modules
    getattr(process,"pfNoMuon"+postfix).verbose = False

    process.p = cms.Path(
    #    #getattr(process,"patPFBRECOSequence"+postfix)
        getattr(process,"patPF2PATSequence"+postfix)
    #  + process.patDefaultSequence
    )

    # Add ntuple production
    process.load("KCMSAnalyses.Configuration.ntuple_template_cff")
    process.goodJets.isMC = runOnMC
    process.event.isMC = runOnMC
    process.p += process.ntupleSequence


