import FWCore.ParameterSet.Config as cms

def initialise(runOnMC, decayMode, doOutModule=False):
    process = cms.Process("PAT")

    process.load("Configuration.StandardSequences.Services_cff")
    process.load("Configuration.Geometry.GeometryDB_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = 10000

    process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    from Configuration.AlCa.autoCond import autoCond
    if runOnMC: process.GlobalTag.globaltag = autoCond['startup']
    else: process.GlobalTag.globaltag = autoCond['com10']

    outputModuleForTriggerMatch = ""
    outputModules = []
    if doOutModule:
        from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
        process.out = cms.OutputModule("PoolOutputModule",
            fileName = cms.untracked.string("out.root"),
            outputCommands = cms.untracked.vstring(
                'drop *',
                'keep recoPFCandidates_particleFlow_*_*',
                *patEventContentNoCleaning
            )
        )
        process.outPath = cms.EndPath(process.out)

        outputModuleForTriggerMatch = "out"
        outputModules.append(process.out)

    ## Load PAT
    process.load("PhysicsTools.PatAlgos.patSequences_cff")

    ## Apply MVA
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

    ## Load trigger matching
    process.load("KCMSAnalyses.Configuration.hltFilters_cff")
    #from PhysicsTools.PatAlgos.tools.trigTools import *
    #switchOnTriggerMatchEmbedding(process, outputModule="")

    ## Apply PF2PAT
    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
    if runOnMC: jecLevels = ['L1FastJet','L2Relative','L3Absolute']
    else: jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

    #usePFBRECO(process,runPFBRECO=True,
    usePF2PAT(process, runPF2PAT=True,
              runOnMC=runOnMC, outputModules = outputModules, postfix="PF",
              jetAlgo="AK5", jetCorrections=("AK5PFchs", jecLevels),
              typeIMetCorrections=True)

    # top projections in PF2PAT:
    process.pfNoPileUpPF.enable = True
    process.pfNoMuonPF.enable = True
    process.pfNoElectronPF.enable = True
    process.pfNoTauPF.enable = False
    process.pfNoJetPF.enable = True

    # verbose flags for the PF2PAT modules
    process.pfNoMuonPF.verbose = False

    # Change DR cone size to 0.3
    process.pfIsolatedMuonsPF.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PF'))
    process.pfIsolatedMuonsPF.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PF')
    process.pfIsolatedMuonsPF.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PF'),
                                                                         cms.InputTag('muPFIsoValueGamma03PF'),)
    process.pfMuonsPF.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PF') )
    process.pfMuonsPF.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PF')
    process.pfMuonsPF.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PF'),
                                                                 cms.InputTag('muPFIsoValueGamma03PF'),)
    process.patMuonsPF.isolationValues.pfNeutralHadrons   = cms.InputTag('muPFIsoValueNeutral03PF')
    process.patMuonsPF.isolationValues.pfChargedAll       = cms.InputTag('muPFIsoValueChargedAll03PF')
    process.patMuonsPF.isolationValues.pfPUChargedHadrons = cms.InputTag('muPFIsoValuePU03PF')
    process.patMuonsPF.isolationValues.pfPhotons          = cms.InputTag('muPFIsoValueGamma03PF')
    process.patMuonsPF.isolationValues.pfChargedHadrons   = cms.InputTag('muPFIsoValueCharged03PF')

    process.pfIsolatedElectronsPF.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPF'))
    process.pfIsolatedElectronsPF.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPF')
    process.pfIsolatedElectronsPF.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPF'),
                                                                             cms.InputTag('elPFIsoValueGamma03PFIdPF'))
    process.pfElectronsPF.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPF'))
    process.pfElectronsPF.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPF')
    process.pfElectronsPF.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPF'),
                                                                     cms.InputTag('elPFIsoValueGamma03PFIdPF'))
    process.patElectronsPF.isolationValues.pfNeutralHadrons   = cms.InputTag('elPFIsoValueNeutral03PFIdPF')
    process.patElectronsPF.isolationValues.pfChargedAll       = cms.InputTag('elPFIsoValueChargedAll03PFIdPF')
    process.patElectronsPF.isolationValues.pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFIdPF')
    process.patElectronsPF.isolationValues.pfPhotons          = cms.InputTag('elPFIsoValueGamma03PFIdPF')
    process.patElectronsPF.isolationValues.pfChargedHadrons   = cms.InputTag('elPFIsoValueCharged03PFIdPF')

    ## Add common filters
    process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
    process.goodOfflinePrimaryVertices.filter = True

    process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
    process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
    if runOnMC: process.eventCleaning += process.eventCleaningMC
    else: process.eventCleaning += process.eventCleaningData

    # event counters
    process.nEventsTotal = cms.EDProducer("EventCountProducer")
    process.nEventsClean = cms.EDProducer("EventCountProducer")
    process.nEventsPAT   = cms.EDProducer("EventCountProducer")
    process.nEventsHLTElEl = cms.EDProducer("EventCountProducer")
    process.nEventsHLTMuMu = cms.EDProducer("EventCountProducer")
    process.nEventsHLTMuEl = cms.EDProducer("EventCountProducer")
    process.nEventsHLTMuJets = cms.EDProducer("EventCountProducer")
    process.nEventsHLTElJets = cms.EDProducer("EventCountProducer")

    process.commonFilterSequence = cms.Sequence(
        process.goodOfflinePrimaryVertices
      * process.eventCleaning
      + process.nEventsClean
    )

    process.patSequenceComplete = cms.Sequence(
    #  + process.patDefaultSequence
    #  + process.patPFBRECOSequencePF
        process.patPF2PATSequencePF
      + process.nEventsPAT
    )

    ## Defile paths
    if decayMode in ("all", "dilepton", "ElEl", "ee"):
        process.pElEl = cms.Path(
            process.nEventsTotal
          + process.commonFilterSequence
          + process.hltElEl + process.nEventsHLTElEl
          + process.patSequenceComplete
        )
    if decayMode in ("all", "dilepton", "MuMu", "mumu"):
        process.pMuMu = cms.Path(
            process.nEventsTotal
          + process.commonFilterSequence
          + process.hltMuMu + process.nEventsHLTMuMu
          + process.patSequenceComplete
        )
    if decayMode in ("all", "dilepton", "MuEl", "emu"):
        process.pMuEl = cms.Path(
            process.nEventsTotal
          + process.commonFilterSequence
          + process.hltMuEl + process.nEventsHLTMuEl
          + process.patSequenceComplete
        )
    if decayMode in ("all", "MuJets"):
        process.pMuJets = cms.Path(
            process.nEventsTotal
          + process.commonFilterSequence
          + process.hltMuJets + process.nEventsHLTMuJets
          + process.patSequenceComplete
        )
        if runOnMC: process.pMuJets.remove(process.hltMuJets)
    if decayMode in ("all", "ElJets"):
        process.pElJets = cms.Path(
            process.nEventsTotal
          + process.commonFilterSequence
          + process.hltElJets + process.nEventsHLTElJets
          + process.patSequenceComplete
        )
        if runOnMC: process.pElJets.remove(process.hltElJets)

    return process

def addNtupleStep(process, runOnMC):
    # Add ntuple production
    process.load("KCMSAnalyses.Configuration.ntuple_template_cff")
    process.goodJets.isMC = runOnMC

    for mode in ('ElEl', 'MuMu', 'MuEl', 'ElJets', 'MuJets'):
        if not hasattr(process, 'p'+mode): continue

        getattr(process, mode).isMC = runOnMC
        getattr(process, 'p'+mode).replace(
            process.nEventsPAT,
            process.nEventsPAT+(getattr(process, 'ntupleSequence'+mode))
        )

