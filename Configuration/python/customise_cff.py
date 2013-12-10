import FWCore.ParameterSet.Config as cms

def customise(process, runOnMC, decayMode="Dilepton"):
    process.load("Configuration.StandardSequences.Services_cff")
    process.load("Configuration.Geometry.GeometryDB_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
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

    #from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
    #process.out = cms.OutputModule("PoolOutputModule",
    #    fileName = cms.untracked.string("out.root"),
    #    outputCommands = cms.untracked.vstring(
    #        'drop *',
    #        'keep recoPFCandidates_particleFlow_*_*',
    #        *patEventContentNoCleaning
    #    )
    #)
    #process.outPath = cms.EndPath(process.out)

    process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
    process.goodOfflinePrimaryVertices.filter = True

    process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
    process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
    if runOnMC:
        process.eventCleaning += process.eventCleaningMC
    else:
        process.eventCleaning += process.eventCleaningData

    process.load("PhysicsTools.PatAlgos.patSequences_cff")
    ## Apply MVA
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

    #from TopQuarkAnalysis.Configuration.patRefSel_triggerMatching_cfi import patMuonTriggerMatch
    #from PhysicsTools.PatAlgos.tools.trigTools import *
    #switchOnTriggerMatchEmbedding(process, outputModule="")

    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT

    postfix = "PFlow"
    if runOnMC: jecLevels = ['L1FastJet','L2Relative','L3Absolute']
    else: jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']
    #usePFBRECO(process,runPFBRECO=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix)
    usePF2PAT(process, jetAlgo="AK5", jetCorrections=("AK5PFchs", jecLevels), runOnMC=runOnMC, typeIMetCorrections=True, outputModules = [], postfix=postfix)

    # top projections in PF2PAT:
    getattr(process,"pfNoPileUp"+postfix).enable = True
    getattr(process,"pfNoMuon"+postfix).enable = True
    getattr(process,"pfNoElectron"+postfix).enable = True
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoJet"+postfix).enable = True

    # verbose flags for the PF2PAT modules
    getattr(process,"pfNoMuon"+postfix).verbose = False

    # Change DR cone size to 0.3
    getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix ))
    getattr( process, 'pfIsolatedMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                              , cms.InputTag( 'muPFIsoValueGamma03' + postfix ))
    getattr( process, 'pfMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix ) )
    getattr( process, 'pfMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    getattr( process, 'pfMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                      , cms.InputTag( 'muPFIsoValueGamma03' + postfix ))
    getattr( process, 'patMuons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'muPFIsoValueChargedAll03' + postfix )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'muPFIsoValueGamma03' + postfix )
    getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'muPFIsoValueCharged03' + postfix )

    getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix ))
    getattr( process, 'pfIsolatedElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                                  , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix ))
    getattr( process, 'pfElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix ))
    getattr( process, 'pfElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    getattr( process, 'pfElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                          , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix ))
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
    getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )

    # event counters
    process.nEventsTotal = cms.EDProducer("EventCountProducer")
    process.nEventsClean = cms.EDProducer("EventCountProducer")
    process.nEventsPAT   = cms.EDProducer("EventCountProducer")

    process.p = cms.Path(
        process.nEventsTotal
      + process.goodOfflinePrimaryVertices * process.eventCleaning
      + process.nEventsClean
    #    #getattr(process,"patPFBRECOSequence"+postfix)
      + getattr(process,"patPF2PATSequence"+postfix)
    #  + process.patDefaultSequence
      + process.nEventsPAT
    )

    # Add ntuple production
    process.load("KCMSAnalyses.Configuration.ntuple_template_cff")
    process.goodJets.isMC = runOnMC
    process.event.isMC = runOnMC
    if decayMode == 'MuJet': process.p += process.ntupleSequenceMuJet
    elif decayMode == 'ElJet': process.p += process.ntupleSequenceElJet
    else: process.p += process.ntupleSequenceDilepton ## Default as dilepton

