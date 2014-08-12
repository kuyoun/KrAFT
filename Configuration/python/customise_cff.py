import FWCore.ParameterSet.Config as cms

def customisePAT(process, runOnMC, outputModules = []):
    ## Load PAT
    process.load("PhysicsTools.PatAlgos.patSequences_cff")

    ## Apply MVA
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

    ## Load trigger matching
    process.load("KrAFT.Configuration.hltFilters_cff")
    #from PhysicsTools.PatAlgos.tools.trigTools import *
    #switchOnTriggerMatchEmbedding(process, outputModule="")

    ## Apply PF2PAT
    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
    if runOnMC: jecLevels = ['L1FastJet','L2Relative','L3Absolute']
    else: jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

    #usePFBRECO(process,runPFBRECO=True,
    usePF2PAT(process, runPF2PAT=True,
              runOnMC=runOnMC, outputModules = outputModules, postfix="PFlow",
              jetAlgo="AK5", jetCorrections=("AK5PFchs", jecLevels),
              typeIMetCorrections=True)

    #put event counter at the end of the seqeuence
    process.nEventsPAT   = cms.EDProducer("EventCountProducer")
    process.patPF2PATSequencePFlow += process.nEventsPAT

    # top projections in PF2PAT:
    process.pfNoPileUpPFlow.enable = True
    process.pfNoMuonPFlow.enable = True
    process.pfNoElectronPFlow.enable = True
    process.pfNoTauPFlow.enable = False
    process.pfNoJetPFlow.enable = True

    # verbose flags for the PF2PAT modules
    process.pfNoMuonPFlow.verbose = False

    # Use non-isolated muons and electrons
    process.patMuonsPFlow.pfMuonSource = "pfMuonsPFlow"
    process.patElectronsPFlow.pfElectronSource = "pfElectronsPFlow"

    # And turn on delta-beta corrections while building pfIsolated*PFlow
    process.pfIsolatedMuonsPFlow.doDeltaBetaCorrection = True
    process.pfIsolatedElectronsPFlow.doDeltaBetaCorrection = True

    # Change DR cone size to 0.3
    process.pfIsolatedMuonsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PFlow'))
    process.pfIsolatedMuonsPFlow.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PFlow')
    process.pfIsolatedMuonsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PFlow'),
                                                                         cms.InputTag('muPFIsoValueGamma03PFlow'),)
    process.pfMuonsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PFlow') )
    process.pfMuonsPFlow.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PFlow')
    process.pfMuonsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PFlow'),
                                                                 cms.InputTag('muPFIsoValueGamma03PFlow'),)
    process.patMuonsPFlow.isolationValues.pfNeutralHadrons   = cms.InputTag('muPFIsoValueNeutral03PFlow')
    process.patMuonsPFlow.isolationValues.pfChargedAll       = cms.InputTag('muPFIsoValueChargedAll03PFlow')
    process.patMuonsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag('muPFIsoValuePU03PFlow')
    process.patMuonsPFlow.isolationValues.pfPhotons          = cms.InputTag('muPFIsoValueGamma03PFlow')
    process.patMuonsPFlow.isolationValues.pfChargedHadrons   = cms.InputTag('muPFIsoValueCharged03PFlow')

    process.pfIsolatedElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFlow'))
    process.pfIsolatedElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
    process.pfIsolatedElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPFlow'),
                                                                             cms.InputTag('elPFIsoValueGamma03PFIdPFlow'))
    process.pfElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFlow'))
    process.pfElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
    process.pfElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPFlow'),
                                                                     cms.InputTag('elPFIsoValueGamma03PFIdPFlow'))
    process.patElectronsPFlow.isolationValues.pfNeutralHadrons   = cms.InputTag('elPFIsoValueNeutral03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfChargedAll       = cms.InputTag('elPFIsoValueChargedAll03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfPhotons          = cms.InputTag('elPFIsoValueGamma03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfChargedHadrons   = cms.InputTag('elPFIsoValueCharged03PFIdPFlow')

