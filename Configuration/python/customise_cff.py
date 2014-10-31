import FWCore.ParameterSet.Config as cms

def customisePAT(process, runOnMC, outputModules = []):
    ## Load PAT
    process.load("PhysicsTools.PatAlgos.patSequences_cff")

    ## Apply MVA
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    #process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    #process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

    ### Customise for 7XY
    process.mvaTrigV0.electronTag = 'gedGsfElectrons'
    process.mvaTrigNoIPV0.electronTag = 'gedGsfElectrons'
    process.mvaNonTrigV0.electronTag = 'gedGsfElectrons'

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
              jetAlgo="AK4", jetCorrections=("AK4PFchs", jecLevels),
              typeIMetCorrections=True)

    #put event counter at the end of the seqeuence
    process.nEventsPAT   = cms.EDProducer("EventCountProducer")
    #process.patPF2PATSequencePFlow += process.nEventsPAT

    # In order to avoid over-subtracting high pT tracks from jets for 2012.
    process.pfPileUpPFlow.checkClosestZVertex = False
    process.pfPileUpPFlow.Vertices = cms.InputTag("goodOfflinePrimaryVertices")

    # top projections in PF2PAT: we are turning off top projection
    process.pfNoPileUpPFlow.enable = True
    process.pfNoMuonPFlow.enable = False #True
    process.pfNoElectronPFlow.enable = False #True
    process.pfNoTauPFlow.enable = False
    process.pfNoJetPFlow.enable = False #True

    # verbose flags for the PF2PAT modules
    process.pfNoMuonPFlow.verbose = False

    # Use non-isolated muons and electrons
    process.patMuonsPFlow.pfMuonSource = "pfMuonsPFlow"
    process.patElectronsPFlow.pfElectronSource = "pfElectronsPFlow"

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

def initialize(runOnMC, processName="KrAFT"):
    process = cms.Process(processName)

    process.load("Configuration.StandardSequences.Services_cff")
    process.load("Configuration.Geometry.GeometryDB_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
    from Configuration.AlCa.autoCond import autoCond
    if runOnMC: process.GlobalTag.globaltag = autoCond['startup']
    else: process.GlobalTag.globaltag = autoCond['com10']

    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring()
    )

    process.out = cms.OutputModule("PoolOutputModule",
        compressionLevel = cms.untracked.int32(4),
        compressionAlgorithm = cms.untracked.string('LZMA'),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        fileName = cms.untracked.string("out.root"),
        outputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_TriggerResults_*_HLT',
            'keep *_TriggerResults_*_%s' % processName,
            'keep edmMergeableCounter_*_*_*',
            'keep *_flat*_*_*',
        ),
        SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring("CANDSEL"),
        ),
    )

    if runOnMC:
        process.out.outputCommands.extend([
            'keep *_partons_*_*',
            #'keep *_pseudoTop_*_*', # recoGenJets/GenParticles from pseudoTop producer
            'keep *_pileupWeight_*_*',
            'keep *_pdfWeight_*_*',
        ])

    return process

