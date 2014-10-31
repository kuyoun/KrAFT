import FWCore.ParameterSet.Config as cms

from KrAFT.GeneratorTools.pileupWeight_cff import *
from KrAFT.GeneratorTools.pdfWeight_cff import *
from KrAFT.GeneratorTools.pseudoTop_cfi import *
from KrAFT.GeneratorTools.partons_cff import *
from KrAFT.RecoSelectorTools.leptonSelector_cfi import *
from KrAFT.RecoSelectorTools.jetSelector_cfi import *
from KrAFT.RecoSelectorTools.jpsiSelector_cfi import *
from KrAFT.GenericNtuple.flatEventInfo_cfi import *
from KrAFT.GenericNtuple.flatCands_cfi import *

for x in jetUncertainties.jecLevels:
    goodJets.jes.extend([cms.InputTag("jetUncertainties", x+"Up"),
                         cms.InputTag("jetUncertainties", x+"Dn")])
goodJets.jes.extend([
    cms.InputTag("jetUncertainties", "res"),
    cms.InputTag("jetUncertainties", "resUp"),
    cms.InputTag("jetUncertainties", "resDn")
])

flatJets.variables.up = cms.InputTag("goodJets", "TotalUp")
flatJets.variables.dn = cms.InputTag("goodJets", "TotalDn")
flatJets.variables.res = cms.InputTag("goodJets", "res")
flatJets.variables.resUp = cms.InputTag("goodJets", "resUp")
flatJets.variables.resDn = cms.InputTag("goodJets", "resDn")

flatMETsUp = flatMETs.clone(src = cms.InputTag("jetUncertainties", "TotalUp"))
flatMETsDn  = flatMETs.clone(src = cms.InputTag("jetUncertainties", "TotalDn"))
flatMETsRes = flatMETs.clone(src = cms.InputTag("jetUncertainties", "res"))
flatMETsResUp = flatMETs.clone(src = cms.InputTag("jetUncertainties", "resDn"))
flatMETsResDn = flatMETs.clone(src = cms.InputTag("jetUncertainties", "resUp"))

analysisObjectSequence = cms.Sequence(
    pileupWeight + pdfWeight
  + goodMuons + goodElectrons
  + jetUncertainties * goodJets
  * jpsiToMuMu + jpsiToElEl

  + flatEventInfo
  * flatMuons + flatElectrons + flatJets
  + flatMETs + flatMETsUp + flatMETsDn
  + flatMETsRes + flatMETsResUp + flatMETsResDn
  + flatJpsiMuMu + flatJpsiElEl
)


