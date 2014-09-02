import FWCore.ParameterSet.Config as cms

from KrAFT.RecoSelectorTools.leptonSelector_cfi import *
from KrAFT.RecoSelectorTools.jetSelector_cfi import *
from KrAFT.RecoSelectorTools.jpsiSelector_cfi import *
from KrAFT.GenericNtuple.flatEventInfo_cfi import *
from KrAFT.GenericNtuple.flatCands_cfi import *

analysisObjectSequence = cms.Sequence(
    goodMuons + goodElectrons
  + jetUncertainties * goodJets
  * jpsiToMuMu + jpsiToElEl

  + flatEventInfo
  * flatMuons + flatElectrons + flatJets
  + flatMETs + flatMETsUp + flatMETsDn
  + flatJpsiMuMu + flatJpsiElEl
)

delattr(flatJets.variables, "res")
delattr(flatJets.variables, "resUp")
delattr(flatJets.variables, "resDn")
