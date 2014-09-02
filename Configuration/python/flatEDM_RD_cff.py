import FWCore.ParameterSet.Config as cms

from KrAFT.RecoSelectorTools.leptonSelector_cfi import *
from KrAFT.RecoSelectorTools.jetSelector_cfi import *
from KrAFT.RecoSelectorTools.jpsiSelector_cfi import *
from KrAFT.GenericNtuple.flatEventInfo_cfi import *
from KrAFT.GenericNtuple.flatCands_cfi import *

analysisObjectSequence = cms.Sequence(
    goodMuons + goodElectrons
  * goodJets
  * jpsiToMuMu + jpsiToElEl

  + flatEventInfo
  * flatMuons + flatElectrons + flatJets
  + flatMETs
  + flatJpsiMuMu + flatJpsiElEl
)

