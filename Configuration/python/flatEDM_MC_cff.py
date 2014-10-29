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

analysisObjectSequence = cms.Sequence(
    pileupWeight + pdfWeight
  + goodMuons + goodElectrons * goodJets
  * jpsiToMuMu + jpsiToElEl

  + flatEventInfo
  * flatMuons + flatElectrons + flatJets
  + flatMETs + flatMETsUp + flatMETsDn
  + flatMETsRes + flatMETsResUp + flatMETsResDn
  + flatJpsiMuMu + flatJpsiElEl
)

goodJets.isMC = True
