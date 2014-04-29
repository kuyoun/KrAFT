#!/usr/bin/env python

def matching( muon, genparticle_list) :
	for genparticle in genparticle_list:
		pass 

from ROOT import *

f = TFile("/pnfs/user/geonmo/kraft/src/test_plot/ntuple.root")

tree = f.Get("MuEl/event")
hJpsiMuPt = TH1F("hJpsiMuPt", "hJpsiMuPt", 100, 0, 100)
hJpsiMuPtEE = TH1F("hJpsiMuPtEE", "hJpsiMuPtEE", 100, 0, 100)
hJpsiMuPtME = TH1F("hJpsiMuPtME", "hJpsiMuPtME", 100, 0, 100)
hJpsiMuPtMM = TH1F("hJpsiMuPtMM", "hJpsiMuPtMM", 100, 0, 100)
for event in tree:
    pdgIds = event.genParticles_pdgId
    jpsis_pt1 = event.jpsis_pt1
    if len(jpsis_pt1) == 0: continue
    hJpsiMuPt.Fill(event.jpsis_pt1[0])

    nGenElectron = 0
    nGenMuon = 0
    for pdgId in pdgIds:
        if ( abs(pdgId) == 13 ): nGenMuon += 1
        elif ( abs(pdgId) == 11 ): nGenElectron += 1
    if nGenElectron == 2:
        hJpsiMuPtEE.Fill(event.jpsis_pt1[0])
    elif nGenMuon == 2:
        hJpsiMuPtMM.Fill(event.jpsis_pt1[0])
    elif nGenMuon == 1 and nGenElectron == 1:
        hJpsiMuPtME.Fill(event.jpsis_pt1[0])
c = TCanvas("c", "c", 500, 500)
hJpsiMuPt.Draw()
cEE = TCanvas("cEE", "cEE", 500, 500)
hJpsiMuPtEE.Draw()
cMM = TCanvas("cMM", "cMM", 500, 500)
hJpsiMuPtMM.Draw()
cME = TCanvas("cME", "cME", 500, 500)
hJpsiMuPtME.Draw()

