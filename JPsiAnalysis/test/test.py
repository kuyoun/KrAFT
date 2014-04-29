#!/usr/bin/env python

def matching( muon, genparticle_list) :
	for genparticle in genparticle_list:
		pass 

from ROOT import *

a = TFile("/pnfs/user/geonmo/kraft/src/test_plot/ntuple.root")

dir = a.Get("MuMu")
tree = dir.Get("event")
h1 = TH1F("mulep1","mulep1",100,0,100)
#h2 = TH1F("ddmulep2","mulep2",100,0,100)
h3 = TH1F("ellep1","ellep1",100,0,100)
#h4 = TH1F("ellep2","ellep2",100,0,100)
for count, event in enumerate(tree) :
	print "Event : %d"%(count)
	pdglist = event.genParticles_pdgId
	muon   = False
	a_muon = False

	electron   = False
	positron   = False
  
	for gen_pdg in pdglist :
		if ( gen_pdg is -13 ) :
			muon = True
			print "found muon"
		if ( gen_pdg is 13 ) :
			a_muon = True
			print "found a_muon"
		if ( gen_pdg is -11) :
			electron = True
			print "found electron"
		if ( gen_pdg is 11 ) :	
			positron = True
			print "found positron"
	if ( muon and a_muon) :
		h1.Fill(event.genParticles_pdgId)		
	if ( electron and positron) :
		h3.Fill(event.genParticles_pdgId)		

c1 = TCanvas("c1","c1",600,600)		
h1.Draw()
c1.SaveAs("c1.png")
c2 = TCanvas("c2","c2",600,600)		
h3.Draw()
c2.SaveAs("c2.png")
	
