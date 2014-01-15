#!/usr/bin/env python

import sys, os, time
from multiprocessing import Pool

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")
from KrAFT.GenericNtuple.NtupleAnalyzer import *

def process(sample, mode):
    ana = NtupleAnalyzer(mode, ["ntuple/%s__%s.root" % (sample, mode)], "hist/%s__%s.root" % (sample, mode))
    ana.setWeightVar("puWeight")

    cut_s1 = "z_m > 12 && lepton1_pt > 20 && lepton2_pt > 20 && lepton1_iso < 0.15 && lepton2_iso < 0.15 && z_Q == 0"
    cut_s2 = "abs(z_m-91.2) > 15"
    cut_s3 = "@jets_pt.size() >= 2"
    cut_s4 = "met_pt > 30"
    cut_s5 = "bjets_n >= 1"

    if mode == "MuEl":
        cut_s2 = "1"
        cut_s4 = "1"

    ana.addH1("zM", "z_m", "MZ;M(l^{+}l^{-}) (GeV/c^{2});Events per 2 GeV/c^{2}", 100, 0, 200)
    ana.addH1("lepton1_pt", "lepton1_pt", "pt1;Leading lepton p_{T} (GeV/c);Events per 2GeV/c", 100, 0, 200)
    ana.addH1("lepton2_pt", "lepton2_pt", "pt2;2nd leading lepton p_{T} (GeV/c);Events per 2GeV/c", 100, 0, 200)
    ana.addH1("njet", "@jets_pt.size()", "njet;Jet multiplicity;Events", 10, 0, 10)
    ana.addH1("nbjet", "bjets_n", "nbjet;B jet multiplicity;Events", 10, 0, 10)
    ana.addH1("met", "met_pt", "met;Missing transverse momentum p_{T} (GeV/c);Events per 2GeV/c", 100, 0, 200)
    ana.addH1("nVertex", "nVertex", "nVertex;Vertex multiplicity;Events", 60, 0, 60)

    ana.addCutStep("S1", cut_s1, "zM,lepton1_pt,lepton2_pt,nVertex")
    ana.addCutStep("S2", cut_s2, "zM,njet,nbjet,met,nVertex")
    ana.addCutStep("S3", cut_s3, "njet,nbjet,met,nVertex")
    ana.addCutStep("S4", cut_s4, "nbjet,met,nVertex")
    ana.addCutStep("S5", cut_s5, "nbjet,met,nVertex")

    ana.process()

if __name__ == '__main__':
    p = Pool(4)
    for mode in ["MuMu", "MuEl", "ElEl"]:
        for f in os.listdir("ntuple"):
            if len(f) < 12: continue
            if f[-11:] != ('__%s.root' % mode): continue
            sample = f.replace("__%s.root" % mode, "")

            p.apply(process, [sample, mode])
    p.join()

