#!/usr/bin/env python

import sys, os, time
import multiprocessing

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")
from KrAFT.GenericNtuple.NtupleAnalyzer import *

def process(sample, mode, files):
    print "Making histogram", sample

    ana = NtupleAnalyzer(mode, files, "hist/%s__%s.root" % (sample, mode))
    ana.setWeightVar("puWeight")

    cut_s1 = "z_m > 12 && lepton1_pt > 20 && lepton2_pt > 20 && lepton1_iso < 0.15 && lepton2_iso < 0.15 && z_Q == 0"
    cut_s2 = "abs(z_m-91.2) > 15"
    cut_s3 = "@jets_pt.size() >= 2"
    cut_s4 = "met_pt > 30"
    cut_s5 = "bjets_n >= 1"

    if mode == "MuEl":
        cut_s2 = "1"
        cut_s4 = "1"
    if "TTJets_MassiveBinDECAY" in sample:
        cut_s1 += "&& decayMode < 1"
    elif "TTJets" in sample or "TTTo" in sample or "TT_" in sample:
        cut_s1 += "&& decayMode >= 1"

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

    print "Done", sample

if __name__ == '__main__':
    samples = {}
    for mode in ["MuMu", "MuEl", "ElEl"]:
        for f in os.listdir("ntuple"):
            if len(f) < 12: continue
            if f[-11:] != ('__%s.root' % mode): continue
            s = f.replace("__%s.root" % mode, "")
            sampleName = ('_'.join(s.split('_')[:-1]), mode)
            if sampleName not in samples: samples[sampleName] = []
            samples[sampleName].append("ntuple/%s" % f)

    p = multiprocessing.Pool(multiprocessing.cpu_count())
    for key in samples.keys():
        sample, mode = key
        files = samples[key]
        p.apply_async(process, [sample, mode, files])
    p.close()
    p.join()

    for sample, mode in samples.keys():
        if 'Run20' in sample: continue
        os.system("hadd -f hist/%s__All.root hist/%s__MuMu.root hist/%s__ElEl.root hist/%s__MuEl.root" % (sample, sample, sample, sample))
    os.system("hadd -f hist/Run2012__ElEl.root hist/DoubleElectron_Run2012*-22Jan2013__ElEl.root")
    os.system("hadd -f hist/Run2012__MuMu.root hist/DoubleMu_Run2012*-22Jan2013__MuMu.root")
    os.system("hadd -f hist/Run2012__MuEl.root hist/MuEG_Run2012*-22Jan2013__MuEl.root")
    os.system("hadd -f hist/Run2012__All.root hist/Run2012__MuMu.root hist/Run2012__ElEl.root hist/Run2012__MuEl.root")
    print "END"

