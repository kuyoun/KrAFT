#!/usr/bin/env python
import sys, os, time
from itertools import product
import multiprocessing

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")

def fed(args):
    mode, file = args
    sample = os.path.basename(file).replace('.root', '')
    print "Running", mode, sample 
    ana = KJpsiTreeAnalyzer(mode, file, "ntuple/%s__%s.root" % (sample, mode))
    ana.run()
    print "Done", mode, sample

if __name__ == '__main__':
    #if os.path.isdir("ntuple"):
    #    print "ntuple directory already exists. rename existing directory and rerun this script"
    #    sys.exit()
    if not os.path.isdir("ntuple"): os.makedirs("ntuple")

    samples = []
    basedir = "MC/"
    samples.extend([basedir+x for x in os.listdir(basedir) if 'root' in x])
    basedir = "RD/"
    samples.extend([basedir+x for x in os.listdir(basedir) if 'root' in x])
    p = multiprocessing.Pool(multiprocessing.cpu_count())
    p.map(fed, [(x, y) for x in ("MuMu", "ElEl", "MuEl") for y in samples])

