#!/usr/bin/env python

import sys, os

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")

if not os.path.isdir("hist"): os.makedirs("hist")

mode = "ElEl"
sample = "DYJetsToLL_M-10To50filter_8TeV-madgraph"

ana = KDileptonTreeReducer("ElEl",
    "/pnfs/user/jhgoh/data/ntuple/20131224_1/%s.root" % sample,
    "hist/%s__%s.root" % (sample, mode))
ana.run()
