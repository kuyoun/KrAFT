#!/usr/bin/env python

from ROOT import *

prefix = "Run2012"
for suffix in ["", "Up", "Dn"]:
    f = TFile("pu%s.root" % suffix)
    h = f.Get("pileup")
    l = []
    for bin in range(1, h.GetNbinsX()+1):
        c = h.GetBinContent(bin)
        l.append(c)
    content = '"%s%s":cms.vdouble(\n' % (prefix, suffix)
    for row in xrange(len(l)/5):
        content += "    "
        content += " ".join("%12e," % l[row*5+col] for col in range(5))
        content += "\n"
    content += '),'

    print content
