#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")
from xml.dom.minidom import parse

lumi = 19.6*1000

## Load plot styles for all samples
plotStyles = []
xml_doc = parse("../data/samples.xml")
for xml_proc in xml_doc.getElementsByTagName("proc"):
    plotColl = []
    title = xml_proc.getAttribute("title")
    color = gROOT.ProcessLine(xml_proc.getAttribute("color")+";")
    for xml_sub in xml_proc.getElementsByTagName("sub"):
        xsec = eval(xml_sub.getAttribute("xsec"))
        sampleName = xml_sub.firstChild.nodeValue
        plotColl.append( (xsec, sampleName) )
    plotStyles.append( (title, color, plotColl) )

## Load histogram structures from an example file
histStr = {}
f = TFile("hist/%s__All.root" % (plotStyles[0][2][0][1]))
for cutStep in [x.GetName() for x in f.GetListOfKeys()]:
    d = f.GetDirectory(cutStep)
    if d == None: continue
    histStr[cutStep] = []
    for histName in [y.GetName() for y in d.GetListOfKeys()]:
        h = d.Get(histName)
        if h == None: continue
        if not h.IsA().InheritsFrom("TH1"): continue
        histStr[cutStep].append(histName)

canvases = []
plotSets = []
for cutStep in sorted(histStr.keys()):
    for histName in histStr[cutStep]:
        c = TCanvas("c%s_%s" % (cutStep, histName), "%s/%s" % (cutStep, histName), 600, 600)
        c.Divide(2,2)
        plotSet = [c]
        ymax = 0.0
        for i, channel in enumerate(("MuMu", "ElEl", "MuEl", "All")):
            hists = []
            hStack = None
            ymaxStack = 0.0
            for title, color, plotColl in plotStyles:
                hSum = None
                for xsec, sampleName in plotColl:
                    f = TFile("hist/%s__%s.root" % (sampleName, channel))
                    scale = lumi*xsec

                    h = f.Get("%s/%s" % (cutStep, histName))
                    if hSum == None:
                        gROOT.cd()
                        hSum = h.Clone()
                        hSum.SetName("%s_%s_%s" % (channel, cutStep, histName))
                        hSum.SetTitle(title)
                        xTitle = hSum.GetXaxis().GetTitle()
                        yTitle = hSum.GetYaxis().GetTitle()
                        hSum.Reset()
                        hSum.SetFillColor(color)
                        hSum.SetOption("hist")
                        hists.append(hSum)
                    hSum.Add(hSum, h, 1, scale)
                if hStack == None:
                    hStack = THStack("hStack_%s_%s_%s" % (cutStep, histName, channel), "%s/%s/%s;%s;%s" % (cutStep, histName, channel, xTitle, yTitle))
                hStack.Add(hSum)

                ymaxStack += max([hSum.GetBinContent(j+1) for j in range(hSum.GetNbinsX())])
            ymax = max(ymax, ymaxStack)

            fData = TFile("hist/Run2012__%s.root" % channel)
            gROOT.cd()
            hData = fData.Get("%s/%s" % (cutStep, histName)).Clone()
            hData.SetName("data_%s_%s_%s" % (channel, cutStep, histName))

            ymax = max(ymax, max([hData.GetBinContent(j+1) for j in range(hData.GetNbinsX())]))

            labels = TLegend(0.65,0.65,0.92,0.92)
            labels.SetBorderSize(0)
            labels.SetFillStyle(0)
            labels.AddEntry(hData, "Data", "lp")
            for h in reversed(hists):
                labels.AddEntry(h, h.GetTitle(), "f")
            c.cd(i+1)
            hStack.Draw("hist")
            hData.Draw("same")
            labels.Draw()
            plotSet.append((hStack, hData, hists, labels))

        for i in range(4):
            pad = c.cd(i+1)
            if histName in ("zM", "njet", "nbjet"):
                pad.SetLogy()
                plotSet[i+1][0].SetMinimum(0.1)
                plotSet[i+1][0].SetMaximum(ymax*10)
            else:
                plotSet[i+1][0].SetMinimum(0)
                plotSet[i+1][0].SetMaximum(ymax*1.2)
        canvases.append(c)

        plotSets.append(plotSet)

for c in canvases:
    c.Update()
