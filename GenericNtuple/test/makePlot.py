#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")
from xml.dom.minidom import parse

## Set lumi value and label
lumi = 19.6*1000 # in pb-1
header = "CMS Work in progress"
#header = "CMS Preliminary"
#header = "CMS"
label  = "#intLdt=%.1f fb^{-1}, #sqrt{s}=8TeV" % (lumi/1000)

## Load plot styles for all samples
plotStyles = []
xml_doc = parse("../data/samples.xml")
for xml_proc in xml_doc.getElementsByTagName("signal") + xml_doc.getElementsByTagName("background"):
    plotColl = []
    title = xml_proc.getAttribute("title")
    color = gROOT.ProcessLine(xml_proc.getAttribute("color")+";")
    for xml_sub in xml_proc.getElementsByTagName("sample"):
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

## Draw everything
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
                if hSum == None: continue
                hSum.AddBinContent(hSum.GetNbinsX(), hSum.GetBinContent(hSum.GetNbinsX()+1))
                if hStack == None:
                    hStack = THStack("hStack_%s_%s_%s" % (cutStep, histName, channel), "%s/%s/%s;%s;%s" % (cutStep, histName, channel, xTitle, yTitle))
                hStack.Add(hSum)

                ymaxStack += max([hSum.GetBinContent(j+1) for j in range(hSum.GetNbinsX()-1)])
            ymax = max(ymax, ymaxStack)

            fData = TFile("hist/Run2012__%s.root" % channel)
            gROOT.cd()
            hData = fData.Get("%s/%s" % (cutStep, histName)).Clone()
            hData.SetName("data_%s_%s_%s" % (channel, cutStep, histName))
            hData.AddBinContent(hData.GetNbinsX(), hData.GetBinContent(hData.GetNbinsX()+1))

            ymax = max(ymax, max([hData.GetBinContent(j+1) for j in range(hData.GetNbinsX()-1)]))

            legends = TLegend(0.70,0.62,0.95,0.90)
            legends.SetBorderSize(0)
            legends.SetFillStyle(0)
            legends.AddEntry(hData, "Data", "lp")
            for h in reversed(hists):
                legends.AddEntry(h, h.GetTitle(), "f")
            c.cd(i+1)
            hStack.Draw("hist")
            hData.Draw("same")
            legends.Draw()
            plotSet.append([hStack, hData, hists, legends])

        for i in range(4):
            pad = c.cd(i+1)
            if histName in ("zM", "njet", "nbjet"):
                pad.SetLogy()
                plotSet[i+1][0].SetMinimum(0.1)
                plotSet[i+1][0].SetMaximum(ymax*100)
            else:
                plotSet[i+1][0].SetMinimum(0)
                plotSet[i+1][0].SetMaximum(ymax*1.5)

            leftMargin, topMargin = pad.GetLeftMargin(), pad.GetTopMargin()
            lh = TLatex(leftMargin+0.04, 1-topMargin-0.08, header)
            ll = TLatex(leftMargin+0.04, 1-topMargin-0.15, label)
            lh.SetNDC()
            ll.SetNDC()
            lh.SetTextSize(0.045)
            ll.SetTextSize(0.035)
            ll.Draw()
            lh.Draw()
            plotSet[i+1].extend([lh, ll])

        canvases.append(c)

        plotSets.append(plotSet)

for c in canvases:
    c.Update()
    c.Print("%s.png" % c.GetName())
