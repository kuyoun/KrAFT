#!/usr/bin/env python
import os

site = "uosaf0007.sscc.uos.ac.kr"
dir = "/cms/store/user/geonmo/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/TTJet_MC_Dilepton_KrAFT_20140505_TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola/fec648b152d83b5e5d7824a66c836bae"

GetListCmd = "xrd "+site+" ls "+dir+" | awk '{ print $5}'"  
list = os.popen(GetListCmd).readlines()
for x in list :
	print x.strip()
