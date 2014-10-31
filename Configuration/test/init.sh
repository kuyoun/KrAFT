#!/bin/bash

eval `scram runtime -sh`

echo "@@@ Staring to calculate PU weights @@@"
JSONDIR=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV
NBIN=60 # 60 bins used in Summer12 S10 sample
MBXSEC=69400 # 69.4mb for 2012 data, 68mb for 2011 data. See twiki:PileupSystematicErrors
MBXSECERR=5 # in percent, +-5% variation according to twiki:PileupSystematicErrors
echo "Number of bins = $NBIN"
echo "MinBias cross section = $MBXSEC"
echo "with uncertainty = += $MBXSECERR %"
MBXSECUP=$((($MBXSEC*(100+$MBXSECERR))/100))
MBXSECDN=$((($MBXSEC*(100-$MBXSECERR))/100))
pileupCalc.py -i $JSONDIR/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt \
              --inputLumiJSON $JSONDIR/PileUp/pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr_v2.txt \
              --calcMode true --minBiasXsec $MBXSEC \
              --maxPileupBin $NBIN --numPileupBin $NBIN ../data/pu.root &
pileupCalc.py -i $JSONDIR/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt \
              --inputLumiJSON $JSONDIR/PileUp/pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr_v2.txt \
              --calcMode true --minBiasXsec $MBXSECUP \
              --maxPileupBin $NBIN --numPileupBin $NBIN ../data/puUp.root &
pileupCalc.py -i $JSONDIR/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt \
              --inputLumiJSON $JSONDIR/PileUp/pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr_v2.txt \
              --calcMode true --minBiasXsec $MBXSECDN \
              --maxPileupBin $NBIN --numPileupBin $NBIN ../data/puDn.root &
wait
echo "Done."

echo "@@@ Extracting numbers from root file @@@"
python analysePU.py
echo "@@@ Put these numbers into KrAFT/GeneratorTools/python/pileupWeight_cff.py"
