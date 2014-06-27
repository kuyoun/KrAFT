KrAFT
=====

Korea CMS Analysis Framwork for Top quark physics
cmsrel CMSSW_5_3_18
cd CMSSW_5_3_18/src
cmsenv
git-cms-addpkg FWCore/Version
git cms-merge-topic -u cms-analysis-tools:5_3_18-updateTopRefSel
git-cms-addpkg EgammaAnalysis/ElectronTools
git-cms-addpkg TopQuarkAnalysis/Configuration
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd -
git clone https://github.com/cms-kr/KrAFT
scram setup lhapdffull # Necessary to speed up PDF weight calculation
scram b clean
scram b -j 20
git remote add geonmo-kraft git@github.com:geonmo/KrAFT.git
git pull geonmo-kraft jpsi-update2
