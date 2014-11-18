KrAFT
=====

Korea CMS Analysis Framwork for Top quark physics
  * Documentation and project overview : https://twiki.cern.ch/twiki/bin/viewauth/CMS/KrAFT
  * PAT recipe : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATRecipes#CMSSW_7_2_X_dev2014
  * TQAF recipe : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTQAFRecipes

## Installation

```sh
# Fresh install your workarea, baseline release is set to 7_X_Y
cmsrel CMSSW_7_2_0
cd CMSSW_7_2_0/src
cmsenv
git-cms-addpkg PhysicsTools/PatAlgos
git-cms-addpkg EgammaAnalysis/ElectronTools
git-cms-addpkg TopQuarkAnalysis/Configuration
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd -

# Install KrAFT
git clone https://github.com/cms-kr/KrAFT

# In case you already forked this project and you will put your contributions...
cd KrAFT
git remote add $USER git@github.com:$(git-config user.github)/KrAFT
git fetch $USER
git checkout 7_X_Y
cd ..

# Patch missing module
git apply KrAFT/missing.patch

# Continue to build whole package
scram setup lhapdffull # Necessary to speed up PDF weight calculation
scram b clean
scram b -j8

```
