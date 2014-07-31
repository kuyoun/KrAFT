KrAFT
=====

Korea CMS Analysis Framwork for Top quark physics
  * Documentation and project overview : https://twiki.cern.ch/twiki/bin/viewauth/CMS/KrAFT
  * PAT recipe : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATReleaseNotes52X
  * TQAF recipe : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTQAFRecipes

## Installation

```sh
# Fresh install your workarea, baseline release is set to 5_3_20
cmsrel CMSSW_5_3_20
cd CMSSW_5_3_20/src
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
git checkout

# Continue to build whole package
cd ..
scram b -j8
```
