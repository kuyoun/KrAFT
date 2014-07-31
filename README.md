KrAFT
=====

Korea CMS Analysis Framwork for Top quark physics


=== Installation ===

cmsrel CMSSW_5_3_20
cd CMSSW_5_3_20/src
cmsenv
git-cms-addpkg PhysicsTools/PatAlgos
git-cms-addpkg EgammaAnalysis/ElectronTools
git-cms-addpkg TopQuarkAnalysis/Configuration
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd -

git clone https://github.com/cms-kr/KrAFT

## In case you already forked this project and wanting to put your contributions...
cd KrAFT
git remote add $USER git@github.com:$(git-config user.github)/KrAFT
git fetch $USER
git checkout

cd ..
scram b -j8
