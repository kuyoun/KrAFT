#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;
using namespace edm;

class FlatCandToNtupleMaker : public edm::EDAnalyzer
{
public:
  FlatCandToNtupleMaker(const edm::ParameterSet& pset);
  ~FlatCandToNtupleMaker() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  typedef edm::ParameterSet PSet;
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;
  typedef std::vector<edm::InputTag> VInputTag;
  typedef StringObjectFunction<reco::Candidate,true> CandFtn;
  typedef StringCutObjectSelector<reco::Candidate,true> CandSel;

  bool skipFailedEvent_;

  std::vector<edm::InputTag> weightLabels_;
  std::vector<edm::InputTag> vWeightLabels_;
  std::vector<edm::InputTag> candLabels_;
  std::vector<std::vector<CandFtn> > exprs_;
  std::vector<std::vector<CandSel> > selectors_;
  std::vector<VInputTag> vmapLabels_;

  TTree* tree_;
  int runNumber_, lumiNumber_, eventNumber_;
  std::vector<double*> weights_;
  std::vector<doubles*> vWeights_;
  std::vector<std::vector<doubles*> > candVars_;

};

FlatCandToNtupleMaker::FlatCandToNtupleMaker(const edm::ParameterSet& pset)
{
  // Output histograms and tree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("event", "event");

  tree_->Branch("run"  , &runNumber_  , "run/I"  );
  tree_->Branch("lumi" , &lumiNumber_ , "lumi/I" );
  tree_->Branch("event", &eventNumber_, "event/I");

  skipFailedEvent_ = pset.getUntrackedParameter<bool>("skipFailedEvent", true);

  edm::ParameterSet weightPSets = pset.getParameter<edm::ParameterSet>("weight");
  const strings weightNames = weightPSets.getParameterNamesForType<edm::ParameterSet>();
  for ( auto& weightName : weightNames )
  {
    edm::ParameterSet weightPSet = weightPSets.getParameter<edm::ParameterSet>(weightName);
    weightLabels_.push_back(weightPSet.getParameter<edm::InputTag>("src"));

    weights_.push_back(new double);
    tree_->Branch(weightName.c_str(), weights_.back(), (weightName+"/D").c_str());
  }

  edm::ParameterSet vWeightPSets = pset.getParameter<edm::ParameterSet>("vWeight");
  const strings vWeightNames = vWeightPSets.getParameterNamesForType<edm::ParameterSet>();
  for ( auto& vWeightName : vWeightNames )
  {
    edm::ParameterSet vWeightPSet = vWeightPSets.getParameter<edm::ParameterSet>(vWeightName);
    vWeightLabels_.push_back(vWeightPSet.getParameter<edm::InputTag>("src"));

    vWeights_.push_back(new doubles);
    tree_->Branch(vWeightName.c_str(), vWeights_.back());
  }

  edm::ParameterSet candPSets = pset.getParameter<edm::ParameterSet>("cands");
  const strings candNames = candPSets.getParameterNamesForType<edm::ParameterSet>();
  for ( auto& candName : candNames )
  {
    edm::ParameterSet candPSet = candPSets.getParameter<edm::ParameterSet>(candName);

    candLabels_.push_back(candPSet.getParameter<edm::InputTag>("src"));

    exprs_.push_back(std::vector<CandFtn>());
    selectors_.push_back(std::vector<CandSel>());
    vmapLabels_.push_back(VInputTag());
    candVars_.push_back(std::vector<doubles*>());
    const string candLabelName = candLabels_.back().label();
    const PSet exprSets = candPSet.getUntrackedParameter<PSet>("exprs", PSet());
    for ( auto& exprName : exprSets.getParameterNamesForType<string>() )
    {
      const string expr = exprSets.getParameter<string>(exprName);
      candVars_.back().push_back(new doubles);
      exprs_.back().push_back(CandFtn(expr));

      tree_->Branch((candName+"_"+exprName).c_str(), candVars_.back().back());
    }
    const PSet selectionSets = candPSet.getUntrackedParameter<PSet>("seletions", PSet());
    for ( auto& selectionName : selectionSets.getParameterNamesForType<string>() )
    {
      const string selection = selectionSets.getParameter<string>(selectionName);
      candVars_.back().push_back(new doubles);
      selectors_.back().push_back(CandSel(selection));

      tree_->Branch((candName+"_"+selectionName).c_str(), candVars_.back().back());
    }
    const strings vmapNames = candPSet.getUntrackedParameter<strings>("vmaps", strings());
    for ( auto& vmapName : vmapNames )
    {
      candVars_.back().push_back(new doubles);
      vmapLabels_.back().push_back(edm::InputTag(candLabelName, vmapName));

      tree_->Branch((candName+"_"+vmapName).c_str(), candVars_.back().back());
    }
  }

}

void FlatCandToNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  typedef edm::View<reco::LeafCandidate> Cands;

  runNumber_   = event.run();
  lumiNumber_  = event.luminosityBlock();
  eventNumber_ = event.id().event();

  for ( size_t i=0, n=weightLabels_.size(); i<n; ++i )
  {
    edm::Handle<double> weightHandle;
    event.getByLabel(weightLabels_[i], weightHandle);
    if ( skipFailedEvent_ and !weightHandle.isValid() ) return;

    *weights_[i] = *weightHandle;
  }

  for ( size_t i=0, n=vWeightLabels_.size(); i<n; ++i )
  {
    edm::Handle<doubles> vWeightHandle;
    event.getByLabel(vWeightLabels_[i], vWeightHandle);
    if ( skipFailedEvent_ and !vWeightHandle.isValid() ) return;

    vWeights_[i]->insert(vWeights_[i]->begin(), vWeightHandle->begin(), vWeightHandle->end());
  }

  const size_t nCand = candLabels_.size();
  for ( size_t iCand=0; iCand < nCand; ++iCand )
  {
    edm::Handle<Cands> srcHandle;
    event.getByLabel(candLabels_[iCand], srcHandle);
    if ( skipFailedEvent_ and !srcHandle.isValid() ) return;

    const vector<CandFtn>& exprs = exprs_[iCand];
    const vector<CandSel>& selectors = selectors_[iCand];
    VInputTag& vmapLabels = vmapLabels_[iCand];
    const size_t nExpr = exprs.size();
    const size_t nSels = selectors.size();
    const size_t nVmap = vmapLabels.size();
    std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVmap);
    for ( size_t iVar=0; iVar<nVmap; ++iVar )
    {
      event.getByLabel(vmapLabels[iVar], vmapHandles[iVar]);
      if ( skipFailedEvent_ and !vmapHandles[iVar].isValid() ) return;
    }

    for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
    {
      edm::Ref<Cands> candRef(srcHandle, i);

      for ( size_t j=0; j<nExpr; ++j )
      {
        const double val = exprs[j](*candRef);
        candVars_[iCand][j]->push_back(val);
      }
      for ( size_t j=0; j<nSels; ++j )
      {
        const double val = selectors[j](*candRef);
        candVars_[iCand][j+nExpr]->push_back(val);
      }
      for ( size_t j=0; j<nVmap; ++j )
      {
        const double val = (*vmapHandles[j])[candRef];
        candVars_[iCand][j+nExpr+nSels]->push_back(val);
      }
    }
  }

  tree_->Fill();
  for ( auto& vWeight : vWeights_ )
  {
    vWeight->clear();
  }

  for ( size_t iCand=0; iCand<nCand; ++iCand )
  {
    const size_t nVar = candVars_[iCand].size();
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      candVars_[iCand][iVar]->clear();
    }
  }
}

DEFINE_FWK_MODULE(FlatCandToNtupleMaker);

