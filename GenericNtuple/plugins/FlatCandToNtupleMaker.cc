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
#include "DataFormats/Candidate/interface/LeafCandidate.h"

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
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;
  typedef std::vector<edm::InputTag> VInputTag;

  bool skipFailedEvent_;

  std::vector<edm::InputTag> weightLabels_;
  std::vector<edm::InputTag> vWeightLabels_;
  std::vector<edm::InputTag> candLabels_;
  std::vector<VInputTag> vmapLabels_;

  TTree* tree_;
  int runNumber_, lumiNumber_, eventNumber_;
  std::vector<double*> weights_;
  std::vector<doubles*> vWeights_;
  std::vector<doubles*> candPt_, candEta_, candPhi_, candM_, candQ_, candPdg_;
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
    const bool fillPt  = candPSet.getUntrackedParameter<bool>("fillPt" , true);
    const bool fillEta = candPSet.getUntrackedParameter<bool>("fillEta", true);
    const bool fillPhi = candPSet.getUntrackedParameter<bool>("fillPhi", true);
    const bool fillM   = candPSet.getUntrackedParameter<bool>("fillM"  , true);
    const bool fillQ   = candPSet.getUntrackedParameter<bool>("fillQ"  , true);
    const bool fillPdg = candPSet.getUntrackedParameter<bool>("fillPdg", true);

    candLabels_.push_back(candPSet.getParameter<edm::InputTag>("src"));

    candPt_ .push_back(new doubles);
    candEta_.push_back(new doubles);
    candPhi_.push_back(new doubles);
    candM_  .push_back(new doubles);
    candQ_  .push_back(new doubles);
    candPdg_.push_back(new doubles);

    if ( fillPt  ) tree_->Branch((candName+"_pt" ).c_str(), candPt_ .back());
    if ( fillEta ) tree_->Branch((candName+"_eta").c_str(), candEta_.back());
    if ( fillPhi ) tree_->Branch((candName+"_phi").c_str(), candPhi_.back());
    if ( fillM   ) tree_->Branch((candName+"_m"  ).c_str(), candM_  .back());
    if ( fillQ   ) tree_->Branch((candName+"_q"  ).c_str(), candQ_  .back());
    if ( fillPdg ) tree_->Branch((candName+"_pdgId").c_str(), candPdg_.back());

    vmapLabels_.push_back(VInputTag());
    candVars_.push_back(std::vector<doubles*>());
    const string candLabelName = candLabels_.back().label();
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

    VInputTag vmapLabels = vmapLabels_[iCand];
    const size_t nVar = vmapLabels.size();
    std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVar);
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      event.getByLabel(vmapLabels[iVar], vmapHandles[iVar]);
      if ( skipFailedEvent_ and !vmapHandles[iVar].isValid() ) return;
    }

    for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
    {
      edm::Ref<Cands> candRef(srcHandle, i);
      candPt_[iCand]->push_back(candRef->pt());
      candEta_[iCand]->push_back(candRef->eta());
      candPhi_[iCand]->push_back(candRef->phi());
      candM_[iCand]->push_back(candRef->mass());
      candQ_[iCand]->push_back(candRef->charge());
      candPdg_[iCand]->push_back(candRef->pdgId());

      for ( size_t iVar=0; iVar<nVar; ++iVar )
      {
        const double var = (*vmapHandles[iVar])[candRef];
        candVars_[iCand][iVar]->push_back(var);
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
    candPt_ [iCand]->clear();
    candEta_[iCand]->clear();
    candPhi_[iCand]->clear();
    candM_  [iCand]->clear();
    candQ_  [iCand]->clear();
    candPdg_[iCand]->clear();
    const size_t nVar = candVars_[iCand].size();
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      candVars_[iCand][iVar]->clear();
    }
  }
}

DEFINE_FWK_MODULE(FlatCandToNtupleMaker);

