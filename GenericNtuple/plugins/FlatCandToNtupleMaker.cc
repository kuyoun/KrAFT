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

#include "KrAFT/GenericNtuple/interface/GenericEvent.h"

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

  std::vector<edm::InputTag> candLabels_;
  std::vector<VInputTag> vmapLabels_;

  TTree* tree_;
  std::vector<doubles*> candPt_, candEta_, candPhi_, candM_;
  std::vector<std::vector<doubles*> > candVars_;

};

FlatCandToNtupleMaker::FlatCandToNtupleMaker(const edm::ParameterSet& pset)
{
  // Output histograms and tree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("event", "event");

  edm::ParameterSet candPSets = pset.getParameter<edm::ParameterSet>("cands");
  const strings candNames = candPSets.getParameterNamesForType<edm::ParameterSet>();
  for ( auto& candName : candNames )
  {
    edm::ParameterSet candPSet = candPSets.getParameter<edm::ParameterSet>(candName);
    candLabels_.push_back(candPSet.getParameter<edm::InputTag>("src"));

    candPt_ .push_back(new doubles);
    candEta_.push_back(new doubles);
    candPhi_.push_back(new doubles);
    candM_  .push_back(new doubles);

    tree_->Branch((candName+"_pt" ).c_str(), candPt_ .back());
    tree_->Branch((candName+"_eta").c_str(), candEta_.back());
    tree_->Branch((candName+"_phi").c_str(), candPhi_.back());
    tree_->Branch((candName+"_m"  ).c_str(), candM_  .back());

    vmapLabels_.push_back(VInputTag());
    candVars_.push_back(std::vector<doubles*>());
    const string candLabelName = candLabels_.back().label();
    const strings vmapNames = candPSet.getParameter<strings>("vmaps");
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

  const size_t nCand = candLabels_.size();
  for ( size_t iCand=0; iCand < nCand; ++iCand )
  {
    edm::Handle<Cands> srcHandle;
    event.getByLabel(candLabels_[iCand], srcHandle);

    VInputTag vmapLabels = vmapLabels_[iCand];
    const size_t nVar = vmapLabels.size();
    std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVar);
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      event.getByLabel(vmapLabels[iVar], vmapHandles[iVar]);
    }

    for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
    {
      edm::Ref<Cands> candRef(srcHandle, i);
      candPt_[iCand]->push_back(candRef->pt());
      candEta_[iCand]->push_back(candRef->eta());
      candPhi_[iCand]->push_back(candRef->phi());
      candM_[iCand]->push_back(candRef->mass());

      for ( size_t iVar=0; iVar<nVar; ++iVar )
      {
        const double var = (*vmapHandles[iVar])[candRef];
        candVars_[iCand][iVar]->push_back(var);
      }
    }
  }

  tree_->Fill();
  for ( size_t iCand=0; iCand<nCand; ++iCand )
  {
    candPt_ [iCand]->clear();
    candEta_[iCand]->clear();
    candPhi_[iCand]->clear();
    candM_  [iCand]->clear();
    const size_t nVar = candVars_[iCand].size();
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      candVars_[iCand][iVar]->clear();
    }
  }
}

DEFINE_FWK_MODULE(FlatCandToNtupleMaker);

