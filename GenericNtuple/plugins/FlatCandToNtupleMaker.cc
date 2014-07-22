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
  typedef edm::View<reco::LeafCandidate> CandView;
  typedef edm::ValueMap<double> Vmap;
  typedef edm::EDGetTokenT<CandView> CandToken;
  typedef edm::EDGetTokenT<Vmap> VmapToken;

  std::vector<CandToken> candTokens_;
  std::vector<std::vector<VmapToken> > vmapTokens_;

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
  for ( size_t i=0, n=candNames.size(); i<n; ++i )
  {
    const string& candName = candNames[i];
    edm::ParameterSet candPSet = candPSets.getParameter<edm::ParameterSet>(candName);
    edm::InputTag candLabel = candPSet.getParameter<edm::InputTag>("src");
    candTokens_.push_back(consumes<CandView>(candLabel));

    candPt_ .push_back(new doubles);
    candEta_.push_back(new doubles);
    candPhi_.push_back(new doubles);
    candM_  .push_back(new doubles);

    tree_->Branch((candName+"_pt" ).c_str(), candPt_ .back());
    tree_->Branch((candName+"_eta").c_str(), candEta_.back());
    tree_->Branch((candName+"_phi").c_str(), candPhi_.back());
    tree_->Branch((candName+"_m"  ).c_str(), candM_  .back());

    vmapTokens_.push_back(std::vector<VmapToken>());
    candVars_.push_back(std::vector<doubles*>());
    const string candLabelName = candLabel.label();
    const strings vmapNames = candPSet.getParameter<strings>("vmaps");
    for ( size_t j=0, m=vmapNames.size(); j<m; ++j )
    {
      const string& vmapName = vmapNames[j];
      candVars_.back().push_back(new doubles);

      edm::InputTag vmapLabel(candLabelName, vmapName);
      vmapTokens_.back().push_back(consumes<Vmap>(vmapLabel));

      tree_->Branch((candName+"_"+vmapName).c_str(), candVars_.back().back());
    }
  }

}

void FlatCandToNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  const size_t nCand = candTokens_.size();
  for ( size_t iCand=0; iCand < nCand; ++iCand )
  {
    edm::Handle<CandView> srcHandle;
    event.getByToken(candTokens_[iCand], srcHandle);

    std::vector<VmapToken>& vmapTokens = vmapTokens_[iCand];
    const size_t nVar = vmapTokens.size();
    std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVar);
    for ( size_t iVar=0; iVar<nVar; ++iVar )
    {
      VmapToken& vmapToken = vmapTokens[iVar];
      event.getByToken(vmapToken, vmapHandles[iVar]);
    }

    for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
    {
      edm::Ref<CandView> candRef(srcHandle, i);
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

