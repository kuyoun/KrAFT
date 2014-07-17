#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/OneToValue.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <memory>
#include <vector>
#include <string>
#include <boost/assign.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost::assign;

class FlatCandProducer : public edm::EDProducer
{
public:
  FlatCandProducer(const edm::ParameterSet& pset);
  virtual ~FlatCandProducer()
  {
    delete loader_;
  };

  void produce(edm::Event& event, const edm::EventSetup& eventSetup);
 
private:
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;

  typedef std::vector<reco::LeafCandidate> Cands;
  typedef edm::ValueMap<double> CandValueMap;
  typedef edm::RefProd<Cands> CandRefProd;
  typedef edm::Ref<Cands> CandRef;

  edm::InputTag srcLabel_;
  strings varNames_;

  struct LoaderBase
  {
    LoaderBase(const int n):N(n)
    {
      v.resize(n);
    };
    virtual ~LoaderBase() {};
    void init()
    {
      for ( std::vector<std::vector<double> >::iterator itr = v.begin();
            itr != v.end(); ++itr ) itr->clear();
      // foreach ( auto& x; v ) x.clear();
    };
    virtual void load(const reco::Candidate&) = 0;
    std::vector<std::vector<double> > v;
    const size_t N;
  } * loader_;

  struct LoadMuon : LoaderBase
  {
    LoadMuon(const int n):LoaderBase(n) {};
    void load(const reco::Candidate& cand)
    {
      const pat::Muon& muon = dynamic_cast<const pat::Muon&>(cand);
      v[0].push_back(muon.pt());
    };
  };

  struct LoadElectron : LoaderBase
  {
    LoadElectron(const int n):LoaderBase(n) {};
    void load(const reco::Candidate& cand)
    {
      const pat::Electron& electron = dynamic_cast<const pat::Electron&>(cand);
      v[0].push_back(electron.pt());
    };
  };

  struct LoadJet : LoaderBase
  {
    LoadJet(const int n):LoaderBase(n) {};
    void load(const reco::Candidate& cand)
    {
      const pat::Jet& jet = dynamic_cast<const pat::Jet&>(cand);
      v[0].push_back(jet.bDiscriminator("CSV"));
    };
  };

};

FlatCandProducer::FlatCandProducer(const edm::ParameterSet& pset)
{
  srcLabel_ = pset.getParameter<edm::InputTag>("src");

  std::string type = pset.getParameter<std::string>("type");
  boost::algorithm::to_lower(type);
  if ( type == "muon" )
  {
    varNames_ += "isTight", "isLoose";
    varNames_ += "relIsoDbeta03", "relIsoDbeta04";
    loader_ = new LoadMuon(varNames_.size());
  }
  else if ( type == "electron" )
  {
    varNames_ += "mva";
    varNames_ += "relIsoDbeta03", "relIsoDbeta04", "relIsoRho03", "relIsoRho04";
    loader_ = new LoadElectron(varNames_.size());
  }
  else if ( type == "jet" )
  {
    varNames_ += "bTagCSV";
    loader_ = new LoadJet(varNames_.size());
  }

  for ( strings::const_iterator itr = varNames_.begin(); itr != varNames_.end(); ++itr )
  {
    produces<CandValueMap>(*itr);
  }
}

void FlatCandProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::Candidate> > srcHandle;
  event.getByLabel(srcLabel_, srcHandle);

  loader_->init();
  std::auto_ptr<Cands> cands(new Cands);
  CandRefProd refProd = event.getRefBeforePut<Cands>();

  // Fill four vector informations
  for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
  {
    const reco::Candidate& srcCand = srcHandle->at(i);
    reco::LeafCandidate cand(srcCand.charge(), srcCand.p4());
    cands->push_back(cand);
    loader_->load(srcCand);
  }

  event.put(cands);
  for ( size_t i=0; i<loader_->N; ++i )
  {
    const std::string& varName = varNames_[i];
    std::auto_ptr<CandValueMap> vmap(new CandValueMap);
    CandValueMap::Filler filler(*vmap);
    filler.insert(refProd, loader_->v[i].begin(), loader_->v[i].end());
    event.put(vmap, varName);
  }
}

DEFINE_FWK_MODULE(FlatCandProducer);

