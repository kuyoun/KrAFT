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

#define MUONVARS "isTight", "isLoose", "relIso", "dxy"
#define ELECTRONVARS "mva", "relIso", "scEta", "dxy", "chargeID"
#define JETVARS "bTagCSV", "JESup", "JESdn", "JER", "JERup", "JERdn"

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
  std::vector<edm::InputTag> vmapLabels_;
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
    virtual void load(const reco::Candidate&, const std::vector<double>&) = 0;
    std::vector<std::vector<double> > v;
    const size_t N;
  } * loader_;

  struct LoadMuon : LoaderBase
  {
    LoadMuon(const int n):LoaderBase(n) {};
    void load(const reco::Candidate& cand, const std::vector<double>& vv)
    {
      const pat::Muon& muon = dynamic_cast<const pat::Muon&>(cand);
      v[0].push_back(0.);
      v[1].push_back(1.*muon.isLooseMuon());
      v[2].push_back(muon.userIso(1));
      v[3].push_back(muon.dB());
    };
  };

  struct LoadElectron : LoaderBase
  {
    LoadElectron(const int n):LoaderBase(n) {};
    void load(const reco::Candidate& cand, const std::vector<double>& vv)
    {
      const pat::Electron& e = dynamic_cast<const pat::Electron&>(cand);
      v[0].push_back(e.electronID("mvaTrigV0"));
      v[1].push_back(e.userIso(2));
      v[2].push_back(e.superCluster()->eta());
      v[3].push_back(e.dB());
      int chargeID = 0;
      if ( e.isGsfCtfScPixChargeConsistent() ) chargeID = 3;
      else if ( e.isGsfScPixChargeConsistent() ) chargeID = 2;
      else if ( e.isGsfCtfChargeConsistent() ) chargeID = 1;
      v[4].push_back(chargeID);
    };
  };

  struct LoadJet : LoaderBase
  {
    LoadJet(const int n):LoaderBase(n) {};
    void load(const reco::Candidate& cand, const std::vector<double>& vv)
    {
      const pat::Jet& jet = dynamic_cast<const pat::Jet&>(cand);
      v[0].push_back(jet.bDiscriminator("combinedSecondaryVertexBJetTags"));
      v[1].push_back(vv[0]);
      v[2].push_back(vv[1]);
      v[3].push_back(vv[2]);
      v[4].push_back(vv[3]);
      v[5].push_back(vv[4]);
    };
  };

};

FlatCandProducer::FlatCandProducer(const edm::ParameterSet& pset)
{
  srcLabel_ = pset.getParameter<edm::InputTag>("src");
  vmapLabels_ = pset.getParameter<std::vector<edm::InputTag> >("vmaps");

  std::string type = pset.getParameter<std::string>("type");
  boost::algorithm::to_lower(type);
  if ( type == "muon" )
  {
    varNames_ += MUONVARS;
    loader_ = new LoadMuon(varNames_.size());
  }
  else if ( type == "electron" )
  {
    varNames_ += ELECTRONVARS;
    loader_ = new LoadElectron(varNames_.size());
  }
  else if ( type == "jet" )
  {
    varNames_ += JETVARS;
    loader_ = new LoadJet(varNames_.size());
  }

  produces<Cands>();
  for ( strings::const_iterator itr = varNames_.begin(); itr != varNames_.end(); ++itr )
  {
    produces<CandValueMap>(*itr);
  }
}

void FlatCandProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::Candidate> > srcHandle;
  event.getByLabel(srcLabel_, srcHandle);

  const size_t nVal = vmapLabels_.size();
  std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVal);
  std::vector<double> values(nVal);
  for ( size_t i=0, n=vmapLabels_.size(); i<n; ++i )
  {
    event.getByLabel(vmapLabels_[i], vmapHandles[i]);
  }

  loader_->init();
  std::auto_ptr<Cands> cands(new Cands);

  // Fill four vector informations
  for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
  {
    edm::Ref<edm::View<reco::Candidate> > candRef(srcHandle, i);
    reco::LeafCandidate cand(candRef->charge(), candRef->p4());
    cands->push_back(cand);
    for ( size_t j=0; j<nVal; ++j )
    {
      values[j] = (*vmapHandles[j])[candRef];
    }
    loader_->load(*candRef, values);
  }

  edm::OrphanHandle<Cands> outHandle = event.put(cands);
  for ( size_t i=0; i<loader_->N; ++i )
  {
    const std::string& varName = varNames_[i];
    std::auto_ptr<CandValueMap> vmap(new CandValueMap);
    CandValueMap::Filler filler(*vmap);
    filler.insert(outHandle, loader_->v[i].begin(), loader_->v[i].end());
    filler.fill();
    event.put(vmap, varName);
  }
}

DEFINE_FWK_MODULE(FlatCandProducer);

