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

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

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
  void produce(edm::Event& event, const edm::EventSetup& eventSetup);
 
private:
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;

  typedef std::vector<reco::LeafCandidate> Cands;
  typedef edm::ValueMap<double> CandValueMap;
  typedef edm::RefProd<Cands> CandRefProd;
  typedef edm::Ref<Cands> CandRef;

  edm::InputTag srcLabel_;
  std::vector<StringObjectFunction<reco::Candidate,true> > exprs_;
  std::vector<StringCutObjectSelector<reco::Candidate,true> > selectors_;
  std::vector<edm::InputTag> vmapLabels_;

  strings varNames_;
  std::vector<doubles> values_;

};

FlatCandProducer::FlatCandProducer(const edm::ParameterSet& pset)
{
  srcLabel_ = pset.getParameter<edm::InputTag>("src");
  edm::ParameterSet vars = pset.getParameter<edm::ParameterSet>("variables");

  const strings strVars = vars.getParameterNamesForType<std::string>();
  for ( auto& varName : strVars )
  {
    const string& varExpr = vars.getParameter<string>(varName);
    exprs_.push_back(StringObjectFunction<reco::Candidate,true>(varExpr));
    varNames_.push_back(varName);
  }
  const strings vmapNames = vars.getParameterNamesForType<edm::InputTag>();
  for ( auto& vmapName : vmapNames )
  {
    edm::InputTag vmapLabel = vars.getParameter<edm::InputTag>(vmapName);
    vmapLabels_.push_back(vmapLabel);
    varNames_.push_back(vmapName);
  }

  edm::ParameterSet sels = pset.getParameter<edm::ParameterSet>("selections");
  const strings strSels = sels.getParameterNamesForType<std::string>();
  for ( auto& selName : strSels )
  {
    const string& selection = sels.getParameter<string>(selName);
    selectors_.push_back(StringCutObjectSelector<reco::Candidate,true>(selection));
    varNames_.push_back(selName);
  }

  values_.resize(varNames_.size());

  produces<Cands>();
  for ( auto& varName : varNames_ )
  {
    produces<CandValueMap>(varName);
  }
}

void FlatCandProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::Candidate> > srcHandle;
  event.getByLabel(srcLabel_, srcHandle);

  const size_t nExpr = exprs_.size();
  const size_t nSele = selectors_.size();
  const size_t nVmap = vmapLabels_.size();
  std::vector<edm::Handle<edm::ValueMap<double> > > vmapHandles(nVmap);
  for ( size_t i=0; i<nVmap; ++i )
  {
    event.getByLabel(vmapLabels_[i], vmapHandles[i]);
  }

  std::auto_ptr<Cands> cands(new Cands);

  // Fill informations
  for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
  {
    edm::Ref<edm::View<reco::Candidate> > candRef(srcHandle, i);
    reco::LeafCandidate cand(candRef->charge(), candRef->p4());
    cand.setPdgId(candRef->pdgId());
    cands->push_back(cand);
    for ( size_t j=0; j<nExpr; ++j )
    {
      values_[j].push_back(exprs_[j](*candRef));
    }
    for ( size_t j=0; j<nVmap; ++j )
    {
      values_[j+nExpr].push_back((*vmapHandles[j])[candRef]);
    }
    for ( size_t j=0; j<nSele; ++j )
    {
      values_[j+nExpr+nVmap].push_back(selectors_[j](*candRef));
    }
  }

  edm::OrphanHandle<Cands> outHandle = event.put(cands);
  for ( size_t i=0, n=nExpr+nVmap+nSele; i<n; ++i )
  {
    std::auto_ptr<CandValueMap> vmap(new CandValueMap);
    CandValueMap::Filler filler(*vmap);
    filler.insert(outHandle, values_[i].begin(), values_[i].end());
    filler.fill();
    values_[i].clear();

    const std::string& varName = varNames_[i];
    event.put(vmap, varName);
  }
}

DEFINE_FWK_MODULE(FlatCandProducer);

