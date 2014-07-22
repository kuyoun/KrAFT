#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

class GenParticleCountFilter : public  edm::EDFilter
{ 
public:
  GenParticleCountFilter(const edm::ParameterSet& pset);
  ~GenParticleCountFilter() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
  
private:
  edm::EDGetTokenT<reco::GenParticleCollection> srcToken_;
  unsigned int minNumber_, maxNumber_;

  StringCutObjectSelector<reco::GenParticle>* select_;
};

GenParticleCountFilter::GenParticleCountFilter(const edm::ParameterSet& pset)
{
  srcToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("src"));
  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");
  std::string cut = pset.getParameter<std::string>("cut");
  select_ = new StringCutObjectSelector<reco::GenParticle>(cut);
}

bool GenParticleCountFilter::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  event.getByToken(srcToken_, genParticleHandle);

  unsigned int nMatch = 0;
  for ( int i=0, n=genParticleHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = genParticleHandle->at(i);

    if ( !(*select_)(p) ) continue;

    ++nMatch;

  }
  if ( nMatch < minNumber_ or nMatch > maxNumber_ ) return false;

  return true;
}

DEFINE_FWK_MODULE(GenParticleCountFilter);
