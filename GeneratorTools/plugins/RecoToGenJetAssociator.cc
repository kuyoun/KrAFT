#ifndef KrAFT_GeneratorTools_RecoToGenJetAssociator_H
#define KrAFT_GeneratorTools_RecoToGenJetAssociator_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "KrAFT/GeneratorTools/interface/Types.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <vector>
#include <string>

class RecoToGenJetAssociator : public edm::EDProducer
{
public:
  RecoToGenJetAssociator(const edm::ParameterSet& pset);
  ~RecoToGenJetAssociator() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::EDGetTokenT<std::vector<pat::Jet> > recoJetToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetToken_;
  double cut_maxDR_, cut_maxDPt_;

};

RecoToGenJetAssociator::RecoToGenJetAssociator(const edm::ParameterSet& pset)
{
  recoJetToken_ = consumes<std::vector<pat::Jet> >(pset.getParameter<edm::InputTag>("recoJets"));
  genJetToken_ = consumes<std::vector<reco::GenJet> >(pset.getParameter<edm::InputTag>("genJets"));

  cut_maxDR_ = cut_maxDPt_ = 1e9;
  edm::ParameterSet cuts = pset.getParameter<edm::ParameterSet>("cuts");
  if ( cuts.exists("maxDR") ) cut_maxDR_ = cuts.getParameter<double>("maxDR");
  if ( cuts.exists("maxDPt") ) cut_maxDPt_ = cuts.getParameter<double>("maxDPt");

  produces<pat::RecoToGenJetMap>();
}

void RecoToGenJetAssociator::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<pat::RecoToGenJetMap> recoToGenJetMap(new pat::RecoToGenJetMap);

  edm::Handle<std::vector<pat::Jet> > recoJetHandle;
  event.getByToken(recoJetToken_, recoJetHandle);

  edm::Handle<std::vector<reco::GenJet> > genJetHandle;
  event.getByToken(genJetToken_, genJetHandle);

  for ( int i=0, n=recoJetHandle->size(); i<n; ++i )
  {
    const pat::Jet& recoJet = recoJetHandle->at(i);
    int matchedJetIndex = -1;
    double maxDR = cut_maxDR_;
    double maxDPt = cut_maxDPt_;

    // Find best pair with dR or dPt cut
    for ( int j=0, m=genJetHandle->size(); j<m; ++j )
    {
      const reco::GenJet& genJet = genJetHandle->at(j);

      const double dR = deltaR(recoJet, genJet);
      const double dPt = abs(recoJet.pt() - genJet.pt());
      if ( dR > maxDR and cut_maxDR_ < 1e9 ) continue;
      if ( dPt > maxDPt and cut_maxDPt_ < 1e9 ) continue;

      maxDR = dR;
      maxDPt = dPt;
      matchedJetIndex = j;
    }

    if ( matchedJetIndex == -1 ) continue;
      
    // Now we have best matching reco->gen jet pair
    edm::Ref<std::vector<pat::Jet> > recoJetRef(recoJetHandle, i);
    edm::Ref<std::vector<reco::GenJet> > genJetRef(genJetHandle, matchedJetIndex);
    recoToGenJetMap->insert(recoJetRef, genJetRef);
  }

  event.put(recoToGenJetMap);
}

DEFINE_FWK_MODULE(RecoToGenJetAssociator);

#endif

