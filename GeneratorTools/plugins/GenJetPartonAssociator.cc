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

using namespace std;

class GenJetPartonAssociator : public edm::EDProducer
{
public:
  GenJetPartonAssociator(const edm::ParameterSet& pset);
  ~GenJetPartonAssociator() {};

  void produce(edm::Event& event, const edm::EventSetup& eventSetup);

  bool hasMother(const reco::Candidate* p, const reco::Candidate* mother);

private:
  edm::InputTag genParticleLabel_;
  edm::InputTag genJetLabel_;

  std::set<unsigned int> pdgIdsToMatch_;
  unsigned int matchAlgo_;
  double cut_maxDR_, cut_maxDPt_;
  int cut_minNConstituent_;
  double cut_minFracConstituent_;

};

GenJetPartonAssociator::GenJetPartonAssociator(const edm::ParameterSet& pset)
{
  genParticleLabel_ = pset.getParameter<edm::InputTag>("genParticles");
  genJetLabel_ = pset.getParameter<edm::InputTag>("genJets");

  // Default values for matching condition - accept everything as a default
  matchAlgo_ = 0;
  cut_maxDR_ = cut_maxDPt_ = 1e9;
  cut_minNConstituent_ = 0; 
  cut_minFracConstituent_ = -999.; 

  edm::ParameterSet cuts = pset.getParameter<edm::ParameterSet>("cuts");
  const std::string matchAlgoName = pset.getParameter<std::string>("matchAlgo");
  if ( matchAlgoName == "deltaR" )
  {
    matchAlgo_ = 1;
    cut_maxDR_ = cuts.getParameter<double>("maxDR");
    if ( cuts.exists("maxDPt") ) cut_maxDPt_ = cuts.getParameter<double>("maxDPt");
  }
  else if ( matchAlgoName == "constituent" )
  {
    matchAlgo_ = 2;
    cut_minNConstituent_ = cuts.getParameter<int>("minNConstituent");
    cut_minFracConstituent_ = cuts.getParameter<double>("minFracConstituent");
  }

  // List of parton's PDG id's to perform matching
  // Give empty list to match all particles with status code = 3
  std::vector<unsigned int> pdgIdsToMatch = pset.getParameter<std::vector<unsigned int> >("pdgIdsToMatch");
  for ( auto& pdgIdToMatch : pdgIdsToMatch ) pdgIdsToMatch_.insert(pdgIdToMatch);

  produces<reco::GenJetToGenParticlesMap>();
}

void GenJetPartonAssociator::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<reco::GenJetToGenParticlesMap> genJetToGenParticlesMap(new reco::GenJetToGenParticlesMap);

  edm::Handle<std::vector<reco::GenJet> > genJetHandle;
  event.getByLabel(genJetLabel_, genJetHandle);

  edm::Handle<std::vector<reco::GenParticle> > genParticleHandle;
  event.getByLabel(genParticleLabel_, genParticleHandle);

  // Collect list of gen particles in hard process
  std::vector<unsigned int> genPartonIndicies;
  for ( unsigned int i=0, n=genParticleHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = genParticleHandle->at(i);
    if ( p.status() != 3 ) continue; // Keep "hard processes" only
    if ( p.numberOfMothers() == 0 ) continue; // ignore initial incident partons. 
    const reco::Candidate* mother = p.mother();
    if ( mother->numberOfDaughters() == 1 and mother->numberOfMothers() == 0 ) continue; // Skip if particle is simple copy of incident parton
    if ( !pdgIdsToMatch_.empty() and pdgIdsToMatch_.find(abs(p.pdgId())) == pdgIdsToMatch_.end() ) continue;

    genPartonIndicies.push_back(i);
  }
  const int nGenParton = genPartonIndicies.size();

  // Find matching between genJet to genParticle
  for ( unsigned int i=0, n=genJetHandle->size(); i<n; ++i )
  {
    auto& genJet = genJetHandle->at(i);
    std::set<unsigned int> matchedPartons;

    if ( matchAlgo_ == 1 ) // Algorithm 1 : deltaR (and deltaPt) matching
    {
      const reco::Candidate::LorentzVector& genJetP4 = genJet.p4();
      const double genJetPt = genJet.pt();

      for ( auto genPartonIndex : genPartonIndicies )
      {
        const reco::GenParticle* genParton = &genParticleHandle->at(genPartonIndex);
        if ( reco::deltaR(genJetP4, genParton->p4()) > cut_maxDR_ ) continue;
        if ( abs(genJetPt-genParton->pt()) > cut_maxDPt_ ) continue;

        matchedPartons.insert(genPartonIndex);
      }
    }
    else if ( matchAlgo_ == 2 ) // Algorithm 2 : Jet constituent overlap check
    {
      std::vector<const reco::GenParticle*> genConstituents = genJet.getGenConstituents();
      const int nConstituent = genConstituents.size();
      if ( nConstituent < cut_minNConstituent_ ) continue;

      for ( int k=0; k<nGenParton; ++k )
      {
        const reco::GenParticle* genParton = &genParticleHandle->at(genPartonIndicies[k]);

        int nMatched = 0;

        for ( auto& pp : genConstituents )
        {
          if ( !hasMother(pp, genParton) ) continue;

          ++nMatched;
        }

        if ( nMatched < cut_minNConstituent_ ) continue;
        if ( 1.*nMatched/nConstituent < cut_minFracConstituent_ ) continue;

        matchedPartons.insert(genPartonIndicies[k]);
      }
    }

    if ( matchedPartons.empty() ) continue;

    // Add mothers into the matched partons list
    std::set<unsigned int> matchedMothers;  
    for ( auto k : matchedPartons )
    {
      const reco::GenParticle* p = &genParticleHandle->at(k);
      const reco::GenParticle* m = p;
      while ( (m = dynamic_cast<const reco::GenParticle*>(m->mother()) ) != 0 )
      {
        if ( m->status() != 3 ) continue;
        if ( m->numberOfMothers() == 0 ) continue; // ignore initial incident partons. 
        const reco::Candidate* gm = m->mother();
        if ( gm->numberOfDaughters() == 1 and gm->numberOfMothers() == 0 ) break; // Skip if particle is simple copy of incident parton

        matchedMothers.insert(m-&genParticleHandle->at(0));
      }
    }
    matchedPartons.insert(matchedMothers.begin(), matchedMothers.end());

    edm::Ref<std::vector<reco::GenJet> > genJetRef(genJetHandle, i);
    for ( auto k : matchedPartons )
    {
      edm::Ref<std::vector<reco::GenParticle> > genParticleRef(genParticleHandle, k);
      genJetToGenParticlesMap->insert(genJetRef, genParticleRef);
    }
  }

  event.put(genJetToGenParticlesMap);
}

bool GenJetPartonAssociator::hasMother(const reco::Candidate* p, const reco::Candidate* mother)
{
  if ( !p or !mother ) return false;
  //if ( p == mother ) return false;
  const reco::Candidate* m = p->mother();
  if ( !m ) return false;

  //if ( m == mother ) return true;
  if ( m->status() == mother->status() and
       m->pdgId() == mother->pdgId() and
       m->p4() == mother->p4() ) return true;

  return hasMother(m, mother);
}

DEFINE_FWK_MODULE(GenJetPartonAssociator);

