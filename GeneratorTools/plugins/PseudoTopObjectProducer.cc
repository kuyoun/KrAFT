#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "KrAFT/GeneratorTools/interface/Types.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include "RecoJets/JetProducers/interface/JetSpecific.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TMath.h"

#include <memory>
#include <vector>
#include <set>

using namespace std;
using namespace edm;
using namespace reco;

class PseudoTopObjectProducer : public edm::EDProducer
{
public:
  PseudoTopObjectProducer(const edm::ParameterSet& pset);
  virtual ~PseudoTopObjectProducer() {};
  virtual void produce(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isFromHadron(const reco::Candidate& p, const int depth=100) const;
  bool isBHadron(const reco::Candidate* p) const;
  bool isBHadron(const unsigned int pdgId) const;

private:
  edm::InputTag srcLabel_;
  double leptonMinPt_, jetMinPt_;

  typedef fastjet::JetDefinition JetDef;
  boost::shared_ptr<JetDef> fjLepDef_, fjJetDef_;
  reco::Particle::Point genVertex_;
};

PseudoTopObjectProducer::PseudoTopObjectProducer(const edm::ParameterSet& pset)
{
  srcLabel_ = pset.getParameter<edm::InputTag>("src");
  leptonMinPt_ = pset.getParameter<double>("leptonMinPt");
  jetMinPt_ = pset.getParameter<double>("jetMinPt");

  const double leptonConeSize = pset.getParameter<double>("leptonConeSize");
  const double jetConeSize = pset.getParameter<double>("jetConeSize");
  fjLepDef_ = boost::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, leptonConeSize));
  fjJetDef_ = boost::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, jetConeSize));

  genVertex_ = reco::Particle::Point(0,0,0);

  produces<reco::GenParticleCollection>("neutrinos");
  produces<reco::GenJetCollection>("leptons");
  produces<reco::GenJetCollection>("jets");
}

void PseudoTopObjectProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::GenParticle> > srcHandle;
  event.getByLabel(srcLabel_, srcHandle);

  std::auto_ptr<reco::GenParticleCollection> neutrinos(new reco::GenParticleCollection);
  std::auto_ptr<reco::GenJetCollection> leptons(new reco::GenJetCollection);
  std::auto_ptr<reco::GenJetCollection> jets(new reco::GenJetCollection);

  // Collect stable leptons and neutrinos
  size_t nStables = 0;
  std::vector<size_t> leptonIdxs;
  std::set<size_t> neutrinoIdxs, bHadronIdxs;
  for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = srcHandle->at(i);
    const int absPdgId = abs(p.pdgId());
    if ( p.status() == 1 )
    {
      ++nStables;
      if ( isFromHadron(p) ) continue;
      switch ( absPdgId )
      {
        case 11: case 13: // Leptons
        case 22: // Photons
          leptonIdxs.push_back(i);
          break;
        case 12: case 14: case 16:
          neutrinoIdxs.insert(i);
          neutrinos->push_back(p);
          break;
      }
    }
    else
    {
      // Collect B-hadrons, to be used in b tagging
      if ( isBHadron(&p) ) bHadronIdxs.insert(i);
    }
  }

  // Sort neutrinos by pT.
  std::sort(neutrinos->begin(), neutrinos->end(), GreaterByPt<reco::GenParticle>());

  // Make dressed leptons with anti-kt(0.1) algorithm
  //// Prepare input particle list
  std::vector<fastjet::PseudoJet> fjLepInputs;
  fjLepInputs.reserve(leptonIdxs.size());
  for ( auto index : leptonIdxs )
  {
    const reco::GenParticle& p = srcHandle->at(index);
    if ( std::isnan(p.pt()) or p.pt() <= 0 ) continue;

    fjLepInputs.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.energy()));
    fjLepInputs.back().set_user_index(index);
  }

  //// Run the jet algorithm
  fastjet::ClusterSequence fjLepClusterSeq(fjLepInputs, *fjLepDef_);
  std::vector<fastjet::PseudoJet> fjLepJets = fastjet::sorted_by_pt(fjLepClusterSeq.inclusive_jets(leptonMinPt_));

  //// Build dressed lepton objects from the FJ output
  leptons->reserve(fjLepJets.size());
  std::set<size_t> lepDauIdxs; // keep lepton constituents to remove from GenJet construction
  for ( auto& fjJet : fjLepJets )
  {
    // Get jet constituents from fastJet
    const std::vector<fastjet::PseudoJet> fjConstituents = fastjet::sorted_by_pt(fjJet.constituents());
    // Convert to CandidatePtr
    std::vector<reco::CandidatePtr> constituents;
    reco::CandidatePtr lepCand;
    for ( auto& fjConstituent : fjConstituents )
    {
      const size_t index = fjConstituent.user_index();
      reco::CandidatePtr cand = srcHandle->ptrAt(index);
      const int absPdgId = abs(cand->pdgId());
      if ( absPdgId == 11 or absPdgId == 13 ) lepCand = cand;;
      constituents.push_back(cand);
    }
    if ( lepCand.isNull() ) continue;

    const reco::Particle::LorentzVector jetP4(fjJet.px(), fjJet.py(), fjJet.pz(), fjJet.E());
    reco::GenJet lepJet;
    reco::writeSpecific(lepJet, jetP4, genVertex_, constituents, eventSetup);

    lepJet.setPdgId(lepCand->pdgId());
    lepJet.setCharge(lepCand->charge());

    const double jetArea = fjJet.has_area() ? fjJet.area() : 0;
    lepJet.setJetArea(jetArea);

    leptons->push_back(lepJet);

    // Keep constituent indices to be used in the next step.
    for ( auto& fjConstituent : fjConstituents )
    {
      lepDauIdxs.insert(fjConstituent.user_index());
    }
  }

  // Now proceed to jets.
  // Jets: anti-kt excluding the e, mu, nu, and photons in selected leptons. 
  //// Prepare input particle list. Remove particles used in lepton clusters, neutrinos
  std::vector<fastjet::PseudoJet> fjJetInputs;
  fjJetInputs.reserve(nStables);
  for ( size_t i=0, n=srcHandle->size(); i<n; ++i )
  {
    const reco::GenParticle& p = srcHandle->at(i);
    if ( p.status() != 1 ) continue;
    if ( std::isnan(p.pt()) or p.pt() <= 0 ) continue;

    if ( neutrinoIdxs.find(i) != neutrinoIdxs.end() ) continue;
    if ( lepDauIdxs.find(i) != lepDauIdxs.end() ) continue;

    fjJetInputs.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.energy()));
    fjJetInputs.back().set_user_index(i);
  }
  //// Also don't forget to put B hadrons
  for ( auto index : bHadronIdxs )
  {
    const reco::GenParticle& p = srcHandle->at(index);
    if ( std::isnan(p.pt()) or p.pt() <= 0 ) continue;

    const double scale = 1e-20/p.p();
    fjJetInputs.push_back(fastjet::PseudoJet(p.px()*scale, p.py()*scale, p.pz()*scale, p.energy()*scale));
    fjJetInputs.back().set_user_index(index);
  }

  //// Run the jet algorithm
  fastjet::ClusterSequence fjJetClusterSeq(fjJetInputs, *fjJetDef_);
  std::vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_pt(fjJetClusterSeq.inclusive_jets(jetMinPt_));

  /// Build jets
  jets->reserve(fjJets.size());
  for ( size_t i=0, n=fjJets.size(); i<n; ++i )
  {
    const fastjet::PseudoJet& fjJet = fjJets[i];

    // Get jet constituents from fastJet
    const std::vector<fastjet::PseudoJet> fjConstituents = fastjet::sorted_by_pt(fjJet.constituents());
    // Convert to CandidatePtr
    std::vector<reco::CandidatePtr> constituents;
    bool hasBHadron = false;
    for ( size_t j=0, m=fjConstituents.size(); j<m; ++j )
    { 
      const size_t index = fjConstituents[j].user_index();
      if ( bHadronIdxs.find(index) != bHadronIdxs.end() ) hasBHadron = true;
      reco::CandidatePtr cand = srcHandle->ptrAt(index);
      constituents.push_back(cand);
    }
    
    const reco::Particle::LorentzVector jetP4(fjJet.px(), fjJet.py(), fjJet.pz(), fjJet.E());
    reco::GenJet genJet;
    reco::writeSpecific(genJet, jetP4, genVertex_, constituents, eventSetup);

    const double jetArea = fjJet.has_area() ? fjJet.area() : 0;
    genJet.setJetArea(jetArea);
    if ( hasBHadron ) genJet.setPdgId(5);

    jets->push_back(genJet);
  }

  event.put(neutrinos, "neutrinos");
  event.put(leptons, "leptons");
  event.put(jets, "jets");
}

bool PseudoTopObjectProducer::isFromHadron(const reco::Candidate& p, const int depth) const
{
  if ( depth <= 0 ) return false;

  for ( int i=0, n=p.numberOfMothers(); i<n; ++i )
  {
    const reco::Candidate& mother = *p.mother(i);
    const int pdgId = abs(mother.pdgId());

    if ( pdgId == 23 or pdgId == 24 or pdgId == 35 ) return false;
    else if ( pdgId == 22 or pdgId < 6 ) return true;
    else if ( pdgId > 100 ) return true;
    else if ( isFromHadron(mother, depth-1) ) return true;
  }
  return false;
}

bool PseudoTopObjectProducer::isBHadron(const reco::Candidate* p) const
{
  const unsigned int absPdgId = abs(p->pdgId());
  if ( !isBHadron(absPdgId) ) return false;

  // Do not consider this particle if it has B hadron daughter
  // For example, B* -> B0 + photon; then we drop B* and take B0 only
  for ( int i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( isBHadron(abs(dau->pdgId())) ) return false;
  }

  return true;
}

bool PseudoTopObjectProducer::isBHadron(const unsigned int absPdgId) const
{
  if ( absPdgId <= 100 ) return false; // Fundamental particles and MC internals
  if ( absPdgId >= 1000000000 ) return false; // Nuclei, +-10LZZZAAAI

  // General form of PDG ID is 7 digit form
  // +- n nr nL nq1 nq2 nq3 nJ
  //const int nJ = absPdgId % 10; // Spin
  const int nq3 = (absPdgId / 10) % 10;
  const int nq2 = (absPdgId / 100) % 10;
  const int nq1 = (absPdgId / 1000) % 10;

  if ( nq3 == 0 ) return false; // Diquarks
  if ( nq1 == 0 and nq2 == 5 ) return true; // B mesons
  if ( nq1 == 5 ) return true; // B baryons

  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PseudoTopObjectProducer);

