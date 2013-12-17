#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"

#include "TMath.h"

#include <memory>
#include <vector>
#include <set>

using namespace std;
using namespace edm;
using namespace reco;

// Produces ghost GenParticles (mainly for B hadrons) to be entered GenJet clustering
// This method is standard in pseudo-top
class GhostGenParticleProducer : public edm::EDProducer
{
public:
  GhostGenParticleProducer(const edm::ParameterSet& pset);
  ~GhostGenParticleProducer() {};
  void produce(edm::Event& event, const edm::EventSetup& eventSetup);

  typedef edm::AssociationMap<edm::OneToOne<reco::GenParticleCollection, reco::GenParticleCollection> > GenParticleToGenParticleMap;

private:
  bool isBHadron(const reco::Candidate* p) const; // B hadron without b hadron decendents
  bool isBHadron(const unsigned int absPdgId) const; // Check |pdgId| to check flavor
  bool isPartonLevel(const reco::Candidate* p) const; // True if particle is parton level and hadronized

private:
  edm::InputTag srcLabel_;
  int userPdgId_;
  bool doPartonLevel_;
  std::set<unsigned int> partonIdsToVeto_;

};

GhostGenParticleProducer::GhostGenParticleProducer(const edm::ParameterSet& pset)
{
  srcLabel_ = pset.getParameter<edm::InputTag>("src");
  userPdgId_ = pset.getUntrackedParameter<int>("userPdgId", 0);
  doPartonLevel_ = pset.getParameter<bool>("doPartonLevel");
  if ( doPartonLevel_ )
  {
    std::vector<unsigned int> partonIdsToVeto = pset.getParameter<std::vector<unsigned int> >("partonIdsToVeto");
    for ( int i=0, n=partonIdsToVeto.size(); i<n; ++i )
    {
      partonIdsToVeto_.insert(partonIdsToVeto[i]);
    }
  }

  produces<reco::GenParticleCollection>();
  produces<std::vector<int> >().setBranchAlias( pset.getParameter<std::string>("@module_label") + "BarCodes" );
  produces<GenParticleToGenParticleMap>(); // GhostParticle to source particle map, to track up original history
}

void GhostGenParticleProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  event.getByLabel(srcLabel_, genParticlesHandle);

  // Prepare output collections
  std::auto_ptr<reco::GenParticleCollection> outGenParticles(new reco::GenParticleCollection());
  std::auto_ptr<std::vector<int> > outGenParticleIndices(new std::vector<int>());
  reco::GenParticleRefProd ref = event.getRefBeforePut<reco::GenParticleCollection>();
  std::auto_ptr<GenParticleToGenParticleMap> ghostToSrcMap(new GenParticleToGenParticleMap);

  // Copy GenParticles and ghost B hadrons to output collection
  const double scale = 1e-20;
  std::map<const reco::Candidate*, int> particleIndexMap;
  std::vector<int> srcParticleIndices;
  const int nParticles = genParticlesHandle->size();
  for ( int i=0, nOutGenParticle=0; i<nParticles; ++i )
  {
    const reco::GenParticle& pSrc = genParticlesHandle->at(i);
    particleIndexMap[&pSrc] = nOutGenParticle;

    reco::GenParticle p = pSrc;
    p.clearMothers();
    p.clearDaughters();
    p.resetMothers(ref.id());
    p.resetDaughters(ref.id());

    outGenParticleIndices->push_back(nOutGenParticle);
    outGenParticles->push_back(p);
    srcParticleIndices.push_back(i);

    ++nOutGenParticle;

    if ( isBHadron(&pSrc) or ( doPartonLevel_ and isPartonLevel(&pSrc) ) )
    {
      reco::GenParticle b = p;

      // Scale down particle momentum by 1e-20 : neglegible contribution to the jet momentum
      // Algorithm 2 : Particle energy is sqrt((1e-20*p)^2+m^2) rather than simple 1e-20*E to reserve particle mass
      // this gives wrong eta distribution for some reason...
      //const reco::Candidate::LorentzVector p4(p.px()*scale, p.py()*scale, p.pz()*scale, std::hypot(p.p()*scale, p.mass()));
      //b.setP4(p4);
      b.setP4(p.p4()*scale);
      b.setStatus(1); // This ghost particle will be considered as a stable particle
      if ( userPdgId_ != 0 ) b.setPdgId(userPdgId_); // Set user defined PDG ID if set

      b.clearDaughters();
      outGenParticleIndices->push_back(nOutGenParticle);
      outGenParticles->push_back(b);
      srcParticleIndices.push_back(i);

      ++nOutGenParticle;

      // Build ghost particle to original particle map
      reco::GenParticleRef ghostRef(ref, nOutGenParticle-1);
      reco::GenParticleRef srcRef(ref, nOutGenParticle-2);
      ghostToSrcMap->insert(ghostRef, srcRef);
    }
  }

  // Make connection between candidate
  for ( int i=0, n=outGenParticles->size(); i<n; ++i )
  {
    reco::GenParticle& p = outGenParticles->at(i);

    const reco::GenParticle& pSrc = genParticlesHandle->at(srcParticleIndices[i]);
    for ( int j=0, m=pSrc.numberOfMothers(); j<m; ++j )
    {
      const reco::Candidate* mother = pSrc.mother(j);
      if ( particleIndexMap.find(mother) == particleIndexMap.end() )
      {
        cout << "Particle index is corrupted\n";
        continue;
      }
      const unsigned int motherIndex = particleIndexMap[mother];
      p.addMother(reco::GenParticleRef(ref, motherIndex));
      reco::GenParticle& outMother = outGenParticles->at(motherIndex);
      outMother.addDaughter(reco::GenParticleRef(ref, i));
    }
  }

  event.put(outGenParticles);
  event.put(outGenParticleIndices);
  event.put(ghostToSrcMap);
}

bool GhostGenParticleProducer::isBHadron(const reco::Candidate* p) const
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

bool GhostGenParticleProducer::isBHadron(const unsigned int absPdgId) const
{
  if ( doPartonLevel_ and absPdgId == 5 ) return true; // select b quarks if "doPartonLevel" is set

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

bool GhostGenParticleProducer::isPartonLevel(const reco::Candidate* p) const
{
  if ( p->status() != 3 ) return false;
  for ( int i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    if ( p->status() == 3 ) return false;
    if ( partonIdsToVeto_.find(p->pdgId()) != partonIdsToVeto_.end() ) return false;
  }
  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GhostGenParticleProducer);
