#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "KrAFT/GeneratorTools/interface/Types.h"

namespace pat {
  edm::RefProd<std::vector<pat::Jet> > dummy00;

  // pat::Jet -> reco::GenJet mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToOne<std::vector<pat::Jet>, std::vector<reco::GenJet> ,unsigned int> > > dummy10;
  RecoToGenJetMap dummy11;
}

namespace reco {
  edm::RefProd<std::vector<reco::GenJet> > dummy01;

  // reco::GenJet -> reco::GenParticle mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>,std::vector<reco::GenParticle>,unsigned int> > > dummy12;

}
