#ifndef KrAFT_GeneratorTools_Types_H
#define KrAFT_GeneratorTools_Types_H

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToOne.h"

namespace pat
{
  typedef edm::AssociationMap<edm::OneToOne<std::vector<pat::Jet>, std::vector<reco::GenJet> > > RecoToGenJetMap;
}

namespace reco
{
  typedef edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>, reco::GenParticleCollection> > GenJetToGenParticlesMap;
  typedef edm::AssociationMap<edm::OneToOne<reco::GenParticleCollection, reco::GenParticleCollection> > GenParticleToGenParticleMap;
}

#endif

