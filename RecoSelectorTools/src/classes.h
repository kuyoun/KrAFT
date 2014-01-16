#include "KrAFT/RecoSelectorTools/interface/Types.h"

namespace pat {
  edm::RefProd<std::vector<pat::Jet> > dummy00;

  // pat::Jet -> double mapping
  edm::Wrapper<edm::AssociationMap<edm::OneToValue<std::vector<pat::Jet>, double, unsigned int> > > dummy10;
  JetToValue dummy11;
}

