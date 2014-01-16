#ifndef KrAFT_RecoSelectorTools_Types_H
#define KrAFT_RecoSelectorTools_Types_H

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToValue.h"

namespace pat
{
  typedef edm::AssociationMap<edm::OneToValue<std::vector<pat::Jet>, double> > JetToValue;
}

#endif

