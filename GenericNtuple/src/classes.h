#include "KrAFT/GenericNtuple/interface/KLeptonJetTreeReducer.h"
#include "KrAFT/GenericNtuple/interface/KDileptonTreeReducer.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
  struct dictionary {
    KDileptonTreeReducer kDileptonTreeReducer;
    KLeptonJetTreeReducer kLeptonJetTreeReducer;
  };

}
