#include "KrAFT/GenericNtuple/interface/KLeptonJetTreeAnalyzer.h"
#include "KrAFT/GenericNtuple/interface/KDileptonTreeAnalyzer.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
  struct dictionary {
    KDileptonTreeAnalyzer kDileptonTreeAnalyzer;
    KLeptonJetTreeAnalyzer kLeptonJetTreeAnalyzer;
  };

}
