#include "KrAFT/FlatTree/interface/FlatEvent.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
  struct dictionary {
    FlatEvent flatEvent;
    edm::Wrapper<FlatEvent> w_flatEvent;
  };

}
