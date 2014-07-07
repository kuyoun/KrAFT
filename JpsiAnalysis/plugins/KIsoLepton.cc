#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include <memory>
#include <vector>
#include <string>

template<typename Lepton>
class KIsoLeptonSelector : public edm::EDFilter
{
public:
  KIsoLeptonSelector(const edm::ParameterSet& pset);
  ~KIsoLeptonSelector() {};

  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
private:
bool isIsolation(const pat::Electron& electron);
bool isIsolation(const pat::Muon& muon);

private:
  edm::InputTag leptonLabel_;
  unsigned int minNumber_, maxNumber_;
	
};

template<typename Lepton>
KIsoLeptonSelector<Lepton>::KIsoLeptonSelector(const edm::ParameterSet& pset)
{
  // isMC_ = true as a default, determine it in the event loop.
  // Changing this default value of this flag should be done carefully,
  // modification can be needed in #ISMC_DEPENDENT_PART#
  leptonLabel_ = pset.getParameter<edm::InputTag>("src");
  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");
  produces<std::vector<Lepton> >();
}

template<typename Lepton>
bool KIsoLeptonSelector<Lepton>::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{

  edm::Handle<edm::View<Lepton> > leptonHandle;
  event.getByLabel(leptonLabel_, leptonHandle);


  std::auto_ptr<std::vector<Lepton> > selectedLeptons(new std::vector<Lepton>());

  for ( int i=0, n=leptonHandle->size(); i<n; ++i )
  {
    const Lepton& srcLepton = leptonHandle->at(i);
    // Build lepton
    Lepton lepton(srcLepton);
	  if ( !isIsolation( lepton) ) continue; 
    selectedLeptons->push_back(lepton);
  }
  const unsigned int nLepton = selectedLeptons->size();
  event.put(selectedLeptons);
  return (nLepton >= minNumber_ and nLepton <= maxNumber_);
}

template<typename Lepton>
bool KIsoLeptonSelector<Lepton>::isIsolation(const pat::Electron& electron)
{
	if ( electron.userIso(2) < 0.15 ) return true;
	else return false; 
}
template<typename Lepton>
bool KIsoLeptonSelector<Lepton>::isIsolation(const pat::Muon& muon)
{
	if ( muon.userIso(1) < 0.15 ) return true;
	else return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"

typedef KIsoLeptonSelector<pat::Muon> KIsoMuonSelector;
typedef KIsoLeptonSelector<pat::Electron> KIsoElectronSelector;

DEFINE_FWK_MODULE(KIsoMuonSelector);
DEFINE_FWK_MODULE(KIsoElectronSelector);
