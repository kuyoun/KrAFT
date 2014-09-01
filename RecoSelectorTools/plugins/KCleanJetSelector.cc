#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/View.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
//#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "KrAFT/RecoSelectorTools/interface/Types.h"

#include <memory>

using namespace edm;
using namespace std;

class KCleanJetSelector : public edm::EDFilter
{
public:
  KCleanJetSelector(const edm::ParameterSet& pset);
  ~KCleanJetSelector() {};

  typedef std::vector<pat::Jet> Jets;
  typedef std::vector<pat::MET> METs;
  typedef std::vector<std::string> strings;

private:
  void beginJob() {};
  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
  void endJob() {};

  bool doClean_;
  double overlapDeltaR_;
  unsigned int minNumber_;
  unsigned int maxNumber_;

  strings jesLevels_;
  int iJER_, iJERDn_, iJERUp_, nJES_;
  edm::InputTag jetLabel_;
  std::vector<edm::InputTag> overlapCandLabels_;

  PFJetIDSelectionFunctor* isGoodJet_;
  double minPt_, maxEta_;

private:
};

KCleanJetSelector::KCleanJetSelector(const edm::ParameterSet& pset)
{
  jetLabel_ = pset.getParameter<edm::InputTag>("jet");
  jesLevels_ = pset.getParameter<strings>("jesLevels");

  // Selection cuts
  edm::ParameterSet selectionPSet = pset.getParameter<edm::ParameterSet>("selection");
  isGoodJet_ = new PFJetIDSelectionFunctor(selectionPSet.getParameter<edm::ParameterSet>("jetId"));
  minPt_ = selectionPSet.getParameter<double>("minPt");
  maxEta_ = selectionPSet.getParameter<double>("maxEta");

  // Cleaning
  edm::ParameterSet cleanPSet = pset.getParameter<edm::ParameterSet>("cleaning");
  doClean_ = cleanPSet.getParameter<bool>("doClean");
  if ( doClean_ )
  {
    overlapDeltaR_ = cleanPSet.getParameter<double>("overlapDeltaR");
    overlapCandLabels_ = cleanPSet.getParameter<std::vector<edm::InputTag> >("overlapCands");
    if ( overlapCandLabels_.empty() ) doClean_ = false;
  }

  minNumber_ = pset.getParameter<unsigned int>("minNumber");
  maxNumber_ = pset.getParameter<unsigned int>("maxNumber");

  produces<std::vector<pat::Jet> >();
  nJES_ = jesLevels_.size();
  for ( int i=0; i<nJES_; ++i )
  {
    auto& jesLevel = jesLevels_[i];
    if ( jesLevel.find("resDn") != string::npos ) iJERDn_ = i;
    else if ( jesLevel.find("resUp") != string::npos ) iJERUp_ = i;
    else if ( jesLevel.find("res") != string::npos ) iJER_ = i;
  }
}

bool KCleanJetSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<pat::Jet> > jetHandle;
  event.getByLabel(jetLabel_, jetHandle);

  // Collect jes factors
  std::vector<edm::Handle<pat::JetToValue> > jesHandles;
  std::vector<std::vector<double> > fJECs;
  for ( auto& jesLevel : jesLevels_ )
  {
    edm::Handle<pat::JetToValue> jesHandle;
    event.getByLabel(edm::InputTag(jetLabel_.label(), jesLevel), jesHandle);
    jesHandles.push_back(jesHandle);
    fJECs.push_back(std::vector<double>());
  }

  // Collect candidates for overlap removal
  std::vector<const reco::Candidate*> overlapCands;
  if ( doClean_ )
  {
    for ( auto& overlapCandLabel : overlapCandLabels_ )
    {
      edm::Handle<edm::View<reco::Candidate> > overlapCandHandle;
      event.getByLabel(overlapCandLabel, overlapCandHandle);

      //if ( !overlapCandHandle.isValid() ) continue;

      for ( auto& overlapCand : *overlapCandHandle )
      {
        overlapCands.push_back(&overlapCand);
      }
    }
  }

  // Now start to build jet collection
  std::auto_ptr<std::vector<pat::Jet> > cleanJets(new std::vector<pat::Jet>());

  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    auto& jet = jetHandle->at(i);
    if ( !(*isGoodJet_)(jet) ) continue;
    const reco::Candidate::LorentzVector jetP4 = jet.p4();
    const double jetPt = jetP4.pt();
    const double jetEta = jetP4.eta();

    if ( abs(jetEta) > 5 ) continue;

    // Check overlap
    if ( doClean_ )
    {
      bool isOverlap = false;
      for ( auto overlapCand : overlapCands )
      {
        if ( deltaR(jet.p4(), overlapCand->p4()) < overlapDeltaR_ )
        {
          isOverlap = true;
        }
      }
      if ( isOverlap ) continue;
    }

    // Copy JES factors
    double fJECDn = 1, fJERDn = 1;
    for ( int j=0; j<nJES_; ++j )
    {
      const double jes = jesHandles[j]->operator[](jetHandle->refAt(i));
      if      ( j == iJERDn_ ) fJERDn = jes;
      else if ( j != iJER_ and j != iJERUp_ ) fJECDn = std::min(fJECDn, jes);
    }

    // Check acceptance
    if ( std::abs(jetEta) > maxEta_ ) continue;
    if ( jetPt*fJECDn*fJERDn < minPt_ ) continue;

    cleanJets->push_back(jet);
    for ( int j=0; j<nJES_; ++j )
    {
      const double jes = jesHandles[j]->operator[](jetHandle->refAt(i));
      fJECs[j].push_back(jes);
    }
  }

  const unsigned int nCleanJet = cleanJets->size();
  edm::OrphanHandle<pat::JetCollection> outColl = event.put(cleanJets);

  for ( int i=0; i<nJES_; ++i )
  {
    std::auto_ptr<pat::JetToValue> fJECMap(new pat::JetToValue);
    pat::JetToValue::Filler fillJEC(*fJECMap);
    fillJEC.insert(outColl, fJECs[i].begin(), fJECs[i].end());
    fillJEC.fill();
    event.put(fJECMap, jesLevels_[i]);
  }

  if ( nCleanJet < minNumber_ or nCleanJet > maxNumber_ ) return false;

  return true;
}

DEFINE_FWK_MODULE(KCleanJetSelector);

