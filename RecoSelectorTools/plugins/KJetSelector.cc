#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/View.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
//#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "KrAFT/GeneratorTools/interface/Types.h"

#include <memory>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH1F.h>
#include <TH2F.h>

using namespace edm;
using namespace std;

class KJetSelector : public edm::EDFilter
{
public:
  KJetSelector(const edm::ParameterSet& pset);
  ~KJetSelector() {};

  typedef std::vector<pat::Jet> Jets;
  typedef edm::AssociationMap<edm::OneToValue<Jets, double> > JetToValueMap;

private:
  void beginJob() {};
  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
  void endJob() {};

  bool doClean_;
  double overlapDeltaR_;
  unsigned int minNumber_;
  unsigned int maxNumber_;

  edm::InputTag jetLabel_;
  edm::InputTag jetUncLabel_;
  std::vector<edm::InputTag> overlapCandLabels_;

  PFJetIDSelectionFunctor* isGoodJet_;
  double minPt_, maxEta_;

private:

};

KJetSelector::KJetSelector(const edm::ParameterSet& pset)
{
  jetLabel_ = pset.getParameter<edm::InputTag>("jet");
  jetUncLabel_ = pset.getParameter<edm::InputTag>("jetUnc");

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

  //produces<std::vector<pat::Jet> >();
  produces<edm::RefVector<Jets> >();
}

bool KJetSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<Jets> jetHandle;
  event.getByLabel(jetLabel_, jetHandle);

  std::vector<const JetToValueMap*> jetScales;
  edm::Handle<JetToValueMap> jetScaleHandle;
  const std::string jetUncLabel = jetUncLabel_.label();
  if ( event.getByLabel(edm::InputTag(jetUncLabel, "up"), jetScaleHandle) ) jetScales.push_back(&*jetScaleHandle);
  if ( event.getByLabel(edm::InputTag(jetUncLabel, "dn"), jetScaleHandle) ) jetScales.push_back(&*jetScaleHandle);
  if ( event.getByLabel(edm::InputTag(jetUncLabel, "res"), jetScaleHandle) ) jetScales.push_back(&*jetScaleHandle);
  if ( event.getByLabel(edm::InputTag(jetUncLabel, "resUp"), jetScaleHandle) ) jetScales.push_back(&*jetScaleHandle);
  if ( event.getByLabel(edm::InputTag(jetUncLabel, "resDn"), jetScaleHandle) ) jetScales.push_back(&*jetScaleHandle);

  //std::auto_ptr<std::vector<pat::Jet> > cleanJets(new std::vector<pat::Jet>());
  std::auto_ptr<edm::RefVector<Jets> > cleanJets(new edm::RefVector<Jets>());

  std::vector<const reco::Candidate*> overlapCands;
  if ( doClean_ )
  {
    for ( int iLabel=0, nLabel=overlapCandLabels_.size(); iLabel<nLabel; ++iLabel )
    {
      edm::Handle<edm::View<reco::Candidate> > overlapCandHandle;
      event.getByLabel(overlapCandLabels_.at(iLabel), overlapCandHandle);

      //if ( !overlapCandHandle.isValid() ) continue;

      for ( int i=0, n=overlapCandHandle->size(); i<n; ++i )
      {
        overlapCands.push_back(&(overlapCandHandle->at(i)));
      }
    }
  }

  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    pat::Jet jet = jetHandle->at(i);
    if ( !(*isGoodJet_)(jet) ) continue;
    const reco::Candidate::LorentzVector jetP4 = jet.p4();

    // Check overlap
    if ( doClean_ )
    {
      bool isOverlap = false;
      for ( int j=0, m=overlapCands.size(); j<m; ++j )
      {
        if ( deltaR(jet.p4(), overlapCands.at(j)->p4()) < overlapDeltaR_ )
        {
          isOverlap = true;
        }
      }
      if ( isOverlap ) continue;
    }

    // Check acceptance
    bool isAccepted = false;
    if ( std::abs(jetP4.eta()) <= maxEta_ ) isAccepted = true;
    if ( jetP4.pt() > minPt_ ) isAccepted = true;

    // Check accepted with JES/JER uncertanty
    edm::Ref<Jets> jetRef(jetHandle, i);
    for ( int j=0, m=jetScales.size(); j<m; ++j )
    {
      JetToValueMap::const_iterator jetScale = jetScales.at(j)->find(jetRef);
      if ( jetScale == jetScales.at(j)->end() ) continue;
      const double scale = jetScale->val;
      if ( jetP4.pt()*scale > minPt_ )
      {
        isAccepted = true;
        break;
      }
    }
    if ( !isAccepted ) continue;

    cleanJets->push_back(jetRef);
  }

  const unsigned int nCleanJet = cleanJets->size();
  event.put(cleanJets);

  if ( nCleanJet < minNumber_ or nCleanJet > maxNumber_ ) return false;

  return true;
}

DEFINE_FWK_MODULE(KJetSelector);

