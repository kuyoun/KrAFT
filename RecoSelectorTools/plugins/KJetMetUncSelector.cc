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
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
//#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "KrAFT/RecoSelectorTools/interface/Types.h"

#include <memory>

using namespace edm;
using namespace std;

class KJetMetUncSelector : public edm::EDFilter
{
public:
  KJetMetUncSelector(const edm::ParameterSet& pset);
  ~KJetMetUncSelector() {};

  typedef std::vector<pat::Jet> Jets;
  typedef std::vector<pat::MET> METs;

private:
  void beginJob() {};
  bool filter(edm::Event& event, const edm::EventSetup& eventSetup);
  void endJob() {};

  bool doClean_;
  double overlapDeltaR_;
  unsigned int minNumber_;
  unsigned int maxNumber_;

  edm::InputTag jetLabel_;
  edm::InputTag metLabel_;
  std::vector<edm::InputTag> overlapCandLabels_;

  PFJetIDSelectionFunctor* isGoodJet_;
  JetCorrectionUncertainty* jecUncCalculator_;
  double minPt_, maxEta_;

  bool isMC_;

private:
  void loadJEC(const double jetPt, const double jetEta, double& fJECUp, double& fJECDn)
  {
    jecUncCalculator_->setJetPt(jetPt);
    jecUncCalculator_->setJetEta(jetEta);
    fJECUp = 1+jecUncCalculator_->getUncertainty(true);
    jecUncCalculator_->setJetPt(jetPt);
    jecUncCalculator_->setJetEta(jetEta);
    fJECDn = 1-jecUncCalculator_->getUncertainty(false);
  };

  void loadJER(const double jetEta, double& cJER, double& cJERUp, double& cJERDn) const
  {
    if      ( jetEta < 0.5 ) { cJER = 1.052; cJERUp = 1.115; cJERDn = 0.990; }
    else if ( jetEta < 1.1 ) { cJER = 1.057; cJERUp = 1.114; cJERDn = 1.001; }
    else if ( jetEta < 1.7 ) { cJER = 1.096; cJERUp = 1.161; cJERDn = 1.032; }
    else if ( jetEta < 2.3 ) { cJER = 1.134; cJERUp = 1.228; cJERDn = 1.042; }
    else if ( jetEta < 5.0 ) { cJER = 1.288; cJERUp = 1.488; cJERDn = 1.089; }
    else { cJER = cJERUp = cJERDn = 1; }
  };

};

KJetMetUncSelector::KJetMetUncSelector(const edm::ParameterSet& pset)
{
  isMC_ = pset.getParameter<bool>("isMC");

  jetLabel_ = pset.getParameter<edm::InputTag>("jet");
  metLabel_ = pset.getParameter<edm::InputTag>("met");

  // JEC correction
  edm::FileInPath jecFilePath(pset.getParameter<string>("jecFile"));
  const std::string jecLevel = pset.getParameter<string>("jecLevel");
  jecUncCalculator_ = new JetCorrectionUncertainty(JetCorrectorParameters(jecFilePath.fullPath(), jecLevel));

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
  produces<pat::JetToValue>("up");
  produces<pat::JetToValue>("dn");
  produces<METs>("up");
  produces<METs>("dn");
  if ( isMC_ )
  {
    produces<pat::JetToValue>("res");
    produces<pat::JetToValue>("resUp");
    produces<pat::JetToValue>("resDn");
    produces<METs>("res");
    produces<METs>("resUp");
    produces<METs>("resDn");
  }

}

bool KJetMetUncSelector::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<Jets> jetHandle;
  event.getByLabel(jetLabel_, jetHandle);

  edm::Handle<METs> metHandle;
  event.getByLabel(metLabel_, metHandle);

  // Prepare JEC and JER factors
  std::vector<double> fJECsUp;
  std::vector<double> fJECsDn;
  std::vector<double> fJERs;
  std::vector<double> fJERsUp;
  std::vector<double> fJERsDn;

  // Prepare MET
  const pat::MET& met = metHandle->at(0);
  const double metX = met.px(), metY = met.py();
  double metUpX = metX, metUpY = metY;
  double metDnX = metX, metDnY = metY;
  double metResX   = metX, metResY   = metY;
  double metResUpX = metX, metResUpY = metY;
  double metResDnX = metX, metResDnY = metY;

  // Collect candidates for overlap removal
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

  // Now start to build jet collection
  std::auto_ptr<std::vector<pat::Jet> > cleanJets(new std::vector<pat::Jet>());

  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    pat::Jet jet = jetHandle->at(i);
    if ( !(*isGoodJet_)(jet) ) continue;
    const reco::Candidate::LorentzVector jetP4 = jet.p4();
    const double jetPt = jetP4.pt();
    const double jetEta = jetP4.eta();

    // Calculate JEC uncertanties
    double fJECUp, fJECDn;
    loadJEC(jetPt, jetEta, fJECUp, fJECDn);
    math::XYZTLorentzVector jetUpP4 = jetP4*fJECUp;
    math::XYZTLorentzVector jetDnP4 = jetP4*fJECDn;

    metUpX += jetP4.px() - jetUpP4.px();
    metUpY += jetP4.py() - jetUpP4.py();
    metDnX += jetP4.px() - jetDnP4.px();
    metDnY += jetP4.py() - jetDnP4.py();

    // JER and uncertainties
    double fJER = 1, fJERUp = 1, fJERDn = 1;
    if ( isMC_ )
    {
      const reco::GenJet* genJet = jet.genJet();
      if ( genJet and genJet->pt() > 10 )
      {
        double cJER, cJERUp, cJERDn;
        loadJER(jetEta, cJER, cJERUp, cJERDn);

        const math::XYZTLorentzVector& rawJetP4 = jet.correctedP4(0);
        const double rawPx = rawJetP4.px();
        const double rawPy = rawJetP4.py();

        const double genJetPt = genJet->pt();
        const double dPt = jetPt-genJetPt;

        fJER   = max(0., (genJetPt+dPt*cJER  )/jetPt);
        fJERUp = max(0., (genJetPt+dPt*cJERUp)/jetPt);
        fJERDn = max(0., (genJetPt+dPt*cJERDn)/jetPt);

        const double metDx   = rawPx*(1-fJER  );
        const double metDxUp = rawPx*(1-fJERUp);
        const double metDxDn = rawPx*(1-fJERDn);

        const double metDy   = rawPy*(1-fJER  );
        const double metDyUp = rawPy*(1-fJERUp);
        const double metDyDn = rawPy*(1-fJERDn);

        // Correct MET
        metDnX    += metDx  ; metDnY    += metDy  ;
        metUpX    += metDx  ; metUpY    += metDy  ;
        metResX   += metDx  ; metResY   += metDy  ;
        metResUpX += metDxUp; metResUpY += metDyUp;
        metResDnX += metDxDn; metResDnY += metDyDn;
      }
    }

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
    if ( std::abs(jetEta) > maxEta_ ) continue;
    if ( jetPt*fJECDn*fJERDn < minPt_ ) continue;

    cleanJets->push_back(jet);
    fJECsUp.push_back(fJECUp);
    fJECsDn.push_back(fJECDn);
    if ( isMC_ )
    {
      fJERs.push_back(fJER);
      fJERsUp.push_back(fJERUp);
      fJERsDn.push_back(fJERDn);
    }
  }

  const unsigned int nCleanJet = cleanJets->size();
  edm::OrphanHandle<pat::JetCollection> outColl = event.put(cleanJets);

  std::auto_ptr<pat::JetToValue> fJECMapUp(new pat::JetToValue);
  std::auto_ptr<pat::JetToValue> fJECMapDn(new pat::JetToValue);
  pat::JetToValue::Filler fillJECUp(*fJECMapUp);
  pat::JetToValue::Filler fillJECDn(*fJECMapDn);
  fillJECUp.insert(outColl, fJECsUp.begin(), fJECsUp.end());
  fillJECDn.insert(outColl, fJECsDn.begin(), fJECsDn.end());
  fillJECUp.fill();
  fillJECDn.fill();
  event.put(fJECMapUp, "up");
  event.put(fJECMapDn, "dn");

  pat::MET metUp, metDn;
  metUp.setP4(reco::Candidate::LorentzVector(metUpX, metUpY, 0, hypot(metUpX, metUpY)));
  metDn.setP4(reco::Candidate::LorentzVector(metDnX, metDnY, 0, hypot(metDnX, metDnY)));
  std::auto_ptr<METs> metsUp(new METs);
  std::auto_ptr<METs> metsDn(new METs);
  metsUp->push_back(metUp);
  metsDn->push_back(metDn);
  event.put(metsUp, "up");
  event.put(metsDn, "dn");

  if ( isMC_ )
  {
    std::auto_ptr<pat::JetToValue> fJERMap(new pat::JetToValue);
    std::auto_ptr<pat::JetToValue> fJERMapUp(new pat::JetToValue);
    std::auto_ptr<pat::JetToValue> fJERMapDn(new pat::JetToValue);
    pat::JetToValue::Filler fillJER(*fJERMap);
    pat::JetToValue::Filler fillJERUp(*fJERMapUp);
    pat::JetToValue::Filler fillJERDn(*fJERMapDn);
    fillJER.insert(outColl, fJERs.begin(), fJERs.end());
    fillJERUp.insert(outColl, fJERsUp.begin(), fJERsUp.end());
    fillJERDn.insert(outColl, fJERsDn.begin(), fJERsDn.end());
    fillJER.fill();
    fillJERUp.fill();
    fillJERDn.fill();
    event.put(fJERMap, "res");
    event.put(fJERMapUp, "resUp");
    event.put(fJERMapDn, "resDn");

    pat::MET metRes, metResUp, metResDn;
    metRes.setP4(reco::Candidate::LorentzVector(metResX, metResY, 0, hypot(metResX, metResY)));
    metResUp.setP4(reco::Candidate::LorentzVector(metResUpX, metResUpY, 0, hypot(metResUpX, metResUpY)));
    metResDn.setP4(reco::Candidate::LorentzVector(metResDnX, metResDnY, 0, hypot(metResDnX, metResDnY)));

    std::auto_ptr<METs> metsRes(new METs);
    std::auto_ptr<METs> metsResUp(new METs);
    std::auto_ptr<METs> metsResDn(new METs);
    metsRes->push_back(metRes);
    metsResUp->push_back(metResUp);
    metsResDn->push_back(metResDn);
    event.put(metsRes, "res");
    event.put(metsResUp, "resUp");
    event.put(metsResDn, "resDn");
  }

  if ( nCleanJet < minNumber_ or nCleanJet > maxNumber_ ) return false;

  return true;
}

DEFINE_FWK_MODULE(KJetMetUncSelector);

