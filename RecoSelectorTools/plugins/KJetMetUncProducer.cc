#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/View.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "KrAFT/RecoSelectorTools/interface/Types.h"

#include <memory>

using namespace edm;
using namespace std;

class KJetMetUncProducer : public edm::EDProducer
{
public:
  KJetMetUncProducer(const edm::ParameterSet& pset);
  ~KJetMetUncProducer() {};

  typedef std::vector<std::string> strings;
  typedef std::vector<pat::Jet> Jets;
  typedef std::vector<pat::MET> METs;

private:
  void beginJob() {};
  void produce(edm::Event& event, const edm::EventSetup& eventSetup);
  void endJob() {};

  edm::InputTag jetLabel_;
  edm::InputTag metLabel_;

  strings jecLevels_;
  std::map<std::string, JetCorrectionUncertainty*> jecCalc_;
  std::map<int, JetCorrectionUncertainty*> jecByFlav_;

  std::map<std::string, std::vector<double> > fJECsUp_, fJECsDn_;
  std::map<std::string, double> metsUpX_, metsUpY_, metsDnX_, metsDnY_;

private:
  void loadJEC(const double jetPt, const double jetEta,
               JetCorrectionUncertainty* jecCalc, double& fJECUp, double& fJECDn)
  {
    if ( !jecCalc ) return;

    jecCalc->setJetPt(jetPt);
    jecCalc->setJetEta(jetEta);
    fJECUp = 1+jecCalc->getUncertainty(true);
    jecCalc->setJetPt(jetPt);
    jecCalc->setJetEta(jetEta);
    fJECDn = 1-jecCalc->getUncertainty(false);
  };

  void loadJER(const double jetEta, double& cJER, double& cJERUp, double& cJERDn) const
  {
    const double absEta = std::abs(jetEta);
    if      ( absEta < 0.5 ) { cJER = 1.079; cJERUp = 1.105; cJERDn = 1.053; }
    else if ( absEta < 1.1 ) { cJER = 1.099; cJERUp = 1.127; cJERDn = 1.071; }
    else if ( absEta < 1.7 ) { cJER = 1.121; cJERUp = 1.150; cJERDn = 1.092; }
    else if ( absEta < 2.3 ) { cJER = 1.208; cJERUp = 1.254; cJERDn = 1.162; }
    else if ( absEta < 2.8 ) { cJER = 1.254; cJERUp = 1.316; cJERDn = 1.192; }
    else if ( absEta < 3.2 ) { cJER = 1.395; cJERUp = 1.458; cJERDn = 1.332; }
    else if ( absEta < 5.0 ) { cJER = 1.056; cJERUp = 1.247; cJERDn = 0.865; }
/* // These values are from 2011
    if      ( absEta < 0.5 ) { cJER = 1.052; cJERUp = 1.115; cJERDn = 0.990; }
    else if ( absEta < 1.1 ) { cJER = 1.057; cJERUp = 1.114; cJERDn = 1.001; }
    else if ( absEta < 1.7 ) { cJER = 1.096; cJERUp = 1.161; cJERDn = 1.032; }
    else if ( absEta < 2.3 ) { cJER = 1.134; cJERUp = 1.228; cJERDn = 1.042; }
    else if ( absEta < 5.0 ) { cJER = 1.288; cJERUp = 1.488; cJERDn = 1.089; }
*/
    else { cJER = cJERUp = cJERDn = 1; }
  };

};

KJetMetUncProducer::KJetMetUncProducer(const edm::ParameterSet& pset)
{
  jetLabel_ = pset.getParameter<edm::InputTag>("jet");
  metLabel_ = pset.getParameter<edm::InputTag>("met");

  // JEC correction
  edm::FileInPath jecFile(pset.getParameter<string>("jecFile"));
  const auto& jecFilePath = jecFile.fullPath();
  for ( auto& jecLevel : pset.getParameter<strings>("jecLevels") )
  {
    if ( jecLevel.find("FlavorPure") != string::npos )
    {
      cerr << "Skipping JEC level " << jecLevel << endl;
      continue;
    }
    jecLevels_.push_back(jecLevel);
    JetCorrectorParameters jecParams(jecFilePath, jecLevel);
    jecCalc_[jecLevel] = new JetCorrectionUncertainty(jecParams);
    fJECsUp_[jecLevel] = std::vector<double>();
    fJECsDn_[jecLevel] = std::vector<double>();
  }

  // Flavour dependent JEC are done manuallly.
  jecByFlav_[22] = new JetCorrectionUncertainty(JetCorrectorParameters(jecFilePath, "FlavorPureGluon"));
  jecByFlav_[ 5] = new JetCorrectionUncertainty(JetCorrectorParameters(jecFilePath, "FlavorPureBottom"));
  jecByFlav_[ 4] = new JetCorrectionUncertainty(JetCorrectorParameters(jecFilePath, "FlavorPureCharm"));
  jecByFlav_[ 0] = new JetCorrectionUncertainty(JetCorrectorParameters(jecFilePath, "FlavourPureQuark"));
  jecByFlav_[1] = jecByFlav_[2] = jecByFlav_[3] = jecByFlav_[0];
  fJECsUp_["Flavor"] = std::vector<double>();
  fJECsDn_["Flavor"] = std::vector<double>();

  for ( const string& jecLevel : jecLevels_ )
  {
    produces<pat::JetToValue>(jecLevel+"Up");
    produces<pat::JetToValue>(jecLevel+"Dn");
    produces<METs>(jecLevel+"Up");
    produces<METs>(jecLevel+"Dn");
  }
  produces<pat::JetToValue>("FlavorUp");
  produces<pat::JetToValue>("FlavorDn");
  produces<METs>("FlavorUp");
  produces<METs>("FlavorDn");

  produces<pat::JetToValue>("res");
  produces<pat::JetToValue>("resUp");
  produces<pat::JetToValue>("resDn");
  produces<METs>("res");
  produces<METs>("resUp");
  produces<METs>("resDn");

}

void KJetMetUncProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<Jets> jetHandle;
  event.getByLabel(jetLabel_, jetHandle);

  edm::Handle<METs> metHandle;
  event.getByLabel(metLabel_, metHandle);

  // Prepare JEC and JER factors, MET
  const pat::MET& met = metHandle->at(0);
  const double metX = met.px(), metY = met.py();
  for ( auto key = fJECsUp_.begin(); key != fJECsUp_.end(); ++key )
  {
    const string& jecLevel = key->first;
    fJECsUp_[jecLevel].clear();
    fJECsDn_[jecLevel].clear();
    metsUpX_[jecLevel] = metsDnX_[jecLevel] = metX;
    metsUpY_[jecLevel] = metsDnY_[jecLevel] = metY;
  }
  std::vector<double> fJERs;
  std::vector<double> fJERsUp;
  std::vector<double> fJERsDn;
  double metResX   = metX, metResY   = metY;
  double metResUpX = metX, metResUpY = metY;
  double metResDnX = metX, metResDnY = metY;

  // Now start to build jet uncertainties
  for ( auto jet : *jetHandle )
  {
    const reco::Candidate::LorentzVector jetP4 = jet.p4();
    const double jetPt = jetP4.pt();
    const double jetEta = jetP4.eta();

    if ( abs(jetEta) > 5 ) continue;

    // Calculate JEC uncertanties
    double fJECUp, fJECDn;
    for ( auto& jecLevel : jecLevels_ )
    {
      loadJEC(jetPt, jetEta, jecCalc_[jecLevel], fJECUp, fJECDn);
      fJECsUp_[jecLevel].push_back(fJECUp);
      fJECsDn_[jecLevel].push_back(fJECDn);

      math::XYZTLorentzVector jetUpP4 = jetP4*fJECUp;
      math::XYZTLorentzVector jetDnP4 = jetP4*fJECDn;

      metsUpX_[jecLevel] += jetP4.px() - jetUpP4.px();
      metsUpY_[jecLevel] += jetP4.py() - jetUpP4.py();
      metsDnX_[jecLevel] += jetP4.px() - jetDnP4.px();
      metsDnY_[jecLevel] += jetP4.py() - jetDnP4.py();
    }

    // Calculate Flavour JEC uncertainties
    const int partonFlav = jet.partonFlavour();
    loadJEC(jetPt, jetEta, jecByFlav_[partonFlav], fJECUp, fJECDn);
    fJECsUp_["Flavor"].push_back(fJECUp);
    fJECsDn_["Flavor"].push_back(fJECDn);

    // JER and uncertainties
    double fJER = 1, fJERUp = 1, fJERDn = 1;
    const reco::GenJet* genJet = jet.genJet();
    if ( genJet and genJet->pt() > 10 )
    {
      double cJER, cJERUp, cJERDn;
      loadJER(jetEta, cJER, cJERUp, cJERDn);
      
      const double genJetPt = genJet->pt();
      const double dPt = jetPt-genJetPt;

      fJER   = max(0., (genJetPt+dPt*cJER  )/jetPt);
      fJERUp = max(0., (genJetPt+dPt*cJERUp)/jetPt);
      fJERDn = max(0., (genJetPt+dPt*cJERDn)/jetPt);

      fJERs.push_back(fJER);
      fJERsUp.push_back(fJERUp);
      fJERsDn.push_back(fJERDn);

      const math::XYZTLorentzVector& rawJetP4 = jet.correctedP4(0);
      const double rawPx = rawJetP4.px();
      const double rawPy = rawJetP4.py();

      const double metDx   = rawPx*(1-fJER  );
      const double metDxUp = rawPx*(1-fJERUp);
      const double metDxDn = rawPx*(1-fJERDn);

      const double metDy   = rawPy*(1-fJER  );
      const double metDyUp = rawPy*(1-fJERUp);
      const double metDyDn = rawPy*(1-fJERDn);

      // Correct MET
      //metDnX    += metDx  ; metDnY    += metDy  ;
      //metUpX    += metDx  ; metUpY    += metDy  ;
      metResX   += metDx  ; metResY   += metDy  ;
      metResUpX += metDxUp; metResUpY += metDyUp;
      metResDnX += metDxDn; metResDnY += metDyDn;

    }
  }

  for ( auto key = fJECsUp_.begin(); key != fJECsUp_.end(); ++key )
  {
    const string& jecLevel = key->first;
    std::auto_ptr<pat::JetToValue> fJECMapUp(new pat::JetToValue);
    std::auto_ptr<pat::JetToValue> fJECMapDn(new pat::JetToValue);
    pat::JetToValue::Filler fillJECUp(*fJECMapUp);
    pat::JetToValue::Filler fillJECDn(*fJECMapDn);
    fillJECUp.insert(jetHandle, fJECsUp_[jecLevel].begin(), fJECsUp_[jecLevel].end());
    fillJECDn.insert(jetHandle, fJECsDn_[jecLevel].begin(), fJECsDn_[jecLevel].end());
    fillJECUp.fill();
    fillJECDn.fill();
    event.put(fJECMapUp, jecLevel+"Up");
    event.put(fJECMapDn, jecLevel+"Dn");

    pat::MET metUp, metDn;
    const double metUpX = metsUpX_[jecLevel];
    const double metUpY = metsUpY_[jecLevel];
    const double metDnX = metsDnX_[jecLevel];
    const double metDnY = metsDnY_[jecLevel];
    metUp.setP4(reco::Candidate::LorentzVector(metUpX, metUpY, 0, hypot(metUpX, metUpY)));
    metDn.setP4(reco::Candidate::LorentzVector(metDnX, metDnY, 0, hypot(metDnX, metDnY)));
    std::auto_ptr<METs> metsUp(new METs);
    std::auto_ptr<METs> metsDn(new METs);
    metsUp->push_back(metUp);
    metsDn->push_back(metDn);
    event.put(metsUp, jecLevel+"Up");
    event.put(metsDn, jecLevel+"Dn");
  }

  std::auto_ptr<pat::JetToValue> fJERMap(new pat::JetToValue);
  std::auto_ptr<pat::JetToValue> fJERMapUp(new pat::JetToValue);
  std::auto_ptr<pat::JetToValue> fJERMapDn(new pat::JetToValue);
  pat::JetToValue::Filler fillJER(*fJERMap);
  pat::JetToValue::Filler fillJERUp(*fJERMapUp);
  pat::JetToValue::Filler fillJERDn(*fJERMapDn);
  fillJER.insert(jetHandle, fJERs.begin(), fJERs.end());
  fillJERUp.insert(jetHandle, fJERsUp.begin(), fJERsUp.end());
  fillJERDn.insert(jetHandle, fJERsDn.begin(), fJERsDn.end());
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

DEFINE_FWK_MODULE(KJetMetUncProducer);

