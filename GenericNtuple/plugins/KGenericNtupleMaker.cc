#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "KrAFT/GenericNtuple/interface/GenericEvent.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class KGenericNtupleMaker : public edm::EDAnalyzer
{
public:
  KGenericNtupleMaker(const edm::ParameterSet& pset);
  ~KGenericNtupleMaker();

  void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool isInAcceptance(const pat::Electron& electron);
  bool isInAcceptance(const pat::Muon& muon);
  bool isInAcceptance(const pat::Jet& jet);

private:
  // Input objects
  edm::InputTag genEventInfoLabel_;
  edm::InputTag genParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag puWeightLabel_;
  edm::InputTag vertexLabel_;

  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag metLabel_;
  edm::InputTag jpsiLabel_;
  std::string bTagType_;

  std::vector<std::string> eventCounterLabels_;

  bool isMC_;

  unsigned int muonMinNumber_, electronMinNumber_;
  unsigned int jetMinNumber_;

  double muonDz_, electronDz_;
  double jetLeptonDeltaR_;

  TH1F* hEventCounter_;
  TNamed* dataType_;

  // Output tree
  TTree* tree_;
  GenericEvent* fevent_;

};

KGenericNtupleMaker::KGenericNtupleMaker(const edm::ParameterSet& pset)
{
  isMC_ = pset.getParameter<bool>("isMC");

  // Input labels
  puWeightLabel_ = pset.getParameter<edm::InputTag>("puWeight");
  vertexLabel_ = pset.getParameter<edm::InputTag>("vertex");

  edm::ParameterSet electronPSet = pset.getParameter<edm::ParameterSet>("electron");
  electronMinNumber_ = electronPSet.getParameter<unsigned int>("minNumber");
  electronLabel_ = electronPSet.getParameter<edm::InputTag>("src");

  edm::ParameterSet muonPSet = pset.getParameter<edm::ParameterSet>("muon");
  muonMinNumber_ = muonPSet.getParameter<unsigned int>("minNumber");
  muonLabel_ = muonPSet.getParameter<edm::InputTag>("src");

  edm::ParameterSet jetPSet = pset.getParameter<edm::ParameterSet>("jet");
  jetMinNumber_ = jetPSet.getParameter<unsigned int>("minNumber");
  jetLeptonDeltaR_ = jetPSet.getParameter<double>("leptonDeltaR");
  jetLabel_ = jetPSet.getParameter<edm::InputTag>("src");
  bTagType_ = jetPSet.getParameter<std::string>("bTagType");

  edm::ParameterSet metPSet = pset.getParameter<edm::ParameterSet>("met");
  metLabel_ = metPSet.getParameter<edm::InputTag>("src");

  edm::ParameterSet jpsiPSet = pset.getParameter<edm::ParameterSet>("jpsi");
  jpsiLabel_ = jpsiPSet.getParameter<edm::InputTag>("src");

  // Event counter
  eventCounterLabels_ = pset.getParameter<std::vector<std::string> >("eventCounters");

  if ( isMC_ )
  {
    genEventInfoLabel_ = pset.getParameter<edm::InputTag>("genEventInfo");
    genParticleLabel_ = pset.getParameter<edm::InputTag>("genParticle");
    genJetLabel_ = pset.getParameter<edm::InputTag>("genJet");
  }

  // Output histograms and tree
  edm::Service<TFileService> fs;
  dataType_ = fs->make<TNamed>("dataType", isMC_ ? "MC" : "Data");
  hEventCounter_ = fs->make<TH1F>("hEventCounter", "Event counter", eventCounterLabels_.size(), 1, eventCounterLabels_.size()+1);
  for ( int i=0, n=eventCounterLabels_.size(); i<n; ++i )
  {
    hEventCounter_->GetXaxis()->SetBinLabel(i+1, eventCounterLabels_.at(i).c_str());
  }

  tree_ = fs->make<TTree>("event", "Mixed event tree");
  fevent_ = new GenericEvent(isMC_);
  fevent_->book(tree_);
}

KGenericNtupleMaker::~KGenericNtupleMaker()
{
}

void KGenericNtupleMaker::endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
{
  for ( int i=0, n=eventCounterLabels_.size(); i<n; ++i )
  {
    edm::Handle<edm::MergeableCounter> eventCounterHandle;
    lumi.getByLabel(edm::InputTag(eventCounterLabels_.at(i)), eventCounterHandle);
    if ( !eventCounterHandle.isValid() ) continue;

    hEventCounter_->Fill(i+1, eventCounterHandle->value);
  }
}

void KGenericNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  fevent_->clear();

  if ( isMC_ and event.isRealData() ) isMC_ = false;

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);
  fevent_->nVertex_ = vertexHandle->size();
  const reco::Vertex& pv = vertexHandle->at(0);

  if ( !isMC_ )
  {
    fevent_->puWeight_   = 1.0;
    fevent_->puWeightUp_ = 1.0;
    fevent_->puWeightDn_ = 1.0;
  }
  else
  {
    edm::Handle<double> puWeightHandle, puWeightUpHandle, puWeightDnHandle;
    event.getByLabel(puWeightLabel_, puWeightHandle);
    event.getByLabel(edm::InputTag(puWeightLabel_.label(), "up"), puWeightUpHandle);
    event.getByLabel(edm::InputTag(puWeightLabel_.label(), "dn"), puWeightDnHandle);
    fevent_->puWeight_   = *(puWeightHandle.product()  );
    fevent_->puWeightUp_ = *(puWeightUpHandle.product());
    fevent_->puWeightDn_ = *(puWeightDnHandle.product());
  }

  edm::Handle<std::vector<pat::Electron> > electronHandle;
  event.getByLabel(electronLabel_, electronHandle);
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    const pat::Electron& e = electronHandle->at(i);
    if ( e.pt() < 10 or std::abs(e.eta()) > 2.5 ) continue;
    if ( !e.isPF() and e.gsfTrack().isNull() ) continue;

    //if ( abs(e.dz(pv.position())) > electronDz_ ) continue;
    const double scEta = e.superCluster()->eta();
    const double dxy = e.dB();
    const double mva = e.electronID("mvaTrigV0");

    int eType = 0;
    // Veto electrons
    if ( 0.0 < mva and mva < 1.0 ) eType += 1;
    if ( (e.gsfTrack().isNonnull() or e.isPF()) and e.passConversionVeto() and
         e.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 and
         mva > 0.5 )
    {
      // Top Single electron ID
      if ( dxy < 0.02 and not (1.4442 < std::abs(scEta) and std::abs(scEta) < 1.5660) ) eType += 10;
      // Top Dilepton electron ID
      if ( dxy < 0.04 ) eType += 100;
    }

    fevent_->electrons_pt_    ->push_back(e.pt()    );
    fevent_->electrons_eta_   ->push_back(e.eta()   );
    fevent_->electrons_phi_   ->push_back(e.phi()   );
    fevent_->electrons_m_     ->push_back(e.mass()  );
    fevent_->electrons_Q_     ->push_back(e.charge());
    fevent_->electrons_type_  ->push_back(eType     );
    fevent_->electrons_relIso_->push_back(e.userIso(2)); // rho corrected isolation

    fevent_->electrons_mva_->push_back(e.mva());
    fevent_->electrons_scEta_->push_back(scEta);
  }
  if ( fevent_->electrons_pt_->size() < electronMinNumber_ ) return;

  edm::Handle<std::vector<pat::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const pat::Muon& mu = muonHandle->at(i);
    if ( mu.pt() < 10 or std::abs(mu.eta()) > 2.5 ) continue;
    //if ( abs(mu.dz(pv.position())) > muonDz_ ) continue;

    int muType = 0;
    if ( mu.isPFMuon() and (mu.isGlobalMuon() or mu.isTrackerMuon()) ) muType += 1;
    if ( muon::isLooseMuon(mu) ) muType += 10;
    if ( muon::isSoftMuon(mu, pv) ) muType += 100;
    if ( muon::isTightMuon(mu, pv) ) muType += 1000;
    if ( muon::isHighPtMuon(mu, pv, reco::improvedTuneP) ) muType += 10000;

    fevent_->muons_pt_->push_back(mu.pt()   );
    fevent_->muons_eta_->push_back(mu.eta() );
    fevent_->muons_phi_->push_back(mu.phi() );
    fevent_->muons_m_->push_back(mu.mass()  );
    fevent_->muons_Q_->push_back(mu.charge());
    fevent_->muons_type_->push_back(muType  );
    fevent_->muons_relIso_->push_back(mu.userIso(1)); // dBeta corrected isolation
  }
  if ( fevent_->muons_pt_->size() < muonMinNumber_ ) return;

  edm::Handle<std::vector<pat::MET> > metHandle, metUpHandle, metDnHandle, metResUpHandle, metResDnHandle;
  event.getByLabel(metLabel_, metHandle);
  event.getByLabel(edm::InputTag(metLabel_.label(), "up"), metUpHandle);
  event.getByLabel(edm::InputTag(metLabel_.label(), "dn"), metDnHandle);
  fevent_->met_pt_ = metHandle->at(0).pt();
  fevent_->metUp_pt_ = metUpHandle->at(0).pt();
  fevent_->metDn_pt_ = metDnHandle->at(0).pt();
  fevent_->met_phi_ = metHandle->at(0).phi();
  fevent_->metUp_phi_ = metUpHandle->at(0).phi();
  fevent_->metDn_phi_ = metDnHandle->at(0).phi();
  if ( isMC_ )
  {
    event.getByLabel(edm::InputTag(metLabel_.label(), "resUp"), metResUpHandle);
    event.getByLabel(edm::InputTag(metLabel_.label(), "resDn"), metResDnHandle);
    fevent_->metResUp_pt_ = metResUpHandle->at(0).pt();
    fevent_->metResDn_pt_ = metResDnHandle->at(0).pt();
    fevent_->metResUp_phi_ = metResUpHandle->at(0).phi();
    fevent_->metResDn_phi_ = metResDnHandle->at(0).phi();
  }

  if ( isMC_ )
  {
    // Event weight and PDF stuffs
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    event.getByLabel(genEventInfoLabel_, genEventInfoHandle);

    const gen::PdfInfo* pdf = genEventInfoHandle->pdf();

    fevent_->pdf_id1_ = pdf->id.first ;
    fevent_->pdf_id2_ = pdf->id.second;
    fevent_->pdf_x1_  = pdf->x.first  ;
    fevent_->pdf_x2_  = pdf->x.second ;
    fevent_->pdf_q_   = pdf->scalePDF ;
    //const double pdf_xPDF1 = pdf->xPDF.first, pdf_xPDF2 = pdf->xPDF.second;

    fevent_->genWeight_ = genEventInfoHandle->weights().at(0);

    // Generator level objects
    edm::Handle<reco::GenParticleCollection> genHandle;
    event.getByLabel(genParticleLabel_, genHandle);
    edm::Handle<reco::GenJetCollection> genJetHandle;
    event.getByLabel(genJetLabel_, genJetHandle);

    // Find top quark from the genParticles
    const reco::GenParticle* firstGen = &genHandle->at(0);
    for ( int i=0, n=genHandle->size(); i<n; ++i )
    {
      const reco::GenParticle& p = genHandle->at(i);
      const int charge = p.charge();
      const int mother1 = dynamic_cast<const reco::GenParticle*>(p.mother(0))-firstGen;
      const int mother2 = mother1+p.numberOfMothers();
      const int daughter1 = dynamic_cast<const reco::GenParticle*>(p.daughter(0))-firstGen;
      const int daughter2 = daughter2+p.numberOfDaughters();

      fevent_->genParticles_pt_ ->push_back(p.pt()  );
      fevent_->genParticles_eta_->push_back(p.eta() );
      fevent_->genParticles_phi_->push_back(p.phi() );
      fevent_->genParticles_m_  ->push_back(p.mass());
      fevent_->genParticles_pdgId_->push_back(p.pdgId());
      fevent_->genParticles_mother1_->push_back(mother1);
      fevent_->genParticles_mother1_->push_back(mother2);
      fevent_->genParticles_mother1_->push_back(daughter1);
      fevent_->genParticles_mother1_->push_back(daughter2);
    }

    for ( int i=0, n=genJetHandle->size(); i<n; ++i )
    {
      const reco::GenJet& p = genJetHandle->at(i);
      if ( p.pt() < 20 or std::abs(p.eta()) > 2.5 ) continue;
      fevent_->genJets_pt_ ->push_back(p.pt()  );
      fevent_->genJets_eta_->push_back(p.eta() );
      fevent_->genJets_phi_->push_back(p.phi() );
      fevent_->genJets_m_  ->push_back(p.mass());
    }
  }

  unsigned int nJet = 0;
  edm::Handle<std::vector<pat::Jet> > jetHandle;
  edm::Handle<std::vector<pat::Jet> > jetUpHandle;
  edm::Handle<std::vector<pat::Jet> > jetDnHandle;
  edm::Handle<std::vector<pat::Jet> > jetResUpHandle;
  edm::Handle<std::vector<pat::Jet> > jetResDnHandle;
  event.getByLabel(jetLabel_, jetHandle);
  event.getByLabel(edm::InputTag(jetLabel_.label(), "up"), jetUpHandle);
  event.getByLabel(edm::InputTag(jetLabel_.label(), "dn"), jetDnHandle);
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetHandle->at(i);
    if ( !isInAcceptance(jet) ) continue;
    fevent_->jets_pt_ ->push_back(jet.pt());
    fevent_->jets_eta_->push_back(jet.eta());
    fevent_->jets_phi_->push_back(jet.phi());
    fevent_->jets_m_  ->push_back(jet.mass());
    fevent_->jets_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  nJet = std::max(nJet, (unsigned int)fevent_->jets_pt_->size());
  for ( int i=0, n=jetUpHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetUpHandle->at(i);
    if ( !isInAcceptance(jet) ) continue;
    fevent_->jetsUp_pt_ ->push_back(jet.pt());
    fevent_->jetsUp_eta_->push_back(jet.eta());
    fevent_->jetsUp_phi_->push_back(jet.phi());
    fevent_->jetsUp_m_  ->push_back(jet.mass());
    fevent_->jetsUp_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  nJet = std::max(nJet, (unsigned int)fevent_->jetsUp_pt_->size());
  for ( int i=0, n=jetDnHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetDnHandle->at(i);
    if ( !isInAcceptance(jet) ) continue;
    fevent_->jetsDn_pt_ ->push_back(jet.pt());
    fevent_->jetsDn_eta_->push_back(jet.eta());
    fevent_->jetsDn_phi_->push_back(jet.phi());
    fevent_->jetsDn_m_  ->push_back(jet.mass());
    fevent_->jetsDn_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  nJet = std::max(nJet, (unsigned int)fevent_->jetsDn_pt_->size());
  if ( isMC_ )
  {
    event.getByLabel(edm::InputTag(jetLabel_.label(), "resUp"), jetResUpHandle);
    event.getByLabel(edm::InputTag(jetLabel_.label(), "resDn"), jetResDnHandle);
    for ( int i=0, n=jetResUpHandle->size(); i<n; ++i )
    {
      const pat::Jet& jet = jetResUpHandle->at(i);
      if ( !isInAcceptance(jet) ) continue;
      fevent_->jetsResUp_pt_  ->push_back(jet.pt()  );
      fevent_->jetsResUp_eta_ ->push_back(jet.eta() );
      fevent_->jetsResUp_phi_ ->push_back(jet.phi() );
      fevent_->jetsResUp_m_   ->push_back(jet.mass());
      fevent_->jetsResUp_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
    }
    nJet = std::max(nJet, (unsigned int)fevent_->jetsResUp_pt_->size());
    for ( int i=0, n=jetResDnHandle->size(); i<n; ++i )
    {
      const pat::Jet& jet = jetResDnHandle->at(i);
      if ( !isInAcceptance(jet) ) continue;
      fevent_->jetsResDn_pt_  ->push_back(jet.pt()  );
      fevent_->jetsResDn_eta_ ->push_back(jet.eta() );
      fevent_->jetsResDn_phi_ ->push_back(jet.phi() );
      fevent_->jetsResDn_m_   ->push_back(jet.mass());
      fevent_->jetsResDn_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
    }
    nJet = std::max(nJet, (unsigned int)fevent_->jetsResDn_pt_->size());
  }
  if ( nJet < jetMinNumber_ ) return;

  edm::Handle<std::vector<reco::VertexCompositeCandidate> > jpsiHandle;
  event.getByLabel(jpsiLabel_, jpsiHandle);
  edm::Handle<std::vector<double> > jpsiLxyHandle;
  event.getByLabel(edm::InputTag(jpsiLabel_.label(), "lxy"), jpsiLxyHandle);
  for ( int i=0, n=jpsiHandle->size(); i<n; ++i )
  {
    const reco::VertexCompositeCandidate& jpsiCand = jpsiHandle->at(i);
    fevent_->jpsis_pt_ ->push_back(jpsiCand.pt()  );
    fevent_->jpsis_eta_->push_back(jpsiCand.eta() );
    fevent_->jpsis_phi_->push_back(jpsiCand.phi() );
    fevent_->jpsis_m_  ->push_back(jpsiCand.mass());
    fevent_->jpsis_lxy_->push_back(jpsiLxyHandle->at(i));
  }

  // Now put jets in current event to the event cache
  fevent_->run_ = event.run();
  fevent_->lumi_ = event.luminosityBlock();
  fevent_->event_ = event.id().event();

  fevent_->tree_->Fill();
}

bool KGenericNtupleMaker::isInAcceptance(const pat::Electron& electron)
{
  if ( electron.pt() < 10 ) return false;
  if ( std::abs(electron.eta()) > 2.5 ) return false;
  return true;
}

bool KGenericNtupleMaker::isInAcceptance(const pat::Muon& muon)
{
  if ( muon.pt() < 10 ) return false;
  if ( std::abs(muon.eta()) > 2.5 ) return false;
  return true;
}

bool KGenericNtupleMaker::isInAcceptance(const pat::Jet& jet)
{
  if ( jet.pt() < 30 ) return false;
  if ( std::abs(jet.eta()) > 2.5 ) return false;
  return true;
}

DEFINE_FWK_MODULE(KGenericNtupleMaker);
