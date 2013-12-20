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

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class KFlatTreeMaker : public edm::EDAnalyzer
{
public:
  KFlatTreeMaker(const edm::ParameterSet& pset);
  ~KFlatTreeMaker();

  void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

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

  unsigned int muonMinNumber_, electronMinNumber_;

  double muonDz_, electronDz_;
  double jetLeptonDeltaR_;

  TH1F* hEventCounter_;
  TNamed* dataType_;

  // Output tree
  TTree* tree_;
  int run_, lumi_, event_;
  double puWeight_, puWeightUp_, puWeightDn_;
  int nVertex_;

  typedef std::vector<int> ints;
  typedef std::vector<unsigned int> uints;
  typedef std::vector<double> doubles;
  typedef ints* intsP;
  typedef uints* uintsP;
  typedef doubles* doublesP;

  doublesP muons_pt_, muons_eta_, muons_phi_, muons_m_;

  intsP    muons_Q_;
  uintsP   muons_type_;
  doublesP muons_relIso_;

  doublesP electrons_pt_, electrons_eta_, electrons_phi_, electrons_m_;
  intsP    electrons_Q_;
  uintsP   electrons_type_;
  doublesP electrons_relIso_;

  doublesP electrons_mva_;
  doublesP electrons_scEta_;

  doublesP jets_pt_, jets_eta_, jets_phi_, jets_m_;
  doublesP jetsUp_pt_, jetsUp_eta_, jetsUp_phi_, jetsUp_m_;
  doublesP jetsDn_pt_, jetsDn_eta_, jetsDn_phi_, jetsDn_m_;
  doublesP jets_bTag_, jetsUp_bTag_, jetsDn_bTag_;
  double met_pt_, met_phi_;
  double metUp_pt_, metUp_phi_;
  double metDn_pt_, metDn_phi_;

  doublesP jpsis_pt_, jpsis_eta_, jpsis_phi_, jpsis_m_;
  doublesP jpsis_lxy_;

  // JER
  doublesP jetsResUp_pt_, jetsResUp_eta_, jetsResUp_phi_, jetsResUp_m_;
  doublesP jetsResDn_pt_, jetsResDn_eta_, jetsResDn_phi_, jetsResDn_m_;
  doublesP jetsResUp_bTag_, jetsResDn_bTag_;
  double metResUp_pt_, metResUp_phi_;
  double metResDn_pt_, metResDn_phi_;

  // Generator level information
  bool isMC_;

  double genWeight_;
  int pdf_id1_, pdf_id2_;
  double pdf_q_, pdf_x1_, pdf_x2_;

  doublesP genMuons_pt_, genMuons_eta_, genMuons_phi_, genMuons_m_;
  doublesP genElectrons_pt_, genElectrons_eta_, genElectrons_phi_, genElectrons_m_;
  doublesP genNeutrinos_pt_, genNeutrinos_eta_, genNeutrinos_phi_;
  intsP genMuons_Q_, genElectrons_Q_;
  doublesP genJets_pt_, genJets_eta_, genJets_phi_, genJets_m_;

};

KFlatTreeMaker::KFlatTreeMaker(const edm::ParameterSet& pset)
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

  muons_pt_ = new doubles; muons_eta_ = new doubles; muons_phi_ = new doubles; muons_m_ = new doubles;
  muons_Q_ = new ints;
  muons_type_ = new uints();
  muons_relIso_ = new doubles;

  electrons_pt_ = new doubles; electrons_eta_ = new doubles; electrons_phi_ = new doubles; electrons_m_ = new doubles;
  electrons_Q_ = new ints;
  electrons_type_ = new uints();
  electrons_relIso_ = new doubles;

  electrons_mva_ = new doubles;
  electrons_scEta_ = new doubles;

  jets_pt_ = new doubles; jets_eta_ = new doubles; jets_phi_ = new doubles; jets_m_ = new doubles;
  jetsUp_pt_ = new doubles; jetsUp_eta_ = new doubles; jetsUp_phi_ = new doubles; jetsUp_m_ = new doubles;
  jetsDn_pt_ = new doubles; jetsDn_eta_ = new doubles; jetsDn_phi_ = new doubles; jetsDn_m_ = new doubles;
  jets_bTag_ = new doubles;
  jetsUp_bTag_ = new doubles;
  jetsDn_bTag_ = new doubles;

  jpsis_pt_ = new doubles; jpsis_eta_ = new doubles; jpsis_phi_ = new doubles; jpsis_m_ = new doubles;
  jpsis_lxy_ = new doubles;

  tree_ = fs->make<TTree>("event", "Mixed event tree");
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/I");

  tree_->Branch("puWeight", &puWeight_, "puWeight/D");
  tree_->Branch("puWeightUp", &puWeightUp_, "puWeightUp/D");
  tree_->Branch("puWeightDn", &puWeightDn_, "puWeightDn/D");
  tree_->Branch("nVertex", &nVertex_, "nVertex/I");

  tree_->Branch("muons_pt"  , muons_pt_  );
  tree_->Branch("muons_eta" , muons_eta_ );
  tree_->Branch("muons_phi" , muons_phi_ );
  tree_->Branch("muons_m"   , muons_m_   );
  tree_->Branch("muons_Q"   , muons_Q_   );
  tree_->Branch("muons_type", muons_type_);
  tree_->Branch("muons_relIso", muons_relIso_);

  tree_->Branch("electrons_pt"  , electrons_pt_  );
  tree_->Branch("electrons_eta" , electrons_eta_ );
  tree_->Branch("electrons_phi" , electrons_phi_ );
  tree_->Branch("electrons_m"   , electrons_m_   );
  tree_->Branch("electrons_Q"   , electrons_Q_   );
  tree_->Branch("electrons_type", electrons_type_);
  tree_->Branch("electrons_relIso", electrons_relIso_);

  tree_->Branch("electrons_mva", electrons_mva_);
  tree_->Branch("electrons_scEta", electrons_scEta_);

  tree_->Branch("jets_pt" , jets_pt_ );
  tree_->Branch("jets_eta", jets_eta_);
  tree_->Branch("jets_phi", jets_phi_);
  tree_->Branch("jets_m"  , jets_m_  );

  tree_->Branch("jetsUp_pt" , jetsUp_pt_ );
  tree_->Branch("jetsUp_eta", jetsUp_eta_);
  tree_->Branch("jetsUp_phi", jetsUp_phi_);
  tree_->Branch("jetsUp_m"  , jetsUp_m_  );

  tree_->Branch("jetsDn_pt" , jetsDn_pt_ );
  tree_->Branch("jetsDn_eta", jetsDn_eta_);
  tree_->Branch("jetsDn_phi", jetsDn_phi_);
  tree_->Branch("jetsDn_m"  , jetsDn_m_  );

  tree_->Branch("jets_bTag"  , jets_bTag_  );
  tree_->Branch("jetsUp_bTag", jetsUp_bTag_);
  tree_->Branch("jetsDn_bTag", jetsDn_bTag_);

  tree_->Branch("met_pt"  , &met_pt_  , "met_pt/D"  );
  tree_->Branch("metUp_pt", &metUp_pt_, "metUp_pt/D");
  tree_->Branch("metDn_pt", &metDn_pt_, "metDn_pt/D");

  tree_->Branch("met_phi"  , &met_phi_  , "met_phi/D"  );
  tree_->Branch("metUp_phi", &metUp_phi_, "metUp_phi/D");
  tree_->Branch("metDn_phi", &metDn_phi_, "metDn_phi/D");

  tree_->Branch("jpsis_pt" , &jpsis_pt_ );
  tree_->Branch("jpsis_eta", &jpsis_eta_);
  tree_->Branch("jpsis_phi", &jpsis_phi_);
  tree_->Branch("jpsis_m"  , &jpsis_m_  );
  tree_->Branch("jpsis_lxy", &jpsis_lxy_);

  if ( isMC_ )
  {
    // JER
    jetsResUp_pt_ = new doubles; jetsResUp_eta_ = new doubles; jetsResUp_phi_ = new doubles; jetsResUp_m_ = new doubles;
    jetsResDn_pt_ = new doubles; jetsResDn_eta_ = new doubles; jetsResDn_phi_ = new doubles; jetsResDn_m_ = new doubles;
    jetsResUp_bTag_ = new doubles;
    jetsResDn_bTag_ = new doubles;

    tree_->Branch("jetsDn_pt" , jetsDn_pt_ );
    tree_->Branch("jetsDn_eta", jetsDn_eta_);
    tree_->Branch("jetsDn_phi", jetsDn_phi_);
    tree_->Branch("jetsDn_m"  , jetsDn_m_  );

    // Generator information
    genMuons_pt_  = new doubles;
    genMuons_eta_ = new doubles;
    genMuons_phi_ = new doubles;
    genMuons_m_   = new doubles;
    genMuons_Q_   = new ints   ;

    genElectrons_pt_  = new doubles;
    genElectrons_eta_ = new doubles;
    genElectrons_phi_ = new doubles;
    genElectrons_m_   = new doubles;
    genElectrons_Q_   = new ints   ;

    genNeutrinos_pt_  = new doubles;
    genNeutrinos_eta_ = new doubles;
    genNeutrinos_phi_ = new doubles;

    genJets_pt_  = new doubles;
    genJets_eta_ = new doubles;
    genJets_phi_ = new doubles;
    genJets_m_   = new doubles;

    tree_->Branch("genWeight", &genWeight_, "genWeight/D");

    tree_->Branch("pdf_id1", &pdf_id1_, "pdf_id1/I");
    tree_->Branch("pdf_id2", &pdf_id2_, "pdf_id2/I");
    tree_->Branch("pdf_x1" , &pdf_x1_ , "pdf_x1/D" );
    tree_->Branch("pdf_x2" , &pdf_x2_ , "pdf_x2/D" );
    tree_->Branch("pdf_q"  , &pdf_q_  , "pdf_q/D"  );

    tree_->Branch("genMuons_pt" , genMuons_pt_ );
    tree_->Branch("genMuons_eta", genMuons_eta_);
    tree_->Branch("genMuons_phi", genMuons_phi_);
    tree_->Branch("genMuons_m"  , genMuons_m_  );
    tree_->Branch("genMuons_Q"  , genMuons_Q_  );

    tree_->Branch("genElectrons_pt ", genElectrons_pt_ );
    tree_->Branch("genElectrons_eta", genElectrons_eta_);
    tree_->Branch("genElectrons_phi", genElectrons_phi_);
    tree_->Branch("genElectrons_m"  , genElectrons_m_  );
    tree_->Branch("genElectrons_Q"  , genElectrons_Q_  );

    tree_->Branch("genNeutrinos_pt" , genNeutrinos_pt_ );
    tree_->Branch("genNeutrinos_eta", genNeutrinos_eta_);
    tree_->Branch("genNeutrinos_phi", genNeutrinos_phi_);

    tree_->Branch("genJets_pt" , genJets_pt_ );
    tree_->Branch("genJets_eta", genJets_eta_);
    tree_->Branch("genJets_phi", genJets_phi_);
    tree_->Branch("genJets_m"  , genJets_m_  );
  }

}

KFlatTreeMaker::~KFlatTreeMaker()
{
}

void KFlatTreeMaker::endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
{
  for ( int i=0, n=eventCounterLabels_.size(); i<n; ++i )
  {
    edm::Handle<edm::MergeableCounter> eventCounterHandle;
    lumi.getByLabel(edm::InputTag(eventCounterLabels_.at(i)), eventCounterHandle);
    if ( !eventCounterHandle.isValid() ) continue;

    hEventCounter_->Fill(i+1, eventCounterHandle->value);
  }
}

void KFlatTreeMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  // Clear up
  electrons_pt_->clear();
  electrons_eta_->clear();
  electrons_phi_->clear();
  electrons_m_->clear();
  electrons_Q_->clear();
  electrons_type_->clear();
  electrons_relIso_->clear();

  electrons_mva_->clear();
  electrons_scEta_->clear();

  muons_pt_->clear();
  muons_eta_->clear();
  muons_phi_->clear();
  muons_m_->clear();
  muons_Q_->clear();
  muons_type_->clear();
  muons_relIso_->clear();

  jets_pt_->clear();
  jets_eta_->clear();
  jets_phi_->clear();
  jets_m_->clear();

  jetsUp_pt_->clear();
  jetsUp_eta_->clear();
  jetsUp_phi_->clear();
  jetsUp_m_->clear();

  jetsDn_pt_->clear();
  jetsDn_eta_->clear();
  jetsDn_phi_->clear();
  jetsDn_m_->clear();

  jets_bTag_->clear(); jetsUp_bTag_->clear(); jetsDn_bTag_->clear();

  jpsis_pt_->clear();
  jpsis_eta_->clear();
  jpsis_phi_->clear();
  jpsis_m_->clear();
  jpsis_lxy_->clear();

  if ( isMC_ and event.isRealData() ) isMC_ = false;

  if ( isMC_ )
  {
    jetsResUp_pt_->clear();
    jetsResUp_eta_->clear();
    jetsResUp_phi_->clear();
    jetsResUp_m_->clear();

    jetsResDn_pt_->clear();
    jetsResDn_eta_->clear();
    jetsResDn_phi_->clear();
    jetsResDn_m_->clear();

    jetsResUp_bTag_->clear();
    jetsResDn_bTag_->clear();

    genMuons_pt_ ->clear();
    genMuons_eta_->clear();
    genMuons_phi_->clear();
    genMuons_m_  ->clear();
    genMuons_Q_  ->clear();

    genElectrons_pt_ ->clear();
    genElectrons_eta_->clear();
    genElectrons_phi_->clear();
    genElectrons_m_  ->clear();
    genElectrons_Q_  ->clear();

    genNeutrinos_pt_ ->clear();
    genNeutrinos_eta_->clear();
    genNeutrinos_phi_->clear();

    genJets_pt_ ->clear();
    genJets_eta_->clear();
    genJets_phi_->clear();
    genJets_m_  ->clear();
  }

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);
  nVertex_ = vertexHandle->size();
  const reco::Vertex& pv = vertexHandle->at(0);

  if ( !isMC_ )
  {
    puWeight_ = puWeightUp_ = puWeightDn_ = 1.0;
  }
  else
  {
    edm::Handle<double> puWeightHandle, puWeightUpHandle, puWeightDnHandle;
    event.getByLabel(puWeightLabel_, puWeightHandle);
    event.getByLabel(edm::InputTag(puWeightLabel_.label(), "up"), puWeightUpHandle);
    event.getByLabel(edm::InputTag(puWeightLabel_.label(), "dn"), puWeightDnHandle);
    puWeight_ = *(puWeightHandle.product());
    puWeightUp_ = *(puWeightUpHandle.product());
    puWeightDn_ = *(puWeightDnHandle.product());
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

    electrons_pt_    ->push_back(e.pt()    );
    electrons_eta_   ->push_back(e.eta()   );
    electrons_phi_   ->push_back(e.phi()   );
    electrons_m_     ->push_back(e.mass()  );
    electrons_Q_     ->push_back(e.charge());
    electrons_type_  ->push_back(eType     );
    electrons_relIso_->push_back(e.userIso(2)); // rho corrected isolation

    electrons_mva_->push_back(e.mva());
    electrons_scEta_->push_back(scEta);
  }
  if ( electrons_pt_->size() < electronMinNumber_ ) return;

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

    muons_pt_->push_back(mu.pt()   );
    muons_eta_->push_back(mu.eta() );
    muons_phi_->push_back(mu.phi() );
    muons_m_->push_back(mu.mass()  );
    muons_Q_->push_back(mu.charge());
    muons_type_->push_back(muType  );
    muons_relIso_->push_back(mu.userIso(1)); // dBeta corrected isolation
  }
  if ( muons_pt_->size() < muonMinNumber_ ) return;

  edm::Handle<std::vector<pat::MET> > metHandle, metUpHandle, metDnHandle, metResUpHandle, metResDnHandle;
  event.getByLabel(metLabel_, metHandle);
  event.getByLabel(edm::InputTag(metLabel_.label(), "up"), metUpHandle);
  event.getByLabel(edm::InputTag(metLabel_.label(), "dn"), metDnHandle);
  met_pt_ = metHandle->at(0).pt();
  metUp_pt_ = metUpHandle->at(0).pt();
  metDn_pt_ = metDnHandle->at(0).pt();
  met_phi_ = metHandle->at(0).phi();
  metUp_phi_ = metUpHandle->at(0).phi();
  metDn_phi_ = metDnHandle->at(0).phi();
  if ( isMC_ )
  {
    event.getByLabel(edm::InputTag(metLabel_.label(), "resUp"), metResUpHandle);
    event.getByLabel(edm::InputTag(metLabel_.label(), "resDn"), metResDnHandle);
    metResUp_pt_ = metResUpHandle->at(0).pt();
    metResDn_pt_ = metResDnHandle->at(0).pt();
    metResUp_phi_ = metResUpHandle->at(0).phi();
    metResDn_phi_ = metResDnHandle->at(0).phi();
  }

  if ( isMC_ )
  {
    // Event weight and PDF stuffs
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    event.getByLabel(genEventInfoLabel_, genEventInfoHandle);

    const gen::PdfInfo* pdf = genEventInfoHandle->pdf();

    pdf_id1_ = pdf->id.first ;
    pdf_id2_ = pdf->id.second;
    pdf_x1_  = pdf->x.first  ;
    pdf_x2_  = pdf->x.second ;
    pdf_q_   = pdf->scalePDF ;
    //const double pdf_xPDF1 = pdf->xPDF.first, pdf_xPDF2 = pdf->xPDF.second;

    genWeight_ = genEventInfoHandle->weights().at(0);

    // Generator level objects
    edm::Handle<reco::GenParticleCollection> genHandle;
    event.getByLabel(genParticleLabel_, genHandle);
    edm::Handle<reco::GenJetCollection> genJetHandle;
    event.getByLabel(genJetLabel_, genJetHandle);

    // Find top quark from the genParticles
    for ( int i=0, n=genHandle->size(); i<n; ++i )
    {
      const reco::GenParticle& p = genHandle->at(i);
      if ( p.status() != 3 ) continue;
      const int charge = p.charge();

      switch(abs(p.pdgId()))
      {
        case 11:
          genElectrons_pt_ ->push_back(p.pt()  );
          genElectrons_eta_->push_back(p.eta() );
          genElectrons_phi_->push_back(p.phi() );
          genElectrons_m_  ->push_back(p.mass());
          genElectrons_Q_->push_back(charge); break;
        case 13:
          genMuons_pt_ ->push_back(p.pt()  );
          genMuons_eta_->push_back(p.eta() );
          genMuons_phi_->push_back(p.phi() );
          genMuons_m_  ->push_back(p.mass());
          genMuons_Q_->push_back(charge); break;
        case 12:
        case 14:
          genNeutrinos_pt_ ->push_back(p.pt() );
          genNeutrinos_eta_->push_back(p.eta());
          genNeutrinos_phi_->push_back(p.phi()); break;
        default: break;
      }
    }

    for ( int i=0, n=genJetHandle->size(); i<n; ++i )
    {
      const reco::GenJet& p = genJetHandle->at(i);
      if ( p.pt() < 20 or std::abs(p.eta()) > 2.5 ) continue;
      genJets_pt_ ->push_back(p.pt()  );
      genJets_eta_->push_back(p.eta() );
      genJets_phi_->push_back(p.phi() );
      genJets_m_  ->push_back(p.mass());
    }
  }

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
    jets_pt_ ->push_back(jet.pt());
    jets_eta_->push_back(jet.eta());
    jets_phi_->push_back(jet.phi());
    jets_m_  ->push_back(jet.mass());
    jets_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  for ( int i=0, n=jetUpHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetUpHandle->at(i);
    jetsUp_pt_ ->push_back(jet.pt());
    jetsUp_eta_->push_back(jet.eta());
    jetsUp_phi_->push_back(jet.phi());
    jetsUp_m_  ->push_back(jet.mass());
    jetsUp_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  for ( int i=0, n=jetDnHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetDnHandle->at(i);
    jetsDn_pt_ ->push_back(jet.pt());
    jetsDn_eta_->push_back(jet.eta());
    jetsDn_phi_->push_back(jet.phi());
    jetsDn_m_  ->push_back(jet.mass());
    jetsDn_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  if ( isMC_ )
  {
    event.getByLabel(edm::InputTag(jetLabel_.label(), "resUp"), jetResUpHandle);
    event.getByLabel(edm::InputTag(jetLabel_.label(), "resDn"), jetResDnHandle);
    for ( int i=0, n=jetResUpHandle->size(); i<n; ++i )
    {
      const pat::Jet& jet = jetResUpHandle->at(i);
      jetsResUp_pt_  ->push_back(jet.pt()  );
      jetsResUp_eta_ ->push_back(jet.eta() );
      jetsResUp_phi_ ->push_back(jet.phi() );
      jetsResUp_m_   ->push_back(jet.mass());
      jetsResUp_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
    }
    for ( int i=0, n=jetResDnHandle->size(); i<n; ++i )
    {
      const pat::Jet& jet = jetResDnHandle->at(i);
      jetsResDn_pt_  ->push_back(jet.pt()  );
      jetsResDn_eta_ ->push_back(jet.eta() );
      jetsResDn_phi_ ->push_back(jet.phi() );
      jetsResDn_m_   ->push_back(jet.mass());
      jetsResDn_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
    }
  }

  edm::Handle<std::vector<reco::VertexCompositeCandidate> > jpsiHandle;
  event.getByLabel(jpsiLabel_, jpsiHandle);
  edm::Handle<std::vector<double> > jpsiLxyHandle;
  event.getByLabel(edm::InputTag(jpsiLabel_.label(), "lxy"), jpsiLxyHandle);
  for ( int i=0, n=jpsiHandle->size(); i<n; ++i )
  {
    const reco::VertexCompositeCandidate& jpsiCand = jpsiHandle->at(i);
    jpsis_pt_ ->push_back(jpsiCand.pt()  );
    jpsis_eta_->push_back(jpsiCand.eta() );
    jpsis_phi_->push_back(jpsiCand.phi() );
    jpsis_m_  ->push_back(jpsiCand.mass());
    jpsis_lxy_->push_back(jpsiLxyHandle->at(i));
  }

  // Now put jets in current event to the event cache
  run_ = event.run();
  lumi_ = event.luminosityBlock();
  event_ = event.id().event();

  tree_->Fill();
}

DEFINE_FWK_MODULE(KFlatTreeMaker);

