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
#include "KrAFT/RecoSelectorTools/interface/Types.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>
#include <TLorentzVector.h>

using namespace std;
class KGenericNtupleMaker : public edm::EDAnalyzer
{
public:
  KGenericNtupleMaker(const edm::ParameterSet& pset);
  virtual ~KGenericNtupleMaker();

  void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  typedef std::vector<pat::Electron> Electrons;
  typedef std::vector<pat::Muon> Muons;
  typedef std::vector<pat::MET> METs;
  typedef std::vector<pat::Jet> Jets;
  typedef edm::RefVector<Jets> JetRefs;
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> VString;

private:
  int findMother(const reco::GenParticle* p, std::vector<const reco::GenParticle*>& genParticles);

private:
  // Input objects
  edm::InputTag genEventInfoLabel_;
  edm::InputTag genParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag puWeightLabel_;
  edm::InputTag vertexLabel_;
  edm::InputTag pdfWeightsLabel_;

  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag metLabel_;
  edm::InputTag uncLabel_;
  edm::InputTag jpsiLabel_;
  std::string bTagType_;

  VString eventCounterLabels_;

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

  edm::ParameterSet jetMETPSet = pset.getParameter<edm::ParameterSet>("jetMET");
  jetLabel_ = jetMETPSet.getParameter<edm::InputTag>("jet");
  metLabel_ = jetMETPSet.getParameter<edm::InputTag>("met");
  uncLabel_ = jetMETPSet.getParameter<edm::InputTag>("unc");
  jetMinNumber_ = jetMETPSet.getParameter<unsigned int>("minNumber");
  jetLeptonDeltaR_ = jetMETPSet.getParameter<double>("leptonDeltaR");
  bTagType_ = jetMETPSet.getParameter<std::string>("bTagType");

  edm::ParameterSet jpsiPSet = pset.getParameter<edm::ParameterSet>("jpsi");
  jpsiLabel_ = jpsiPSet.getParameter<edm::InputTag>("src");

  // Event counter
  eventCounterLabels_ = pset.getParameter<VString>("eventCounters");

  if ( isMC_ )
  {
    genEventInfoLabel_ = pset.getParameter<edm::InputTag>("genEventInfo");
    pdfWeightsLabel_ = pset.getParameter<edm::InputTag>("pdfWeights");
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
  //std::cout<<"Event : "<<event.id()<<std::endl;
  using namespace std;

  fevent_->clear();

  if ( isMC_ and event.isRealData() ) isMC_ = false;

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);
  const reco::Vertex& pv = vertexHandle->at(0);

//  fevent_->nVertex_ = vertexHandle->size();
  fevent_->nVertex_ = 0;
  for ( int i=0, n=vertexHandle->size(); i<n; ++i )
  {
    const reco::Vertex& v = vertexHandle->at(i);
    if ( !v.isFake() and v.ndof()>4 and std::abs(v.z())<=24.0 and v.position().Rho()<=2.0 ) ++(fevent_->nVertex_);
  }

  if ( !isMC_ )
  {
    fevent_->puWeight_   = 1.0;
    fevent_->puWeightUp_ = 1.0;
    fevent_->puWeightDn_ = 1.0;
    fevent_->nPileup_    = -1;
  }
  else
  {
    const std::string& puWeightName = puWeightLabel_.label();

    edm::Handle<double> puWeightHandle, puWeightUpHandle, puWeightDnHandle;
    event.getByLabel(puWeightLabel_, puWeightHandle);
    event.getByLabel(edm::InputTag(puWeightName, "up"), puWeightUpHandle);
    event.getByLabel(edm::InputTag(puWeightName, "dn"), puWeightDnHandle);
    fevent_->puWeight_   = *(puWeightHandle.product()  );
    fevent_->puWeightUp_ = *(puWeightUpHandle.product());
    fevent_->puWeightDn_ = *(puWeightDnHandle.product());

    edm::Handle<int> nPileupHandle;
    event.getByLabel(edm::InputTag(puWeightName, "nTrueInteraction"), nPileupHandle);

    fevent_->nPileup_ = *nPileupHandle;

    edm::Handle<std::vector<double> > pdfWeightsHandle;
    event.getByLabel(pdfWeightsLabel_, pdfWeightsHandle);
    std::copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), std::back_inserter(*fevent_->pdfWeights_));
  }

  edm::Handle<Electrons> electronHandle;
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

    fevent_->electrons_mva_->push_back(mva);
    fevent_->electrons_scEta_->push_back(scEta);

    int qConsistent = 0;
    if ( e.isGsfCtfScPixChargeConsistent() ) qConsistent = 3;
    else if ( e.isGsfScPixChargeConsistent() ) qConsistent = 2;
    else if ( e.isGsfCtfChargeConsistent() ) qConsistent = 1;
    fevent_->electrons_qConsistent_->push_back(qConsistent);
  }
	//std::cout<<"MinNumber"<<electronMinNumber_<<"    electron size : "<<fevent_->electrons_pt_->size()<<std::endl;
  if ( fevent_->electrons_pt_->size() < electronMinNumber_ ) 
  { 
	//std::cout<<"Wow! electron minimum out!!"<<std::endl; 
	return; 
  }

  edm::Handle<Muons> muonHandle;
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
	//std::cout<<"MinNumber"<<muonMinNumber_<<"    muon size : "<<fevent_->muons_pt_->size()<<std::endl;
  if ( fevent_->muons_pt_->size() < muonMinNumber_ ) 
  { 
	//std::cout<<"Wow! muon minimum out!!"<<std::endl; 
	return; 
  }

  edm::Handle<METs> metHandle, metJESUpHandle, metJESDnHandle;
  edm::Handle<METs> metJERHandle, metJERUpHandle, metJERDnHandle;
  event.getByLabel(metLabel_, metHandle);
  event.getByLabel(edm::InputTag(uncLabel_.label(), "up"), metJESUpHandle);
  event.getByLabel(edm::InputTag(uncLabel_.label(), "dn"), metJESDnHandle);
  fevent_->met_pt_ = metHandle->at(0).pt();
  fevent_->metJESUp_pt_ = metJESUpHandle->at(0).pt();
  fevent_->metJESDn_pt_ = metJESDnHandle->at(0).pt();
  fevent_->met_phi_ = metHandle->at(0).phi();
  fevent_->metJESUp_phi_ = metJESUpHandle->at(0).phi();
  fevent_->metJESDn_phi_ = metJESDnHandle->at(0).phi();
  if ( isMC_ )
  {
    event.getByLabel(edm::InputTag(uncLabel_.label(), "res"), metJERHandle);
    event.getByLabel(edm::InputTag(uncLabel_.label(), "resUp"), metJERUpHandle);
    event.getByLabel(edm::InputTag(uncLabel_.label(), "resDn"), metJERDnHandle);
    fevent_->metJER_pt_ = metJERHandle->at(0).pt();
    fevent_->metJERUp_pt_ = metJERUpHandle->at(0).pt();
    fevent_->metJERDn_pt_ = metJERDnHandle->at(0).pt();
    fevent_->metJER_phi_ = metJERHandle->at(0).phi();
    fevent_->metJERUp_phi_ = metJERUpHandle->at(0).phi();
    fevent_->metJERDn_phi_ = metJERDnHandle->at(0).phi();
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
    std::vector<const reco::GenParticle*> genParticlesToStore;
    for ( int i=0, n=genHandle->size(); i<n; ++i )
    {
      const reco::GenParticle& p = genHandle->at(i);
      if ( p.status() != 3 ) continue;
      if ( p.pt() == 0 ) continue;
      genParticlesToStore.push_back(&p);
    }
    for ( int i=0, n=genParticlesToStore.size(); i<n; ++i )
    {
      const reco::GenParticle* p = genParticlesToStore[i];
      //const int charge = p.charge();
      const int mother = findMother(p, genParticlesToStore);

      fevent_->genParticles_pt_ ->push_back(p->pt()  );
      fevent_->genParticles_eta_->push_back(p->eta() );
      fevent_->genParticles_phi_->push_back(p->phi() );
      fevent_->genParticles_m_  ->push_back(p->mass());
      fevent_->genParticles_pdgId_->push_back(p->pdgId());
      fevent_->genParticles_mother_->push_back(mother);
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

  // Do jets
  edm::Handle<JetRefs> jetHandle;
  edm::Handle<pat::JetToValue> fJESUpHandle, fJESDnHandle;
  edm::Handle<pat::JetToValue> fJERHandle;
  edm::Handle<pat::JetToValue> fJERUpHandle, fJERDnHandle;
  event.getByLabel(jetLabel_, jetHandle);
  event.getByLabel(edm::InputTag(uncLabel_.label(), "up"), fJESUpHandle);
  event.getByLabel(edm::InputTag(uncLabel_.label(), "dn"), fJESDnHandle);
  if ( isMC_ )
  {
    event.getByLabel(edm::InputTag(uncLabel_.label(), "res"), fJERHandle);
    event.getByLabel(edm::InputTag(uncLabel_.label(), "resUp"), fJERUpHandle);
    event.getByLabel(edm::InputTag(uncLabel_.label(), "resDn"), fJERDnHandle);
  }
  double fJER = 1, fJERUp = 1, fJERDn = 1;
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    edm::Ref<Jets> jetRef = jetHandle->at(i);
    const pat::Jet& jet = *jetRef;
    const double jetEta = jet.eta();
    if ( std::abs(jetEta) > 2.5 ) continue;

    const double fJESUp = (*fJESUpHandle)[jetRef];
    const double fJESDn = (*fJESDnHandle)[jetRef];
    double maxPtScale = max(max(1., fJESUp), fJESDn);
    if ( isMC_ )
    {
      fJER   = (*fJERHandle)[jetRef];
      fJERUp = (*fJERUpHandle)[jetRef];
      fJERDn = (*fJERDnHandle)[jetRef];
      maxPtScale = max(max(max(fJER, fJERUp), fJERDn), maxPtScale);
    }
    const double jetPt = jet.pt();
    if ( jetPt*maxPtScale < 30 ) continue;

    fevent_->jets_pt_ ->push_back(jetPt);
    fevent_->jets_eta_->push_back(jetEta);
    fevent_->jets_phi_->push_back(jet.phi());
    fevent_->jets_m_  ->push_back(jet.mass());
    fevent_->jets_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
    fevent_->jets_partonflavor_->push_back(jet.partonFlavour());

    fevent_->jets_JESUp_->push_back(fJESUp);
    fevent_->jets_JESDn_->push_back(fJESDn);
    if ( isMC_ )
    {
      fevent_->jets_JER_->push_back(fJER);
      fevent_->jets_JERUp_->push_back(fJERUp);
      fevent_->jets_JERDn_->push_back(fJERDn);
    }
  }
	//std::cout<<"MinNumber"<<jetMinNumber_<<"    jet size : "<<fevent_->jets_pt_->size()<<std::endl;
  if ( fevent_->jets_pt_->size() < jetMinNumber_ ) 
  { 
	//std::cout<<"Wow! jet minimum out!!"<<std::endl; 
	return; 
  }

  // Do Jpsi
  edm::Handle<std::vector<reco::VertexCompositeCandidate> > jpsiHandle;
  event.getByLabel(jpsiLabel_, jpsiHandle);
  edm::Handle<doubles> jpsiLxyHandle;
  edm::Handle<doubles> jpsiL3DHandle;
  event.getByLabel(edm::InputTag(jpsiLabel_.label(), "lxy"), jpsiLxyHandle);
  event.getByLabel(edm::InputTag(jpsiLabel_.label(), "l3D"), jpsiL3DHandle);
  for ( int i=0, n=jpsiHandle->size(); i<n; ++i )
  {
    //std::cout<<"Start Loop "<<std::endl;
    const reco::VertexCompositeCandidate& jpsiCand = jpsiHandle->at(i);
    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(jpsiCand.daughter(0));
    const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(jpsiCand.daughter(1));

		TLorentzVector mu1(muon1->px(), muon1->py(), muon1->pz(), muon1->energy());
		TLorentzVector mu2(muon2->px(), muon2->py(), muon2->pz(), muon2->energy());
		TLorentzVector sumMu = mu1+mu2;
		TVector3 sumBoost( sumMu.BoostVector());
		mu1.Boost(-1*sumBoost);
		mu2.Boost(-1*sumBoost);

    fevent_->jpsis_pt_ ->push_back(jpsiCand.pt()  );
    fevent_->jpsis_eta_->push_back(jpsiCand.eta() );
    fevent_->jpsis_phi_->push_back(jpsiCand.phi() );
    fevent_->jpsis_m_  ->push_back(jpsiCand.mass());
    fevent_->jpsis_lxy_->push_back(jpsiLxyHandle->at(i));
    fevent_->jpsis_l3D_->push_back(jpsiL3DHandle->at(i));
		fevent_->jpsis_vProb_->push_back(TMath::Prob(  jpsiCand.vertexChi2(),(int)jpsiCand.vertexNdof()));
    fevent_->jpsis_cos_->push_back( TMath::Cos( mu1.Angle(mu2.Vect() )) );
		fevent_->jpsis_vProb_->push_back(TMath::Prob(  jpsiCand.vertexChi2(),(int)jpsiCand.vertexNdof()));

    fevent_->jpsis_pt1_ ->push_back(muon1->pt() );
    fevent_->jpsis_eta1_->push_back(muon1->eta());
    fevent_->jpsis_phi1_->push_back(muon1->phi());
		fevent_->jpsis_cos1_->push_back( TMath::Cos(mu1.Theta()));

    fevent_->jpsis_pt2_ ->push_back(muon2->pt() );
    fevent_->jpsis_eta2_->push_back(muon2->eta());
    fevent_->jpsis_phi2_->push_back(muon2->phi());
		fevent_->jpsis_cos2_->push_back( TMath::Cos(mu2.Theta()));

    reco::TrackRef muonTrack1 = muon1->improvedMuonBestTrack();
    reco::TrackRef muonTrack2 = muon2->improvedMuonBestTrack();
    fevent_->jpsis_nPixHits1_->push_back(muonTrack1->hitPattern().numberOfValidPixelHits());
    fevent_->jpsis_nPixHits2_->push_back(muonTrack2->hitPattern().numberOfValidPixelHits());
  }

  // Now put jets in current event to the event cache
  fevent_->run_ = event.run();
  fevent_->lumi_ = event.luminosityBlock();
  fevent_->event_ = event.id().event();

  fevent_->tree_->Fill();
}

int KGenericNtupleMaker::findMother(const reco::GenParticle* p, std::vector<const reco::GenParticle*>& genParticles)
{
  if ( !p ) return -1;
  for ( int i=0, n=genParticles.size(); i<n; ++i )
  {
    if ( p == genParticles[i] ) return i;
  }
  return findMother(dynamic_cast<const reco::GenParticle*>(p->mother()), genParticles);
}

DEFINE_FWK_MODULE(KGenericNtupleMaker);
