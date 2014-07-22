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
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

template<typename Lepton1, typename Lepton2>
class DoubleLeptonAnalyzer : public edm::EDAnalyzer
{
public:
  DoubleLeptonAnalyzer(const edm::ParameterSet& pset);
  ~DoubleLeptonAnalyzer();

  void beginJob();
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  void endJob();
  void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);

  double getId(const pat::Muon& mu, const std::string& idName);
  double getId(const pat::Electron& ele, const std::string& idName);

  bool isMatched(const pat::Muon& mu, std::vector<const reco::GenParticle*> genLeptons);
  bool isMatched(const pat::Electron& ele, std::vector<const reco::GenParticle*> genLeptons);

private:
  edm::EDGetTokenT<edm::View<Lepton1> > lepton1Token_;
  edm::EDGetTokenT<edm::View<Lepton2> > lepton2Token_;
  edm::EDGetTokenT<edm::View<pat::MET> > metToken_;
  //edm::InputTag beamSpotLabel_;
  edm::EDGetTokenT<double> weightToken_;
  edm::EDGetTokenT<edm::MergeableCounter> eventCounterToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

private:
  bool isMC_;
  TH1F* hEvents_;
  TTree* tree_;
  int run_, lumi_, event_;
  double weight_;

  double pt1_, eta1_, phi1_;
  double pt2_, eta2_, phi2_;
  double iso1_, isoRho1_, isoDbeta1_;
  double iso2_, isoRho2_, isoDbeta2_;
  //double chIso1_, nhIso1_, phIso1_;
  //double chIso2_, nhIso2_, phIso2_;
  int q_;
  double m_;
  double met_;
  int genMatch1_;
  int genMatch2_;

  std::vector<std::string> idNames1_, idNames2_;
  std::vector<double> ids1_, ids2_;

  bool isSameLepton_;
};

template<typename Lepton1, typename Lepton2>
DoubleLeptonAnalyzer<Lepton1, Lepton2>::DoubleLeptonAnalyzer(const edm::ParameterSet& pset)
{
  edm::InputTag lepton1Label = pset.getParameter<edm::InputTag>("lepton1");
  edm::InputTag lepton2Label = pset.getParameter<edm::InputTag>("lepton2");
  isSameLepton_ = (lepton1Label == lepton2Label);
  lepton1Token_ = consumes<edm::View<Lepton1> >(lepton1Label);
  lepton2Token_ = consumes<edm::View<Lepton2> >(lepton2Label);
  metToken_ = consumes<edm::View<pat::MET> >(pset.getParameter<edm::InputTag>("met"));
  //beamSpotLabel_ = pset.getParameter<edm::InputTag>("beamSpot");
  weightToken_ = consumes<double>(pset.getParameter<edm::InputTag>("weight"));
  eventCounterToken_ = consumes<edm::MergeableCounter>(pset.getParameter<edm::InputTag>("eventCounter"));
  genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genPariclesLabel"));

  idNames1_ = pset.getParameter<std::vector<std::string> >("idNames1");
  idNames2_ = pset.getParameter<std::vector<std::string> >("idNames2");
  ids1_.resize(idNames1_.size());
  ids2_.resize(idNames2_.size());
}

template<typename Lepton1, typename Lepton2>
DoubleLeptonAnalyzer<Lepton1, Lepton2>::~DoubleLeptonAnalyzer()
{
}

template<typename Lepton1, typename Lepton2>
void DoubleLeptonAnalyzer<Lepton1, Lepton2>::endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::MergeableCounter> eventCounterHandle;
  if ( lumi.getByToken(eventCounterToken_, eventCounterHandle) )
  {
    hEvents_->Fill(0., double(eventCounterHandle->value));
  }
}

template<typename Lepton1, typename Lepton2>
void DoubleLeptonAnalyzer<Lepton1, Lepton2>::beginJob()
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "Lepton tree");

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/I");
  tree_->Branch("weight", &weight_, "weight/D");

  tree_->Branch("pt1" , &pt1_ , "pt1/D" );
  tree_->Branch("eta1", &eta1_, "eta1/D");
  tree_->Branch("phi1", &phi1_, "phi1/D");

  tree_->Branch("pt2" , &pt2_ , "pt2/D" );
  tree_->Branch("eta2", &eta2_, "eta2/D");
  tree_->Branch("phi2", &phi2_, "phi2/D");

  tree_->Branch("iso1", &iso1_, "iso1/D");
  tree_->Branch("iso2", &iso2_, "iso2/D");
  tree_->Branch("isoRho1", &isoRho1_, "isoRho1/D");
  tree_->Branch("isoRho2", &isoRho2_, "isoRho2/D");
  tree_->Branch("isoDbeta1", &isoDbeta1_, "isoDbeta1/D");
  tree_->Branch("isoDbeta2", &isoDbeta2_, "isoDbeta2/D");

  //tree_->Branch("chIso1", &chIso1_, "chIso1/D");
  //tree_->Branch("nhIso1", &nhIso1_, "nhIso1/D");
  //tree_->Branch("phIso1", &phIso1_, "phIso1/D");
  //tree_->Branch("chIso2", &chIso2_, "chIso2/D");
  //tree_->Branch("nhIso2", &nhIso2_, "nhIso2/D");
  //tree_->Branch("phIso2", &phIso2_, "phIso2/D");

  tree_->Branch("q", &q_, "q/I");
  tree_->Branch("m", &m_, "m/D");
  tree_->Branch("met", &met_, "met/D");
  tree_->Branch("genMatch1", &genMatch1_, "genMatch1/I");
  tree_->Branch("genMatch2", &genMatch2_, "genMatch2/I");

  for ( int i=0, n=idNames1_.size(); i<n; ++i )
  {
    tree_->Branch(Form("id1_%s", idNames1_[i].c_str()), &ids1_[i], Form("id1_%s/D", idNames1_[i].c_str()));
  }
  for ( int i=0, n=idNames2_.size(); i<n; ++i )
  {
    tree_->Branch(Form("id2_%s", idNames2_[i].c_str()), &ids2_[i], Form("id2_%s/D", idNames2_[i].c_str()));
  }

  hEvents_ = fs->make<TH1F>("hEvents", "Number of events", 1, 0, 1);
  hEvents_->GetXaxis()->SetBinLabel(1, "Total");
}

template<typename Lepton1, typename Lepton2>
void DoubleLeptonAnalyzer<Lepton1, Lepton2>::endJob()
{
}

template<typename Lepton1, typename Lepton2>
void DoubleLeptonAnalyzer<Lepton1, Lepton2>::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;
  using namespace edm;

  isMC_ = !event.isRealData();

  event_ = event.id().event();
  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();

  // Event weight
  edm::Handle<double> weightHandle;
  event.getByToken(weightToken_, weightHandle);
  weight_ = 1;
  if ( weightHandle.isValid() ) weight_ = *weightHandle;

  edm::Handle<edm::View<Lepton1> > lepton1Handle;
  event.getByToken(lepton1Token_, lepton1Handle);

  edm::Handle<edm::View<Lepton2> > lepton2Handle;
  event.getByToken(lepton2Token_, lepton2Handle);

  edm::Handle<edm::View<pat::MET> > metHandle;
  event.getByToken(metToken_, metHandle);
  met_ = metHandle->at(0).pt();

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(genParticlesToken_, genParticles);

  std::vector<const reco::GenParticle*> genLeptons;
  if ( genParticles.isValid() ) {
    for ( int i = 0, ngen = genParticles->size(); i < ngen; ++i ) {
      const reco::GenParticle& p = genParticles->at(i);
      if ( p.status() != 3 ) continue;
//      const int absPdgId = abs(p.pdgId());
	  if ( abs(p.pdgId()) == 11 or abs(p.pdgId()) == 13 ) genLeptons.push_back(&p);
    }
  }

  //edm::Handle<reco::BeamSpot> beamSpotHandle_;
  //event.getByToken(beamSpotLabel_, beamSpotHandle_); 

  const int nLep1 = lepton1Handle->size();
  const int nLep2 = lepton2Handle->size();
  for ( int iLep1 = 0; iLep1 < nLep1; ++iLep1 )
  {
    const Lepton1& lep1 = lepton1Handle->at(iLep1);
    if ( lep1.pt() < 20 ) continue;
    const double dxy1 = fabs(lep1.dB());
    if ( dxy1 > 0.04 ) continue;

    for ( int iLep2 = (isSameLepton_ ? iLep1+1 : 0); iLep2 < nLep2; ++iLep2 )
    {
      const Lepton2& lep2 = lepton2Handle->at(iLep2);
      if ( lep2.pt() < 20 ) continue;
      const double dxy2 = fabs(lep2.dB());
      if ( dxy2 > 0.04 ) continue;

      if ( lep1.p4() == lep2.p4() ) continue;

      pt1_ = lep1.pt(); eta1_ = lep1.eta(); phi1_ = lep1.phi();
      pt2_ = lep2.pt(); eta2_ = lep2.eta(); phi2_ = lep2.phi();
      iso1_ = lep1.userIsolation("User1Iso");
      iso2_ = lep2.userIsolation("User1Iso");
      isoDbeta1_ = lep1.userIsolation("User2Iso");
      isoDbeta2_ = lep2.userIsolation("User2Iso");
      isoRho1_ = lep1.userIsolation("User3Iso");
      isoRho2_ = lep2.userIsolation("User3Iso");

      q_ = lep1.charge() + lep2.charge();
      m_ = (lep1.p4()+lep2.p4()).M();

      for ( int iId=0, nId=idNames1_.size(); iId<nId; ++iId ) ids1_[iId] = getId(lep1, idNames1_[iId]);
      for ( int iId=0, nId=idNames2_.size(); iId<nId; ++iId ) ids2_[iId] = getId(lep2, idNames2_[iId]);

// doMCMatching
      genMatch1_ = isMatched(lep1, genLeptons);
      genMatch2_ = isMatched(lep2, genLeptons);

      tree_->Fill();
    }
  }
}

template<typename Lepton1, typename Lepton2>
double DoubleLeptonAnalyzer<Lepton1, Lepton2>::getId(const pat::Electron& ele, const std::string& idName)
{
  return ele.electronID(idName);
}

template<typename Lepton1, typename Lepton2>
double DoubleLeptonAnalyzer<Lepton1, Lepton2>::getId(const pat::Muon& mu, const std::string& idName)
{
  return mu.muonID(idName);
}

template<typename Lepton1, typename Lepton2>
bool DoubleLeptonAnalyzer<Lepton1, Lepton2>::isMatched(const pat::Muon& mu, std::vector<const reco::GenParticle*> genLeptons)
{
  double minDR = 0.3;
  const reco::GenParticle* matchedGenParticle = 0;
  for ( int i=0, n=genLeptons.size(); i<n; ++i )
  {
    const reco::GenParticle* p = genLeptons.at(i);

    if ( abs(p->pdgId()) != 13 ) continue;

    const double dR = deltaR(mu, *p);

    if ( dR > minDR ) continue;

    minDR = dR;
    matchedGenParticle = p;
  }

  if ( matchedGenParticle ) return true;
  return false;
}
template<typename Lepton1, typename Lepton2>
bool DoubleLeptonAnalyzer<Lepton1, Lepton2>::isMatched(const pat::Electron& ele, std::vector<const reco::GenParticle*> genLeptons)
{
  double minDR = 0.3;
  const reco::GenParticle* matchedGenParticle = 0;
  for ( int i=0, n=genLeptons.size(); i<n; ++i )
  {
    const reco::GenParticle* p = genLeptons.at(i);

    if ( abs(p->pdgId()) != 11 ) continue;

    const double dR = deltaR(ele, *p);

    if ( dR > minDR ) continue;

    minDR = dR;
    matchedGenParticle = p;
  }

  if ( matchedGenParticle ) return true;
  return false;
}

typedef DoubleLeptonAnalyzer<pat::Electron, pat::Electron> DoubleElectronAnalyzer;
typedef DoubleLeptonAnalyzer<pat::Muon, pat::Muon> DoubleMuonAnalyzer;
typedef DoubleLeptonAnalyzer<pat::Muon, pat::Electron> MuEGAnalyzer;

DEFINE_FWK_MODULE(DoubleElectronAnalyzer);
DEFINE_FWK_MODULE(DoubleMuonAnalyzer);
DEFINE_FWK_MODULE(MuEGAnalyzer);

