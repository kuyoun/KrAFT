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

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/OneToOne.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

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

  //void beginJob();
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool doMCMatch_;

  // Input objects
  edm::InputTag genLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag recoToGenJetMapLabel_;
  edm::InputTag genJetToPartonMapLabel_;
  std::string weightLabelStr_;
  edm::InputTag vertexLabel_;

  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  std::string jetLabelStr_;
  std::string metLabelStr_;
  std::string bTagType_;

  std::vector<std::string> eventCounterLabels_;

  unsigned int muonMinNumber_, electronMinNumber_;
  unsigned int muonMaxNumber_, electronMaxNumber_;

  double muonDz_, electronDz_;
  double jetLeptonDeltaR_;

  TH1F* hEventCounter_;

  // Output tree
  TTree* eventTree_;
  int run_, lumi_, event_;
  double weight_, weightUp_, weightDn_;
  int nVertex_;

  typedef math::XYZTLorentzVector LVec;
  typedef std::vector<LVec> LVecs;
  typedef std::vector<int> ints;
  typedef std::vector<double> doubles;

  LVecs*   muons_    , * electrons_    ;
  ints*    muons_Q_  , * electrons_Q_  ;
  doubles* muons_Iso_, * electrons_Iso_;
  LVecs* jets_, * jetsUp_, * jetsDn_;
  doubles* jets_bTag_, * jetsUp_bTag_, * jetsDn_bTag_;
  LVec* met_, * metUp_, * metDn_;

  // Generator level information
  //TTree* genTree_;
  // jet MC matching
  //std::vector<int> jets_motherId_;
  //std::vector<math::XYZTLorentzVector> genJet_;
  //std::vector<int> genJetMotherId_;
  //LVecs* genMuons_, * genElectrons_;
  //LVecs* genMuonNus_, * genElectronNus_;

};

KFlatTreeMaker::KFlatTreeMaker(const edm::ParameterSet& pset)
{
  //doMCMatch_ = pset.getParameter<bool>("doMCMatch");

  // Input labels
  //genLabel_ = pset.getParameter<edm::InputTag>("gen");
  //recoToGenJetMapLabel_ = pset.getParameter<edm::InputTag>("recoToGenJetMap");
  //genJetToPartonMapLabel_ = pset.getParameter<edm::InputTag>("genJetToPartonsMap");

  weightLabelStr_ = pset.getParameter<std::string>("weight");
  vertexLabel_ = pset.getParameter<edm::InputTag>("vertex");

  edm::ParameterSet electronPSet = pset.getParameter<edm::ParameterSet>("electron");
  electronMinNumber_ = electronPSet.getParameter<unsigned int>("minNumber");
  electronMaxNumber_ = electronPSet.getParameter<unsigned int>("maxNumber");
  electronLabel_ = electronPSet.getParameter<edm::InputTag>("src");

  edm::ParameterSet muonPSet = pset.getParameter<edm::ParameterSet>("muon");
  muonMinNumber_ = muonPSet.getParameter<unsigned int>("minNumber");
  muonMaxNumber_ = muonPSet.getParameter<unsigned int>("maxNumber");
  muonLabel_ = muonPSet.getParameter<edm::InputTag>("src");

  edm::ParameterSet jetPSet = pset.getParameter<edm::ParameterSet>("jet");
  jetLeptonDeltaR_ = jetPSet.getParameter<double>("leptonDeltaR");
  jetLabelStr_ = jetPSet.getParameter<std::string>("src");
  bTagType_ = jetPSet.getParameter<std::string>("bTagType");

  edm::ParameterSet metPSet = pset.getParameter<edm::ParameterSet>("met");
  metLabelStr_ = metPSet.getParameter<std::string>("src");

  // Event counter
  eventCounterLabels_ = pset.getParameter<std::vector<std::string> >("eventCounters");

  // Output histograms and tree
  edm::Service<TFileService> fs;
  hEventCounter_ = fs->make<TH1F>("hEventCounter", "Event counter", eventCounterLabels_.size(), 1, eventCounterLabels_.size()+1);
  for ( int i=0, n=eventCounterLabels_.size(); i<n; ++i )
  {
    hEventCounter_->GetXaxis()->SetBinLabel(i+1, eventCounterLabels_.at(i).c_str());
  }

  LVecs* jets_, * jetsUp_, * jetsDn_;
  doubles* jets_bTag_;
  LVec* met_, * metUp_, * metDn_;

  muons_     = new LVecs()  ; electrons_     = new LVecs()  ;
  muons_Q_   = new ints()   ; electrons_Q_   = new ints()   ;
  muons_Iso_ = new doubles(); electrons_Iso_ = new doubles();
  jets_   = new LVecs(); jets_bTag_   = new doubles(); met_   = new LVec();
  jetsUp_ = new LVecs(); jetsUp_bTag_ = new doubles(); metUp_ = new LVec();
  jetsDn_ = new LVecs(); jetsDn_bTag_ = new doubles(); metDn_ = new LVec();

  eventTree_ = fs->make<TTree>("event", "Mixed event tree");
  eventTree_->Branch("run", &run_, "run/I");
  eventTree_->Branch("lumi", &lumi_, "lumi/I");
  eventTree_->Branch("event", &event_, "event/I");

  eventTree_->Branch("weight", &weight_, "weight/D");
  eventTree_->Branch("weightUp", &weightUp_, "weightUp/D");
  eventTree_->Branch("weightDn", &weightDn_, "weightDn/D");
  eventTree_->Branch("nVertex", &nVertex_, "nVertex/I"); 

  const char* lvecsTypeName = "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >";
  eventTree_->Branch("electrons", lvecsTypeName, electrons_);
  eventTree_->Branch("electrons_Q", electrons_Q_);
  eventTree_->Branch("electrons_Iso", electrons_Iso_);

  eventTree_->Branch("muons", lvecsTypeName, muons_);
  eventTree_->Branch("muons_Q", muons_Q_);
  eventTree_->Branch("muons_Iso", muons_Iso_);

  eventTree_->Branch("jets", lvecsTypeName, jets_);
  eventTree_->Branch("jetsUp", lvecsTypeName, jetsUp_);
  eventTree_->Branch("jetsDn", lvecsTypeName, jetsDn_);
  eventTree_->Branch("jets_bTag", jets_bTag_);
  eventTree_->Branch("jetsUp_bTag", jetsUp_bTag_);
  eventTree_->Branch("jetsDn_bTag", jetsDn_bTag_);

  eventTree_->Branch("met"  , "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", met_  );
  eventTree_->Branch("metUp", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", metUp_);
  eventTree_->Branch("metDn", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", metDn_);

  //genTree_ = fs->make<TTree>("gen", "Gen level event tree");
  //genTree_->Branch("muons", &genMuons_);
  //genTree_->Branch("electrons", &genElectrons_);
  //genTree_->Branch("muonNus", &genMuonNus_);
  //genTree_->Branch("electronNus", &genElectronNus_);

  //genTree_->Branch("jets_motherId", &jets_motherId_);

}

KFlatTreeMaker::~KFlatTreeMaker()
{
}

void KFlatTreeMaker::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup)
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
  electrons_->clear(); electrons_Q_->clear(); electrons_Iso_->clear();
  muons_    ->clear(); muons_Q_    ->clear(); muons_Iso_    ->clear();
  jets_     ->clear(); jetsUp_     ->clear(); jetsDn_     ->clear();
  jets_bTag_->clear(); jetsUp_bTag_->clear(); jetsDn_bTag_->clear();
  //jets_motherId_->clear();

  //genMuons_->clear();
  //genElectrons_->clear();
  //genMuonNus_->clear();
  //genElectronNus_->clear();

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);
  nVertex_ = vertexHandle->size();
  const reco::Vertex& pv = vertexHandle->at(0);

  if ( event.isRealData() )
  {
    weight_ = weightUp_ = weightDn_ = 1.0;
  }
  else
  {
    edm::Handle<double> weightHandle, weightUpHandle, weightDnHandle;
    event.getByLabel(edm::InputTag(weightLabelStr_), weightHandle);
    event.getByLabel(edm::InputTag(weightLabelStr_, "up"), weightUpHandle);
    event.getByLabel(edm::InputTag(weightLabelStr_, "dn"), weightDnHandle);
    weight_ = *(weightHandle.product());
    weightUp_ = *(weightUpHandle.product());
    weightDn_ = *(weightDnHandle.product());
  }

  edm::Handle<std::vector<pat::Electron> > electronHandle;
  event.getByLabel(electronLabel_, electronHandle);
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    const pat::Electron& e = electronHandle->at(i);
    //if ( abs(e.dz(pv.position())) > electronDz_ ) continue;

    electrons_    ->push_back(e.p4());
    electrons_Q_  ->push_back(e.charge());
    electrons_Iso_->push_back(e.userIso(2)); // rho corrected isolation
  }
  if ( electrons_->size() < electronMinNumber_ ) return;
  if ( electrons_->size() > electronMaxNumber_ ) return;

  edm::Handle<std::vector<pat::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const pat::Muon& mu = muonHandle->at(i);
    //if ( abs(mu.dz(pv.position())) > muonDz_ ) continue;

    muons_    ->push_back(mu.p4());
    muons_Q_  ->push_back(mu.charge());
    muons_Iso_->push_back(mu.userIso(1)); // dBeta corrected isolation
  }
  if ( muons_->size() < muonMinNumber_ ) return;
  if ( muons_->size() > muonMaxNumber_ ) return;

  edm::Handle<std::vector<pat::MET> > metHandle, metUpHandle, metDnHandle;
  event.getByLabel(edm::InputTag(metLabelStr_), metHandle);
  event.getByLabel(edm::InputTag(metLabelStr_, "up"), metUpHandle);
  event.getByLabel(edm::InputTag(metLabelStr_, "dn"), metDnHandle);
  *met_ = metHandle->at(0).p4();
  *metUp_ = metUpHandle->at(0).p4();
  *metDn_ = metDnHandle->at(0).p4();

/*
  typedef edm::AssociationMap<edm::OneToMany<std::vector<reco::GenJet>, reco::GenParticleCollection> > GenJetToGenParticlesMap;
  typedef edm::AssociationMap<edm::OneToOne<std::vector<pat::Jet>, std::vector<reco::GenJet> > > RecoToGenJetMap;
  edm::Handle<GenJetToGenParticlesMap> genJetToPartonMapHandle;
  edm::Handle<RecoToGenJetMap> recoToGenJetMapHandle;

  // This while loop runs just for once, a "break" statement must be kept in the end of loop
  // It reduces nested loop
  while ( doMCMatch_ ) 
  {
    if ( event.isRealData() ) doMCMatch_ = false;

    edm::Handle<reco::GenParticleCollection> genHandle;
    event.getByLabel(genLabel_, genHandle);
    event.getByLabel(genJetToPartonMapLabel_, genJetToPartonMapHandle);
    event.getByLabel(recoToGenJetMapLabel_, recoToGenJetMapHandle);
    if ( !genHandle.isValid() or
         !genJetToPartonMapHandle.isValid() or 
         !recoToGenJetMapHandle.isValid() )
    {
      doMCMatch_ = false;
      break;
    }

    // Find top quark from the genParticles
    for ( int i=0, n=genHandle->size(); i<n; ++i )
    {
      const reco::GenParticle& p = genHandle->at(i);
      if ( p.status() != 3 ) continue;

      switch(abs(p.pdgId()))
      {
        case 11: genElectrons_.push_back(p.p4())  ; break;
        case 13: genMuons_.push_back(p.p4())      ; break;
        case 12: genElectronNus_.push_back(p.p4()); break;
        case 14: genMuonNus_.push_back(p.p4())    ; break;
        default: break;
      }
    }

    break;
  }
*/

  edm::Handle<std::vector<pat::Jet> > jetHandle;
  edm::Handle<std::vector<pat::Jet> > jetUpHandle;
  edm::Handle<std::vector<pat::Jet> > jetDnHandle;
  event.getByLabel(edm::InputTag(jetLabelStr_), jetHandle);
  event.getByLabel(edm::InputTag(jetLabelStr_, "up"), jetUpHandle);
  event.getByLabel(edm::InputTag(jetLabelStr_, "dn"), jetDnHandle);
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetHandle->at(i);

    jets_->push_back(jet.p4());
    jets_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  for ( int i=0, n=jetUpHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetUpHandle->at(i);
    jetsUp_->push_back(jet.p4());
    jetsUp_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }
  for ( int i=0, n=jetDnHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetDnHandle->at(i);
    jetsDn_->push_back(jet.p4());
    jetsDn_bTag_->push_back(jet.bDiscriminator(bTagType_.c_str()));
  }

/*
    int jetMotherId = 0;

    while ( doMCMatch_ )
    {
      edm::Ref<std::vector<pat::Jet> > jetRef(jetHandle, i);
      RecoToGenJetMap::const_iterator recoToGenJet = recoToGenJetMapHandle->find(jetRef);
      if ( recoToGenJet == recoToGenJetMapHandle->end() ) break;

      const edm::Ref<std::vector<reco::GenJet> >& genJet = recoToGenJet->val;
      GenJetToGenParticlesMap::const_iterator genJetToParton = genJetToPartonMapHandle->find(genJet);
      if ( genJetToParton == genJetToPartonMapHandle->end() ) break;

      const edm::RefVector<reco::GenParticleCollection>& genPartons = genJetToParton->val;
      for ( int j=0, m=genPartons.size(); j<m; ++j )
      {
        const int partonId = genPartons.at(j)->pdgId();
        // NOTE : Maybe there can be better way to set mother parton's id
        if ( abs(partonId) == 24 or partonId == 23 ) jetMotherId = partonId;
        else if ( jetMotherId == 0 and abs(partonId) == 6 ) jetMotherId = partonId;
        else if ( jetMotherId == 0 and partonId == 25 ) jetMotherId = partonId;
      }

      break;
    }

    jets_motherId_.push_back(jetMotherId);
  }
*/

  // Now put jets in current event to the event cache
  run_ = event.run();
  lumi_ = event.luminosityBlock();
  event_ = event.id().event();

  eventTree_->Fill();
//  genTree_->Fill();
}

DEFINE_FWK_MODULE(KFlatTreeMaker);

