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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "KrAFT/GeneratorTools/interface/Types.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Common/interface/View.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <memory>
#include <vector>
#include <string>

using namespace std;

class PseudoTopAnalyzer : public edm::EDAnalyzer
{
public:
  PseudoTopAnalyzer(const edm::ParameterSet& pset);
  ~PseudoTopAnalyzer() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  const reco::GenParticle* match(const math::XYZTLorentzVector& p4, std::vector<const reco::GenParticle*>& v)
  {
    const reco::GenParticle* matched = 0;
    double matchedDR = 0.5;
    for ( int i=0, n=v.size(); i<n; ++i )
    {
      const reco::GenParticle* cand = v.at(i);
      const double dR = deltaR(p4, cand->p4());
      if ( dR < matchedDR )
      {
        matched = cand;
        matchedDR = dR;
      }
    }
    return matched;
  };

  const reco::GenJet* match(const math::XYZTLorentzVector& p4, const reco::GenJetCollection& v)
  {
    const reco::GenJet* matched = 0;
    double matchedDR = 0.5;
    for ( int i=0, n=v.size(); i<n; ++i )
    {
      const reco::GenJet* cand = &v.at(i);
      const double dR = deltaR(p4, cand->p4());
      if ( dR < matchedDR )
      {
        matched = cand;
        matchedDR = dR;
      }
    }
    return matched;
  };

private:
  edm::InputTag genParticlesLabel_;
  edm::InputTag genJetsLabel_;
  edm::InputTag pseudoTopLabel_;

  typedef TH1F* H1P;
  typedef TH2F* H2P;
  H1P hPlepN_, hDlepN_, hBlepN_;
  H2P hPlepPt_DlepPt_, hPlepPt_BlepPt_, hDlepPt_BlepPt_;
  H1P hDlepPlepDPt_, hBlepPlepDPt_, hDlepBlepDPt_;
  H1P hDlepPlepDR_, hBlepPlepDR_, hDlepBlepDR_;

};

PseudoTopAnalyzer::PseudoTopAnalyzer(const edm::ParameterSet& pset)
{
  genParticlesLabel_ = pset.getParameter<edm::InputTag>("genParticles");
  genJetsLabel_ = pset.getParameter<edm::InputTag>("genJets");
  pseudoTopLabel_ = pset.getParameter<edm::InputTag>("pseudoTop");

  edm::Service<TFileService> fs;
  // P* : Parton*, D* : Dressed*, B* : Bare*
  hPlepN_ = fs->make<TH1F>("hPLepN", "Parton lepton multiplicity;Lepton multiplicity;Events", 10, 0, 10);
  hDlepN_ = fs->make<TH1F>("hDLepN", "Dressed lepton multiplicity;Lepton multiplicity;Events", 10, 0, 10);
  hBlepN_ = fs->make<TH1F>("hBLepN", "Bare lepton multiplicity;Lepton multiplicity;Events", 10, 0, 10);

  hPlepPt_DlepPt_ = fs->make<TH2F>("hPlepPt_DlepPt", "Parton lepton p_{T} vs Dressed lepton p_{T};Parton lepton p_{T} (GeV/c);Dressed lepton p_{T}", 100, 0, 100, 100, 0, 100);
  hPlepPt_BlepPt_ = fs->make<TH2F>("hPlepPt_BlepPt", "Parton lepton p_{T} vs Bare lepton p_{T};Parton lepton p_{T} (GeV/c);Bare lepton p_{T}", 100, 0, 100, 100, 0, 100);
  hDlepPt_BlepPt_ = fs->make<TH2F>("hDlepPt_BlepPt", "Dressed lepton p_{T} vs Bare lepton p_{T};Dressed lepton p_{T} (GeV/c);Bare lepton p_{T}", 100, 0, 100, 100, 0, 100);

  hDlepPlepDPt_ = fs->make<TH1F>("hPLepDLepDPt", "lepton#Delta p_{T}(Dressed, Parton);#Delta p_{T} (GeV/c);Entries", 100, -5, 5);
  hBlepPlepDPt_ = fs->make<TH1F>("hPLepBLepDPt", "lepton#Delta p_{T}(Bare, Parton);#Delta p_{T} (GeV/c);Entries", 100, -5, 5);
  hDlepBlepDPt_ = fs->make<TH1F>("hDLepBLepDPt", "lepton#Delta p_{T}(Bare, Dressed);#Delta p_{T} (GeV/c);Entries", 100, -5, 5);

  hDlepPlepDR_ = fs->make<TH1F>("hDlepPlepDR", "lepton matching #Delta R(Dressed, Parton);#Delta R;Entries", 100, 0, 0.1);
  hBlepPlepDR_ = fs->make<TH1F>("hBlepPlepDR", "lepton matching #Delta R(Bare, Parton);#Delta R;Entries", 100, 0, 0.1);
  hDlepBlepDR_ = fs->make<TH1F>("hDlepBlepDR", "lepton matching #Delta R(Bare, Dressed);#Delta R;Entries", 100, 0, 0.1);


}

void PseudoTopAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  event.getByLabel(genParticlesLabel_, genParticlesHandle);

  edm::Handle<reco::GenJetCollection> genJetsHandle;
  event.getByLabel(genJetsLabel_, genJetsHandle);

  edm::Handle<reco::GenParticleCollection> ptopNeutrinosHandle;
  event.getByLabel(edm::InputTag(pseudoTopLabel_.label(), "neutrinos"), ptopNeutrinosHandle);

  edm::Handle<reco::GenJetCollection> ptopLeptonsHandle;
  event.getByLabel(edm::InputTag(pseudoTopLabel_.label(), "leptons"), ptopLeptonsHandle);

  edm::Handle<reco::GenJetCollection> ptopJetsHandle;
  event.getByLabel(edm::InputTag(pseudoTopLabel_.label(), "jets"), ptopJetsHandle);

  // Collect parton level objects
  std::vector<const reco::GenParticle*> genLeptons;
  std::vector<const reco::GenParticle*> genNeutrinos;
  std::vector<const reco::GenParticle*> genBquarks;
  std::vector<const reco::GenParticle*> bareLeptons, bareNeutrinos;
  for ( int i=0, n=genParticlesHandle->size(); i<n; ++i )
  {
    const reco::GenParticle* p = &genParticlesHandle->at(i);
    const int absPdgId = abs(p->pdgId());

    if ( p->status() == 1 )
    {
      switch ( absPdgId )
      {
        case 11: case 13:
          bareLeptons.push_back(p);
          break;
        case 12: case 14: case 16:
          bareNeutrinos.push_back(p);
          break;
      }
    }
    if ( p->status() != 3 ) continue;

    switch ( absPdgId )
    {
      case 11: case 13:
        genLeptons.push_back(p);
        break;
      case 12: case 14: case 16:
        genNeutrinos.push_back(p);
        break;
      case 5:
        genBquarks.push_back(p);
        break;
    }
  }

  // Fill leptons
  hPlepN_->Fill(genLeptons.size());
  hDlepN_->Fill(ptopLeptonsHandle->size());
  hBlepN_->Fill(bareLeptons.size());
  for ( int i=0, n=genLeptons.size(); i<n; ++i )
  {
    const reco::GenParticle* partonLep = genLeptons[i];
    const reco::GenJet* matchedDlep = match(partonLep->p4(), *(ptopLeptonsHandle.product()));
    const reco::GenParticle* matchedBlep = match(partonLep->p4(), bareLeptons);

    if ( matchedDlep )
    {
      hDlepPlepDR_->Fill(deltaR(partonLep->p4(), matchedDlep->p4()));
      hDlepPlepDPt_->Fill(matchedDlep->pt()-partonLep->pt());
      hPlepPt_DlepPt_->Fill(partonLep->pt(), matchedDlep->pt());
    }
    if ( matchedBlep )
    {
      hBlepPlepDR_->Fill(deltaR(partonLep->p4(), matchedBlep->p4()));
      hBlepPlepDPt_->Fill(matchedBlep->pt()-partonLep->pt());
      hPlepPt_BlepPt_->Fill(partonLep->pt(), matchedBlep->pt());
    }
  }
  for ( int i=0, n=ptopLeptonsHandle->size(); i<n; ++i )
  {
    const reco::GenJet* dressedLep = &ptopLeptonsHandle->at(i);
    const reco::GenParticle* matchedBlep = match(dressedLep->p4(), bareLeptons);

    if ( matchedBlep )
    {
      hDlepBlepDR_->Fill(deltaR(dressedLep->p4(), matchedBlep->p4()));
      hDlepBlepDPt_->Fill(dressedLep->pt()-matchedBlep->pt());
      hDlepPt_BlepPt_->Fill(dressedLep->pt(), matchedBlep->pt());
    }
  }

  // Fill Neutrinos
}

DEFINE_FWK_MODULE(PseudoTopAnalyzer);
