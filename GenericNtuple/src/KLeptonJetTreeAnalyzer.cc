#include "KrAFT/GenericNtuple/interface/KLeptonJetTreeAnalyzer.h"

using namespace std;

KLeptonJetTreeAnalyzer::KLeptonJetTreeAnalyzer(const std::string modeName,
                                           const std::string inputFileName,
                                           const std::string outputFileName):
  KFlatTreeAnalyzerBase(modeName, inputFileName, outputFileName)
{

  outTree_->Branch("lepton_pt" , &lepton_pt_ , "lepton_pt/D" );
  //outTree_->Branch("lepton_eta", &lepton_eta_, "lepton_eta/D");
  //outTree_->Branch("lepton_phi", &lepton_phi_, "lepton_phi/D");
  outTree_->Branch("lepton_iso", &lepton_iso_, "lepton_iso/D");

  jets_pt_   = new doubles;
  jetsUp_pt_ = new doubles;
  jetsDn_pt_ = new doubles;
  outTree_->Branch("jets_pt" , &jets_pt_ );
  outTree_->Branch("jetsUp_pt", &jetsUp_pt_ );
  outTree_->Branch("jetsDn_pt", &jetsDn_pt_ );
  if ( isMC_ )
  {
    jetsResUp_pt_ = new doubles;
    jetsResDn_pt_ = new doubles;
    outTree_->Branch("jetsResUp_pt", &jetsUp_pt_ );
    outTree_->Branch("jetsResDn_pt", &jetsUp_pt_ );
  }
}

bool KLeptonJetTreeAnalyzer::analyze()
{
  if ( event_->iVars_["muons_type"]->size() == 0 ) return false;
  if ( event_->iVars_["muons_type"]->at(0) == 0 ) return false;
  // Select leptons
  LorentzVector leptonP4;
  const double muon_pt   = event_->fVars_["muons_pt" ]->at(0);
  const double muon_eta  = event_->fVars_["muons_eta"]->at(0);
  const double muon_phi  = event_->fVars_["muons_phi"]->at(0);
  const double muon_mass = event_->fVars_["muons_m"  ]->at(0);
  
  if ( muon_pt < 26 or abs(muon_eta) > 2.1 ) return false;
  if ( event_->fVars_["muons_relIso"]->at(0) >= 0.12 ) return false;

  leptonP4.SetPtEtaPhiM(muon_pt, muon_eta, muon_phi, muon_mass);

/*
  for ( int i=0, n=event_->jets_pt_->size(); i<n; ++i )
  {
    const double jetPt = event_->jets_pt_->at(i);
    jets_pt_->push_back(jetPt);

    LorentzVector jetP4;
    jetP4.SetPtEtaPhiM(jetPt, event_->jets_eta_->at(i), event_->jets_phi_->at(i), event_->jets_m_->at(i));
  }
  bTagCut_ = 0.898; // Tight cut
  bTagCut_ = 0.679; // Medium cut
  bTagCut_ = 0.244; // Loose cut

*/

  return true;
}

