#include "KrAFT/GenericNtuple/interface/KDileptonTreeReducer.h"

using namespace std;

KDileptonTreeReducer::KDileptonTreeReducer(const std::string modeName,
                                           const std::string inputFileName,
                                           const std::string outputFileName):
  KFlatTreeReducerBase(modeName, inputFileName, outputFileName)
{
  if ( modeName_ == "MuMu" )
  {
    mode_ = 1;
    leptons1_pt_  = leptons2_pt_  = event_->muons_pt_ ;
    leptons1_eta_ = leptons2_eta_ = event_->muons_eta_;
    leptons1_phi_ = leptons2_phi_ = event_->muons_phi_;
    leptons1_m_   = leptons2_m_   = event_->muons_m_  ;
    leptons1_Q_   = leptons2_Q_   = event_->muons_Q_  ;
    leptons1_iso_ = leptons2_iso_ = event_->muons_relIso_;
  }
  else if ( modeName_ == "ElEl")
  {
    mode_ = 2;
    leptons1_pt_  = leptons2_pt_  = event_->electrons_pt_ ;
    leptons1_eta_ = leptons2_eta_ = event_->electrons_eta_;
    leptons1_phi_ = leptons2_phi_ = event_->electrons_phi_;
    leptons1_m_   = leptons2_m_   = event_->electrons_m_  ;
    leptons1_Q_   = leptons2_Q_   = event_->electrons_Q_  ;
    leptons1_iso_ = leptons2_iso_ = event_->electrons_relIso_;
  }
  else
  {
    mode_ = 3;

    leptons1_pt_  = event_->muons_pt_ ;
    leptons1_eta_ = event_->muons_eta_;
    leptons1_phi_ = event_->muons_phi_;
    leptons1_m_   = event_->muons_m_  ;
    leptons1_Q_   = event_->muons_Q_  ;
    leptons1_iso_ = event_->muons_relIso_;

    leptons2_pt_  = event_->electrons_pt_ ;
    leptons2_eta_ = event_->electrons_eta_;
    leptons2_phi_ = event_->electrons_phi_;
    leptons2_m_   = event_->electrons_m_  ;
    leptons2_Q_   = event_->electrons_Q_  ;
    leptons2_iso_ = event_->electrons_relIso_;
  }

  outTree_->Branch("lepton1_pt" , &lepton1_pt_ , "lepton1_pt/D" );
  //outTree_->Branch("lepton1_eta", &lepton1_eta_, "lepton1_eta/D");
  //outTree_->Branch("lepton1_phi", &lepton1_phi_, "lepton1_phi/D");
  outTree_->Branch("lepton1_iso", &lepton1_iso_, "lepton1_iso/D");

  outTree_->Branch("lepton2_pt" , &lepton2_pt_ , "lepton2_pt/D" );
  //outTree_->Branch("lepton2_eta", &lepton2_eta_, "lepton2_eta/D");
  //outTree_->Branch("lepton2_phi", &lepton2_phi_, "lepton2_phi/D");
  outTree_->Branch("lepton2_iso", &lepton2_iso_, "lepton2_iso/D");

  outTree_->Branch("z_m" , &z_m_ , "z_m/D" );
  outTree_->Branch("z_pt", &z_pt_, "z_pt/D");
  outTree_->Branch("z_Q" , &z_Q_ , "z_Q/I" );

  outTree_->Branch("met_pt" , &met_pt_ , "met_pt/D" );
  outTree_->Branch("met_phi", &met_phi_, "met_phi/D");
  outTree_->Branch("metUp_pt" , &metUp_pt_ , "metUp_pt/D" );
  outTree_->Branch("metUp_phi", &metUp_phi_, "metUp_phi/D");
  outTree_->Branch("metDn_pt" , &metDn_pt_ , "metDn_pt/D" );
  outTree_->Branch("metDn_phi", &metDn_phi_, "metDn_phi/D");

  jets_pt_   = new doubles;
  jetsUp_pt_ = new doubles;
  jetsDn_pt_ = new doubles;
  outTree_->Branch("jets_pt" , &jets_pt_ );
  outTree_->Branch("jetsUp_pt", &jetsUp_pt_ );
  outTree_->Branch("jetsDn_pt", &jetsDn_pt_ );
  if ( isMC_ )
  {
    outTree_->Branch("metResUp_pt" , &metResUp_pt_ , "metResUp_pt/D" );
    outTree_->Branch("metResUp_phi", &metResUp_phi_, "metResUp_phi/D");
    outTree_->Branch("metResDn_pt" , &metResDn_pt_ , "metResDn_pt/D" );
    outTree_->Branch("metResDn_phi", &metResDn_phi_, "metResDn_phi/D");

    jetsResUp_pt_ = new doubles;
    jetsResDn_pt_ = new doubles;
    outTree_->Branch("jetsResUp_pt", &jetsUp_pt_ );
    outTree_->Branch("jetsResDn_pt", &jetsUp_pt_ );
  }
}

bool KDileptonTreeReducer::analyze()
{
  // Select leptons
  LorentzVector lepton1P4, lepton2P4, zP4;
  z_Q_ = -999;
  for ( int i=0, n=leptons1_pt_->size(); i<n; ++i )
  {
    if ( mode_ != 2 )
    {
      if ( event_->muons_type_->at(i) == 0 ) continue; // MuMu or MuEl
    }
    else if ( event_->electrons_type_->at(i) < 100 ) continue; // ElEl mode

    lepton1_pt_  = leptons1_pt_ ->at(i);
    lepton1_eta_ = leptons1_eta_->at(i);
    if ( lepton1_pt_ < 20 or std::abs(lepton1_eta_) > 2.5 ) continue;

    lepton1_phi_ = leptons1_phi_->at(i);
    const double lepton1M = leptons1_m_->at(i);
    const int lepton1Q = leptons1_Q_->at(i);

    lepton1P4.SetPtEtaPhiM(lepton1_pt_, lepton1_eta_, lepton1_phi_, lepton1M);
    lepton1_iso_ = leptons1_iso_->at(i);

    for ( int j=((mode_ == 3) ? 0 : i+1); j<n; ++j )
    {
      if ( mode_ == 1 )
      {
        if ( event_->muons_type_->at(j) == 0 ) continue; // MuMu mode
      }
      else if ( event_->electrons_type_->at(j) < 100 ) continue; // ElEl and MuEl

      lepton2_pt_  = leptons2_pt_->at(i);
      lepton2_eta_ = leptons2_eta_->at(i);
      if ( lepton2_pt_ < 20 or std::abs(lepton2_eta_) > 2.5 ) continue;

      lepton2_phi_ = leptons2_phi_->at(j);
      const double lepton2M = leptons2_m_->at(j);
      const int lepton2Q = leptons2_Q_->at(j);

      lepton2P4.SetPtEtaPhiM(lepton2_pt_, lepton2_eta_, lepton2_phi_, lepton2M);
      lepton2_iso_ = leptons2_iso_->at(j);

      zP4 = lepton1P4 + lepton2P4;
      z_Q_ = lepton1Q + lepton2Q;

      break;
    }
    break;
  }
  if ( z_Q_ == -999 ) return false;
  z_pt_ = zP4.Pt();
  z_m_  = zP4.M();

  met_pt_  = event_->met_pt_ ;
  met_phi_ = event_->met_phi_;
  metUp_pt_  = event_->metUp_pt_ ;
  metUp_phi_ = event_->metUp_phi_;
  metDn_pt_  = event_->metDn_pt_ ;
  metDn_phi_ = event_->metDn_phi_;
  if ( isMC_ )
  {
    metResUp_pt_  = event_->metResUp_pt_ ;
    metResUp_phi_ = event_->metResUp_phi_;
    metResDn_pt_  = event_->metResDn_pt_ ;
    metResDn_phi_ = event_->metResDn_phi_;
  }

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

