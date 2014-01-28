#include "KrAFT/GenericNtuple/interface/KDileptonTreeReducer.h"

using namespace std;

KDileptonTreeReducer::KDileptonTreeReducer(const std::string modeName,
                                           const std::string inputFileName,
                                           const std::string outputFileName):
  KFlatTreeReducerBase(modeName, inputFileName, outputFileName)
{
  if ( !event_ ) return;

  outTree_->Branch("nVertex", &nVertex_, "nVertex/I");

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
  outTree_->Branch("z_dphi", &z_dphi_, "z_dphi/D");

  outTree_->Branch("met_pt" , &met_pt_ , "met_pt/D" );
  outTree_->Branch("met_phi", &met_phi_, "met_phi/D");
  outTree_->Branch("metJESUp_pt" , &metJESUp_pt_ , "metJESUp_pt/D" );
  outTree_->Branch("metJESUp_phi", &metJESUp_phi_, "metJESUp_phi/D");
  outTree_->Branch("metJESDn_pt" , &metJESDn_pt_ , "metJESDn_pt/D" );
  outTree_->Branch("metJESDn_phi", &metJESDn_phi_, "metJESDn_phi/D");

  jets_pt_   = new doubles;
  jetsUp_pt_ = new doubles;
  jetsDn_pt_ = new doubles;
  outTree_->Branch("jets_pt"  , &jets_pt_   );
  outTree_->Branch("jetsUp_pt", &jetsUp_pt_ );
  outTree_->Branch("jetsDn_pt", &jetsDn_pt_ );

  jets_bTag_   = new doubles;
  jetsUp_bTag_ = new doubles;
  jetsDn_bTag_ = new doubles;
  outTree_->Branch("jets_bTag"  , &jets_bTag_  );
  outTree_->Branch("jetsUp_bTag", &jetsUp_bTag_);
  outTree_->Branch("jetsDn_bTag", &jetsDn_bTag_);

  outTree_->Branch("bjets_n"  , &bjets_n_  , "bjets_n/i"  );
  outTree_->Branch("bjetsUp_n", &bjetsUp_n_, "bjestUp_n/i");
  outTree_->Branch("bjetsDn_n", &bjetsDn_n_, "bjetsDn_n/i");

  outTree_->Branch("lb1_m", &lb1_m_, "lb1_m/D");
  outTree_->Branch("lb2_m", &lb2_m_, "lb2_m/D");
  outTree_->Branch("ttbar_vsumM", &ttbar_vsumM_, "ttbar_vsumM/D");

  if ( isMC_ )
  {
    outTree_->Branch("decayMode", &decayMode_, "decayMode/I");

    outTree_->Branch("puWeight", &puWeight_, "puWeight/D");
    outTree_->Branch("puWeightUp", &puWeightUp_, "puWeightUp/D");
    outTree_->Branch("puWeightDn", &puWeightDn_, "puWeightDn/D");

    outTree_->Branch("metJER_pt" , &metJER_pt_ , "metJER_pt/D" );
    outTree_->Branch("metJER_phi", &metJER_phi_, "metJER_phi/D");
    outTree_->Branch("metJERUp_pt" , &metJERUp_pt_ , "metJERUp_pt/D" );
    outTree_->Branch("metJERUp_phi", &metJERUp_phi_, "metJERUp_phi/D");
    outTree_->Branch("metJERDn_pt" , &metJERDn_pt_ , "metJERDn_pt/D" );
    outTree_->Branch("metJERDn_phi", &metJERDn_phi_, "metJERDn_phi/D");

    jetsResUp_pt_ = new doubles;
    jetsResDn_pt_ = new doubles;
    outTree_->Branch("jetsResUp_pt", &jetsUp_pt_);
    outTree_->Branch("jetsResDn_pt", &jetsUp_pt_);

    jetsResUp_bTag_ = new doubles;
    jetsResDn_bTag_ = new doubles;
    outTree_->Branch("jetsResUp_bTag", &jetsUp_bTag_);
    outTree_->Branch("jetsResDn_bTag", &jetsUp_bTag_);

    outTree_->Branch("bjetsResUp_n", &bjetsResUp_n_, "bjetsResUp_n/i");
    outTree_->Branch("bjetsResDn_n", &bjetsResDn_n_, "bjetsResDn_n/i");
  }
}

bool KDileptonTreeReducer::analyze()
{
  // Initialize tree
  z_Q_ = -999;
  jets_pt_  ->clear(); jets_bTag_  ->clear();
  jetsUp_pt_->clear(); jetsUp_bTag_->clear();
  jetsDn_pt_->clear(); jetsDn_bTag_->clear();
  bjets_n_ = 0;
  bjetsUp_n_ = 0;
  bjetsDn_n_ = 0;
  if ( isMC_ )
  {
    decayMode_ = 0;
    jetsResUp_pt_->clear(); jetsResUp_bTag_->clear();
    jetsResDn_pt_->clear(); jetsResDn_bTag_->clear();
    bjetsResUp_n_ = 0;
    bjetsResDn_n_ = 0;
  }

  // Fill common variables
  nVertex_ = event_->nVertex_;

  // Select leptons
  LorentzVector lepton1P4, lepton2P4, zP4;
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

    for ( int j=((mode_ == 3) ? 0 : i+1), m=leptons2_pt_->size(); j<m; ++j )
    {
      if ( mode_ == 1 )
      {
        if ( event_->muons_type_->at(j) == 0 ) continue; // MuMu mode
      }
      else if ( event_->electrons_type_->at(j) < 100 ) continue; // ElEl and MuEl

      lepton2_pt_  = leptons2_pt_->at(j);
      lepton2_eta_ = leptons2_eta_->at(j);
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
  z_dphi_ = lepton1P4.DeltaPhi(lepton2P4);

  met_pt_  = event_->met_pt_ ;
  met_phi_ = event_->met_phi_;
  metJESUp_pt_  = event_->metJESUp_pt_ ;
  metJESUp_phi_ = event_->metJESUp_phi_;
  metJESDn_pt_  = event_->metJESDn_pt_ ;
  metJESDn_phi_ = event_->metJESDn_phi_;
  if ( isMC_ )
  {
    metJER_pt_  = event_->metJER_pt_ ;
    metJER_phi_ = event_->metJER_phi_;
    metJERUp_pt_  = event_->metJERUp_pt_ ;
    metJERUp_phi_ = event_->metJERUp_phi_;
    metJERDn_pt_  = event_->metJERDn_pt_ ;
    metJERDn_phi_ = event_->metJERDn_phi_;
  }
  TLorentzVector metP4;
  metP4.SetPtEtaPhiM(met_pt_, 0, met_phi_, 0);

  if ( isMC_ )
  {
    // Decay mode
    unsigned int nTop = 0, nMuon = 0, nElectron = 0, nTau = 0;
    for ( int i=0, n=event_->genParticles_pdgId_->size(); i<n; ++i )
    {
      const int pdgId = abs(event_->genParticles_pdgId_->at(i));
      if ( pdgId == 6 ) ++nTop;
      else if ( pdgId == 11 ) ++nElectron;
      else if ( pdgId == 13 ) ++nMuon;
      else if ( pdgId == 15 ) ++nTau;
    }
    if ( nTop == 2 )
    {
      if ( nMuon == 2 ) decayMode_ = 1;
      else if ( nElectron == 2 ) decayMode_ = 2;
      else if ( nMuon == 1 and nElectron == 1 ) decayMode_ = 3;
      else decayMode_ = 4;
    }

    // Weights
    puWeight_ = event_->puWeight_;
    puWeightUp_ = event_->puWeightUp_;
    puWeightDn_ = event_->puWeightDn_;
  }

  std::vector<TLorentzVector> jets;
  for ( int i=0, n=event_->jets_pt_->size(); i<n; ++i )
  {
    const double jetPt = event_->jets_pt_->at(i);
    if ( jetPt < 30 ) continue;
    jets_pt_->push_back(jetPt);

    const double jetBtag = event_->jets_bTag_->at(i);
    //if ( jetBtag > 0.244 ) ++bjets_n_;
    if ( jetBtag > 0.679 ) ++bjets_n_;
    //if ( jetBtag > 0.898 ) ++bjets_n_;
    jets_bTag_->push_back(jetBtag);

    LorentzVector jetP4;
    jetP4.SetPtEtaPhiM(jetPt, event_->jets_eta_->at(i), event_->jets_phi_->at(i), event_->jets_m_->at(i));
    jets.push_back(jetP4);
  }

  ttbar_vsumM_ = -1;
  if ( jets.size() >= 2 )
  {
    const double lb1122_dphi = abs(lepton1P4.DeltaPhi(jets[0]))+abs(lepton2P4.DeltaPhi(jets[1]));
    const double lb1221_dphi = abs(lepton1P4.DeltaPhi(jets[1]))+abs(lepton2P4.DeltaPhi(jets[0]));
    if ( lb1122_dphi < lb1221_dphi )
    {
      lb1_m_ = (lepton1P4+jets[0]).M();
      lb2_m_ = (lepton2P4+jets[1]).M();
    }
    else
    {
      lb1_m_ = (lepton1P4+jets[1]).M();
      lb2_m_ = (lepton2P4+jets[0]).M();
    }
    ttbar_vsumM_ = (zP4+metP4+jets[0]+jets[1]).M();
  }

  return true;
}

