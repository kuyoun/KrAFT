#include "KrAFT/GenericNtuple/interface/KJpsiTreeAnalyzer.h"
#include<boost/iterator/counting_iterator.hpp>
#include<boost/range/algorithm/sort.hpp>

using namespace std;

struct GreaterByBtag
{
  GreaterByBtag(std::vector<double>* x):x_(x) {};
  bool operator()(size_t i, size_t j) { return x_->at(i) > x_->at(j); }
  std::vector<double>* x_;
};

KJpsiTreeAnalyzer::KJpsiTreeAnalyzer(const std::string modeName,
                                           const std::string inputFileName,
                                           const std::string outputFileName):
  KFlatTreeAnalyzerBase(modeName, inputFileName, outputFileName)
{
  if ( !event_ ) return;

  outTree_->Branch("nVertex", &nVertex_, "nVertex/I");
  outTree_->Branch("nEvent", &nEvent_, "nEvent/I");

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

  jpsis_eta1_ = new doubles;
  jpsis_eta2_ = new doubles;
  jpsis_pt1_ = new doubles;
  jpsis_pt2_ = new doubles;


  outTree_->Branch("bjets_n"  , &bjets_n_  , "bjets_n/i"  );
  outTree_->Branch("bjetsUp_n", &bjetsUp_n_, "bjestUp_n/i");
  outTree_->Branch("bjetsDn_n", &bjetsDn_n_, "bjetsDn_n/i");

  outTree_->Branch("st", &st_, "st/D");
  outTree_->Branch("lb1_m", &lb1_m_, "lb1_m/D");
  outTree_->Branch("lb2_m", &lb2_m_, "lb2_m/D");
  outTree_->Branch("ttbar_vsumM", &ttbar_vsumM_, "ttbar_vsumM/D");

  jpsis_pt_ = new doubles;
  jpsis_eta_ = new doubles;
  jpsis_phi_ = new doubles;
  jpsis_m_ = new doubles;
  jpsis_l3D_ = new doubles;
  jpsis_vProb_ = new doubles;
  jpsis_jetdR_ = new doubles;
  
  outTree_->Branch("jpsis_pt",&jpsis_pt_);
  outTree_->Branch("jpsis_eta",&jpsis_eta_);
  outTree_->Branch("jpsis_phi",&jpsis_phi_);
  outTree_->Branch("jpsis_m",&jpsis_m_);
  outTree_->Branch("jpsis_l3D",&jpsis_l3D_);
  outTree_->Branch("jpsis_vProb",&jpsis_vProb_);
  outTree_->Branch("jpsis_jetdR",&jpsis_jetdR_);

  l1jpsi_m_ = new doubles;
  l2jpsi_m_ = new doubles;
  l1jpsi_angle_ = new doubles;
  l2jpsi_angle_ = new doubles;

  outTree_->Branch("l1jpsi_m",&l1jpsi_m_);
  outTree_->Branch("l2jpsi_m",&l2jpsi_m_);
  outTree_->Branch("l1jpsi_angle",&l1jpsi_angle_);
  outTree_->Branch("l2jpsi_angle",&l2jpsi_angle_);
  
  jpsis_eta1_ = new doubles; 
  jpsis_eta2_ = new doubles; 
  jpsis_phi1_ = new doubles; 
  jpsis_phi2_ = new doubles; 
  jpsis_pt1_ = new doubles; 
  jpsis_pt2_ = new doubles; 

  outTree_->Branch("jpsis_eta1",&jpsis_eta1_);
  outTree_->Branch("jpsis_eta2",&jpsis_eta2_);
  outTree_->Branch("jpsis_pt1",&jpsis_pt1_);
  outTree_->Branch("jpsis_pt2",&jpsis_pt2_);

  

 
  if ( isMC_ )
  {
    outTree_->Branch("decayMode", &decayMode_, "decayMode/I");

    outTree_->Branch("puWeight", &puWeight_, "puWeight/D");
    outTree_->Branch("puWeightUp", &puWeightUp_, "puWeightUp/D");
    outTree_->Branch("puWeightDn", &puWeightDn_, "puWeightDn/D");

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

bool KJpsiTreeAnalyzer::analyze()
{
  // Initialize tree
  z_Q_ = -999;
  jets_pt_  ->clear(); jets_bTag_  ->clear();
  jetsUp_pt_->clear(); jetsUp_bTag_->clear();
  jetsDn_pt_->clear(); jetsDn_bTag_->clear();
  bjets_n_ = 0;
  bjetsUp_n_ = 0;
  bjetsDn_n_ = 0;
  st_ = 0;
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
  nEvent_ = event_->event_;

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

  // Run lepton loop again for st calculation
  for ( int i=0, n=event_->muons_pt_->size(); i<n; ++i )
  {
    const double muonPt = event_->muons_pt_->at(i);
    const double muonEta = abs(event_->muons_eta_->at(i));
    if ( muonPt < 20 or muonEta > 2.5 ) continue;
    st_ += muonPt;
  }
  for ( int i=0, n=event_->electrons_pt_->size(); i<n; ++i )
  {
    const double electronPt = event_->electrons_pt_->at(i);
    const double electronEta = abs(event_->electrons_eta_->at(i));
    if ( electronPt < 20 or electronEta > 2.5 ) continue;
    st_ += event_->electrons_pt_->at(i);
  }

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

  // Get jet indices by bTag
  const int nJets = event_->jets_pt_->size();
  std::vector<int> jetIndices(nJets);
  std::copy(boost::counting_iterator<int>(0),
            boost::counting_iterator<int>(nJets), jetIndices.begin());
  GreaterByBtag greaterByBtag(event_->jets_bTag_);
  boost::sort(jetIndices, greaterByBtag);

  // Make jets four vector, insert in bTag-order
  std::vector<TLorentzVector> jets;
  for ( int i=0; i<nJets; ++i )
  {
    const int j = jetIndices[i];

    const double pt   = event_->jets_pt_->at(j);
    const double eta  = event_->jets_eta_->at(j);
    const double phi  = event_->jets_phi_->at(j);
    const double m    = event_->jets_m_->at(j);
    const double bTag = event_->jets_bTag_->at(j);

    jets_pt_->push_back(pt);
    jets_bTag_->push_back(bTag);

    st_ += pt;

    jets.push_back(LorentzVector());
    jets.back().SetPtEtaPhiM(pt, eta, phi, m);

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
  }
  l1jpsi_m_->clear();
  l2jpsi_m_->clear(); 
  l1jpsi_angle_->clear();
  l2jpsi_angle_->clear(); 
  jpsis_pt_->clear();
  jpsis_eta_->clear();
  jpsis_phi_->clear();
  jpsis_m_->clear();
  jpsis_l3D_->clear();
  jpsis_vProb_->clear();
  jpsis_jetdR_->clear();
  //st_ += pt;
  jpsis_eta1_->clear();
  jpsis_eta2_->clear();
  jpsis_pt1_->clear();
  jpsis_pt2_->clear();



  const int nJpsis = event_->jpsis_pt_->size();
  for ( int i=0; i<nJpsis; ++i )
  {
    const double pt   = event_->jpsis_pt_->at(i);
    const double eta  = event_->jpsis_eta_->at(i);
    const double phi  = event_->jpsis_phi_->at(i);
    const double m    = event_->jpsis_m_->at(i);

    const double eta1 = event_->jpsis_eta1_->at(i);
    const double eta2 = event_->jpsis_eta2_->at(i);
    const double pt1 = event_->jpsis_pt1_->at(i);
    const double pt2 = event_->jpsis_pt2_->at(i);

    const double l3D    = event_->jpsis_l3D_->at(i);
    const double vProb    = event_->jpsis_vProb_->at(i);
    TLorentzVector jpsi;
    jpsi.SetPtEtaPhiM( pt,eta,phi,m);
		double lep1_plus_jpsi_mass = (lepton1P4+jpsi).M();
		double lep1_plus_jpsi_angle = lepton1P4.Angle( jpsi.Vect());
		double lep2_plus_jpsi_mass = (lepton2P4+jpsi).M();
		double lep2_plus_jpsi_angle = lepton2P4.Angle( jpsi.Vect());
		l1jpsi_m_->push_back( lep1_plus_jpsi_mass) ;
		l2jpsi_m_->push_back( lep2_plus_jpsi_mass) ;
		l1jpsi_angle_->push_back( lep1_plus_jpsi_angle) ;
		l2jpsi_angle_->push_back( lep2_plus_jpsi_angle) ;

    double minDR = 999.;
    for ( unsigned int k=0 ; k< jets.size() ; k++) {
			double jet_dR = jpsi.DeltaR( jets[k] );
			if ( jet_dR < minDR ) minDR = jet_dR ; 
    }

    jpsis_pt_->push_back(pt);
    jpsis_eta_->push_back(eta);
    jpsis_phi_->push_back(phi);
    jpsis_m_->push_back(m);
		jpsis_l3D_->push_back(l3D);
		jpsis_vProb_->push_back(vProb);
		jpsis_jetdR_->push_back(minDR);
    //st_ += pt;
		jpsis_eta1_->push_back(eta1);
		jpsis_eta2_->push_back(eta2);
		jpsis_pt1_->push_back(pt1);
		jpsis_pt2_->push_back(pt2);
		

    //jets.push_back(LorentzVector());
    //jets.back().SetPtEtaPhiM(pt, eta, phi, m);

  }

  return true;
}

