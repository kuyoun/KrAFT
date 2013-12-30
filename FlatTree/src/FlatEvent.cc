#include "KrAFT/FlatTree/interface/FlatEvent.h"

FlatEvent::FlatEvent(bool isMC)
{
  isMC_ = isMC;

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

  if ( isMC_ )
  {
    // JER
    jetsResUp_pt_ = new doubles; jetsResUp_eta_ = new doubles; jetsResUp_phi_ = new doubles; jetsResUp_m_ = new doubles;
    jetsResDn_pt_ = new doubles; jetsResDn_eta_ = new doubles; jetsResDn_phi_ = new doubles; jetsResDn_m_ = new doubles;
    jetsResUp_bTag_ = new doubles;
    jetsResDn_bTag_ = new doubles;

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

    genParticles_pt_  = new doubles;
    genParticles_eta_ = new doubles;
    genParticles_phi_ = new doubles;
    genParticles_m_   = new doubles;
    genParticles_pdgId_ = new ints;
  }
}

void FlatEvent::book(TTree* tree)
{
  tree_ = tree;

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
    tree_->Branch("jetsResDn_pt" , jetsResDn_pt_ );
    tree_->Branch("jetsResDn_eta", jetsResDn_eta_);
    tree_->Branch("jetsResDn_phi", jetsResDn_phi_);
    tree_->Branch("jetsResDn_m"  , jetsResDn_m_  );

    tree_->Branch("jetsResUp_pt" , jetsResUp_pt_ );
    tree_->Branch("jetsResUp_eta", jetsResUp_eta_);
    tree_->Branch("jetsResUp_phi", jetsResUp_phi_);
    tree_->Branch("jetsResUp_m"  , jetsResUp_m_  );

    tree_->Branch("jetsResUp_bTag", jetsResUp_bTag_);
    tree_->Branch("jetsResDn_bTag", jetsResDn_bTag_);

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

    tree_->Branch("genElectrons_pt" , genElectrons_pt_ );
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

    tree_->Branch("genParticles_pt" , genParticles_pt_ );
    tree_->Branch("genParticles_eta", genParticles_eta_);
    tree_->Branch("genParticles_phi", genParticles_phi_);
    tree_->Branch("genParticles_m"  , genParticles_m_  );
    tree_->Branch("genParticles_pdgId", genParticles_pdgId_);
  }
}

void FlatEvent::clear()
{
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

    genParticles_pt_ ->clear();
    genParticles_eta_->clear();
    genParticles_phi_->clear();
    genParticles_m_  ->clear();
    genParticles_pdgId_->clear();
  }
}

void FlatEvent::setBranch(TTree* tree)
{
  tree_ = tree;

  tree_->SetBranchAddress("run", &run_);
  tree_->SetBranchAddress("lumi", &lumi_);
  tree_->SetBranchAddress("event", &event_);

  tree_->SetBranchAddress("puWeight", &puWeight_);
  tree_->SetBranchAddress("puWeightUp", &puWeightUp_);
  tree_->SetBranchAddress("puWeightDn", &puWeightDn_);
  tree_->SetBranchAddress("nVertex", &nVertex_);

  tree_->SetBranchAddress("muons_pt"  , &muons_pt_  );
  tree_->SetBranchAddress("muons_eta" , &muons_eta_ );
  tree_->SetBranchAddress("muons_phi" , &muons_phi_ );
  tree_->SetBranchAddress("muons_m"   , &muons_m_   );
  tree_->SetBranchAddress("muons_Q"   , &muons_Q_   );
  tree_->SetBranchAddress("muons_type", &muons_type_);
  tree_->SetBranchAddress("muons_relIso", &muons_relIso_);

  tree_->SetBranchAddress("electrons_pt"  , &electrons_pt_  );
  tree_->SetBranchAddress("electrons_eta" , &electrons_eta_ );
  tree_->SetBranchAddress("electrons_phi" , &electrons_phi_ );
  tree_->SetBranchAddress("electrons_m"   , &electrons_m_   );
  tree_->SetBranchAddress("electrons_Q"   , &electrons_Q_   );
  tree_->SetBranchAddress("electrons_type", &electrons_type_);
  tree_->SetBranchAddress("electrons_relIso", &electrons_relIso_);

  tree_->SetBranchAddress("electrons_mva", &electrons_mva_);
  tree_->SetBranchAddress("electrons_scEta", &electrons_scEta_);

  tree_->SetBranchAddress("jets_pt" , &jets_pt_ );
  tree_->SetBranchAddress("jets_eta", &jets_eta_);
  tree_->SetBranchAddress("jets_phi", &jets_phi_);
  tree_->SetBranchAddress("jets_m"  , &jets_m_  );

  tree_->SetBranchAddress("jetsUp_pt" , &jetsUp_pt_ );
  tree_->SetBranchAddress("jetsUp_eta", &jetsUp_eta_);
  tree_->SetBranchAddress("jetsUp_phi", &jetsUp_phi_);
  tree_->SetBranchAddress("jetsUp_m"  , &jetsUp_m_  );

  tree_->SetBranchAddress("jetsDn_pt" , &jetsDn_pt_ );
  tree_->SetBranchAddress("jetsDn_eta", &jetsDn_eta_);
  tree_->SetBranchAddress("jetsDn_phi", &jetsDn_phi_);
  tree_->SetBranchAddress("jetsDn_m"  , &jetsDn_m_  );

  tree_->SetBranchAddress("jets_bTag"  , &jets_bTag_  );
  tree_->SetBranchAddress("jetsUp_bTag", &jetsUp_bTag_);
  tree_->SetBranchAddress("jetsDn_bTag", &jetsDn_bTag_);

  tree_->SetBranchAddress("met_pt"  , &met_pt_  );
  tree_->SetBranchAddress("metUp_pt", &metUp_pt_);
  tree_->SetBranchAddress("metDn_pt", &metDn_pt_);

  tree_->SetBranchAddress("met_phi"  , &met_phi_  );
  tree_->SetBranchAddress("metUp_phi", &metUp_phi_);
  tree_->SetBranchAddress("metDn_phi", &metDn_phi_);

  tree_->SetBranchAddress("jpsis_pt" , &jpsis_pt_ );
  tree_->SetBranchAddress("jpsis_eta", &jpsis_eta_);
  tree_->SetBranchAddress("jpsis_phi", &jpsis_phi_);
  tree_->SetBranchAddress("jpsis_m"  , &jpsis_m_  );
  tree_->SetBranchAddress("jpsis_lxy", &jpsis_lxy_);

  if ( isMC_ )
  {
    tree_->SetBranchAddress("jetsResDn_pt" , &jetsResDn_pt_ );
    tree_->SetBranchAddress("jetsResDn_eta", &jetsResDn_eta_);
    tree_->SetBranchAddress("jetsResDn_phi", &jetsResDn_phi_);
    tree_->SetBranchAddress("jetsResDn_m"  , &jetsResDn_m_  );

    tree_->SetBranchAddress("jetsResUp_pt" , &jetsResUp_pt_ );
    tree_->SetBranchAddress("jetsResUp_eta", &jetsResUp_eta_);
    tree_->SetBranchAddress("jetsResUp_phi", &jetsResUp_phi_);
    tree_->SetBranchAddress("jetsResUp_m"  , &jetsResUp_m_  );

    tree_->SetBranchAddress("jetsResUp_bTag", &jetsResUp_bTag_);
    tree_->SetBranchAddress("jetsResDn_bTag", &jetsResDn_bTag_);

    tree_->SetBranchAddress("genWeight", &genWeight_);

    tree_->SetBranchAddress("pdf_id1", &pdf_id1_);
    tree_->SetBranchAddress("pdf_id2", &pdf_id2_);
    tree_->SetBranchAddress("pdf_x1" , &pdf_x1_ );
    tree_->SetBranchAddress("pdf_x2" , &pdf_x2_ );
    tree_->SetBranchAddress("pdf_q"  , &pdf_q_  );

    tree_->SetBranchAddress("genMuons_pt" , &genMuons_pt_ );
    tree_->SetBranchAddress("genMuons_eta", &genMuons_eta_);
    tree_->SetBranchAddress("genMuons_phi", &genMuons_phi_);
    tree_->SetBranchAddress("genMuons_m"  , &genMuons_m_  );
    tree_->SetBranchAddress("genMuons_Q"  , &genMuons_Q_  );

    tree_->SetBranchAddress("genElectrons_pt" , &genElectrons_pt_ );
    tree_->SetBranchAddress("genElectrons_eta", &genElectrons_eta_);
    tree_->SetBranchAddress("genElectrons_phi", &genElectrons_phi_);
    tree_->SetBranchAddress("genElectrons_m"  , &genElectrons_m_  );
    tree_->SetBranchAddress("genElectrons_Q"  , &genElectrons_Q_  );

    tree_->SetBranchAddress("genNeutrinos_pt" , &genNeutrinos_pt_ );
    tree_->SetBranchAddress("genNeutrinos_eta", &genNeutrinos_eta_);
    tree_->SetBranchAddress("genNeutrinos_phi", &genNeutrinos_phi_);

    tree_->SetBranchAddress("genJets_pt" , &genJets_pt_ );
    tree_->SetBranchAddress("genJets_eta", &genJets_eta_);
    tree_->SetBranchAddress("genJets_phi", &genJets_phi_);
    tree_->SetBranchAddress("genJets_m"  , &genJets_m_  );

    tree_->SetBranchAddress("genParticles_pt" , &genParticles_pt_ );
    tree_->SetBranchAddress("genParticles_eta", &genParticles_eta_);
    tree_->SetBranchAddress("genParticles_phi", &genParticles_phi_);
    tree_->SetBranchAddress("genParticles_m"  , &genParticles_m_  );
    tree_->SetBranchAddress("genParticles_pdgId", &genParticles_pdgId_);
  }
}


