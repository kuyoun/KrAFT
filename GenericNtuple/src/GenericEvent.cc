#include "KrAFT/GenericNtuple/interface/GenericEvent.h"
#include <boost/assign/std/vector.hpp>

GenericEvent::GenericEvent(bool isMC)
{
  isMC_ = isMC;

  // Initialize list of variable names
  using namespace boost::assign;
  fVarNames_ += "muons_pt", "muons_eta", "muons_phi", "muons_m";
  fVarNames_ += "muons_relIso";
  iVarNames_ += "muons_Q", "muons_type";

  fVarNames_ += "electrons_pt", "electrons_eta", "electrons_phi", "electrons_m";
  fVarNames_ += "electrons_relIso";
  fVarNames_ += "electrons_mva", "electrons_scEta";
  fVarNames_ += "electrons_Q", "electrons_type", "electrons_qConsistent";

  fVarNames_ += "jets_pt", "jets_eta", "jets_phi", "jets_m";
  fVarNames_ += "jets_bTag", "jets_JESUp", "jets_JESDn";
  if ( isMC_ )
  {
    fVarNames_ += "jets_JER", "jets_JERUp", "jets_JERDn";
    iVarNames_ += "jets_partonflavor";
  }

  if ( isMC_ )
  {
    fVarNames_ += "pdfWeights";

    // GenJets
    fVarNames_ += "genJets_pt", "genJets_eta", "genJets_phi", "genJets_m";
    iVarNames_ += "genJets_pdgId";

    // Generator information
    fVarNames_ += "genParticles_pt", "genParticles_eta", "genParticles_phi", "genParitcles_m";
    iVarNames_ += "genParticles_pdgId", "genParticles_mother";
  }

  fVarNames_ += "jpsis_pt", "jpsis_eta", "jpsis_phi", "jpsis_m";
  fVarNames_ += "jpsis_lxy", "jpsis_l3D";
  fVarNames_ += "jpsis_pt1", "jpsis_eta1", "jpsis_phi1";
  fVarNames_ += "jpsis_pt2", "jpsis_eta2", "jpsis_phi2";
  iVarNames_ += "jpsis_nPixHits1", "jpsis_nPixHits2";
}

void GenericEvent::book(TTree* tree)
{
  tree_ = tree;

  for ( strings::const_iterator itr = iVarNames_.begin(); itr != iVarNames_.end(); ++itr )
  {
    const char* varName = itr->c_str();
    ints* v = iVars_[varName] = new ints;
    tree_->Branch(varName, v);
  }

  for ( strings::const_iterator itr = fVarNames_.begin(); itr != iVarNames_.end(); ++itr )
  {
    const char* varName = itr->c_str();
    doubles* v = fVars_[varName] = new doubles;
    tree_->Branch(varName, v);
  }

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/I");

  tree_->Branch("puWeight", &puWeight_, "puWeight/D");
  tree_->Branch("puWeightUp", &puWeightUp_, "puWeightUp/D");
  tree_->Branch("puWeightDn", &puWeightDn_, "puWeightDn/D");
  tree_->Branch("nVertex", &nVertex_, "nVertex/I");
  tree_->Branch("nPileup", &nPileup_, "nPileup/I");

  tree_->Branch("met_pt"     , &met_pt_     , "met_pt/D"     );
  tree_->Branch("metJESUp_pt", &metJESUp_pt_, "metJESUp_pt/D");
  tree_->Branch("metJESDn_pt", &metJESDn_pt_, "metJESDn_pt/D");

  tree_->Branch("met_phi"     , &met_phi_     , "met_phi/D"     );
  tree_->Branch("metJESUp_phi", &metJESUp_phi_, "metJESUp_phi/D");
  tree_->Branch("metJESDn_phi", &metJESDn_phi_, "metJESDn_phi/D");

  if ( isMC_ )
  {
    tree_->Branch("metJER_pt"  , &metJER_pt_  , "metJER_pt/D"  );
    tree_->Branch("metJERUp_pt", &metJERUp_pt_, "metJERUp_pt/D");
    tree_->Branch("metJERDn_pt", &metJERDn_pt_, "metJERDn_pt/D");

    tree_->Branch("metJER_phi"  , &metJER_phi_  , "metJER_phi/D"  );
    tree_->Branch("metJERUp_phi", &metJERUp_phi_, "metJERUp_phi/D");
    tree_->Branch("metJERDn_phi", &metJERDn_phi_, "metJERDn_phi/D");

    tree_->Branch("genWeight", &genWeight_, "genWeight/D");

    tree_->Branch("pdf_id1", &pdf_id1_, "pdf_id1/I");
    tree_->Branch("pdf_id2", &pdf_id2_, "pdf_id2/I");
    tree_->Branch("pdf_x1" , &pdf_x1_ , "pdf_x1/D" );
    tree_->Branch("pdf_x2" , &pdf_x2_ , "pdf_x2/D" );
    tree_->Branch("pdf_q"  , &pdf_q_  , "pdf_q/D"  );
  }
}

void GenericEvent::clear()
{
  // Clear up
  for ( FVars::const_iterator itr = fVars_.begin(); itr != fVars_.end(); ++itr )
  {
    itr->second->clear();
  }

  for ( IVars::const_iterator itr = iVars_.begin(); itr != iVars_.end(); ++itr )
  {
    itr->second->clear();
  }
}

void GenericEvent::setBranch(TTree* tree)
{
  tree_ = tree;

  tree_->SetBranchAddress("run", &run_);
  tree_->SetBranchAddress("lumi", &lumi_);
  tree_->SetBranchAddress("event", &event_);

  tree_->SetBranchAddress("puWeight", &puWeight_);
  tree_->SetBranchAddress("puWeightUp", &puWeightUp_);
  tree_->SetBranchAddress("puWeightDn", &puWeightDn_);
  tree_->SetBranchAddress("nVertex", &nVertex_);
  tree_->SetBranchAddress("nPileup", &nPileup_);

  for ( strings::const_iterator itr = iVarNames_.begin(); itr != iVarNames_.end(); ++itr )
  {
    const char* varName = itr->c_str();
    ints*& v = iVars_[varName] = new ints;
    tree_->SetBranchAddress(varName, &v);
  }

  for ( strings::const_iterator itr = fVarNames_.begin(); itr != fVarNames_.end(); ++itr )
  {
    const char* varName = itr->c_str();
    doubles*& v = fVars_[varName] = new doubles;
    tree_->SetBranchAddress(varName, &v);
  }

  tree_->SetBranchAddress("met_pt"     , &met_pt_     );
  tree_->SetBranchAddress("metJESUp_pt", &metJESUp_pt_);
  tree_->SetBranchAddress("metJESDn_pt", &metJESDn_pt_);

  tree_->SetBranchAddress("met_phi"     , &met_phi_     );
  tree_->SetBranchAddress("metJESUp_phi", &metJESUp_phi_);
  tree_->SetBranchAddress("metJESDn_phi", &metJESDn_phi_);

  if ( isMC_ )
  {
    tree_->SetBranchAddress("metJER_pt"  , &metJER_pt_  );
    tree_->SetBranchAddress("metJERUp_pt", &metJERUp_pt_);
    tree_->SetBranchAddress("metJERDn_pt", &metJERDn_pt_);

    tree_->SetBranchAddress("metJER_phi"  , &metJER_phi_  );
    tree_->SetBranchAddress("metJERUp_phi", &metJERUp_phi_);
    tree_->SetBranchAddress("metJERDn_phi", &metJERDn_phi_);

    tree_->SetBranchAddress("genWeight", &genWeight_);

    tree_->SetBranchAddress("pdf_id1", &pdf_id1_);
    tree_->SetBranchAddress("pdf_id2", &pdf_id2_);
    tree_->SetBranchAddress("pdf_x1" , &pdf_x1_ );
    tree_->SetBranchAddress("pdf_x2" , &pdf_x2_ );
    tree_->SetBranchAddress("pdf_q"  , &pdf_q_  );
  }
}

GenericEvent::~GenericEvent()
{
  for ( FVars::const_iterator itr = fVars_.begin();
        itr != fVars_.end(); ++itr )
  {
    delete itr->second;
  }

  for ( IVars::const_iterator itr = iVars_.begin();
        itr != iVars_.end(); ++itr )
  {
    delete itr->second;
  }
}
