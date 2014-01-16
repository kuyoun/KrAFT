#ifndef KrAFT_GenericNtuple_GenericEvent_H
#define KrAFT_GenericNtuple_GenericEvent_H

#include <vector>
#include "TTree.h"

struct GenericEvent
{
public:
  GenericEvent(bool isMC=false);
  void clear();
  void book(TTree* tree); // book leaves to fill the tree
  void setBranch(TTree* tree);

public:
  TTree* tree_;

  int run_, lumi_, event_;
  double puWeight_, puWeightUp_, puWeightDn_;
  int nVertex_;

  typedef std::vector<int> ints;
  typedef std::vector<unsigned int> uints;
  typedef std::vector<double> doubles;
  typedef ints* intsP;
  typedef uints* uintsP;
  typedef doubles* doublesP;

  doublesP muons_pt_, muons_eta_, muons_phi_, muons_m_;

  intsP    muons_Q_;
  uintsP   muons_type_;
  doublesP muons_relIso_;

  doublesP electrons_pt_, electrons_eta_, electrons_phi_, electrons_m_;
  intsP    electrons_Q_;
  uintsP   electrons_type_;
  doublesP electrons_relIso_;

  doublesP electrons_mva_;
  doublesP electrons_scEta_;

  doublesP jets_pt_, jets_eta_, jets_phi_, jets_m_;
  doublesP jetsUp_pt_, jetsUp_eta_, jetsUp_phi_, jetsUp_m_;
  doublesP jetsDn_pt_, jetsDn_eta_, jetsDn_phi_, jetsDn_m_;
  doublesP jets_bTag_, jetsUp_bTag_, jetsDn_bTag_;
  double met_pt_, met_phi_;
  double metUp_pt_, metUp_phi_;
  double metDn_pt_, metDn_phi_;

  doublesP jpsis_pt_, jpsis_eta_, jpsis_phi_, jpsis_m_;
  doublesP jpsis_lxy_;

  // JER
  doublesP jetsResUp_pt_, jetsResUp_eta_, jetsResUp_phi_, jetsResUp_m_;
  doublesP jetsResDn_pt_, jetsResDn_eta_, jetsResDn_phi_, jetsResDn_m_;
  doublesP jetsResUp_bTag_, jetsResDn_bTag_;
  double metResUp_pt_, metResUp_phi_;
  double metResDn_pt_, metResDn_phi_;

  // Generator level information
  bool isMC_;

  double genWeight_;
  int pdf_id1_, pdf_id2_;
  double pdf_q_, pdf_x1_, pdf_x2_;

  doublesP genJets_pt_, genJets_eta_, genJets_phi_, genJets_m_;
  doublesP genParticles_pt_, genParticles_eta_, genParticles_phi_, genParticles_m_;
  intsP genParticles_pdgId_;
  intsP genParticles_mother1_, genParticles_mother2_;
  intsP genParticles_daughter1_, genParticles_daughter2_;

};

#endif

