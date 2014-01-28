#ifndef KrAFT_GenericNtuple_KDileptonTreeReducer_H
#define KrAFT_GenericNtuple_KDileptonTreeReducer_H

#include "KrAFT/GenericNtuple/interface/KFlatTreeReducerBase.h"

class KDileptonTreeReducer : public KFlatTreeReducerBase
{
public:
  KDileptonTreeReducer(const std::string modeName,
                       const std::string inputFileName,
                       const std::string outputFileName);
  void run() { KFlatTreeReducerBase::run(); }

private:
  virtual bool analyze();

private:
  // Cache for Leptons, depend on decay modes
  typedef doubles* doublesP;
  typedef ints* intsP;
  doublesP leptons1_pt_, leptons1_eta_, leptons1_phi_, leptons1_m_, leptons1_iso_;
  doublesP leptons2_pt_, leptons2_eta_, leptons2_phi_, leptons2_m_, leptons2_iso_;
  intsP leptons1_Q_, leptons2_Q_;
  int mode_;

private:
  // Branch items for output tree
  int decayMode_;
  int nVertex_;
  double puWeight_, puWeightUp_, puWeightDn_;

  double lepton1_pt_, lepton1_eta_, lepton1_phi_, lepton1_iso_;
  double lepton2_pt_, lepton2_eta_, lepton2_phi_, lepton2_iso_;
  double z_m_, z_pt_, z_dphi_;
  int z_Q_;

  unsigned int bjets_n_, bjetsUp_n_, bjetsDn_n_, bjetsResUp_n_, bjetsResDn_n_;
  doublesP jets_pt_, jets_bTag_;
  doublesP jetsUp_pt_, jetsUp_bTag_;
  doublesP jetsDn_pt_, jetsDn_bTag_;
  doublesP jetsResUp_pt_, jetsResUp_bTag_;
  doublesP jetsResDn_pt_, jetsResDn_bTag_;

  double met_pt_, met_phi_;
  double metJESUp_pt_, metJESUp_phi_;
  double metJESDn_pt_, metJESDn_phi_;
  double metJER_pt_, metJER_phi_;
  double metJERUp_pt_, metJERUp_phi_;
  double metJERDn_pt_, metJERDn_phi_;

  double lb1_m_, lb2_m_;
  double ttbar_vsumM_;
};

#endif
