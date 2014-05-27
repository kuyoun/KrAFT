#ifndef KrAFT_GenericNtuple_KJpsiTreeAnalyzer_H
#define KrAFT_GenericNtuple_KJpsiTreeAnalyzer_H

#include "KrAFT/GenericNtuple/interface/KFlatTreeAnalyzerBase.h"

class KJpsiTreeAnalyzer : public KFlatTreeAnalyzerBase
{
public:
  KJpsiTreeAnalyzer(const std::string modeName,
                       const std::string inputFileName,
                       const std::string outputFileName);
  void run() { KFlatTreeAnalyzerBase::run(); }

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
  int nVertex_, nEvent_;
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


  double st_;
  double lb1_m_, lb2_m_;
  double ttbar_vsumM_;

  doublesP jpsis_pt_, jpsis_eta_ , jpsis_phi_, jpsis_m_,jpsis_jetdR_;
  doublesP jpsis_eta1_, jpsis_phi1_, jpsis_pt1_; 
  doublesP jpsis_eta2_, jpsis_phi2_, jpsis_pt2_;
  doublesP jpsis_vProb_;
  doublesP jpsis_l3D_;

  doublesP l1jpsi_m_;
  doublesP l2jpsi_m_;
  doublesP l1jpsi_angle_;
  doublesP l2jpsi_angle_;


};

#endif
