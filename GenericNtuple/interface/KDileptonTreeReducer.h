#ifndef KrAFT_GenericNtuple_KDileptonTreeReducer_H
#define KrAFT_GenericNtuple_KDileptonTreeReducer_H

#include "KrAFT/GenericNtuple/interface/KFlatTreeReducerBase.h"

class KDileptonTreeReducer : public KFlatTreeReducerBase
{
public:
  KDileptonTreeReducer(const std::string modeName,
                       const std::string inputFileName,
                       const std::string outputFileName)
    :KFlatTreeReducerBase(modeName, inputFileName, outputFileName) {}
  void run() { KFlatTreeReducerBase::run(); }

private:
  virtual void init();
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
  double lepton1_pt_, lepton1_eta_, lepton1_phi_, lepton1_iso_;
  double lepton2_pt_, lepton2_eta_, lepton2_phi_, lepton2_iso_;
  double z_m_, z_pt_;
  int z_Q_;
};

#endif
