#ifndef KFlatTreeReducer_h
#define KFlatTreeReducer_h

#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"

#include "KrAFT/GenericNtuple/interface/GenericEvent.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

void printEntryFraction(int i, int n);

class KFlatTreeReducer
{
public:
  KFlatTreeReducer(const std::string modeName, const std::string inputFileName, const std::string outputFileName);
  virtual ~KFlatTreeReducer();
  virtual void analyze();

  void setEventScale(const double scale);
  void setCrossSection(const double crossSection);

private:
  TFile* inputFile_;
  TFile* outputFile_;
  GenericEvent* event_;
  TDirectory* outDir_;
  std::vector<TDirectory*> cutStepDirs_;
  TH1* hEvent_;
  TTree* outTree_;

  bool isMC_;
  double eventScale_;
};

#endif
