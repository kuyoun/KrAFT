#ifndef KFlatTreeAnalyzer_h
#define KFlatTreeAnalyzer_h

#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"

#include "KrAFT/FlatTree/interface/FlatEvent.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

using namespace ROOT::Math::VectorUtil;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

void printEntryFraction(int i, int n);

class KFlatTreeAnalyzer
{
public:
  KFlatTreeAnalyzer(const std::string modeName, const std::string inputFileName, const std::string outputFileName);
  virtual ~KFlatTreeAnalyzer();
  virtual void analyze();

private:
  TFile* inputFile_;
  TFile* outputFile_;
  FlatEvent* event_;
};

#endif
