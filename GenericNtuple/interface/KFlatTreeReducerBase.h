#ifndef KFlatTreeReducerBase_h
#define KFlatTreeReducerBase_h

#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"

#include "KrAFT/GenericNtuple/interface/GenericEvent.h"

//#include "Math/GenVector/LorentzVector.h"
//#include "Math/GenVector/VectorUtil.h"
//using namespace ROOT::Math::VectorUtil;
//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

#include <TLorentzVector.h>
typedef TLorentzVector LorentzVector;
typedef std::vector<LorentzVector> LorentzVectors;
typedef std::vector<double> doubles;
typedef std::vector<int> ints;

void printEntryFraction(int i, int n);

class KFlatTreeReducerBase
{
public:
  KFlatTreeReducerBase(const std::string modeName,
                       const std::string inputFileName, 
                       const std::string outputFileName);
  virtual ~KFlatTreeReducerBase();

  void run();

protected:
  virtual void init() = 0;
  virtual bool analyze() = 0;
  void getP4(const doubles* pt, const doubles* eta, const doubles* phi, const doubles* m, LorentzVectors& p4);
  std::string modeName_;
  TFile* inputFile_;
  TFile* outputFile_;
  GenericEvent* event_;
  TDirectory* outDir_;
  TH1* hEvent_;
  TTree* outTree_;

  bool isMC_;
};

#endif
