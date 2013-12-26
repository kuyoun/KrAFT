#include "KrAFT/FlatTree/interface/KFlatTreeAnalyzer.h"
#include "TString.h"

using namespace std;

void printEntryFraction(int i, int n)
{
  static int digit = 1;
  if ( i%max(1, digit/10) == 0 )
  {
    cout << Form(TString("% ")+Form("%d", int(log10(n))+2)+"d/%d processed\r", i, n);
    cout.flush();
  }
  if ( i%digit == 0 ) digit *= 10;

  if ( i==(n-1) )
  {
    cout << Form(TString("% ")+Form("%d", int(log10(n))+2)+"d/%d processed", i, n) << endl;
    digit = 1;
  }
}

KFlatTreeAnalyzer::KFlatTreeAnalyzer(const string inputFileName, const string outputFileName)
{
  inputFile_ = TFile::Open(inputFileName.c_str());
  outputFile_ = TFile::Open(outputFileName.c_str(), "RECREATE");
}

KFlatTreeAnalyzer::~KFlatTreeAnalyzer()
{
  inputFile_->Close();
  outputFile_->Write();
  outputFile_->Close();
}

