#include "KrAFT/FlatTree/interface/KFlatTreeAnalyzer.h"
#include "TString.h"
#include "TNamed.h"

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

KFlatTreeAnalyzer::KFlatTreeAnalyzer(const string modeName, const string inputFileName, const string outputFileName)
{
  inputFile_ = TFile::Open(inputFileName.c_str());
  TTree* tree = dynamic_cast<TTree*>(inputFile_->Get(Form("%s/event", modeName.c_str())));
  const string dataType = inputFile_->Get(Form("%s/dataType", modeName.c_str()))->GetTitle();
  const bool isMC = (dataType == "MC");
  event_ = new FlatEvent(isMC);
  event_->setBranch(tree);

  outputFile_ = TFile::Open(outputFileName.c_str(), "RECREATE");
}

KFlatTreeAnalyzer::~KFlatTreeAnalyzer()
{
  outputFile_->Write();
  outputFile_->Close();
}

void KFlatTreeAnalyzer::analyze()
{
  const int nEvent = event_->tree_->GetEntries();
  for ( int i=0; i<nEvent; ++i )
  {
    printEntryFraction(i, nEvent);
    event_->tree_->GetEntry(i);

    cout << event_->event_;
  }
}

