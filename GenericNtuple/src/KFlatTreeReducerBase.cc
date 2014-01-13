#include "KrAFT/GenericNtuple/interface/KFlatTreeReducerBase.h"
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

KFlatTreeReducerBase::KFlatTreeReducerBase(const string modeName, const string inputFileName, const string outputFileName)
{
  event_ = 0;
  outputFile_ = 0;

  modeName_ = modeName;
  inputFile_ = TFile::Open(inputFileName.c_str());
  TTree* tree = dynamic_cast<TTree*>(inputFile_->Get(Form("%s/event", modeName.c_str())));
  if ( !tree ) return;
  const string dataType = inputFile_->Get(Form("%s/dataType", modeName.c_str()))->GetTitle();
  isMC_ = (dataType == "MC");
  event_ = new GenericEvent(isMC_);
  event_->setBranch(tree);

  if ( isMC_ )
  {
    hEvent_ = (TH1F*)inputFile_->Get(Form("%s/hEvent", modeName.c_str()));
  }

  outputFile_ = TFile::Open(outputFileName.c_str(), "UPDATE");
  outDir_ = outputFile_->mkdir(modeName.c_str());
  outDir_->cd();
  outTree_ = new TTree("ntuple", "ntuple");
}

KFlatTreeReducerBase::~KFlatTreeReducerBase()
{
  if ( outputFile_ and outputFile_->IsOpen() ) outputFile_->Close();
}

void KFlatTreeReducerBase::run()
{
  if ( !event_ ) return;

  const int nEvent = event_->tree_->GetEntries();
  for ( int i=0; i<nEvent; ++i )
  {
    printEntryFraction(i, nEvent);
    event_->tree_->GetEntry(i);
    const bool isAccepted = analyze();
    if ( !isAccepted ) continue;

    outTree_->Fill();
  }
  outDir_->cd();
  outTree_->Write();
}

void KFlatTreeReducerBase::getP4(const doubles* pts, const doubles* etas, const doubles* phis, const doubles* ms, LorentzVectors& p4s)
{
  p4s.clear();
  LorentzVector p4;
  for ( int i=0, n=pts->size(); i<n; ++i )
  {
    p4.SetPtEtaPhiM(pts->at(i), etas->at(i), phis->at(i), ms->at(i));
    p4s.push_back(p4);
  }
}
