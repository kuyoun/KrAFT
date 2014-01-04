#include "KrAFT/GenericNtuple/interface/KFlatTreeReducer.h"
#include "TString.h"
#include "TNamed.h"

using namespace std;

using namespace ROOT::Math::VectorUtil;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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

KFlatTreeReducer::KFlatTreeReducer(const string modeName, const string inputFileName, const string outputFileName)
{
  eventScale_ = 1;

  inputFile_ = TFile::Open(inputFileName.c_str());
  TTree* tree = dynamic_cast<TTree*>(inputFile_->Get(Form("%s/event", modeName.c_str())));
  const string dataType = inputFile_->Get(Form("%s/dataType", modeName.c_str()))->GetTitle();
  isMC_ = (dataType == "MC");
  event_ = new GenericEvent(isMC_);
  event_->setBranch(tree);

  if ( isMC_ )
  {
    hEvent_ = (TH1F*)inputFile_->Get(Form("%s/hEvent", modeName.c_str()));
  }

  outputFile_ = TFile::Open(outputFileName.c_str(), "RECREATE");
  outputFile_->cd();
  TDirectory* outDir = outputFile_->mkdir(modeName.c_str());
  outTree_ = new TTree("ntuple", "ntuple");
}

void KFlatTreeReducer::setEventScale(const double scale)
{
  eventScale_ = scale;
}

void KFlatTreeReducer::setCrossSection(const double crossSection)
{
  if ( not isMC_ ) return;
  const double nEventGen = hEvent_->GetBinContent(1);
  eventScale_ = crossSection/nEventGen;
}

KFlatTreeReducer::~KFlatTreeReducer()
{
  //outputFile_->cd();
  outTree_->Write();
  //outputFile_->Write();
  //outputFile_->Close();
}

void KFlatTreeReducer::analyze()
{
  const int nEvent = event_->tree_->GetEntries();
  for ( int i=0; i<nEvent; ++i )
  {
    printEntryFraction(i, nEvent);
    event_->tree_->GetEntry(i);
  }
}

