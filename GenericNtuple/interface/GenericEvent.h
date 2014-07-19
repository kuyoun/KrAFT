#ifndef KrAFT_GenericNtuple_GenericEvent_H
#define KrAFT_GenericNtuple_GenericEvent_H

#include <string>
#include <vector>
#include <map>
#include "TTree.h"

struct GenericEvent
{
public:
  GenericEvent(bool isMC=false);
  ~GenericEvent();
  void clear();
  void book(TTree* tree); // book leaves to fill the tree
  void setBranch(TTree* tree);

  void append(const std::string varName, const double value)
  {
    fVars_[varName]->push_back(value);
  };
  void append(const std::string varName, const int value)
  {
    iVars_[varName]->push_back(value);
  };

public:
  TTree* tree_;

  typedef std::vector<int> ints;
  typedef std::vector<unsigned int> uints;
  typedef std::vector<double> doubles;
  typedef std::vector<std::string> strings;
  typedef std::map<std::string, doubles*> FVars;
  typedef std::map<std::string, ints*> IVars;

  strings iVarNames_, fVarNames_;
  IVars iVars_;
  FVars fVars_;

  int run_, lumi_, event_;
  double puWeight_, puWeightUp_, puWeightDn_;
  int nVertex_, nPileup_;

  double met_pt_, met_phi_;
  double metJESUp_pt_, metJESUp_phi_;
  double metJESDn_pt_, metJESDn_phi_;
  double metJER_pt_  , metJER_phi_  ;
  double metJERUp_pt_, metJERUp_phi_;
  double metJERDn_pt_, metJERDn_phi_;

  // Generator level information
  bool isMC_;

  double genWeight_;
  int pdf_id1_, pdf_id2_;
  double pdf_q_, pdf_x1_, pdf_x2_;
};

#endif

