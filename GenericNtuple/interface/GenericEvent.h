#ifndef KrAFT_GenericNtuple_GenericEvent_H
#define KrAFT_GenericNtuple_GenericEvent_H

#include <vector>
#include "TTree.h"

struct GenericEvent
{
public:
  GenericEvent(bool isMC=false);
  ~GenericEvent();
  void clear();
  void book(TTree* tree); // book leaves to fill the tree
  void setBranch(TTree* tree);

public:
  TTree* tree_;

  int run_, lumi_, event_;
  double puWeight_, puWeightUp_, puWeightDn_;
  int nVertex_, nPileup_;

  typedef std::vector<int> ints;
  typedef std::vector<unsigned int> uints;
  typedef std::vector<double> doubles;
  typedef ints* intsP;
  typedef uints* uintsP;
  typedef doubles* doublesP;

  doublesP muons_pt_, muons_eta_, muons_phi_, muons_m_;

  intsP    muons_Q_;
  uintsP   muons_type_;
  doublesP muons_relIso_;

  doublesP electrons_pt_, electrons_eta_, electrons_phi_, electrons_m_;
  intsP    electrons_Q_;
  uintsP   electrons_type_;
  doublesP electrons_relIso_;

  doublesP electrons_mva_;
  doublesP electrons_scEta_;
  uintsP   electrons_qConsistent_;

  doublesP jets_pt_, jets_eta_, jets_phi_, jets_m_;
  doublesP jets_bTag_, jets_partonflavor_;
  doublesP jets_JESUp_, jets_JESDn_;
  doublesP jets_JER_, jets_JERUp_, jets_JERDn_;

  double met_pt_, met_phi_;
  double metJESUp_pt_, metJESUp_phi_;
  double metJESDn_pt_, metJESDn_phi_;
  double metJER_pt_  , metJER_phi_  ;
  double metJERUp_pt_, metJERUp_phi_;
  double metJERDn_pt_, metJERDn_phi_;

  doublesP jpsis_pt_, jpsis_eta_, jpsis_phi_, jpsis_m_;
  doublesP jpsis_cos_;
  doublesP jpsis_lxy_,jpsis_l3D_, jpsis_vProb_;
  doublesP jpsis_pt1_, jpsis_eta1_, jpsis_phi1_;
  doublesP jpsis_pt2_, jpsis_eta2_, jpsis_phi2_;
  intsP jpsis_nPixHits1_, jpsis_nPixHits2_;

  // Generator level information
  bool isMC_;

  double genWeight_;
  int pdf_id1_, pdf_id2_;
  double pdf_q_, pdf_x1_, pdf_x2_;
  doublesP pdfWeights_;

  doublesP genJets_pt_, genJets_eta_, genJets_phi_, genJets_m_;
  doublesP genParticles_pt_, genParticles_eta_, genParticles_phi_, genParticles_m_;
  intsP genParticles_pdgId_;
  intsP genParticles_mother_;

};

#endif

