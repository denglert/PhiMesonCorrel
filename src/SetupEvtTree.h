#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

class Evts {
public :
   Evts(){};
   ~Evts(){};

   // Declaration of leaf types
//   Int_t           run;
//   Int_t           evt;
//   Int_t           lumi;
//   Float_t         vx;
//   Float_t         vy;
     Float_t         vz;
//   Int_t           hiBin;
//   Float_t         hiHF;
//   Float_t         hiHFplus;
//   Float_t         hiHFminus;
//   Float_t         hiHFplusEta4;
//   Float_t         hiHFminusEta4;
//   Float_t         hiZDC;
//   Float_t         hiZDCplus;
//   Float_t         hiZDCminus;
//   Float_t         hiHFhit;
//   Float_t         hiHFhitPlus;
//   Float_t         hiHFhitMinus;
//   Float_t         hiET;
//   Float_t         hiEE;
//   Float_t         hiEB;
//   Float_t         hiEEplus;
//   Float_t         hiEEminus;
//   Int_t           hiNpix;
//   Int_t           hiNpixelTracks;
     Int_t           hiNtracks;
//   Int_t           hiNtracksPtCut;
//   Int_t           hiNtracksEtaCut;
//   Int_t           hiNtracksEtaPtCut;
//   Int_t           hiNevtPlane;
//   Float_t         hiEvtPlanes[126]; 

   // List of branches
//   TBranch        *b_run;   //!
//   TBranch        *b_evt;   //!
//   TBranch        *b_lumi;   //!
//   TBranch        *b_vx;   //!
//   TBranch        *b_vy;   //!
     TBranch        *b_vz;   //!
//   TBranch        *b_hiBin;   //!
//   TBranch        *b_hiHF;   //!
//   TBranch        *b_hiHFplus;   //!
//   TBranch        *b_hiHFminus;   //!
//   TBranch        *b_hiHFplusEta4;   //!
//   TBranch        *b_hiHFminusEta4;   //!
//   TBranch        *b_hiZDC;   //!
//   TBranch        *b_hiZDCplus;   //!
//   TBranch        *b_hiZDCminus;   //!
//   TBranch        *b_hiHFhit;   //!
//   TBranch        *b_hiHFhitPlus;   //!
//   TBranch        *b_hiHFhitMinus;   //!
//   TBranch        *b_hiET;   //!
//   TBranch        *b_hiEE;   //!
//   TBranch        *b_hiEB;   //!
//   TBranch        *b_hiEEplus;   //!
//   TBranch        *b_hiEEminus;   //!
//   TBranch        *b_hiNpix;   //!
//   TBranch        *b_hiNpixelTracks;   //!
     TBranch        *b_hiNtracks;   //!
//   TBranch        *b_hiNtracksPtCut;   //!
//   TBranch        *b_hiNtracksEtaCut;   //!
//   TBranch        *b_hiNtracksEtaPtCut;   //!
//   TBranch        *b_hiNevtPlane;   //!
//   TBranch        *b_hiEvtPlanes;   //!
};

void setupEvtTree(TFile *f, TTree *t, Evts &tEvts, bool doCheck = 1);
