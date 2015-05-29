#include <TTree.h>
#include <TBranch.h>

class Skims {
public :
   Skims(){};
   ~Skims(){};

   // Declaration of leaf types
//   Int_t           superFilterPath;
//   Int_t           reco_extra;
//   Int_t           reco_extra_jet;
//   Int_t           pat_step;
//   Int_t           ana_step;
//   Int_t           pcollisionEventSelection;
//   Int_t           pHBHENoiseFilter;
//   Int_t           phiEcalRecHitSpikeFilter;
     Int_t           pPAcollisionEventSelectionPA;
//   Int_t           phfPosFilter3;
//   Int_t           phfNegFilter3;
//   Int_t           phfPosFilter2;
//   Int_t           phfNegFilter2;
//   Int_t           phfPosFilter1;
//   Int_t           phfNegFilter1;
//   Int_t           phltPixelClusterShapeFilter;
//   Int_t           pprimaryvertexFilter;
//   Int_t           pPAprimaryVertexFilter;
//   Int_t           phfCoincFilter;
//   Int_t           ppurityFractionFilter;
//   Int_t           pBeamScrapingFilter;
     Int_t           pVertexFilterCutG;
//   Int_t           pVertexFilterCutGloose;
//   Int_t           pVertexFilterCutGtight;
//   Int_t           pVertexFilterCutE;
//   Int_t           pVertexFilterCutEandG;
//   Int_t           pVertexFilterCutGplus;
//   Int_t           hltAna;

   // List of branches
//   TBranch        *b_superFilterPath;   //!
//   TBranch        *b_reco_extra;   //!
//   TBranch        *b_reco_extra_jet;   //!
//   TBranch        *b_pat_step;   //!
//   TBranch        *b_ana_step;   //!
//   TBranch        *b_pcollisionEventSelection;   //!
//   TBranch        *b_pHBHENoiseFilter;   //!
//   TBranch        *b_phiEcalRecHitSpikeFilter;   //!
     TBranch        *b_pPAcollisionEventSelectionPA;   //!
//   TBranch        *b_phfPosFilter3;   //!
//   TBranch        *b_phfNegFilter3;   //!
//   TBranch        *b_phfPosFilter2;   //!
//   TBranch        *b_phfNegFilter2;   //!
//   TBranch        *b_phfPosFilter1;   //!
//   TBranch        *b_phfNegFilter1;   //!
//   TBranch        *b_phltPixelClusterShapeFilter;   //!
//   TBranch        *b_pprimaryvertexFilter;   //!
//   TBranch        *b_pPAprimaryVertexFilter;   //!
//   TBranch        *b_phfCoincFilter;   //!
//   TBranch        *b_ppurityFractionFilter;   //!
//   TBranch        *b_pBeamScrapingFilter;   //!
     TBranch        *b_pVertexFilterCutG;   //!
//   TBranch        *b_pVertexFilterCutGloose;   //!
//   TBranch        *b_pVertexFilterCutGtight;   //!
//   TBranch        *b_pVertexFilterCutE;   //!
//   TBranch        *b_pVertexFilterCutEandG;   //!
//   TBranch        *b_pVertexFilterCutGplus;
//   TBranch        *b_hltAna;   //!

};

void setupSkimTree(TFile *f, TTree *t,Skims &tSkims,bool doCheck = 1);
