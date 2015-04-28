#include "SetupSkimTree.h"

bool isEventGood(TTree *SkimAna, int iEv)
{
	SkimAna->GetEntry(iEv);
	
}

void setupSkimTree(TFile *f, TTree *t,Skims &tSkims,bool doCheck = 1)
{
     t = (TTree*)f->Get("ppTrack/trackTree");
   // Set branch addresses and branch pointers
//   t->SetBranchAddress("superFilterPath", &tSkims.superFilterPath, &tSkims.b_superFilterPath);
//   t->SetBranchAddress("reco_extra", &tSkims.reco_extra, &tSkims.b_reco_extra);
//   t->SetBranchAddress("reco_extra_jet", &tSkims.reco_extra_jet, &tSkims.b_reco_extra_jet);
//   t->SetBranchAddress("pat_step", &tSkims.pat_step, &tSkims.b_pat_step);
//   t->SetBranchAddress("ana_step", &tSkims.ana_step, &tSkims.b_ana_step);
//   t->SetBranchAddress("pcollisionEventSelection", &tSkims.pcollisionEventSelection, &tSkims.b_pcollisionEventSelection);
//   t->SetBranchAddress("pHBHENoiseFilter", &tSkims.pHBHENoiseFilter, &tSkims.b_pHBHENoiseFilter);
//   t->SetBranchAddress("phiEcalRecHitSpikeFilter", &tSkims.phiEcalRecHitSpikeFilter, &tSkims.b_phiEcalRecHitSpikeFilter);
     t->SetBranchAddress("pPAcollisionEventSelectionPA", &tSkims.pPAcollisionEventSelectionPA, &tSkims.b_pPAcollisionEventSelectionPA);
//   t->SetBranchAddress("phfPosFilter3", &tSkims.phfPosFilter3, &tSkims.b_phfPosFilter3);
//   t->SetBranchAddress("phfNegFilter3", &tSkims.phfNegFilter3, &tSkims.b_phfNegFilter3);
//   t->SetBranchAddress("phfPosFilter2", &tSkims.phfPosFilter2, &tSkims.b_phfPosFilter2);
//   t->SetBranchAddress("phfNegFilter2", &tSkims.phfNegFilter2, &tSkims.b_phfNegFilter2);
//   t->SetBranchAddress("phfPosFilter1", &tSkims.phfPosFilter1, &tSkims.b_phfPosFilter1);
//   t->SetBranchAddress("phfNegFilter1", &tSkims.phfNegFilter1, &tSkims.b_phfNegFilter1);
//   t->SetBranchAddress("phltPixelClusterShapeFilter", &tSkims.phltPixelClusterShapeFilter, &tSkims.b_phltPixelClusterShapeFilter);
//   t->SetBranchAddress("pprimaryvertexFilter", &tSkims.pprimaryvertexFilter, &tSkims.b_pprimaryvertexFilter);
//   t->SetBranchAddress("pPAprimaryVertexFilter", &tSkims.pPAprimaryVertexFilter, &tSkims.b_pPAprimaryVertexFilter);
//   t->SetBranchAddress("pBeamScrapingFilter", &tSkims.pBeamScrapingFilter, &tSkims.b_pBeamScrapingFilter);
//   t->SetBranchAddress("phfCoincFilter", &tSkims.phfCoincFilter, &tSkims.b_phfCoincFilter);
//   t->SetBranchAddress("ppurityFractionFilter", &tSkims.ppurityFractionFilter, &tSkims.b_ppurityFractionFilter);
     t->SetBranchAddress("pVertexFilterCutG", &tSkims.pVertexFilterCutG, &tSkims.b_pVertexFilterCutG);
//   t->SetBranchAddress("pVertexFilterCutGloose", &tSkims.pVertexFilterCutGloose, &tSkims.b_pVertexFilterCutGloose);
//   t->SetBranchAddress("pVertexFilterCutGtight", &tSkims.pVertexFilterCutGtight, &tSkims.b_pVertexFilterCutGtight);
//   t->SetBranchAddress("pVertexFilterCutE", &tSkims.pVertexFilterCutE, &tSkims.b_pVertexFilterCutE);
//   t->SetBranchAddress("pVertexFilterCutEandG", &tSkims.pVertexFilterCutEandG, &tSkims.b_pVertexFilterCutEandG);
//   t->SetBranchAddress("pVertexFilterCutGplus", &tSkims.pVertexFilterCutGplus, &tSkims.b_pVertexFilterCutGplus);
//
//   t->SetBranchAddress("hltAna", &tSkims.hltAna, &tSkims.b_hltAna);
   if (doCheck) {
   }
}
