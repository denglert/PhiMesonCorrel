#include "SetupEvtTree.h"

void setupEvtTree(TFile *f, TTree *t, Evts &tEvts, bool doCheck = 1)
{

     t = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
   // Set branch addresses and branch pointers
//   t->SetBranchAddress("run", &tEvts.run, &tEvts.b_run);
//   t->SetBranchAddress("evt", &tEvts.evt, &tEvts.b_evt);
//   t->SetBranchAddress("lumi", &tEvts.lumi, &tEvts.b_lumi);
//   t->SetBranchAddress("vx", &tEvts.vx, &tEvts.b_vx);
//   t->SetBranchAddress("vy", &tEvts.vy, &tEvts.b_vy);
     t->SetBranchAddress("vz", &tEvts.vz, &tEvts.b_vz);
//   t->SetBranchAddress("hiBin", &tEvts.hiBin, &tEvts.b_hiBin);
//   t->SetBranchAddress("hiHF", &tEvts.hiHF, &tEvts.b_hiHF);
//   t->SetBranchAddress("hiHFplus", &tEvts.hiHFplus, &tEvts.b_hiHFplus);
//   t->SetBranchAddress("hiHFminus", &tEvts.hiHFminus, &tEvts.b_hiHFminus);
//   t->SetBranchAddress("hiHFplusEta4", &tEvts.hiHFplusEta4, &tEvts.b_hiHFplusEta4);
//   t->SetBranchAddress("hiHFminusEta4", &tEvts.hiHFminusEta4, &tEvts.b_hiHFminusEta4);
//   t->SetBranchAddress("hiZDC", &tEvts.hiZDC, &tEvts.b_hiZDC);
//   t->SetBranchAddress("hiZDCplus", &tEvts.hiZDCplus, &tEvts.b_hiZDCplus);
//   t->SetBranchAddress("hiZDCminus", &tEvts.hiZDCminus, &tEvts.b_hiZDCminus);
//   t->SetBranchAddress("hiHFhit", &tEvts.hiHFhit, &tEvts.b_hiHFhit);
//   t->SetBranchAddress("hiHFhitPlus", &tEvts.hiHFhitPlus, &tEvts.b_hiHFhitPlus);
//   t->SetBranchAddress("hiHFhitMinus", &tEvts.hiHFhitMinus, &tEvts.b_hiHFhitMinus);
//   t->SetBranchAddress("hiET", &tEvts.hiET, &tEvts.b_hiET);
//   t->SetBranchAddress("hiEE", &tEvts.hiEE, &tEvts.b_hiEE);
//   t->SetBranchAddress("hiEB", &tEvts.hiEB, &tEvts.b_hiEB);
//   t->SetBranchAddress("hiEEplus", &tEvts.hiEEplus, &tEvts.b_hiEEplus);
//   t->SetBranchAddress("hiEEminus", &tEvts.hiEEminus, &tEvts.b_hiEEminus);
//   t->SetBranchAddress("hiNpix", &tEvts.hiNpix, &tEvts.b_hiNpix);
//   t->SetBranchAddress("hiNpixelTracks", &tEvts.hiNpixelTracks, &tEvts.b_hiNpixelTracks);
     t->SetBranchAddress("hiNtracks", &tEvts.hiNtracks, &tEvts.b_hiNtracks);
//   t->SetBranchAddress("hiNtracksPtCut", &tEvts.hiNtracksPtCut, &tEvts.b_hiNtracksPtCut);
//   t->SetBranchAddress("hiNtracksEtaCut", &tEvts.hiNtracksEtaCut, &tEvts.b_hiNtracksEtaCut);
//   t->SetBranchAddress("hiNtracksEtaPtCut", &tEvts.hiNtracksEtaPtCut, &tEvts.b_hiNtracksEtaPtCut);
//   if (t->GetBranch("hiNevtPlane")) t->SetBranchAddress("hiNevtPlane", &tEvts.hiNevtPlane, &tEvts.b_hiNevtPlane);
//   if (t->GetBranch("hiEvtPlanes")) t->SetBranchAddress("hiEvtPlanes", tEvts.hiEvtPlanes, &tEvts.b_hiEvtPlanes);
//   if (doCheck) {
//     if (t->GetMaximum("hiNevtPlane")>126) { cout <<"FATAL ERROR: Arrary size of hiNevtPlane too small!!!  "<<t->GetMaximum("hiNevtPlane")<<endl; exit(0);
//     }   }
}




