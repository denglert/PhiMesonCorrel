#ifndef SETUPCUSTOMTRACKTREE_H
#define SETUPCUSTOMTRACKTREE_H
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

class Tracks {
public :
   Tracks(){};
   ~Tracks(){};

   // Declaration of leaf types
   Int_t           nTrk;
   Float_t 			 dedx[maxTracks];
   Float_t         trkPt[maxTracks];   //[nTrk]
   Float_t         trkPtError[maxTracks];   //[nTrk]
   Float_t         trkEta[maxTracks];   //[nTrk]
   Float_t         trkPhi[maxTracks];   //[nTrk]
   Int_t           trkCharge[maxTracks];   //[nTrk]
   Bool_t          highPurity[maxTracks];   //[nTrk]
   Float_t         trkDxy1[maxTracks];   //[nTrk]
   Float_t         trkDxyError1[maxTracks];   //[nTrk]
   Float_t         trkDz1[maxTracks];   //[nTrk]
   Float_t         trkDzError1[maxTracks];   //[nTrk]

   // List of branches
   TBranch        *b_nTrk;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
};

class Particles {
public :
   Particles(){};
   ~Particles(){};

   // Declaration of leaf types
	Int_t           nParticle;
   Int_t           pStatus[maxTracks];   //[nParticle]
   Int_t           pPId[maxTracks];   //[nParticle]
   Float_t         pEta[maxTracks];   //[nParticle]
   Float_t         pPhi[maxTracks];   //[nParticle]
   Float_t         pPt[maxTracks];   //[nParticle]

   // List of branches
   TBranch        *b_nParticle;   //!
   TBranch        *b_pStatus;   //!
   TBranch        *b_pPId;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pPt;   //!
};


// Tracks class for track correction calculation
class Tracks_c {
public :
   Tracks_c(){};
   ~Tracks_c(){};

   // Declaration of leaf types
   Int_t           nEv;
   Int_t           nLumi;
   Int_t           nBX;
   Int_t           nRun;
   Int_t           nv;
   Float_t         vx[MAXVTX];   //[nv]
   Float_t         vy[MAXVTX];   //[nv]
   Float_t         vz[MAXVTX];   //[nv]
   Float_t         vxErr[MAXVTX];   //[nv]
   Float_t         vyErr[MAXVTX];   //[nv]
   Float_t         vzErr[MAXVTX];   //[nv]
   Int_t           nDaugher[MAXVTX];   //[nv]
   Int_t           nVtx;
   Int_t           maxVtx;
   Int_t           nTrkVtx[MAXVTX];   //[nVtx]
   Float_t         xVtx[MAXVTX];   //[nVtx]
   Float_t         yVtx[MAXVTX];   //[nVtx]
   Float_t         zVtx[MAXVTX];   //[nVtx]
   Float_t         xVtxErr[MAXVTX];   //[nVtx]
   Float_t         yVtxErr[MAXVTX];   //[nVtx]
   Float_t         zVtxErr[MAXVTX];   //[nVtx]
   Int_t           nTrk;
   Float_t 			 dedx[maxTracks];
   Float_t         trkPt[maxTracks];   //[nTrk]
   Float_t         trkPtError[maxTracks];   //[nTrk]
   Int_t           trkNHit[maxTracks];   //[nTrk]
   Int_t           trkNlayer[maxTracks];   //[nTrk]
   Float_t         trkEta[maxTracks];   //[nTrk]
   Float_t         trkPhi[maxTracks];   //[nTrk]
   Int_t           trkCharge[maxTracks];   //[nTrk]
   Bool_t          highPurity[maxTracks];   //[nTrk]
   Float_t         trkChi2[maxTracks];   //[nTrk]
   Float_t         trkNdof[maxTracks];   //[nTrk]
   Float_t         trkDxy1[maxTracks];   //[nTrk]
   Float_t         trkDxyError1[maxTracks];   //[nTrk]
   Float_t         trkDz1[maxTracks];   //[nTrk]
   Float_t         trkDzError1[maxTracks];   //[nTrk]
   Bool_t          trkFake[maxTracks];   //[nTrk]
   Float_t         trkAlgo[maxTracks];   //[nTrk]
   Float_t         trkStatus[maxTracks];   //[nTrk]
	Int_t           nParticle;
   Int_t         pStatus[maxTracks];   //[nParticle]
   Int_t         pPId[maxTracks];   //[nParticle]
   Float_t         pEta[maxTracks];   //[nParticle]
   Float_t         pPhi[maxTracks];   //[nParticle]
   Float_t         pPt[maxTracks];   //[nParticle]
   Float_t         pAcc[maxTracks];   //[nParticle]
   Float_t         pAccPair[maxTracks];   //[nParticle]
   Int_t         pNRec[maxTracks];   //[nParticle]
   Int_t           pNHit[maxTracks];   //[nParticle]
   Float_t         mtrkdedx[maxTracks];   //[nParticle]
   Float_t         mtrkPt[maxTracks];   //[nParticle]
   Float_t         mtrkPtError[maxTracks];   //[nParticle]
   Int_t           mtrkNHit[maxTracks];   //[nParticle]
   Int_t           mtrkQual[maxTracks];   //[nParticle]
   Float_t         mtrkChi2[maxTracks];   //[nParticle]
   Float_t         mtrkNdof[maxTracks];   //[nParticle]
   Float_t         mtrkDz1[maxTracks];   //[nParticle]
   Float_t         mtrkDzError1[maxTracks];   //[nParticle]
   Float_t         mtrkDxy1[maxTracks];   //[nParticle]
   Float_t         mtrkDxyError1[maxTracks];   //[nParticle]
   Float_t         mtrkDz2[maxTracks];   //[nParticle]
   Float_t         mtrkDzError2[maxTracks];   //[nParticle]
   Float_t         mtrkDxy2[maxTracks];   //[nParticle]
   Float_t         mtrkDxyError2[maxTracks];   //[nParticle]


   // List of branches
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_nv;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vxErr;   //!
   TBranch        *b_vyErr;   //!
   TBranch        *b_vzErr;   //!
   TBranch        *b_nDaugher;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_maxVtx;   //!
   TBranch        *b_nTrkVtx;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_xVtxErr;   //!
   TBranch        *b_yVtxErr;   //!
   TBranch        *b_zVtxErr;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!
   TBranch        *b_trkStatus;   //!
   TBranch        *b_nParticle;   //!
   TBranch        *b_pStatus;   //!
   TBranch        *b_pPId;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pPt;   //!
   TBranch        *b_pAcc;   //!
   TBranch        *b_pAccPair;   //!
   TBranch        *b_pNRec;   //!
   TBranch        *b_pNHit;   //!
   TBranch        *b_mtrkdedx;   //!
   TBranch        *b_mtrkPt;   //!
   TBranch        *b_mtrkPtError;   //!
   TBranch        *b_mtrkNHit;   //!
   TBranch        *b_mtrkNlayer;   //!
   TBranch        *b_mtrkNlayer3D;   //!
   TBranch        *b_mtrkQual;   //!
   TBranch        *b_mtrkChi2;   //!
   TBranch        *b_mtrkNdof;   //!
   TBranch        *b_mtrkDz1;   //!
   TBranch        *b_mtrkDzError1;   //!
   TBranch        *b_mtrkDxy1;   //!
   TBranch        *b_mtrkDxyError1;   //!
   TBranch        *b_mtrkDz2;   //!
   TBranch        *b_mtrkDzError2;   //!
   TBranch        *b_mtrkDxy2;   //!
   TBranch        *b_mtrkDxyError2;   //!


};


void setupTrackTree(TTree *t, Tracks &tTracks );
void setupTrackTree_c(TTree *t,Tracks_c &tTracks,bool doCheck);
void setupParticleTree(TTree *t, Particles &tTracks);
#endif
