#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>

class EvtSelection {
public :
   EvtSelection(){};
   ~EvtSelection(){};

	TTree *SkimAna;

   Int_t    pPAcollisionEventSelectionPA;
   Int_t    pVertexFilterCutGplus;

   TBranch *b_pPAcollisionEventSelectionPA;
   TBranch *b_pVertexFilterCutGplus;

	void setupSkimTree_pPb( TFile *f , bool isOLD);
	bool isGoodEv_pPb( int iEv );
	int GetEntries();
};

#endif
