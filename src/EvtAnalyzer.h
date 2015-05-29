#ifndef EVENTANALYZER_H
#define EVENTANALYZER_H

#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <iostream>

class EvtAnalyzer {
public :
   EvtAnalyzer();
   ~EvtAnalyzer(){};

	// TTree
	TTree *EvtAna;

	// Variables
   Int_t    hiNtracks;
   Float_t    vz;

	float vz_min;
	float vz_max;
	int Ntrk_min;
	int Ntrk_max;

	// Branches
   TBranch *b_hiNtracks;  
   TBranch *b_vz;  

	// Functions
	void setupEvtAnaTree( TFile *f );
	void GetEntry( int iEv );
	int gethiNtracks();
	double getvz();
	bool isInside_vz();
	bool isEvPass();

};

#endif
