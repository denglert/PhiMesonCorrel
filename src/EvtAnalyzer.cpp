#include "EvtAnalyzer.h"

EvtAnalyzer::EvtAnalyzer()
{
	// Default
	vz_min   = -13;
	vz_max   =  13;
	Ntrk_min =   0;
	Ntrk_max = 120;
};



void EvtAnalyzer::setupEvtAnaTree( TFile *f )
{
	  EvtAna = (TTree*)f->Get("hiEvtAnalyzer/HiTree");

     EvtAna->SetBranchAddress("hiNtracks", &hiNtracks, &b_hiNtracks);
     EvtAna->SetBranchAddress("vz", &vz, &b_vz);
}


void EvtAnalyzer::GetEntry( int iEv )
{
	EvtAna->GetEntry(iEv);
}

int EvtAnalyzer::gethiNtracks()
{
	return hiNtracks;
}

double EvtAnalyzer::getvz()
{
	return vz;
}

bool EvtAnalyzer::isEvPass()
{
	if (  ((  vz_min < vz        ) && (  vz_max < vz        )) &&
		   ((Ntrk_min <= hiNtracks) && (Ntrk_max <= hiNtracks))    )
	{return true;}
	else false;
};

bool EvtAnalyzer::isInside_vz( )
{
	if ( (vz_min < vz) || (vz_max < vz) ) {return true;}
	else false;
};
