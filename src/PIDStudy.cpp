#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include "PIDUtils.h"
#include "AnalysisFW.h"
#include "EvtAnalyzer.h"
#include "EvtSelection.h"

int main( int argc, const char *argv[] )
{

  if(argc != 5)
  {
    std::cerr << "Usage: PIDStudy <sample source> <isMC> <PIDconfig> <nEvents> " << std::endl;
	 exit(1);
  }

 std::string inpFilename = argv[1];
 std::string isMC 		 = argv[2];
 int nEvMax 	  		    = atoi( argv[3] );
 std::string PIDconfig	 = argv[4];

 // Open file
 TFile *f = TFile::Open( inpFilename.c_str() );
 if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}

 //////////////////
 // Initializing //
 //////////////////
 
 PIDUtil *pidutil = new PIDUtil;
 pidutil->ReadInConfig(PIDconfig);

 LogFile log("log");
 log.repeat = 1000;

 TTree *trackTree;

 // trackTree
 if ( isMC == "no" )
 { trackTree = (TTree*)f->Get("pptracks/trackTree"); }
 else if ( isMC == "yes" )
 { trackTree = (TTree*)f->Get("ppTrack/trackTree");  }


 Tracks tTracks;

 bool doMC;
 if (isMC == "no")
 { doMC = false; }
 else if ( isMC == "yes")
 { doMC = true;  }
 
 setupTrackTree(trackTree, tTracks, doMC);

 // EventAnalyzer
 EvtAnalyzer EvAna;
 EvAna.setupEvtAnaTree( f );

 log.wr(Form("EvtAna.vz:  [ %.2f - %.2f ] ", EvAna.vz_min, EvAna.vz_max));
 log.wr(Form("EvtAna.Ntr: [ %d - %d ] ",     EvAna.Ntrk_min, EvAna.Ntrk_max));

 // Event Selection (SkimAnalysis)
 EvtSelection EvSel;
 EvSel.setupSkimTree_pPb( f, doMC);

 
 // Analysis specific //
 
 double pBins[npBinslog+1];
 double dEdxBins[ndEdxBinslog+1];

 double l10 = TMath::Log(10);
 
 double dplog   = ( TMath::Log(pmaxlog)   -TMath::Log(pminlog)    )/npBinslog/l10;
 double ddEdxlog = ( TMath::Log(dEdxmaxlog)-TMath::Log(dEdxminlog) )/ndEdxBinslog/l10;

 for (int i=0; i<=npBinslog; i++)
 { pBins[i] = TMath::Exp(l10*(i*dplog + TMath::Log(pminlog)/l10)); }

 for (int i=0; i<=ndEdxBinslog; i++)
 { dEdxBins[i] = TMath::Exp(l10*(i*ddEdxlog+ TMath::Log(dEdxminlog)/l10)); }

 TH2D *dEdxvsPMapsLin[npt];
 TH2D *dEdxvsPMapsLog[npt];

 for( int ptBin = 0; ptBin < npt; ptBin++ )
 {

   double pt1 = (ptMin + (ptBin  ) * ptbw);
	double pt2 = (ptMin + (ptBin+1) * ptbw);

	dEdxvsPMapsLin[ptBin] = new TH2D(Form("dEdxvsPMapsLin_pt_%.2f-%.2f", pt1, pt2), ";p [GeV/c];dE/dx [MeV/cm]", npBins, pminlin, pmaxlin, ndEdxBins, dEdxminlin, dEdxmaxlin);
	dEdxvsPMapsLog[ptBin] = new TH2D(Form("dEdxvsPMapsLog_pt_%.2f-%.2f", pt1, pt2), ";p [GeV/c];dE/dx [MeV/cm]", npBinslog, pBins, ndEdxBinslog, dEdxBins);

 }

 ////////////////
 // Event loop //
 ////////////////
 
 if (nEvMax == -1) {nEvMax = trackTree->GetEntries();}

 for (int iEvA = 0; iEvA < nEvMax; iEvA++)
 {

 	// EventCounter
 	log.EventCounter(iEvA);
 	EvAna.GetEntry  (iEvA);
 	
 	// Event Selection
 	if ( !EvSel.isGoodEv_pPb( iEvA ) ) continue;
 	if ( !EvAna.isEvPass( )          ) continue;

 	
 	// Tracks & particles
 	trackTree->GetEntry(iEvA);
 	
 	int nTrk = tTracks.nTrk;
 	
 	// Track loop
 	for (int iTrk = 0; iTrk < nTrk; iTrk++)
 	{

 		// *** Track selection *** //
 		if ( !TrackSelection(tTracks, iTrk ) ) continue;

 	}

 	
 }


}
