#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include "SetupCustomTrackTree.h"
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <iostream>
#include <TLorentzVector.h>
#include "CorrelationUtils.h"
#include "AnalysisFW.h"
#include "PhyConst.h"
#include "PIDUtils.h"
#include "EvtSelection.h"

const bool DoMtrk = true;
const bool doMC = true;

int main(int argc, const char *argv[])
{ 

  if(argc != 4)
  {
    std::cerr << "Usage: preprocess <.root file to be processed> <job-id> <nEvents>" << std::endl;
	 exit(1);
  }


 TString inpFilename   = 		 argv[1];
 std::string jobid 	  = 		 argv[2];
 int nEvMax 	  		  = atoi( argv[3] );

 ///////////////////////////////////////////
 // Opening input file, setting up TTrees //
 ///////////////////////////////////////////
 
 // Open file
 TFile *f = TFile::Open(inpFilename);
 if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << "Successfully opened file " << inpFilename << std::endl ;}

 // trackTree
 TTree *trackTree = (TTree*)f->Get("ppTrack/trackTree");
 Tracks tTracks;
 setupTrackTree(trackTree, tTracks, doMC);

 EvtSelection EvSel; const bool isOLD = false;
 EvSel.setupSkimTree_pPb( f, isOLD);

 LogFile log("log");
 log.repeat = 1000;

 ////////////////////
 // Setting output //
 ////////////////////
 
 TFile *output = new TFile(Form("./dedxmaps_%s.root", jobid.c_str() ),"RECREATE");
 output->cd();

 //////////////////
 // Initializing //
 //////////////////
 
 dEdxMaps *dEdxmaps_trk  = new dEdxMaps("trk");
 dEdxMaps *dEdxmaps_mtrk = new dEdxMaps("mtrk");

 TH1D *pDistr = new TH1D("pDistr", ";p;nTrks", 30, 0., 3.);

 // Histograms with errors
 TH1::SetDefaultSumw2( );
 TH2::SetDefaultSumw2( );

 /////////////////////////////
 // ----- EVENT LOOP ------ //
 /////////////////////////////
 
 if (nEvMax == -1) {nEvMax = trackTree->GetEntries();}
 for (int iEv = 0; iEv < nEvMax; iEv++)
 {
	
	// Event counter info
	log.EventCounter(iEv);

	// EventSelection
	if ( !EvSel.isGoodEv_pPb( iEv ) ) continue;

	// Load in tracks
	trackTree->GetEntry(iEv);
	int nTrkA = tTracks.nTrk;

	////////////////////////
	// --- TRACK LOOP --- //
	////////////////////////
 	for (int iTrk = 0; iTrk < nTrkA; iTrk++)
	{
		// *** Track selection *** //
		if ( !TrackSelection(tTracks, iTrk ) ) continue;

		double pt  = tTracks.trkPt[iTrk];
		double eta = tTracks.trkEta[iTrk];
		float p    = pt * cosh(eta);
		float dedx = tTracks.dedx[iTrk];

		int PID = GetPID(p, dedx, eta);

		dEdxmaps_trk->Fill(PID, p, dedx);

	}

	if ( DoMtrk )
	{

	// === Particle loop === //
	for(int iPart = 0; iPart < tTracks.nParticle; iPart++)
	{
		double meta = tTracks.pEta[iPart]; 

		//// matched track selection
		if( !( mTrackSelection( tTracks, iPart) )) continue;
		// particle selection
		if( 2.4 < fabs(meta) ) continue;

		double mpt   = tTracks.mtrkPt [iPart];
		double mdedx = tTracks.mtrkdedx [iPart];
		float  mp    = mpt * cosh(meta);
		int    mPID  = GetPID(mp, mdedx, meta);

		pDistr->Fill(mp);

		dEdxmaps_mtrk->Fill(mPID, mp, mdedx);

	}


 	}

 }

 // Save dEdx value map
 //viewdEdxvsP(dEdxvsp, Form("./PIDtest/dedxvsp_%s.png", jobid.c_str()), delta );
 
 //for(int ptBin = 0; ptBin < nPtBins; ptBin++)
 //{ viewdEdxvsP(dEdxvspPt[ptBin], Form("./PIDtest/dedxvsp_%s_pt_%.2f-%.2f.png", jobid.c_str(), trigptbins[ptBin][0], trigptbins[ptBin][1]), delta );}


 dEdxmaps_trk->PlotFigs("lel");

  log.Close();

 // Output .root file
 output->Write();
 output->Close();

}
