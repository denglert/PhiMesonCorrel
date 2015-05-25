#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <iostream>
#include <deque>
#include <TLorentzVector.h>
#include "AnalysisFW.h"
#include "AnalysisBinning.h"
#include "PIDUtils.h"
#include "CorrelationUtils.h"
#include "SetupCustomTrackTree.h"
#include "EvtSelection.h"

int main( int argc, const char *argv[] )
{ 

  if(argc != 8)
  {
    std::cerr << "Usage: preprocess <.root file to be preprocessed> <dotrkCorr> <tag> <trkCorrFileName> <PIDconfig>  <tag> <nEvents>" << std::endl;
	 exit(1);
  }

 std::string inpFilename = argv[1];
 std::string isMC        = argv[2];
 TString dotrkCorr_str 	 = argv[3];
 std::string trkCorrFilename = argv[4];
 std::string PIDconfig   = argv[5];
 std::string tag		    = argv[6];
 int nEvMax 	  		    = atoi( argv[7] );

 // Binning
 int nCorrTyp			  = nCorrTyp_; 
 int *nPtBins = new int[nCorrTyp_];

 for(int TypBin = 0; TypBin < nCorrTyp; TypBin++)
 { nPtBins[TypBin] = nPtBins_[TypBin]; }

 int nMultiplicityBins_Ana = nMultiplicityBins_Ana_HDR;
 int nMultiplicityBins_EvM = nMultiplicityBins_EvM_HDR;
 int nZvtxBins 		      = nZvtxBins_; 

 PIDUtil *pidutil = new PIDUtil();
 pidutil->ReadInConfig( PIDconfig );

 // Log
 LogFile *log = new LogFile(Form("log_%s", tag.c_str()));
 log->repeat = 1000;
 log->label = "DATA";

 //////////////////////////////////////
 //                                  //
 // ****** Opening input file ****** //
 //                                  //
 //////////////////////////////////////
 
 // Open file
 TFile *f = NULL;
 std::cout << "\nBefore open f: " << f << std::endl;
 f = TFile::Open(inpFilename.c_str());
 std::cout << "\nAfter open f: " << f << std::endl;
 if ( f->IsZombie() || (f == NULL) ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << Form("TFile %s seem to have loaded.\n", inpFilename.c_str()); }

 bool doMC;
 if (isMC == "no")
 { doMC = false; std::cout << "doMC: false" << std::endl; }
 else if ( isMC == "yes")
 { doMC = true; std::cout << "doMC: false" << std::endl; }

 TTree *trackTree;

 // trackTree
 if ( doMC == false )
 { trackTree = (TTree*)f->Get("pptracks/trackTree"); std::cout << "No MC - pptracks/trackTree" << std::endl; }
 else if ( doMC == true )
 { trackTree = (TTree*)f->Get("ppTrack/trackTree"); std::cout << "Yes MC - ppTrack/trackTree" << std::endl; }

 Tracks tTracks;
 
 setupTrackTree(trackTree, tTracks);

 // hiEvtAnalyzer
 TTree *EvtAna = NULL; 
 std::cout << "EvtAnA: (before load)" << EvtAna << std::endl;
 EvtAna = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
 int hiNtracks; EvtAna->SetBranchAddress("hiNtracks", &hiNtracks);
 float vz; EvtAna->SetBranchAddress("vz", &vz);

 std::cout << "EvtAnA: (after load)" << EvtAna << std::endl;

 // Event Selection (SkimAnalysis)
 EvtSelection EvSel;
 EvSel.setupSkimTree_pPb( f, doMC);

 ////////////////////////////////
 //                            //
 // ***** Setting output ***** //
 //                            //
 ////////////////////////////////
 
 TFile *output = new TFile(Form("./correl_analysis_%s.root", tag.c_str() ),"RECREATE");
 output->cd();

 //////////////////////////////
 //                          //
 // ***** Initializing ***** //
 //                          //
 //////////////////////////////
 
 TH1::SetDefaultSumw2( );
 TH2::SetDefaultSumw2( );

 // EventCache
 std::deque< EventData > **EventCache;
 Setup_EventCache(EventCache, nMultiplicityBins_EvM, nZvtxBins);

 std::deque< EventData > *** EventCache_ptr = &EventCache;

 // dEdxvsp map
 TH2D *dEdxvsp  = new TH2D ("dEdxVsP",";p(GeV/c);dE/dx", npBins, pMin, pMax, ndEdxBins, dEdxMin, dEdxMax);

 // Number of tracks distribution
 TH1D *nTrkDistr_signal = new TH1D("nTrkDistr_signal","Track distribution;Multiplicity", 350, 0, 350);


 // Correlation Framework
 CorrelationFramework CFW(nCorrTyp, nPtBins, nMultiplicityBins_Ana, nMultiplicityBins_EvM);
 std::cout << "Correlation Analysis Framework loaded." << std::endl;

 CFW.DoSelfCorrelation = false;
 if ( CFW.DoSelfCorrelation ) { std::cout << "Analysis includes self correlation computation." << std::endl;}

 CFW.SetupForPreprocess();

 // TrackCorrection file
 TFile *f_trkCorr = NULL; 
 f_trkCorr = new TFile(trkCorrFilename.c_str(), "READ");
 if ( (f_trkCorr->IsZombie()) || (f_trkCorr == NULL) ) {std::cerr << "Error opening file: " << trkCorrFilename << std::endl; exit(-1);}
 else {std::cout << Form("trkCorr file %s successfully opened.", trkCorrFilename.c_str()) << std::endl;}

 TrackCorr *trkCorr = new TrackCorr;
 trkCorr->table = Read_TH3D_1Darray(f_trkCorr, "hcorr3D typ", nCorrTyp);

 std::cout << "dotrkCorr: " << dotrkCorr_str << std::endl;

 bool dotrkCorr;
      if( dotrkCorr_str == "yes" ) { dotrkCorr =  true; }
 else if( dotrkCorr_str == "no" ) { dotrkCorr = false; }
 else {std::cerr << "dotrkCorr not defined." << std::endl; exit(-1);}

 trkCorr->DoTrackWeight = dotrkCorr;

 std::cerr << "detaching trkCorr" << std::endl;
 for (int i = 0; i < nCorrTyp ; i++)
 { trkCorr->table[i]->SetDirectory(0); }

 f_trkCorr->Close();
 std::cerr << "trkCorr file closed." << std::endl;

 // EventData
 EventData *ev;
 ev = new EventData;

 ev->Setup_nTriggerParticles(nCorrTyp, nPtBins);
 std::cerr << "TriggerParticles set up " << std::endl;

 ///////////////////////////////////////////
 //                                       //
 // **** PRELOAD MIXEVENTS IN MEMORY **** //
 //                                       //
 ///////////////////////////////////////////
 
 log->wr(Form("trackTree entries: %d", trackTree->GetEntries()));
 log->wr(Form("EventSelection (SkimAna) entries: %d", EvSel.GetEntries()));
 log->wr(Form("nEvMax: %d", nEvMax));

 std::cout << "Preloading events in memory..." << std::endl;
 for(int multBin = 0; multBin < nMultiplicityBins_EvM; multBin++)
 for(int zvtxBin = 0; zvtxBin < nZvtxBins_; zvtxBin++)
 {
	
	int count = 0;
	int nev = 10;
	
	int iEv = 0;


 	while ( (count < (nev+1)) && (iEv < trackTree->GetEntries() )) 
	{

		// Get current event info
		EvtAna->GetEntry(iEv);

		// Event Selection //

		if ( EvSel.isGoodEv_pPb( iEv ) )
		if (     zvtxbin(vz, nZvtxBins_)    == zvtxBin )
		if ( multiplicitybin_EvM(hiNtracks) == multBin )
		{

			tTracks.b_nTrk->GetEntry(iEv);
			if ( maxTracks < tTracks.nTrk) continue;
			trackTree->GetEntry(iEv);
			int nTrk = tTracks.nTrk;
			if ( maxTracks < tTracks.nTrk) continue;

	 		for (int iTrk = 0; iTrk < nTrk; iTrk++)
			{

				// TRACK SELECTION //
				if ( !TrackSelection(tTracks, iTrk ) ) continue;
		
				float pt  = tTracks.trkPt [iTrk];
				float eta = tTracks.trkEta[iTrk];
				float phi = tTracks.trkPhi[iTrk];
		
		  		track trk;
		
				trk.IsInsideReferencePtRange = ( (ptref1 < pt) && ( pt < ptref2 ));
				if ( !trk.IsInsideReferencePtRange) continue;
		
				trk.w0 = trkCorr->trackWeight( 0 , pt, eta, phi); 
				
		  		// Track fill up
		  		trk.charge  = tTracks.trkCharge[iTrk];
		  		trk.phi     = phi;
		  		trk.eta     = eta;

				ev->AddTrack(trk);
			}
					
			ev->EventID = iEv;
			(*EventCache_ptr)[multBin][zvtxBin].push_back( (*ev) );
			count++;

			ev->Clear(nCorrTyp, nPtBins);

		}

		// Event counter info
		if ( ((iEv % 10000) == 0) || (count == 11) )
		{ std::cout << Form("multBin: %02d, zvtxBin: %02d, event: %05d, found: %02d/10", multBin, zvtxBin, iEv, count) << std::endl; }


		iEv++;
	 }

 }

 std::cout << "Preloading completed." << std::endl << std::endl;
 std::cout << "EventCache statistics:" << std::endl;

 // Deleting first elements of deque
 for(int multBin = 0; multBin < nMultiplicityBins_EvM; multBin++)
 for(int zvtxBin = 0; zvtxBin < nZvtxBins_; zvtxBin++)
 { 
	EventCache[multBin][zvtxBin].pop_front();
   std::cout << Form("multBin: %3d, zvtxBin: %3d, found: %2d/10", multBin, zvtxBin, EventCache[multBin][zvtxBin].size()) << std::endl;
 }

 log->wr(Form("trackTree entries: %d", trackTree->GetEntries()));
 log->wr(Form("EventSelection (SkimAna) entries: %d", EvSel.GetEntries()));
 log->wr(Form("nEvMax: %d", nEvMax));

 ///////////////////////////
 //                       //
 // ***** ANALYISIS ***** //
 //                       //
 ///////////////////////////
 

 if (nEvMax == -1) {nEvMax = trackTree->GetEntries();}

 log->wr(Form("trackTree entries: %d", trackTree->GetEntries()));
 log->wr(Form("SkimAna entries: %d", EvSel.GetEntries()));
 log->wr(Form("nEvMax: %d", nEvMax));

 for (int iEvA = 0; iEvA < nEvMax; iEvA++)
 {
	
	log->EventCounter(iEvA);

	// Get current event info
	EvtAna->GetEntry(iEvA);

	// Event Selection
	if ( !EvSel.isGoodEv_pPb( iEvA ) ) continue;
	if ( zvtxbin(vz, nZvtxBins) == -1 ) continue;
	CFW.nEvents_Processed_signal_total->Fill(0.);
	if ( multiplicitybin_Ana(hiNtracks, nMultiplicityBins_Ana) == -1) continue;

 	ev->Clear(nCorrTyp, nPtBins);

	ev->EventID = iEvA;
	ev->SetnTrk(hiNtracks);
	ev->SetzVtx(vz);

	// Statistics
	CFW.nEvents_Processed_signal_total->Fill(1.);
	nTrkDistr_signal->Fill( hiNtracks );
	
	CFW.ResetCurrentEventCorrelation();

	tTracks.b_nTrk->GetEntry(iEvA);
	if ( maxTracks < tTracks.nTrk) continue;
	trackTree->GetEntry(iEvA);

	// Read in event
	ev->ReadInDATA(tTracks, pidutil, trkCorr, CFW.spectrum);

	CFW.SignalCorrelation(ev);
	CFW.MixedCorrelation(ev, EventCache_ptr);
	CFW.AddCurrentEventCorrelation(ev);

 }
 
 //////////////////////
 //                  //
 // **** OUTPUT **** //
 //                  //
 //////////////////////
 
 output->Write();
 output->Close();


}
