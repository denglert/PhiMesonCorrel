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

  if(argc != 5)
  {
    std::cerr << "Usage: preprocess <.root file to be preprocessed> <tag> <nEvents>" << std::endl;
	 exit(1);
  }

 std::cout << "preprocess_genlevel started" << std::endl;


 std::string inpFilename   = argv[1];
 std::string tag		  = argv[2];
 std::string PIDconfig = argv[3];
 int nEvMax 	  		  = atoi( argv[4] );

 // Binning
 int nCorrTyp			  = nCorrTyp_; 
 int *nPtBins = new int[nCorrTyp_];

 for(int TypBin = 0; TypBin < nCorrTyp; TypBin++)
 { nPtBins[TypBin] = nPtBins_[TypBin]; }

 int nMultiplicityBins_Ana = nMultiplicityBins_Ana_HDR;
 int nMultiplicityBins_EvM = nMultiplicityBins_EvM_HDR;
 int nZvtxBins 		      = nZvtxBins_; 

 PIDUtil *pidutil = new PIDUtil;
 pidutil->ReadInConfig(PIDconfig);

 // Log
 LogFile *log = new LogFile("log");
 log->repeat = 1000;
 log->label = "DATA";


 //////////////////////////////////////
 //                                  //
 // ****** Opening input file ****** //
 //                                  //
 //////////////////////////////////////

 // Open file
 TFile *f;
 f = NULL;
 f = TFile::Open(inpFilename.c_str());
 if ( f->IsZombie() || (f == NULL) ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << Form("TFile %s seem to have loaded.\n", inpFilename.c_str()); }

 // trackTree
 TTree *trackTree = (TTree*)f->Get("ppTrack/trackTree");
 Particles tTracks;
 bool isMC = true;
 setupParticleTree(trackTree, tTracks);

 // hiEvtAnalyzer
 TTree *EvtAna= (TTree*)f->Get("hiEvtAnalyzer/HiTree");
 int hiNtracks; EvtAna->SetBranchAddress("hiNtracks", &hiNtracks);
 float vz; EvtAna->SetBranchAddress("vz", &vz);

 // Event Selection (SkimAnalysis)
 EvtSelection EvSel;
 EvSel.setupSkimTree_pPb( f, isMC);

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

 CFW.DoTrackWeight = false;
 CFW.SetupForPreprocess();
 if ( CFW.DoSelfCorrelation ) { std::cout << "Analysis includes self correlation computation." << std::endl;}



 // EventData
 EventData *ev;
 ev = new EventData;

 ev->Setup_nTriggerParticles(nCorrTyp, nPtBins);

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

			trackTree->GetEntry(iEv);
			int nPar = tTracks.nParticle;

	 		for (int iPar = 0; iPar < nPar; iPar++)
			{

				// Particle selection //
				// Same kinematical cuts as at RECO level!
				if ( 2.4 < abs(tTracks.pEta[iPar]) )  continue;
				bool isOutsideReferencePartPtRange = ( ( tTracks.pPt[iPar] < ptref1 ) || ( ptref2 < tTracks.pPt[iPar] ) );
				if ( isOutsideReferencePartPtRange ) continue;	
				// Particle selection //


				// Particle fill up
				track particle;
				// Warning: no charge fill!;
				particle.charge  = 0;

				particle.phi     = tTracks.pPhi[iPar];
				particle.eta     = tTracks.pEta[iPar];
	
				ev->AddTrack(particle);
			}

					
			ev->EventID = iEv;
			EventCache[multBin][zvtxBin].push_back( (*ev) );
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

	// Load in tracks
	trackTree->GetEntry(iEvA);

	// Read in event
	ev->ReadInMC( tTracks, pidutil );

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
