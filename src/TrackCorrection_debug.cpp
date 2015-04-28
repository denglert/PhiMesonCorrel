#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <iostream>
#include <deque>
#include <TLorentzVector.h>
//#include "AnalysisFW.h"
//#include "../HiForestAnalysis/hiForest.h"
#include "AnalysisBinning.h"
#include "PIDUtils.h"
#include "AnalysisFW.h"
#include "SetupCustomTrackTree.h"
#include "dataType.h"
#include "EvtSelection.h"
#include "EvtAnalyzer.h"

// zvtx distribution parameters
const int nZvtxDistrBins  =  13;
const int nCorrTyp = 5;
const double zVtxDistrMin = -13;
const double zVtxDistrMax =  13;

const int npt[nCorrTyp]      = { 15,     8,    7,   14,     6 };
const double ptMin[nCorrTyp] = { 0.3, 0.20, 0.20, 0.20,  1.00 };
const double ptMax[nCorrTyp] = { 3.0, 1.00, 0.90, 1.60,  1.60 };

const double ptminlog[nCorrTyp] = { 0.3, 0.20, 0.20, 0.20, 1.00 };
const double ptmaxlog[nCorrTyp] = { 3.0, 1.00, 0.90, 1.60, 1.60 };
const int      nptlog[nCorrTyp] = { 15,     8,    7,   14,    6 };



// const int npt[4]      = {1,    1,   1,  1};
// const double ptMin[4] = {0.3, 0.15, 0.15, 0.15};
// const double ptMax[4] = {0.5, 0.40, 0.40, 0.4 };

const double ptGlobalMin = 0.1;

//const int neta[4]      = { 48 , 16 , 16  , 16 };
const int neta[nCorrTyp]      = {  12 ,    4,    4,    4,    4 };
const double etaMin[nCorrTyp] = { -2.4, -0.8, -0.8, -0.8, -0.8 };
const double etaMax[nCorrTyp] = {  2.4,  0.8,  0.8,  0.8,  0.8 };


//const int nphi = 70;
const int nphi = 20;
const double phiMin = -TMath::Pi();
const double phiMax =  TMath::Pi();

TH3D** Setup_TH3D_nCorrTyp(const char histoname[], const char histolabel[], const int nXBins[], const double xmin[], const double xmax[], const int nYBins[], const double ymin[], const double ymax[], int nZBins, double zmin, double zmax )
{
	TH3D **histos = new TH3D*[nCorrTyp];

	for( int i = 0; i < nCorrTyp; i++)
	{
		 histos[i] = new TH3D( Form("%s typ %d", histoname, i), histolabel, nXBins[i], xmin[i], xmax[i], nYBins[i], ymin[i], ymax[i], nZBins, zmin, zmax);
	}

	return histos;
};

TH2D** Setup_TH2D_nCorrTyp(const char histoname[], const char histolabel[], const int nXBins[], const double xmin[], const double xmax[], const int nYBins[], const double ymin[], const double ymax[])
{
	TH2D **histos = new TH2D*[nCorrTyp];

	for( int i = 0; i < nCorrTyp; i++)
	{
		 histos[i] = new TH2D( Form("%s typ %d", histoname, i), histolabel, nXBins[i], xmin[i], xmax[i], nYBins[i], ymin[i], ymax[i]);
	}

	return histos;
};

TH1D** Setup_TH1D_nCorrTyp(const char histoname[], const char histolabel[], const int nXBins[], const double xmin[], const double xmax[])
{
	TH1D **histos = new TH1D*[nCorrTyp];

	for( int i = 0; i < nCorrTyp; i++)
	{
		 histos[i] = new TH1D( Form("%s typ %d", histoname, i), histolabel, nXBins[i], xmin[i], xmax[i]);
	}

	return histos;
};

bool isInsidePt(int PID, double pt)
{
	if (PID == 0) {};
}


int main( int argc, const char *argv[] )
{


  if(argc != 7)
  {
    std::cerr << "Usage: TrackCorrection <MC sample> <DATA sample> <tag> <nEventsDATA> <nEventsMC> <PIDconfig>" << std::endl;
	 exit(1);
  }

 std::string inpFilenameDATA = argv[1];
 std::string inpFilenameMC   = argv[2];
 std::string tag		        = argv[3];
 int nEvMaxDATA 	  		     = atoi( argv[4] );
 int nEvMaxMC 	  		        = atoi( argv[5] );
 std::string PIDconfig		  = argv[6];
 
 PIDUtil *pidutil = new PIDUtil;
 pidutil->ReadInConfig(PIDconfig);
 
 // Binning
 int *nPtBins = new int[nCorrTyp_];

 for(int TypBin = 0; TypBin < nCorrTyp; TypBin++)
 { nPtBins[TypBin] = nPtBins_[TypBin]; }

 int nMultiplicityBins_Ana = nMultiplicityBins_Ana_HDR;
 int nMultiplicityBins_EvM = nMultiplicityBins_EvM_HDR;
 int nZvtxBins 		      = nZvtxBins_; 

 LogFile log("log");
 log.repeat = 1000;

 ///////////////////////////
 // Global initialization //
 ///////////////////////////

 TH1D *zvtxDistrMC   = new TH1D("zvtxDistrMC  ",";zvtx;nEvents", nZvtxDistrBins, zVtxDistrMin, zVtxDistrMax);
 TH1D *zvtxDistrDATA = new TH1D("zvtxDistrDATA",";zvtx;nEvents", nZvtxDistrBins, zVtxDistrMin, zVtxDistrMax);
 TH1D *ratiozvtx;

 // ptbins
 const int nptlogmax = 15;

 const double l10 = TMath::Log(10);
 double ptbinslog[nCorrTyp][nptlogmax];

 double dptlog[nCorrTyp];
 
 for(int i = 0; i < nCorrTyp; i++)
 {
 	dptlog[i] = ( TMath::Log(ptmaxlog[i]) - TMath::Log(ptminlog[i]) )/nptlog[i]/l10;
 }


/////////////////////////////////////
//                                 //
// ====== zvtx distribution ====== //
//                                 //
/////////////////////////////////////

 ////////////////////////
 // === DATA z-vtx === //
 ////////////////////////

 ///////////////
 // Open file //
 ///////////////
 
 TFile *fdt = NULL;
 fdt = TFile::Open(inpFilenameDATA.c_str());
 if ( fdt->IsZombie() || (fdt == NULL) ) {std::cerr << "Error opening file: " << inpFilenameDATA.c_str() << std::endl; exit(-1);}
 else {std::cout << Form("TFile %s seem to have loaded.\n", inpFilenameDATA.c_str()); }

 //////////////////
 // Initializing //
 //////////////////
 
 EvtAnalyzer EvAnaDATA;
 EvAnaDATA.setupEvtAnaTree( fdt );

 EvtSelection EvSelDATA;
 EvSelDATA.setupSkimTree_pPb( fdt, false);

 ////////////////
 // Event loop //
 ////////////////
 
 log.wr( "Starting to process DATA for zvtx distribution." );

 if (nEvMaxDATA == -1)
 { nEvMaxDATA = EvSelDATA.SkimAna->GetEntries(); }

 log.wr( Form("nEvMaxDATA: %d", nEvMaxDATA) );

 for (int iEv = 0; iEv < nEvMaxDATA; iEv++)
 {

	// EventCounter
	log.EventCounter(iEv);

	// EventSelection
	if ( !EvSelDATA.isGoodEv_pPb( iEv ) ) continue;
	
	EvAnaDATA.GetEntry( iEv );
   zvtxDistrDATA->Fill( EvAnaDATA.getvz() );

 } 

 log.wr( "DATA processed." );
 log.wr( "zVtx distribution." );
 fdt->Close();
 delete fdt;

 //////////////////////
 // === MC z-vtx === //
 //////////////////////
 
  ///////////////
 // Open file //
 ///////////////
 
 TFile *fmc = NULL;
 fmc = TFile::Open(inpFilenameMC.c_str());
 if ( (fmc->IsZombie()) || (fmc == NULL) ) {std::cerr << "Error opening file: " << inpFilenameMC << std::endl; exit(-1);}
 else {std::cout << Form("TFile %s seem to have loaded.\n", inpFilenameMC.c_str()); }

 //////////////////
 // Initializing //
 //////////////////
 
 EvtAnalyzer EvAnaMC;
 EvAnaMC.setupEvtAnaTree( fmc );

 EvtSelection EvSelMC;
 EvSelMC.setupSkimTree_pPb( fmc, true);

 ////////////////
 // Event loop //
 ////////////////
 
 log.wr( "Starting to process MC for zvtx distribution." );
 
 if (nEvMaxMC == -1)
 { nEvMaxMC = EvSelMC.SkimAna->GetEntries(); }

 log.wr( Form("nEvMaxMC: %d", nEvMaxMC) );


 for (int iEv = 0; iEv < nEvMaxMC; iEv++)
 {

	// Event counter info
	log.EventCounter(iEv);

	// EventSelection
	if ( !EvSelMC.isGoodEv_pPb( iEv ) ) continue;
	
	EvAnaMC.GetEntry( iEv );
   zvtxDistrMC->Fill( EvAnaMC.getvz() );

 }

 ////////////////////
 // Setting output //
 ////////////////////

 TFile *output = new TFile(Form("./trkCorrections_%s.root", tag.c_str() ),"RECREATE");
 output->cd();

 ratiozvtx = (TH1D*)zvtxDistrDATA->Clone("zvtxratio");
 ratiozvtx->Divide( zvtxDistrMC );

//////////////////////////////////////
//                                  //
// ====== MC eff calculation ====== //
//                                  //
//////////////////////////////////////

 ///////////////
 // Open file //
 ///////////////

 TTree *trackTree = (TTree*)fmc->Get("ppTrack/trackTree");
 Tracks_c tTracks;
 setupTrackTree_c(trackTree, tTracks, true);

 //////////////////////////////
 //                          //
 // ***** Initializing ***** //
 //                          //
 //////////////////////////////
 
 std::cout << "Initializing..." << std::endl;

  TH3D **hmatched3D   = Setup_TH3D_nCorrTyp("hmatched3D",   ";p_{T};#eta;#phi", npt,ptMin,ptMax,neta,etaMin,etaMax,nphi,phiMin,phiMax);
  TH3D **hgen3D       = Setup_TH3D_nCorrTyp("hgen3D",       ";p_{T};#eta;#phi", npt,ptMin,ptMax,neta,etaMin,etaMax,nphi,phiMin,phiMax);
  TH3D *heff3D[nCorrTyp];

  std::cout << "hmatched3D, hgen3D, heff3D initialized." << std::endl;

  TH2D **hmatched2D   = Setup_TH2D_nCorrTyp("hmatched2D",   ";p_{T};#eta", npt,ptMin,ptMax,neta,etaMin,etaMax);
  TH2D **hgen2D       = Setup_TH2D_nCorrTyp("hgen2D",       ";p_{T};#eta", npt,ptMin,ptMax,neta,etaMin,etaMax);
  TH2D **hfake2D      = Setup_TH2D_nCorrTyp("hfake2D",      ";p_{T};#eta", npt,ptMin,ptMax,neta,etaMin,etaMax);
  TH2D **hreco2D      = Setup_TH2D_nCorrTyp("hreco2D",      ";p_{T};#eta", npt,ptMin,ptMax,neta,etaMin,etaMax);
  TH2D **hsecondary2D = Setup_TH2D_nCorrTyp("hsecondary2D", ";p_{T};#eta", npt,ptMin,ptMax,neta,etaMin,etaMax);
  TH2D **hreal2D      = Setup_TH2D_nCorrTyp("hreal2D",      ";p_{T};#eta", npt,ptMin,ptMax,neta,etaMin,etaMax);
  TH2D **hmultrec2D   = Setup_TH2D_nCorrTyp("hmultrec2D",   ";p_{T};#eta", npt,ptMin,ptMax,neta,etaMin,etaMax);
  TH2D *hmultrecrate2D[nCorrTyp];

  std::cout << "hmatched2D, hgen2D, hfake2D, hreco2D, hsecondary2D, hreco2D, hmultrecrate2D, hmultrecrate2D initialized." << std::endl;

  TH1D **hmatched1D   = Setup_TH1D_nCorrTyp("hmatched1D",   ";#eta", neta,etaMin,etaMax);
  TH1D **hgen1D       = Setup_TH1D_nCorrTyp("hgen1D",       ";#eta", neta,etaMin,etaMax);
  TH1D **hfake1D      = Setup_TH1D_nCorrTyp("hfake1D",      ";#eta", neta,etaMin,etaMax);
  TH1D **hreco1D      = Setup_TH1D_nCorrTyp("hreco1D",      ";#eta", neta,etaMin,etaMax);
  TH1D **hsecondary1D = Setup_TH1D_nCorrTyp("hsecondary1D", ";#eta", neta,etaMin,etaMax);
  TH1D **hreal1D      = Setup_TH1D_nCorrTyp("hreal1D",      ";#eta", neta,etaMin,etaMax);
  TH1D **hmultrec1D   = Setup_TH1D_nCorrTyp("hmultrec1D",   ";#eta", neta,etaMin,etaMax);
  TH1D *hmultrecrate1D[nCorrTyp];

  std::cout << "hmatched1D, hgen1D, hfakse1D, hreco1D, hsecondary1D, hreal1D, hmultrec1D, hmultrecrate1D initialized." << std::endl;

  TH3D **hcorr3D      = Setup_TH3D_nCorrTyp("hcorr3D",      ";p_{T};#eta;#phi", npt,ptMin,ptMax,neta,etaMin,etaMax,nphi,phiMin,phiMax);


  std::cout << "hcorr3D initialized." << std::endl;

 ///////////////////////////////////////
 //                                   //
 // ***** Calculate corrections ***** //
 //                                   //
 ///////////////////////////////////////
 
 log.wr( "Starting to process MC to calclulate track corrrections..." );

 TH1D *hist1 = new TH1D("hist1", ";eta;", 20, -1.0, 1.0);
 TH1D *hist2 = new TH1D("hist2", ";eta;", 20, -1.0, 1.0);
 TH1D *hist3 = new TH1D("hist3", ";eta;", 20, -1.0, 1.0);

 TH2D* dEdxvsplinlinall1 = new TH2D ("dEdxVsP lin-lin1 " ,";p(GeV/c);dE/dx [MeV/cm]", npBins, pminlin, pmaxlin, ndEdxBins, dEdxminlin, dEdxmaxlin);
 TH2D* dEdxvsplinlinall2 = new TH2D ("dEdxVsP lin-lin2 " ,";p(GeV/c);dE/dx [MeV/cm]", npBins, pminlin, pmaxlin, ndEdxBins, dEdxminlin, dEdxmaxlin);
 TH2D* dEdxvsplinlinall3 = new TH2D ("dEdxVsP lin-lin3 " ,";p(GeV/c);dE/dx [MeV/cm]", npBins, pminlin, pmaxlin, ndEdxBins, dEdxminlin, dEdxmaxlin);

 TH2D* dEdxvsEta1 = new TH2D ("dEdxVsEta1" ,";#eta;dE/dx [MeV/cm]", 40, -2.0, 2.0, ndEdxBins, dEdxminlin, 8.0);
 TH2D* dEdxvsEta2 = new TH2D ("dEdxVsEta2" ,";#eta;dE/dx [MeV/cm]", 40, -2.0, 2.0, ndEdxBins, dEdxminlin, 8.0);
 TH2D* dEdxvsEta3 = new TH2D ("dEdxVsEta3" ,";#eta;dE/dx [MeV/cm]", 40, -2.0, 2.0, ndEdxBins, dEdxminlin, 8.0);

 double ptlim1 = 0.3;
 double ptlim2 = 0.4;

  // Event loop 
  if ( nEvMaxMC == -1 )
  { nEvMaxMC = trackTree->GetEntries(); }

  log.wr( Form("nEvMaxMC: %d", nEvMaxMC));

  for (int iEvA = 0; iEvA < nEvMaxMC; iEvA++)
  {

		// EventCounter
		log.EventCounter(iEvA);
		
		// Event Selection
		if ( !EvSelMC.isGoodEv_pPb( iEvA ) ) continue;;
		
		// Event zvtx
		EvAnaMC.GetEntry( iEvA );
		float vz = EvAnaMC.getvz();
		int vzB = ratiozvtx->FindBin(vz) ;
		int wzvtx = ratiozvtx->GetBinContent(vzB);

		if ( (vzB == -1) || ((nZvtxDistrBins) == vzB)) continue;
	
		// Tracks & particles
		trackTree->GetEntry(iEvA);
		
		int nTrk = tTracks.nTrk;
		
		// Track loop
		for (int iTrk = 0; iTrk < nTrk; iTrk++)
		{

			// *** Track selection *** //
			if ( !TrackSelection(tTracks, iTrk ) ) continue;

			double pt  = tTracks.trkPt [iTrk];
			double eta = tTracks.trkEta[iTrk]; 
			double phi = tTracks.trkPhi[iTrk];
			float p    = pt * cosh(eta);
			int PID 	= pidutil->GetID( tTracks, iTrk);
			bool isPID = (PID != 99);
			
			bool isInsideChrPt = ( ptbin(  0, pt) != -1 );
			bool isInsidePIDPt = ( ptbin(PID, pt) != -1 );

			// Fake tracks
			if( tTracks.trkFake[iTrk] )
			{ 

				           hfake2D[ 0 ]->Fill(pt,eta, wzvtx);
							  if (isInsideChrPt)
							  { hfake1D[ 0 ]->Fill( eta, wzvtx);}
				if (isPID) {hfake2D[PID]->Fill(pt,eta, wzvtx);
						  		if (isInsidePIDPt)
								{ hfake1D[PID]->Fill(eta, wzvtx);}};
			}
			else 
			// Real track
			{ 
			             hreal2D[ 0 ]->Fill(pt,eta, wzvtx);
							 if (isInsideChrPt)
							 {hreal1D[ 0 ]->Fill(   eta, wzvtx);}

			  if (isPID) {hreal2D[PID]->Fill(pt,eta, wzvtx);

						  	  if (isInsidePIDPt)
							  {hreal1D[PID]->Fill(   eta, wzvtx);}}
			  // Real secondary track
			  if( tTracks.trkStatus[iTrk] < 0 )
			  {
			               hsecondary2D[ 0 ]->Fill(pt,eta,wzvtx);

							   if (isInsideChrPt)
								{hsecondary1D[ 0 ]->Fill(   eta,wzvtx);}
			    if (isPID) { hsecondary2D[PID]->Fill(pt,eta,wzvtx);

						  	  		if (isInsidePIDPt)
									{hsecondary1D[PID]->Fill(   eta,wzvtx);}}
			  };
		   }
		
			// Reconstructed tracks
		               hreco2D[ 0 ]->Fill(pt,eta,wzvtx);
							if (isInsideChrPt)
		               {hreco1D[ 0 ]->Fill(   eta,wzvtx);}
			if (isPID) {hreco2D[PID]->Fill(pt,eta,wzvtx);
						  	if (isInsidePIDPt)
							{hreco1D[PID]->Fill(   eta,wzvtx);}}
		}
		
		// === Particle loop === //
		for(int iPart = 0; iPart < tTracks.nParticle; iPart++)
		{

			double pt  = tTracks.pPt [iPart];
			double eta = tTracks.pEta[iPart];
			double phi = tTracks.pPhi[iPart];
			double p   = pt * cosh(eta);
			double dEdx = tTracks.mtrkdedx[iPart];

			int PID 	= pidutil->GetIDgenPart_trkCorr( tTracks, iPart);
			bool isPID = (PID != 99);

			int    mPID = pidutil->GetIDmTrk_trkCorr(tTracks, iPart);
			bool   PID_match = (PID == mPID);

		   bool isInsideChrPt = ( ptbin(  0, pt) != -1 );
			bool isInsidePIDPt = ( ptbin(PID, pt) != -1 );

			// Particle selection
			if( fabs(eta) > 2.4 ) continue;
			if( pt < ptGlobalMin ) continue;
			//
		
								hgen3D[ 0 ]->Fill( pt,eta,phi,wzvtx );
								hgen2D[ 0 ]->Fill( pt,eta,    wzvtx );

						  	  	if (isInsideChrPt)
								{hgen1D[ 0 ]->Fill(    eta,    wzvtx );}


			if ( isPID )
			{

				// Debuggg 
			//	if ( 0.6 < pt)
			//	{
			//	std::cerr << Form("PID: %d p: %.2f eta: %.2f", PID, p, eta) << std::endl;
			//	}

				hgen3D[PID]->Fill( pt,eta,phi,wzvtx ); 
				hgen2D[PID]->Fill( pt,eta    ,wzvtx );

				if (isInsidePIDPt)
				{hgen1D[PID]->Fill(    eta,    wzvtx );}

				if ( (PID == 2) && (ptlim1 < pt) && (pt < ptlim2))
				{ 
					hist1->Fill(eta);
					dEdxvsplinlinall1->Fill(p,dEdx);
					dEdxvsEta1->Fill(eta,dEdx);
				}

				if ( (PID == 2) && (ptlim1 < pt) && (pt < ptlim2) && !PID_match )
				{ 
					hist2->Fill(eta);
					dEdxvsplinlinall2->Fill(p,dEdx);
					dEdxvsEta2->Fill(eta,dEdx);
				}

					if ( PID_match )
					{ 
						if ( (PID == 2) && (ptlim1 < pt) && (pt < ptlim2) )
						{ 
							hist3->Fill(eta); 
							dEdxvsplinlinall3->Fill(p,dEdx);
							dEdxvsEta3->Fill(eta,dEdx);
						}

						hmatched3D[PID]->Fill(pt,eta,phi,wzvtx); 
						hmatched2D[PID]->Fill(pt,eta,    wzvtx); 
						hmultrec2D[PID]->Fill(pt,eta,wzvtx*tTracks.pNRec[iPart]);

						if (isInsidePIDPt)
						{
							hmatched1D[PID]->Fill(   eta,    wzvtx); 
							hmultrec1D[PID]->Fill(   eta,wzvtx*tTracks.pNRec[iPart]);
						}
		
					}

			};

			// matched track selection
			// WARNING
			// mTrackSelection is OFF
//			if( !( mTrackSelection_c( tTracks, iPart) )) continue;


				hmatched3D[ 0 ]->Fill(pt,eta,phi,wzvtx); 
				hmatched2D[ 0 ]->Fill(pt,eta,    wzvtx); 
				hmultrec2D[ 0 ]->Fill(pt,eta,wzvtx*tTracks.pNRec[iPart]);
				if (isInsideChrPt)
				{		  
				hmultrec1D[ 0 ]->Fill(   eta,wzvtx*tTracks.pNRec[iPart]);
				hmatched1D[ 0 ]->Fill(   eta,    wzvtx);
				}


		   }
  }

  for( int iPar = 0; iPar < nCorrTyp; iPar++)
  {
   	hfake2D     [iPar] -> Divide( hreco2D   [iPar] );
   	hfake1D     [iPar] -> Divide( hreco1D   [iPar] );
   	hsecondary2D[iPar] -> Divide( hreal2D   [iPar] );
   	hsecondary1D[iPar] -> Divide( hreal1D   [iPar] );
   	heff3D      [iPar]  = (TH3D*)hmatched3D [iPar]->Clone(Form("heff3D part %d", iPar));
   	heff3D      [iPar] -> Divide( hgen3D    [iPar] );

   	hmultrecrate2D [iPar] = (TH2D*)hmultrec2D [iPar]->Clone(Form("hmultrecrate2D part %d", iPar));
   	hmultrecrate1D [iPar] = (TH1D*)hmultrec1D [iPar]->Clone(Form("hmultrecrate1D part %d", iPar));

   	hmultrecrate2D [iPar] -> Divide( hmatched2D[iPar] );
   	hmultrecrate1D [iPar] -> Divide( hmatched1D[iPar] );
  }

  // Filling track correction table
  double ptbw[4];
  for(int i = 0; i < nCorrTyp; i++)
  {
	ptbw[i] = (ptMax[i]-ptMin[i])/npt[i];
  }

  double etabw[4]; 
  for(int i = 0; i < nCorrTyp; i++)
  { etabw[i] = (etaMax[i]-etaMin[i])/neta[i];}

  const double phibw = (phiMax-phiMin)/nphi;

  for(int i = 0; i < nCorrTyp; i++)
  for(int x = 1; x < npt[i]  +1; x++)
  for(int y = 1; y < neta[i] +1; y++)
  for(int z = 1; z < nphi +1; z++)
  {
	  double pt  = ptMin[i]+x*ptbw[i]; 
	  double eta = etaMin[i]+y*etabw[i]; 
	  double phi = phiMin+z*phibw; 

     double value = (1.0-hfake2D[i]->GetBinContent(x,y))*(1.0-hsecondary2D[i]->GetBinContent(x,y)) / ((heff3D[i]->GetBinContent(x,y,z)) * (hmultrecrate2D[i]->GetBinContent(x,y)));
	  hcorr3D[i]->SetBinContent(x,y,z,value);
 	  log.wr(Form("%d %.3f %.2f %.2f : %.4f", i, pt, eta, phi, value));
  }

 //////////////////////
 //                  //
 // **** OUTPUT **** //
 //                  //
 //////////////////////
  
  output->cd();
  output->Write();
  
  //heff->wr();
  //hfake->wr();
  //hmultrec->wr();
  //hsecondary->wr();
  //
  
  log.Close();
  
  output->Close();
  delete output;

  printf("Done.\n");
}
