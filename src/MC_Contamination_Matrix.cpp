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
const double zVtxDistrMin = -13;
const double zVtxDistrMax =  13;

const int npt = 14;

const double ptMin = 0.2;
const double ptMax = 1.6;
const double ptbw  = 0.1;

int GetPtBin_sl (double pt)
{
   if ( (pt < ptMin) || (ptMax < pt))
	{return -1;}
	else
	{return floor( (pt-ptMin)/ptbw );}
}

// const int npt_regions = 3;
// const int npt[npt_regions] = {7, 1, 6};
// const int ptMin[npt_regions] = {0.2, 0.9, 1.0};
// const int ptMax[npt_regions] = {0.9, 1.0, 1.6};
// const double ptbw[npt_regions] = {0.1, 0.1, 0.1};
// 
// const int nmxBins[npt_regions]  = {  4,   3,  3 };
// const double mxMin[npt_regions] = {0.0, 0.0, 0.0};
// const double mxMax[npt_regions] = {4.0, 3.0, 3.0};

// int GetPtBin_sl (double pt)
// {
//    if ( (pt < ptMin) || (ptMax < pt))
// 	{return -1;}
// 	else
// 	{return floor( (pt-ptMin)/ptbw );}
// }
// 
// int GetPtRegionBin ( double pt )
// {  
// 	if ( (pt < ptMin[0]) || (ptMax[npt_regions-1] < pt))
// 	{return -1;}
// 
// 	for (int i = 0; i < npt_regions, i++)
// 	{
// 		if ( (ptMin[i] < pt) && (pt < ptMax[i]) )
// 		{ return i; }
// 	}
// }
// 
// int GetPtBin ( int PtRegionBin, double pt )
// {   return ( floor( (pt-ptMin[PtRegionBin])/ptbw[PtRegionBin]) );  }


int main( int argc, const char *argv[] )
{

  if(argc != 7)
  {
    std::cerr << "Usage: MC_PID_Contamination_matrix <MC sample> <DATA sample> <tag> <nEventsDATA> <nEventsMC> <PIDconfig>" << std::endl;
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
 int nCorrTyp = nCorrTyp_; 
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
 
 TFile *fdt = TFile::Open(inpFilenameDATA.c_str());
 if ( fdt->IsZombie() ) {std::cerr << "Error opening file: " << inpFilenameDATA << std::endl; exit(-1);}

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
 
 TFile *fmc = TFile::Open(inpFilenameMC.c_str());
 if ( fmc->IsZombie() ) {std::cerr << "Error opening file: " << inpFilenameMC << std::endl; exit(-1);}

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

 TFile *output = new TFile(Form("./PID_Contamination_matrix_%s.root", tag.c_str() ),"RECREATE");
 output->cd();

 ratiozvtx = (TH1D*)zvtxDistrDATA->Clone("zvtxratio");
 ratiozvtx->Divide( zvtxDistrMC );

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

 TH2D *ptres = new TH2D("ptres",";p_{T,GEN};p_{T,RECO}", 40, 0.0, 2.0, 40, 0.0, 2.0 );

 // dEdxvsPMaps
 for( int ptBin = 0; ptBin < npt; ptBin++ )
 {
   double pt1 = (ptMin + (ptBin  ) * ptbw);
	double pt2 = (ptMin + (ptBin+1) * ptbw);
	dEdxvsPMapsLin[ptBin] = new TH2D(Form("dEdxvsPMapsLin_pt_%.2f-%.2f", pt1, pt2), ";p [GeV/c];dE/dx [MeV/cm]", npBins, pminlin, pmaxlin, ndEdxBins, dEdxminlin, dEdxmaxlin);
	dEdxvsPMapsLog[ptBin] = new TH2D(Form("dEdxvsPMapsLog_pt_%.2f-%.2f", pt1, pt2), ";p [GeV/c];dE/dx [MeV/cm]", npBinslog, pBins, ndEdxBinslog, dEdxBins);
 }

 // Contamination Matrices
// TH2D ***matrix;
//
// matrix = new TH2D**[npt_regions];
//
// for( int ptRegion = 0; ptRegion < npt_regions;   ptRegion++ )
// for( int ptBin    = 0; ptBin    < npt[ptRegion]; ptBin++   )
// {
//   double pt1 = (ptMin[ptRegion] + (ptBin  ) * ptbw[ptRegion]);
//	double pt2 = (ptMin[ptRegion] + (ptBin+1) * ptbw[ptRegion]);
//
//	matrix[ptRegion][ptBin] = new TH2D (Form("contmatrix_pt_%.2f-%.2f", pt1, pt2),";RECO;MC", nmxBins[ptRegion], mxMin[ptRegion], mxMax[ptRegion], nmxBins[ptRegion], mxMin[ptRegion], mxMax[ptRegion]);
// }

 TH2D *matrix[npt];

 for( int ptBin = 0; ptBin < npt; ptBin++ )
 {
   double pt1 = (ptMin + (ptBin  ) * ptbw);
	double pt2 = (ptMin + (ptBin+1) * ptbw);

	matrix[ptBin]      = new TH2D(Form("contmatrix_pt_%.2f-%.2f", pt1, pt2),";RECO;MC", 4, 0.0, 4.0, 4, 0.0, 4.0);

 }


 ///////////////////////////////////////
 //                                   //
 // ***** Calculate corrections ***** //
 //                                   //
 ///////////////////////////////////////
 
 log.wr( "Starting to process MC to calclulate contamination matrix..." );

  // Event loop 
  if ( nEvMaxMC == -1 )
  { nEvMaxMC = trackTree->GetEntries(); }

  log.wr( Form("nEvMaxMC: %d", nEvMaxMC));

  for (int iEvA = 0; iEvA < nEvMaxMC; iEvA++)
  {

		// EventCounter
		log.EventCounter(iEvA);
		
		// Event Selection
		if ( !EvSelMC.isGoodEv_pPb( iEvA ) ) continue;
		
		// Event zvtx
		EvAnaMC.GetEntry( iEvA );
		float vz = EvAnaMC.getvz();
		int vzB = ratiozvtx->FindBin(vz) ;
		int wzvtx = ratiozvtx->GetBinContent(vzB);

		if ( (vzB == -1) || ((nZvtxDistrBins) == vzB)) continue;
	
		// Tracks & particles
		trackTree->GetEntry(iEvA);
		int nTrk = tTracks.nTrk; 
		// === Particle loop === //
		for(int iPart = 0; iPart < tTracks.nParticle; iPart++)
		{

			// matched track selection
			// particle selection
			if( !( mTrackSelection_c( tTracks, iPart) )) continue;

			double mpt = tTracks.mtrkPt [iPart];
			double eta = tTracks.pEta[iPart]; 
			float  mp  = mpt * cosh( eta);
			double pt  = tTracks.pPt [iPart];
			float dEdx = tTracks.mtrkdedx[iPart];

			//int ptBin_RECO_sl = GetPtBin_sl( mpt );
			//int ptBin_RECO_sl = GetPtBin( PtRegionBin_RECO, mpt );

			ptres->Fill(pt, mpt);

			double PID_RECO = pidutil->GetID_cm( tTracks, iPart );
			if (   PID_RECO == 99 ) continue; // un-ID or 0.8 < |eta| particle

			int  ptBin_RECO =  ptbin( floor(PID_RECO+1.0), mpt);
			if ( ptBin_RECO == -1 ) continue; // particle outside pt range

			double PID_GENE = McPID2AnaPID_cm  ( tTracks, iPart );

			int ptBin_GENE_sl = GetPtBin_sl(  pt );
			int ptBin_RECO_sl = GetPtBin_sl( mpt );

			if ( (ptBin_GENE_sl == -1) || (ptBin_RECO_sl == -1) ) continue;


//			std::cerr << "PID_GEN: " << PID_GEN << std::endl;
//			std::cerr << "PID_RECO: " << PID_RECO << std::endl;

			dEdxvsPMapsLin[ptBin_RECO_sl]->Fill(mp, dEdx);
			dEdxvsPMapsLog[ptBin_RECO_sl]->Fill(mp, dEdx);

			matrix[ptBin_RECO_sl]->Fill(PID_RECO, PID_GENE);
		
		 }
		
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
