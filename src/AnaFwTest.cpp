#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrix.h>
#include <deque>
#include <vector>
#include <TNtuple.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TVector.h>
#include "AnalysisFW.h"
#include "TGraphErrors.h"
#include "PIDUtils.h"
#include "CorrelationUtils.h"
#include "ContMatrix.h"


//void SetupCorr( TH2D ***&correl_currev_signal )
//{
//
// int nCorrTyp = nCorrTyp_;
//
// correl_currev_signal = new TH2D**[nCorrTyp];
//
// //int nPtBins = nPtBins_;
//
// for( int TypBin=0; TypBin < nCorrTyp; TypBin++)
// {
// 	correl_currev_signal[TypBin] = new TH2D*[nPtBins];
// }
//
// for( int TypBin=0; TypBin < nCorrTyp; TypBin++)
// for( int ptBin=0 ;  ptBin < nPtBins ; ptBin++)
// {
//	double pt1 = pt(ptBin, 0, nPtBins);
//	double pt2 = pt(ptBin, 1, nPtBins);
//
// 	correl_currev_signal[TypBin][ptBin] = new TH2D( Form("correl_currev_signal_typ_%1d_pt_%.2f-%.2f", TypBin, pt1, pt2),
//															"2D Correlation function;#Delta #eta; #Delta #Phi",
//	                                          ndEtaBins,dEtaMin,dEtaMax,ndPhiBins,dPhiMin,dPhiMax);
//	
// 
//   
// }
//}
//
//void charargtest( const char arg[] )
//{
//	std::cerr << arg << std::endl;
//};
//
//class AnaConfig {
//public :
//   AnaConfig(){};
//   ~AnaConfig(){};
//
//	// Variables
//
//	// Functions
//};
//
//AnaConfig::AnaConfig()
//{
//}


int main( int argc, const char *argv[] )
{ 

//	std::cerr << "hello" << std::endl;
//	std::deque<int> deku[3];
//	deku[0].push_back(1);
//	deku[0].push_back(2);
//	deku[0].push_back(3);
//	deku[1].push_back(4);
//	deku[1].push_back(4);
//	deku[1].push_back(4);
//	deku[2].push_back(5);
//	deku[2].push_back(4);
//	deku[2].push_back(4);
//
//	// Debuggg 
//	std::cerr << "deku[0][1]: " << deku[0][1] << std::endl;
//	std::cerr << "deku[0][2]: " << deku[0][2] << std::endl;
//
//	track trk;
//	trk.pid = 2;
//	trk.charge = 1;
//	trk.pt = 1.323;
//	trk.phi = 0.231;
//	trk.eta = 1.21;
//
//	EventData *ev;
//	ev = new EventData;
//
//	ev->AddTrack(trk);
//
//	std::deque< EventData > events[15][10];
//	events[4][5].push_back( (*ev));
//
//	std::cerr << "pid" << events[4][5][0].tracks[0].pid << std::endl; 
//
//	// Debuggg 
//	std::cerr << "ArrayTest" << std::endl;
//
//	double array[3];
//
//	// Debuggg 
//	std::cerr << "array: " << array << std::endl;


 //TFile *output = new TFile("./AnaFwTest.root" ,"RECREATE");


 //TH2D ***correl_currev_signal;
 // Debuggg 
 //std::cerr << "cucc " << std::endl;
 //SetupCorr(correl_currev_signal);
 //correl_currev_signal[0][0]->Fill(0.1,0.1);

 //int *array = new int[3];
 //array[0] = 0;
 //array[1] = 0;
 //array[3] = 0;

 //// Debuggg 
 //std::cerr << "hel " << std::endl;
 //std::cerr << "array[0]: " << array[0] << std::endl;

 //array[0]++;

 //std::cerr << "array[0]: " << array[0] << std::endl;

 //std::deque<int>** hagyma;
 //hagyma = new std::deque<int>*[3];


 //float *x = new float;
 //double *BB_Pion_low_parr = new double[3];
 //BB_Pion_low_parr[0] = BB_Pion_low_par[0];
 //BB_Pion_low_parr[1] = BB_Pion_low_par[1];
 //BB_Pion_low_parr[2] = BB_Pion_low_par[2];

 //double arr[3];
 //arr[0] = BB_Pion_low_par[0];
 //arr[1] = BB_Pion_low_par[1];
 //arr[2] = BB_Pion_low_par[2];

 //(*x) = 2;
 //double golya = 10.;

 //// Debuggg 
 //std::cerr << "bbcurve1: " << BBcurve1c(x, BB_Pion_hig_par) << std::endl;

 //TH1D *correl_1D_back  = new TH1D("correl_1D_back",  "1D Correlation", ndPhiBins, dPhiMin, dPhiMax);
	
// int TypBin = 99;
// double pt = 2.9;
// std::cerr << Form("ptbin(%d, %f): ", TypBin ,pt ) << ptbin(TypBin, pt) << std::endl;
//
// std::vector< double >  xvalues(1); 
// std::vector< double >  yvalues(1); 
// std::vector< double > xEvalues(1); 
// std::vector< double > yEvalues(1); 
// xvalues[0]  = 0;
// yvalues[0]  = 0;
// xEvalues[0] = 0;
// yEvalues[0] = 0;
// TGraphErrors *g1 = new TGraphErrors(xvalues, yvalues, xEvalues, yEvalues);
//
//
 
// const int nCorrTyp = 4;
// int nPtBins[4] = {1, 4, 4, 7};
// int nMultiplicityBins = 4;
//
// CorrelationResults Results;
// Results.Setup(nCorrTyp, nPtBins, nMultiplicityBins);
//
// Results.V2[1][0][0] = 4.4;
// Results.V2[1][1][0] = 5.4;
// Results.V2[1][2][0] = 6.4;
// Results.V2[1][3][0] = 7.4;
//
// std::vector<double > vecc(4);
// vecc = Results.V2vec(1,0);
//
// // Debuggg 
// std::cerr << "vec: " << vecc[3]  << std::endl;
//
// int mult1 = multiplicity(3, 0, nMultiplicityBins);
// std::string str = Form("heheheh %d", mult1);
// // Debuggg 
// std::cerr << "mult1: " << str << std::endl;
//
// double var = 23423;
// // Debuggg 
// std::cerr << "var: " << var << std::endl;
// // Debuggg 
// std::cerr << Form("%03d", 1) << std::endl;
//
// int a1 = multiplicity_Ana(0, 0, 1);
// int a2 = multiplicity_Ana(0, 1, 1);
//
// // Debuggg 
// std::cerr << "multiplicit_Ana(0, 0, 1): " << a1 << std::endl;
// std::cerr << "multiplicit_Ana(0, 0, 1): " << a2 << std::endl;

// Tracks tTracks;
//
// tTracks.trkEta[0] = 2.39999;
// tTracks.highPurity[0] = true;
//
// tTracks.trkDz1[0] = 3.001;
// tTracks.trkDzError1[0] = 1.0;
//
// tTracks.trkDxy1[0] = 3.000;
// tTracks.trkDxyError1[0] = 1.0;
//
// tTracks.trkPtError[0] = 0.0999;
// tTracks.trkPt[0] = 1.00;
//
// bool result = TrackSelection(tTracks, 0);
// std::cerr << "result " << result << std::endl;

// === Histogram projection test ==== ///
//  TFile *output = new TFile("file.root","RECREATE");
//
//	const double dEtaMin = -4.8;
//	const double dEtaMax =  4.8;
//	
//	const double dEtaMin_plot = -2.8;
//	const double dEtaMax_plot =  2.8;
//	
//	const double dPhiMin = -TMath::Pi()/2;
//	const double dPhiMax = 3*TMath::Pi()/2;
//	
//	const int ndEtaBins = 96;
//	const int ndPhiBins = 96;
//
//  TH2D histo2D = TH2D("histo2D","title; #DELTA #eta; #Delta #phi",ndEtaBins,dEtaMin,dEtaMin, ndPhiBins,dPhiMin, dPhiMax);
//  TH1D histo1D = TH1D("histo1D","title; #DELTA #phi; y axis",ndPhiBins, dPhiMin, dPhiMax);
//
//  histo2D.Fill(3.95,0.);
//
//  int bin1 = 87;
//  int bin2 = 88;
//
//  histo2D.ProjectionY("histo1D",bin1,bin2);
//
//  histo2D.Draw();
//  histo1D.Draw();
//
//  charargtest("test");
//
//
// output->Write();
// output->Close();
	

//	std::cout << "zvtxbin: " << zvtxbin(-12.5, nZvtxBins_) << std::endl;
//	std::cout << "zvtxbin: " << zvtxbin(-13.5, nZvtxBins_) << std::endl;
//	std::cout << "zvtxbin: " << zvtxbin(12.5, nZvtxBins_) << std::endl;

// std::string configfile = "./PIDUtils/config/config_default";
// PIDUtil * pidutil = new PIDUtil;
// pidutil->ReadInConfig( configfile );

//TH2D *matrix = new TH2D("matrix",";RECO;GEN", 4, 0.0, 4.0, 4, 0.0, 4.0 );
//TH2D *matrix2;
//
//matrix->Fill(0.5,0.5);
//matrix->Fill(0.5,0.5);
//matrix->Fill(0.5,0.5);
//
//matrix->Fill(0.5,1.5);
//matrix->Fill(0.5,1.5);
//matrix->Fill(0.5,1.5);
//
//matrix->Fill(0.5,2.5);
//matrix->Fill(0.5,2.5);
//matrix->Fill(0.5,2.5);
//matrix->Fill(0.5,2.5);
//matrix->Fill(0.5,2.5);
//matrix->Fill(0.5,2.5);
//
//matrix->Fill(1.5,0.5);
//
//matrix->Fill(2.5,0.5);
//
//matrix->Fill(3.5,0.5);
//matrix->Fill(3.5,0.5);
//
//matrix->Fill(3.5,1.5);
//matrix->Fill(3.5,1.5);
//matrix->Fill(3.5,1.5);
//
//matrix->Fill(3.5,3.5);
//
//CM::displayMatrix_TH2D(matrix);
//
//CM::CopyTH2DtoTH2D( matrix, matrix2, "new", 123);
//
//CM::displayMatrix_TH2D(matrix2);
//
//CM::normalizeColoumn_TH2D(matrix2, 1);
//
//CM::displayMatrix_TH2D(matrix2);
//
//CM::plotContMatrix(matrix2, 234, "lolmatrix" );


	std::string filename = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/preprocessed/EPOS_RECO_level_trkCorr_no_nEv_-1/correl_analysis_0.root";

	TFile *f = new TFile( filename.c_str(), "READ");
	if ( f->IsZombie() ) {std::cerr << "Error opening file: " << filename << std::endl; }
	else{ std::cout << Form("TFile %s seems to be loaded.", filename.c_str()) << std::endl; };

	TH2D *corr1 = (TH2D*)f->Get("correl2D_signal_Typ_1_pt_0.20-0.30_nTrk_000-120");
	TH2D *corr2 = (TH2D*)f->Get("correl2D_signal_Typ_2_pt_0.20-0.30_nTrk_000-120");
	TH2D *corr3 = (TH2D*)f->Get("correl2D_signal_Typ_3_pt_0.20-0.30_nTrk_000-120");

	TH2D *corr = new TH2D( "lol", "2D Correlation function;#Delta #eta; #Delta #Phi",
		                                            ndEtaBins,dEtaMin,dEtaMax,ndPhiBins,dPhiMin,dPhiMax);
;
	corr->Add(corr1, -0.03);
	corr->Add(corr2, -0.14);
	corr->Add(corr3,  1.07);


	std::string figurename = "correlation";
	std::string title = "lol";
	std::string zaxistitle = "lol";
	std::string leftlabel  = "leftlabel";
	std::string rightlabel = "rightlabel";

	plot2DCorrelation(corr , figurename, title, zaxistitle, leftlabel, rightlabel);

}
