#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <deque>
#include "AnalysisFW.h"
#include "AnalysisBinning.h"
#include "PIDUtils.h"
#include "SetupCustomTrackTree.h"
#include "TLatex.h"
#include "GraphTools.h"
#include "ContMatrix.h"


const int    npt   = 14;
const double ptMin = 0.2;
const double ptMax = 1.6;
const double ptbw  = 0.1;

int main( int argc, const char *argv[] )
{ 


 if(argc != 3)
 {
   std::cerr << "Usage: MC_Contamination_viewer <trkCorr.root file to be displayed> <PIDConfig>" << std::endl;
   exit(1);
 }

 TString inpFilename     = argv[1];
 std::string PIDconfig	 = argv[2];

 // Open file
 TFile *f = TFile::Open(inpFilename);
 if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << "File successfully opened." << std::endl;}

 ////////////////////////
 // Read In histograms //
 ////////////////////////
 
 TH2D *dEdxvsPMapsLin[npt];
 TH2D *dEdxvsPMapsLog[npt];

 TH2D *matrix[npt];

 TH2D *ptres = (TH2D*)f->Get("ptres");

 for( int ptBin = 0; ptBin < npt; ptBin++ )
 {

   double pt1 = (ptMin + (ptBin  ) * ptbw);
	double pt2 = (ptMin + (ptBin+1) * ptbw);

	matrix[ptBin] = (TH2D*)f->Get(Form("contmatrix_pt_%.2f-%.2f", pt1, pt2));

	dEdxvsPMapsLin[ptBin] = (TH2D*)f->Get(Form("dEdxvsPMapsLin_pt_%.2f-%.2f", pt1, pt2) );
	dEdxvsPMapsLog[ptBin] = (TH2D*)f->Get(Form("dEdxvsPMapsLog_pt_%.2f-%.2f", pt1, pt2) );

 }

 ////////////////////////////////
 //                            //
 // ***** Setting output ***** //
 //                            //
 ////////////////////////////////
 
 TFile *output = new TFile( "./Cont_Matrix.root", "RECREATE");
 output->cd();

 /////////////////////////
 //                     //
 // ***** Prepare ***** //
 //                     //
 /////////////////////////

 TH2D *cont_matrix_nor_TH2D[CM::nPt];
 TH2D *cont_matrix_inv_TH2D[CM::nPt];
 TMatrix *cont_matrix_nor_TMatrix[CM::nPt];
 TMatrix *cont_matrix_inv_TMatrix[CM::nPt];

 for( int RgnBin = 0; RgnBin < CM::nRegions; RgnBin++ )
 for( int i = 0; i < CM::nPtBins[RgnBin]; i++ )
 {

		int ptBin = CM::PtBins[1][RgnBin][i];

   	double pt1 = (CM::ptMin + (ptBin  ) * CM::PtBw);
		double pt2 = (CM::ptMin + (ptBin+1) * CM::PtBw);

		std::string matrixname = Form("cont_matrix_pt_%.2f-%.2f", pt1, pt2);

		CM::CopyTH2DtoTH2D   ( matrix[ptBin], cont_matrix_nor_TH2D[ptBin],    matrixname.c_str(), CM::Indices[RgnBin] );
		CM::normalizeMatrix_TH2D( cont_matrix_nor_TH2D[ptBin] );

		CM::CopyTH2DtoTMatrix( cont_matrix_nor_TH2D[ptBin], cont_matrix_nor_TMatrix[ptBin],   CM::Indices_red[RgnBin] );
		cont_matrix_inv_TMatrix[ptBin] = (TMatrix*)cont_matrix_nor_TMatrix[ptBin]->Clone();
		cont_matrix_inv_TMatrix[ptBin]->Invert();

		CM::plotContMatrix( cont_matrix_nor_TH2D[ptBin], CM::Indices[RgnBin], matrixname.c_str() );

		CM::displayMatrix_TH2D( matrix[ptBin] );
		CM::displayMatrix_TH2D( cont_matrix_nor_TH2D[ptBin] );
		CM::displayMatrix_TMatrix( cont_matrix_nor_TMatrix[ptBin] );
		CM::displayMatrix_TMatrix( cont_matrix_inv_TMatrix[ptBin] );
 }


 TH2D *th2D_matrix[CM::nPt];

 ////////////////////////

// for( int ptBin = 0; ptBin < npt; ptBin++ )
// {
//   double pt1 = (ptMin + (ptBin  ) * ptbw);
//	double pt2 = (ptMin + (ptBin+1) * ptbw);
//
//	TCanvas canvas_matrix ("contmatrix", ";RECO;GEN", 800, 600);
//
//	gStyle->SetOptStat(0);
//	gStyle->SetPaintTextFormat("2.2f");
//	matrix[ptBin]->SetMarkerSize(2);
//
//	matrix[ptBin]->Draw("TEXT");
//
// 	matrix[ptBin]->GetXaxis()->SetRange(1, 3);
// 	matrix[ptBin]->GetYaxis()->SetRange(1, 3);
//
// 	matrix[ptBin]->GetXaxis()->SetLabelSize(0.05);
// 	matrix[ptBin]->GetYaxis()->SetLabelSize(0.05);
//
//	matrix[ptBin]->GetXaxis()->SetBinLabel(1, "#pi");
//	matrix[ptBin]->GetXaxis()->SetBinLabel(2, "K");
//	matrix[ptBin]->GetXaxis()->SetBinLabel(3, "p");
//	matrix[ptBin]->GetYaxis()->SetBinLabel(1, "#pi");
//	matrix[ptBin]->GetYaxis()->SetBinLabel(2, "K");
//	matrix[ptBin]->GetYaxis()->SetBinLabel(3, "p");
//
//	std::string filebasename = Form("matrix_pt_%.2f-%.2f", pt1, pt2);
//	std::string outPDF = filebasename+".pdf"; 
//	std::string outPNG = filebasename+".png";
//
//	canvas_matrix.SaveAs(outPDF.c_str());
//	canvas_matrix.SaveAs(outPNG.c_str());
//
// 	std::string dEdxvsplinFigBase = Form( "dEdxvspLin_typ_%.2f-%.2f", pt1, pt2 	);
// 	std::string dEdxvsplogFigBase = Form( "dEdxvspLog_typ_%.2f-%.2f", pt1, pt2 	);
//
// 	std::string dEdxvsplinFigPNG = dEdxvsplinFigBase+".png";
// 	std::string dEdxvsplinFigPDF = dEdxvsplinFigBase+".pdf";
// 	std::string dEdxvsplogFigPNG = dEdxvsplogFigBase+".png";
// 	std::string dEdxvsplogFigPDF = dEdxvsplogFigBase+".pdf";
//
//	makedEdxvspFigloglog( dEdxvsPMapsLog[ptBin], PIDconfig, dEdxvsplogFigPNG);
//	makedEdxvspFigloglog( dEdxvsPMapsLog[ptBin], PIDconfig ,dEdxvsplogFigPDF);
//	makedEdxvspFiglinlin( dEdxvsPMapsLin[ptBin], PIDconfig, dEdxvsplinFigPNG);
//	makedEdxvspFiglinlin( dEdxvsPMapsLin[ptBin], PIDconfig, dEdxvsplinFigPDF);
//
// }


  for( int ptBin = 0; ptBin < npt; ptBin++ )
  {

  	double pt1 = (ptMin + (ptBin  ) * ptbw);
 	double pt2 = (ptMin + (ptBin+1) * ptbw);


 	std::string dEdxvsplinFigBase = Form( "dEdxvspLin_typ_%.2f-%.2f", pt1, pt2 	);
 	std::string dEdxvsplogFigBase = Form( "dEdxvspLog_typ_%.2f-%.2f", pt1, pt2 	);

 	std::string dEdxvsplinFigPNG = dEdxvsplinFigBase+".png";
 	std::string dEdxvsplinFigPDF = dEdxvsplinFigBase+".pdf";
 	std::string dEdxvsplogFigPNG = dEdxvsplogFigBase+".png";
 	std::string dEdxvsplogFigPDF = dEdxvsplogFigBase+".pdf";

	makedEdxvspFigloglog( dEdxvsPMapsLog[ptBin], PIDconfig, dEdxvsplogFigPNG);
	makedEdxvspFigloglog( dEdxvsPMapsLog[ptBin], PIDconfig ,dEdxvsplogFigPDF);
	makedEdxvspFiglinlin( dEdxvsPMapsLin[ptBin], PIDconfig, dEdxvsplinFigPNG);
	makedEdxvspFiglinlin( dEdxvsPMapsLin[ptBin], PIDconfig, dEdxvsplinFigPDF);

  }

 plotTH2D(ptres, ";p_{T,gen} [GeV/c];p_{T,reco} [GeV/c]", "ptgenreco", "COLZ");

 output->cd();

 for( int RgnBin = 0; RgnBin < CM::nRegions; RgnBin++ )
 for( int i = 0; i < CM::nPtBins[RgnBin]; i++ )
 {
	int ptBin = CM::PtBins[1][RgnBin][i];

   double pt1 = (CM::ptMin + (ptBin  ) * CM::PtBw);
	double pt2 = (CM::ptMin + (ptBin+1) * CM::PtBw);

	std::string str_cont_matrix_nor_TH2D    = Form("cont_matrix_nor_TH2D_pt_%.2f-%2.2f", pt1, pt2);
	std::string str_cont_matrix_nor_TMatrix = Form("cont_matrix_nor_TMatrix_pt_%.2f-%2.2f", pt1, pt2);
	std::string str_cont_matrix_inv_TMatrix = Form("cont_matrix_inv_TMatrix_pt_%.2f-%2.2f", pt1, pt2);

	cont_matrix_nor_TH2D[ptBin]->Write(str_cont_matrix_nor_TH2D.c_str());
	cont_matrix_nor_TMatrix[ptBin]->Write(str_cont_matrix_nor_TMatrix.c_str());
	cont_matrix_inv_TMatrix[ptBin]->Write(str_cont_matrix_inv_TMatrix.c_str());

 }

 output->Close();

 
}
