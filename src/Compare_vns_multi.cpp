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
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSpline.h"
#include "GraphStyles.h"


int main( int argc, const char *argv[] )
{ 

 std::cout << "Compare_vns_multi binary started." << std::endl;


 const int nFiles = 5;
 int centerbin = 2;

 std::string inpFilenames[nFiles];
 std::string       labels[nFiles];

 inpFilenames[0] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/PIDscan_HighMult_Full_KaondEdxminWweep_Full_1nw_config_0_trkCorr_no_temp_fit_limunc/dump.root";
 inpFilenames[1] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/PIDscan_HighMult_Full_KaondEdxminSweep_Full_1nw_config_1_trkCorr_no_temp_fit_limunc/dump.root";
 inpFilenames[2] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/PIDscan_HighMult_Full_KaondEdxminSweep_Full_1nw_config_2_trkCorr_no_temp_fit_limunc/dump.root";
 inpFilenames[3] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/PIDscan_HighMult_Full_KaondEdxminSweep_Full_1nw_config_3_trkCorr_no_temp_fit_limunc/dump.root";
 inpFilenames[4] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/PIDscan_HighMult_Full_KaondEdxminSweep_Full_1nw_config_4_trkCorr_no_temp_fit_limunc/dump.root";

 labels[0] = "PID_0";
 labels[1] = "PID_1";
 labels[2] = "PID_2";
 labels[3] = "PID_3";
 labels[4] = "PID_4";

// 0 - 120
 const int nMultiplicityBins = 1;
 int multbins[nMultiplicityBins][2] = {
	  											  {  120,  180 }
 												  };



 std::string CMSsystemlabel = Form("#splitline{CMS (work in progress)}{HighMult}");

 // Plot style
 
 double MarkerSizes[7] = {       2,      2,     2,     2,      2,        2,    2 };
 int   MarkerStyles[7] = {      32,     27,     5,    25,     30,       28,   26 };
 int     LineStyles[7] = {       9,      9,     9,     9,      9,        9,    9 };
//int         Colors[7] = { kYellow, kGreen, kCyan, kBlue, kBlack, kMagenta, kRed };
 int         Colors[7] = { kCyan, kBlue, kBlack, kMagenta, kRed,  kYellow, kGreen };
 double       shift[7] = {   -0.02,   -0.01,    0,  0.01,   0.02,     0.03, 0.04 };

 // Binning
 int *nPtBins;

 ///////////////////////////////////////
 //                                   //
 // ****** Opening input files ****** //
 //                                   //
 ///////////////////////////////////////
 
 TFile *inpFiles[nFiles];

 for(int i = 0; i < nFiles; i++)
 { 
	inpFiles[i] = new TFile( inpFilenames[i].c_str(), "READ");
	if ( inpFiles[i]->IsZombie() ) {std::cerr << "Error opening file : " << inpFilenames[i] << std::endl; exit(-1); }
 }
 
 //////////////////////////////
 //                          //
 // ***** Initializing ***** //
 //                          //
 //////////////////////////////
 
 // Declare data points holders (TGraphErrors)
 TGraphErrors *data[nFiles];

 // Declare fit functions (TF1)
 TF1 *fit[nFiles];

 // MultBin loop
 for(int multBin = 0; multBin < nMultiplicityBins; multBin++)
 {

 int mult1 = multbins[multBin][0];
 int mult2 = multbins[multBin][1];

 std::string multlabel = Form("%3d #leq N_{trk}^{offline} #leq %3d", mult1, mult2);

 // Read in data points
 for(int j = 0; j < nFiles; j++)
 {
   std::string dataname = Form("phimv2vspt_avg_Ntrk_%03d-%03d",   mult1, mult2 );
 	data[j] = (TGraphErrors*)inpFiles[j]->Get( dataname.c_str() );
	data[j] -> RemovePoint(0);
	if ( data[j]->IsZombie() ){ std::cerr << Form("TGraphErrors %s not", dataname.c_str()) << std::endl; return -1;}
 }

 int nPtBins = data[0]->GetN();

 // Data placeholders
 double       *xpts[nFiles];
 double       *ypts[nFiles];
 double   *yreldiff[nFiles];
 double *yreldiff_E[nFiles];

 	for(int j = 0; j < nFiles; j++)
	{
		      xpts[j] = data[j]->GetX();
		      ypts[j] = data[j]->GetY();
	}

	// shift
 	for(int j = 0; j < nFiles; j++)
 	for(int ptBin = 0; ptBin < nPtBins; ptBin++)
	{
		    xpts[j][ptBin] = xpts[j][ptBin] + shift[j];
	}

 // Calculate relative difference
 for(int j = 0; j < nFiles; j++)
 {
		      yreldiff[j] = data[j]->GetX();
 }

 //////////////////////
 //                  //
 // ***** Plot ***** //
 //                  //
 //////////////////////
 
 // Debuggg 
 std::cout << "Plotting section." << std::endl;

 for(int j = 0; j < nFiles; j++)
 {

 	// Cosmetics
	data[j]->SetMarkerStyle( MarkerStyles[j] );
	data[j]->SetMarkerSize( MarkerSizes[j] );
	data[j]->SetMarkerColor( Colors[j] );
	data[j]->SetLineColor(   Colors[j] );

 }

 // *** Plotting the graphs *** //
 gStyle->SetPadTickY(1);
 gStyle->SetPadTickX(1);
 
 	TCanvas canvas_v2_vs_pT ("v2 vs pT", "v_{2} values as a function of p_{T}; p_{T} [GeV/c];v_{2}", 1024, 768);
 	
 	canvas_v2_vs_pT.SetLeftMargin(0.10);
 	canvas_v2_vs_pT.SetBottomMargin(0.12);
 	canvas_v2_vs_pT.SetRightMargin(0.05);
 	canvas_v2_vs_pT.SetTopMargin(0.05);
 	
 	double v2vspt_ptmin = 0.0;
 	double v2vspt_ptmax = 2.5;
 	double v2vspt_v2min = 0.0;
 	double v2vspt_v2max = 0.16;
 	
 	data[0]->SetTitle("");
 	data[0]->GetXaxis()->SetLimits(v2vspt_ptmin,v2vspt_ptmax);
 	data[0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
 	data[0]->GetXaxis()->CenterTitle(1);
 	data[0]->GetXaxis()->SetTitleOffset(1.2);
 	data[0]->GetXaxis()->SetTitleSize(figuretextsize);
 	data[0]->GetYaxis()->SetRangeUser(v2vspt_v2min,v2vspt_v2max);
 	data[0]->GetYaxis()->SetTitle("v_{2}");
 	data[0]->GetYaxis()->CenterTitle(1);
 	data[0]->GetYaxis()->SetTitleOffset(1.2);
 	data[0]->GetYaxis()->SetTitleSize(figuretextsize);
 
 	////////////////////
 	// Making figures //
	
	data[0]->Draw("AP");
	for(int j = 1; j < nFiles; j++)
	{data[j]->Draw("P"); }

	double gr_legend_x1=0.6;
	double gr_legend_y1=0.4;
	double gr_legend_x2=gr_legend_x1+.20;
	double gr_legend_y2=gr_legend_y1+.20;

	double CMSsystemlabelposx = 0.14;
	double CMSsystemlabelposy = 0.84;

	double multlabelposx = 0.62;
	double multlabelposy = 0.24;

	double corrtypelabelposx = 0.62;
	double corrtypelabelposy = 0.80;

	TLegend gr_legend (gr_legend_x1, gr_legend_y1, gr_legend_x2, gr_legend_y2);
	gr_legend.SetFillStyle(0);
	gr_legend.SetBorderSize(0);
	for(int j = 0; j < nFiles; j++)
	{ gr_legend.AddEntry(data[j], labels[j].c_str(), "P");}
	gr_legend.SetTextSize(figuretextsize);
	gr_legend.Draw("SAME");

	TLatex tCMSsystemlabel( CMSsystemlabelposx,CMSsystemlabelposy, CMSsystemlabel.c_str()); 
	tCMSsystemlabel.SetTextSize(figuretextsize);
	tCMSsystemlabel.SetNDC(kTRUE);
	tCMSsystemlabel.Draw();

	TLatex tmultlabel( multlabelposx,multlabelposy, multlabel.c_str()); 
	tmultlabel.SetTextSize(figuretextsize);
	tmultlabel.SetNDC(kTRUE);
	tmultlabel.Draw();

//	TLatex tcorrtypelabel( corrtypelabelposx,corrtypelabelposy, corrtypelabel.c_str()); 
//	tcorrtypelabel.SetTextSize(figuretextsize);
//	tcorrtypelabel.SetNDC(kTRUE);
//	tcorrtypelabel.Draw();

	
 	std::string label = Form("v2vspt_nTrk_%03d-%03d",  mult1, mult2); 
	std::string	pngfigure = label+".png";
	std::string	pdffigure = label+".pdf";
	canvas_v2_vs_pT.SaveAs( pngfigure.c_str() );
	canvas_v2_vs_pT.SaveAs( pdffigure.c_str() );

 	}


}
