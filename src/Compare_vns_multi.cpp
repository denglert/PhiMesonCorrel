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


 const int nFiles = 7;
 int centerbin = 2;

 std::string inpFilenames[nFiles];
 std::string       labels[nFiles];

 //inpFilenames[0] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_trueHighMultiplicity_config_looser_2_trkCorr_no/dump.root";
 //inpFilenames[1] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_trueHighMultiplicity_config_looser_1_trkCorr_no/dump.root";
 //inpFilenames[2] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_trueHighMultiplicity_config_default_trkCorr_no/dump.root";
 //inpFilenames[3] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_trueHighMultiplicity_config_strict_1_trkCorr_no/dump.root";
 //inpFilenames[4] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_trueHighMultiplicity_config_strict_2_trkCorr_no/dump.root";

// inpFilenames[0] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_config_looser_2_trkCorr_no/dump.root";
// inpFilenames[1] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_config_looser_1_trkCorr_no/dump.root";
// inpFilenames[2] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_config_default_trkCorr_no/dump.root";
// inpFilenames[3] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_config_strict_1_trkCorr_no/dump.root";
// inpFilenames[4] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_config_strict_2_trkCorr_no/dump.root";

 inpFilenames[0] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_full_config_0_trkCorr_no/dump.root";
 inpFilenames[1] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_full_config_1_trkCorr_no/dump.root";
 inpFilenames[2] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_full_config_2_trkCorr_no/dump.root";
 inpFilenames[3] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_full_config_3_trkCorr_no/dump.root";
 inpFilenames[4] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_full_config_4_trkCorr_no/dump.root";
 inpFilenames[5] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_full_config_5_trkCorr_no/dump.root";
 inpFilenames[6] = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/PIDscan_MinBias_dEdxminsweep_full_config_6_trkCorr_no/dump.root";


 labels[0] = "config_0";
 labels[1] = "config_1";
 labels[2] = "config_2";
 labels[3] = "config_3";
 labels[4] = "config_4";
 labels[5] = "config_5";
 labels[6] = "config_6";

// const int nMultiplicityBins = 4;
// int multbins[nMultiplicityBins][2] = {
//												  { 120, 150 },
//												  { 150, 185 },
//												  { 185, 220 },
//												  { 220, 260 }
// 												  };

// 0 - 120
 const int nMultiplicityBins = 1;
 int multbins[nMultiplicityBins][2] = {
	  											  {   0,  120 }
 												  };



 std::string CMSsystemlabel = Form("#splitline{CMS (work in progress)}{MinBias}");


 double fit_interval[4][2] = {
									  {0.2, 3.0},
									  {0.2, 0.9},
									  {0.2, 0.9},
									  {0.2, 1.4}
 									  };

 // Plot style
 
 double MarkerSizes[7] = {       2,      2,     2,     2,      2,        2,    2 };
 int   MarkerStyles[7] = {      32,     27,     5,    25,     30,       28,   26 };
 int     LineStyles[7] = {       9,      9,     9,     9,      9,        9,    9 };
 int         Colors[7] = { kYellow, kGreen, kCyan, kBlue, kBlack, kMagenta, kRed };

 // Binning
 int nCorrTyp = 4; 
 int *nPtBins = new int[nCorrTyp];

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
 TGraphErrors *data[nCorrTyp][nFiles];

 // Declare fit functions (TF1)
 TF1 *fit[nCorrTyp][nFiles];

 // MultBin loop
 for(int multBin = 0; multBin < nMultiplicityBins; multBin++)
 {

 int mult1 = multbins[multBin][0];
 int mult2 = multbins[multBin][1];

 std::string multlabel = Form("%3d #leq N_{trk}^{offline} #leq %3d", mult1, mult2);

 // Read in data points
 for(int i = 0; i < nCorrTyp; i++)
 for(int j = 0; j < nFiles; j++)
 {
   std::string dataname = Form("%s_Ntrk_%03d-%03d", particletypelabel(i).c_str(),  mult1, mult2 );
 	data[i][j] = (TGraphErrors*)inpFiles[j]->Get( dataname.c_str() );
	if ( data[i][j]->IsZombie() ){ std::cerr << Form("TGraphErrors %s not", dataname.c_str()) << std::endl; return -1;}
 }

 for(int i = 0; i < nCorrTyp; i++)
 {  nPtBins[i] = data[i][0]->GetN(); }

 // Data placeholders
 double       *xpts[nCorrTyp][nFiles];
 double       *ypts[nCorrTyp][nFiles];
 double   *yreldiff[nCorrTyp][nFiles];
 double *yreldiff_E[nCorrTyp][nFiles];

 // Fill in arrays with 'x' and 'y' values
 for(int i = 0; i < nCorrTyp; i++)
 { 
	nPtBins[i] = data[i][0]->GetN();
 	for(int j = 0; j < nFiles; j++)
	{
		      xpts[i][j] = data[i][j]->GetX();
		      ypts[i][j] = data[i][j]->GetX();
	}
 };

 // Calculate relative difference
 for(int i = 0; i < nCorrTyp; i++)
 for(int j = 0; j < nFiles; j++)
 {
		      yreldiff[i][j] = data[i][j]->GetX();
 }

 //////////////////////
 //                  //
 // ***** Plot ***** //
 //                  //
 //////////////////////
 
 // Debuggg 
 std::cout << "Plotting section." << std::endl;

 for(int i = 0; i < nCorrTyp; i++)
 for(int j = 0; j < nFiles; j++)
 {

 	// Discard fit that is in the file
   data[i][j] -> GetFunction(Form("%sv2_fit", particletypelabel(i).c_str() ))->SetBit(TF1::kNotDraw);

	// Create new fit functions
	fit[i][j] = new TF1 (Form("%s_fit_%d", particletypelabel(i).c_str(), j), "[0]+[1]*x+[2]*x*x", fit_interval[i][0], fit_interval[i][1]);


 	// Cosmetics
	data[i][j]->SetMarkerStyle( MarkerStyles[j] );
	data[i][j]->SetMarkerSize( MarkerSizes[j] );
	data[i][j]->SetMarkerColor( Colors[j] );
	data[i][j]->SetLineColor(   Colors[j] );

	fit[i][j]->SetLineColor( Colors[j] );
	fit[i][j]->SetLineStyle( LineStyles[j] );

	// Fit the data
//   data[i][j]->Fit( Form("%s_fit_%d", particletypelabel(i).c_str(), j), "R");

 }


 // *** Plotting the graphs *** //
 gStyle->SetPadTickY(1);
 gStyle->SetPadTickX(1);
 
 for(int i = 0; i < nCorrTyp; i++)
 {

 	TCanvas canvas_v2_vs_pT ("v2 vs pT", "v_{2} values as a function of p_{T}; p_{T} [GeV/c];v_{2}", 1024, 768);
 	
 	canvas_v2_vs_pT.SetLeftMargin(0.10);
 	canvas_v2_vs_pT.SetBottomMargin(0.12);
 	canvas_v2_vs_pT.SetRightMargin(0.05);
 	canvas_v2_vs_pT.SetTopMargin(0.05);
 	
 	double v2vspt_ptmin = 0.0;
 	double v2vspt_ptmax = 2.5;
 	double v2vspt_v2min = 0.0;
 	double v2vspt_v2max = 0.16;
 	
 	data[i][0]->SetTitle("");
 	data[i][0]->GetXaxis()->SetLimits(v2vspt_ptmin,v2vspt_ptmax);
 	data[i][0]->GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
 	data[i][0]->GetXaxis()->CenterTitle(1);
 	data[i][0]->GetXaxis()->SetTitleOffset(1.2);
 	data[i][0]->GetXaxis()->SetTitleSize(figuretextsize);
 	data[i][0]->GetYaxis()->SetRangeUser(v2vspt_v2min,v2vspt_v2max);
 	data[i][0]->GetYaxis()->SetTitle("v_{2}");
 	data[i][0]->GetYaxis()->CenterTitle(1);
 	data[i][0]->GetYaxis()->SetTitleOffset(1.2);
 	data[i][0]->GetYaxis()->SetTitleSize(figuretextsize);
 
 	////////////////////
 	// Making figures //
	
	data[i][0]->Draw("AP");
	for(int j = 1; j < nFiles; j++)
	{data[i][j]->Draw("P"); }

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

	std::string corrtypelabel = particletype(i);

	TLegend gr_legend (gr_legend_x1, gr_legend_y1, gr_legend_x2, gr_legend_y2);
	gr_legend.SetFillStyle(0);
	gr_legend.SetBorderSize(0);
	for(int j = 0; j < nFiles; j++)
	{ gr_legend.AddEntry(data[i][j], labels[j].c_str(), "P");}
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

	TLatex tcorrtypelabel( corrtypelabelposx,corrtypelabelposy, corrtypelabel.c_str()); 
	tcorrtypelabel.SetTextSize(figuretextsize);
	tcorrtypelabel.SetNDC(kTRUE);
	tcorrtypelabel.Draw();

	
 	std::string label = Form("v2vspt_%s_nTrk_%03d-%03d", particletypelabel(i).c_str(), mult1, mult2); 
	std::string	pngfigure = label+".png";
	std::string	pdffigure = label+".pdf";
	canvas_v2_vs_pT.SaveAs( pngfigure.c_str() );
	canvas_v2_vs_pT.SaveAs( pdffigure.c_str() );

	}


   ///////////////////////////////////////////////////////////
	// Plot relative difference compared to centerbin figure //
	///////////////////////////////////////////////////////////
	
//	// Declare data points holders (TGraphErrors)
//	 TGraphErrors *[nCorrTyp][nFiles];
//	
//	 // Declare fit functions (TF1)
//	 TF1 *fit[nCorrTyp][nFiles];
//	
//	 for(int multBin = 0; multBin < nMultiplicityBins; multBin++)
//	 {
//	
//	 int mult1 = multbins[multBin][0];
//	 int mult2 = multbins[multBin][1];
//	
//	 std::string multlabel = Form("%3d #leq N_{trk}^{offline} #leq %3d", mult1, mult2);
//	
//	 // Read in data points
//	 for(int i = 0; i < nCorrTyp; i++)
//	 for(int j = 0; j < nFiles; j++)
//	 {
//	   std::string dataname = Form("%s_Ntrk_%03d-%03d", particletypelabel(i).c_str(),  mult1, mult2 );
//	 	data[i][j] = (TGraphErrors*)inpFiles[j]->Get( dataname.c_str() );
//		if ( data[i][j]->IsZombie() ){ std::cerr << Form("TGraphErrors %s not", dataname.c_str()) << std::endl; return -1;}
//	 }
//	
//	
//		TCanvas canvas_reldiff ("rv2 vs pT", "v_{2} values as a function of p_{T}; p_{T} [GeV/c];v_{2}", 1024, 768);
//	
//		canvas_reldiff.SetLeftMargin (0.10);
//	   canvas_reldiff.SetBottomMargin(0.12);
//	   canvas_reldiff.SetRightMargin(0.05);
//	   canvas_reldiff.SetTopMargin  (0.05);
//		
//		double reldiffvspt_ptmin = 0.0;
//		double reldiffvspt_ptmax = 2.5;
//		double reldiffvspt_rmin  = -0.5;
//		double reldiffvspt_rmax  =  0.5;
//	
//		TGraphErrors cpar_reldiffgraph = cpar_v2( nPtBins[0], cpar_2_xpts, cpar_reldiff, 0, cpar_reldiff_E );
//		TGraphErrors pion_reldiffgraph = pion_v2( nPtBins[1], pion_2_xpts, pion_reldiff, 0, pion_reldiff_E );
//		TGraphErrors kaon_reldiffgraph = kaon_v2( nPtBins[2], kaon_2_xpts, kaon_reldiff, 0, kaon_reldiff_E );
//		TGraphErrors prot_reldiffgraph = prot_v2( nPtBins[3], prot_2_xpts, prot_reldiff, 0, prot_reldiff_E );
//	
//		cpar_reldiffgraph.SetTitle("");
//	   cpar_reldiffgraph.GetXaxis()->SetLimits(reldiffvspt_ptmin, reldiffvspt_ptmax);
//		cpar_reldiffgraph.GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
//		cpar_reldiffgraph.GetXaxis()->CenterTitle(1);
//		cpar_reldiffgraph.GetXaxis()->SetTitleOffset(1.2);
//		cpar_reldiffgraph.GetXaxis()->SetTitleSize(figuretextsize);
//		cpar_reldiffgraph.GetYaxis()->SetRangeUser(reldiffvspt_rmin,reldiffvspt_rmax);
//		cpar_reldiffgraph.GetYaxis()->SetLimits(reldiffvspt_rmin,reldiffvspt_rmax);
//		cpar_reldiffgraph.GetYaxis()->SetTitle("relative difference");
//		cpar_reldiffgraph.GetYaxis()->CenterTitle(1);
//		cpar_reldiffgraph.GetYaxis()->SetTitleOffset(1.2);
//		cpar_reldiffgraph.GetYaxis()->SetTitleSize(figuretextsize);
//	
//		std::string label_reldiff = Form("reldiff_MC_comprarison_nTrk_%03d-%03d", mult1, mult2); 
//		std::string pngfigure_reldiff = label_reldiff+".png";
//		std::string pdffigure_reldiff = label_reldiff+".pdf";
//	
//		cpar_reldiffgraph.SetMarkerStyle(25);
//		pion_reldiffgraph.SetMarkerStyle(26);
//		kaon_reldiffgraph.SetMarkerStyle(30);
//		prot_reldiffgraph.SetMarkerStyle(24);
//	
//		cpar_reldiffgraph.Draw("AP");
//		pion_reldiffgraph.Draw("P");
//		kaon_reldiffgraph.Draw("P");
//		prot_reldiffgraph.Draw("P");
//	
//	   double reldiff_min = -0.05;
//	   double reldiff_max =  0.05;
//	
//		TLine *line_upp  = new TLine(reldiffvspt_ptmin,reldiff_min,reldiffvspt_ptmax,reldiff_min); 
//		TLine *line_mid  = new TLine(reldiffvspt_ptmin,0,reldiffvspt_ptmax,0);
//		TLine *line_low  = new TLine(reldiffvspt_ptmin,reldiff_max,reldiffvspt_ptmax,reldiff_max); 
//	
//		line_upp->SetLineStyle(7);
//		line_mid->SetLineStyle(7);
//		line_low->SetLineStyle(7);
//	
//		line_upp->Draw("SAME");
//		line_upp->Draw("SAME");
//		line_low->Draw("SAME");
//	
//		canvas_reldiff.SaveAs( pngfigure_reldiff.c_str() );
//		canvas_reldiff.SaveAs( pdffigure_reldiff.c_str() );
//
//
//
  }



}
