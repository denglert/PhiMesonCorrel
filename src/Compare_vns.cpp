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

const int mult1 = 0;
const int mult2 = 120;

int main( int argc, const char *argv[] )
{ 

  if(argc != 5)
  {
    std::cerr << "Usage: Compare_vns" << std::endl;
	 exit(1);
  }

  std::string label1 = argv[1];
  std::string label2 = argv[2];
 TString inpFilename_1 = argv[3];
 TString inpFilename_2 = argv[4];

 // Binning
 int nCorrTyp			  = nCorrTyp_; 
 int *nPtBins = new int[nCorrTyp_];

 for(int TypBin = 0; TypBin < nCorrTyp; TypBin++)
 { nPtBins[TypBin] = nPtBins_[TypBin]; }

 int nMultiplicityBins_Ana = nMultiplicityBins_Ana_HDR;
 int nMultiplicityBins_EvM = nMultiplicityBins_EvM_HDR;
 int nZvtxBins 		      = nZvtxBins_; 

 std::cout << "Compare_vns running." << std::endl;


 ///////////////////////////////////////
 //                                   //
 // ****** Opening input files ****** //
 //                                   //
 ///////////////////////////////////////
 
 TFile *f1 = new TFile(inpFilename_1, "READ");
 if ( f1->IsZombie() ) {std::cerr << "Error opening file 1: " << inpFilename_1 << std::endl; exit(-1);}

 TFile *f2 = new TFile(inpFilename_2, "READ");
 if ( f2->IsZombie() ) {std::cerr << "Error opening file 2: " << inpFilename_2 << std::endl; exit(-1);}

 //////////////////////////////
 //                          //
 // ***** Initializing ***** //
 //                          //
 //////////////////////////////

 TGraphErrors *cpar_1 = (TGraphErrors*)f1->Get( Form("cpar_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *pion_1 = (TGraphErrors*)f1->Get( Form("pion_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *kaon_1 = (TGraphErrors*)f1->Get( Form("kaon_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *prot_1 = (TGraphErrors*)f1->Get( Form("prot_Ntrk_%03d-%03d", mult1, mult2 ) );

 TGraphErrors *cpar_2 = (TGraphErrors*)f2->Get( Form("cpar_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *pion_2 = (TGraphErrors*)f2->Get( Form("pion_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *kaon_2 = (TGraphErrors*)f2->Get( Form("kaon_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *prot_2 = (TGraphErrors*)f2->Get( Form("prot_Ntrk_%03d-%03d", mult1, mult2 ) );

 cpar_1->SetPoint(0, -50, -5);
 pion_1->SetPoint(0, -50, -5);
 kaon_1->SetPoint(0, -50, -5);
 prot_1->SetPoint(0, -50, -5);

 cpar_2->SetPoint(0, -50, -5);
 pion_2->SetPoint(0, -50, -5);
 kaon_2->SetPoint(0, -50, -5);
 prot_2->SetPoint(0, -50, -5);

 //////////////////////
 //                  //
 // ***** Plot ***** //
 //                  //
 //////////////////////
 
	cpar_1->SetMarkerStyle(21);
	pion_1->SetMarkerStyle(22);
	kaon_1->SetMarkerStyle(29);
	prot_1->SetMarkerStyle(20);
	
	cpar_2->SetMarkerStyle(25);
	pion_2->SetMarkerStyle(26);
	kaon_2->SetMarkerStyle(30);
	prot_2->SetMarkerStyle(4);

//	cpar_1->GetFunction("cparv2_fit")->SetBit(TF1::kNotDraw);
//	pion_1->GetFunction("pionv2_fit")->SetBit(TF1::kNotDraw);
//	kaon_1->GetFunction("kaonv2_fit")->SetBit(TF1::kNotDraw);
//	prot_1->GetFunction("protv2_fit")->SetBit(TF1::kNotDraw);
	//cpar_2->GetFunction("cparv2_fit")->SetBit(TF1::kNotDraw);
	//pion_2->GetFunction("pionv2_fit")->SetBit(TF1::kNotDraw);
	//kaon_2->GetFunction("kaonv2_fit")->SetBit(TF1::kNotDraw);
	//prot_2->GetFunction("protv2_fit")->SetBit(TF1::kNotDraw);

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
	double v2vspt_v2max = 0.6;

	cpar_1->SetTitle("");
   cpar_1->GetXaxis()->SetLimits(v2vspt_ptmin,v2vspt_ptmax);
	cpar_1->GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
	cpar_1->GetXaxis()->CenterTitle(1);
	cpar_1->GetXaxis()->SetTitleOffset(1.2);
	cpar_1->GetXaxis()->SetTitleSize(figuretextsize);
	cpar_1->GetYaxis()->SetRangeUser(v2vspt_v2min,v2vspt_v2max);
	cpar_1->GetYaxis()->SetTitle("v_{2}");
	cpar_1->GetYaxis()->CenterTitle(1);
	cpar_1->GetYaxis()->SetTitleOffset(1.2);
	cpar_1->GetYaxis()->SetTitleSize(figuretextsize);

	TSpline3 *cpar_spl = new TSpline3("cpar_spl", cpar_1);
	TSpline3 *pion_spl = new TSpline3("pion_spl", pion_1);
	TSpline3 *kaon_spl = new TSpline3("kaon_spl", kaon_1);
	TSpline3 *prot_spl = new TSpline3("prot_spl", prot_1);

	cpar_spl->SetLineColor(kBlack);
	pion_spl->SetLineColor(kRed);
	kaon_spl->SetLineColor(kGreen);
	prot_spl->SetLineColor(kBlue);

	TF1 *cparv2_fit = new TF1 ("cparv2_fit", "[0]+[1]*x+[2]*x*x",0.2,3);
	TF1 *pionv2_fit = new TF1 ("pionv2_fit", "[0]+[1]*x+[2]*x*x",0.2,0.9);
	TF1 *kaonv2_fit = new TF1 ("kaonv2_fit", "[0]+[1]*x+[2]*x*x",0.2,0.9);
	TF1 *protv2_fit = new TF1 ("protv2_fit", "[0]+[1]*x+[2]*x*x",0.2,1.4);

	cparv2_fit->SetLineColor(kBlack);
	pionv2_fit->SetLineColor(kRed);
	kaonv2_fit->SetLineColor(kGreen);
	protv2_fit->SetLineColor(kBlue);

	cparv2_fit->SetLineStyle(9);
	pionv2_fit->SetLineStyle(9);
	kaonv2_fit->SetLineStyle(9);
	protv2_fit->SetLineStyle(9);

   cpar_1->Fit("cparv2_fit", "0R");
   pion_1->Fit("pionv2_fit", "0R");
   kaon_1->Fit("kaonv2_fit", "0R");
   prot_1->Fit("protv2_fit", "0R");

	// Draw
	cpar_1->Draw("AP");
	pion_1->Draw("P");
	kaon_1->Draw("P");
	prot_1->Draw("P");

	cpar_2->Draw("P");
	pion_2->Draw("P");
	kaon_2->Draw("P");
	prot_2->Draw("P");

//	cpar_spl->Draw("SAME");
//   pion_spl->Draw("SAME");
//   kaon_spl->Draw("SAME");
//   prot_spl->Draw("SAME");

	// Legends, texts
	double legend_x1=.14;
	double legend_y1=0.56;
	double legend_x2=legend_x1+.20;
	double legend_y2=legend_y1+.20;

	double gr_legend_x1=0.50;
	double gr_legend_y1=0.30;
	double gr_legend_x2=gr_legend_x1+.10;
	double gr_legend_y2=gr_legend_y1+.10;


	double CMSsystemlabelposx = 0.14;
	double CMSsystemlabelposy = 0.84;
	double multlabelposx = 0.58;
	double multlabelposy = 0.20;


	TLegend v2vsptlegend (legend_x1, legend_y1, legend_x2, legend_y2);
	v2vsptlegend.SetFillStyle(0);
	v2vsptlegend.SetBorderSize(0);
	v2vsptlegend.AddEntry(cpar_1,"charged", "P");
	v2vsptlegend.AddEntry(pion_1,"#pi", "P");
	v2vsptlegend.AddEntry(kaon_1,"K", "P");
	v2vsptlegend.AddEntry(prot_1,"p", "P");
	v2vsptlegend.SetTextSize(figuretextsize);
	v2vsptlegend.Draw("SAME");

	TF1* func_1 = new TF1("func_1", "1", 0, 2);
	TF1* func_2 = new TF1("func_2", "1", 0, 2);

	func_1->SetLineColor(kBlack);
	func_2->SetLineColor(kBlack);
	func_1->SetMarkerStyle(21);
	func_2->SetMarkerStyle(25);

	TLegend gr_legend (gr_legend_x1, gr_legend_y1, gr_legend_x2, gr_legend_y2);
	gr_legend.SetFillStyle(0);
	gr_legend.SetBorderSize(0);
	gr_legend.AddEntry(func_1,Form("%s", label1.c_str()) , "P");
   gr_legend.AddEntry(func_2,Form("%s", label2.c_str()) , "P");
	gr_legend.SetTextSize(figuretextsize);
	gr_legend.Draw("SAME");


	std::string CMSsystemlabel = Form("#splitline{CMS (work in progress)}{pPb EPOS LHC MinBias}", mult1, mult2);
	std::string multlabel = Form("0.3 < p_{T}^{assoc} < 3.0 GeV/c", mult1, mult2);

	TLatex tCMSsystemlabel( CMSsystemlabelposx,CMSsystemlabelposy, CMSsystemlabel.c_str()); 
	tCMSsystemlabel.SetTextSize(figuretextsize);
	tCMSsystemlabel.SetNDC(kTRUE);
	tCMSsystemlabel.Draw();

	TLatex tmultlabel( multlabelposx,multlabelposy, multlabel.c_str()); 
	tmultlabel.SetTextSize(figuretextsize);
	tmultlabel.SetNDC(kTRUE);
	tmultlabel.Draw();
	
	std::string label = Form("v2vspt_MC_comprarison_nTrk_%03d-%03d", mult1, mult2); 
	std::string	pngfigure = label+".png";
	std::string	pdffigure = label+".pdf";
	canvas_v2_vs_pT.SaveAs( pngfigure.c_str() );
	canvas_v2_vs_pT.SaveAs( pdffigure.c_str() );

	double *cpar_1_Eypts = cpar_1->GetEY();
	double *pion_1_Eypts = pion_1->GetEY();
	double *kaon_1_Eypts = kaon_1->GetEY();
	double *prot_1_Eypts = prot_1->GetEY();

	double *cpar_2_Eypts = cpar_2->GetEY();
	double *pion_2_Eypts = pion_2->GetEY();
	double *kaon_2_Eypts = kaon_2->GetEY();
	double *prot_2_Eypts = prot_2->GetEY();

	double *cpar_1_ypts = cpar_1->GetY();
	double *pion_1_ypts = pion_1->GetY();
	double *kaon_1_ypts = kaon_1->GetY();
	double *prot_1_ypts = prot_1->GetY();

	double *cpar_2_ypts = cpar_2->GetY();
	double *pion_2_ypts = pion_2->GetY();
	double *kaon_2_ypts = kaon_2->GetY();
	double *prot_2_ypts = prot_2->GetY();

	double *cpar_2_xpts = cpar_2->GetX();
	double *pion_2_xpts = pion_2->GetX();
	double *kaon_2_xpts = kaon_2->GetX();
	double *prot_2_xpts = prot_2->GetX();

	double *cpar_1_relE = new double[nPtBins[0]];
	double *pion_1_relE = new double[nPtBins[1]];
	double *kaon_1_relE = new double[nPtBins[2]];
	double *prot_1_relE = new double[nPtBins[3]];

	double *cpar_2_relE = new double[nPtBins[0]];
	double *pion_2_relE = new double[nPtBins[1]];
	double *kaon_2_relE = new double[nPtBins[2]];
	double *prot_2_relE = new double[nPtBins[3]];

//	for(int ptBin = 0; ptBin < nPtBins[0]; ptBin++)
//	{
//		cpar_1_ypts[ptBin]  = cpar_1_ypts[ptBin]  / cparv2_fit->Eval( cpar_2_xpts[ptBin] );
//		cpar_2_ypts[ptBin]  = cpar_2_ypts[ptBin]  / cparv2_fit->Eval( cpar_2_xpts[ptBin] );
//		cpar_1_Eypts[ptBin] = cpar_1_Eypts[ptBin] / cparv2_fit->Eval( cpar_2_xpts[ptBin] );
//		cpar_2_Eypts[ptBin] = cpar_2_Eypts[ptBin] / cparv2_fit->Eval( cpar_2_xpts[ptBin] );
//	}
//
//	for(int ptBin = 0; ptBin < nPtBins[1]; ptBin++)
//	{
//		pion_1_ypts[ptBin]  = pion_1_ypts[ptBin]  / pionv2_fit->Eval( pion_2_xpts[ptBin] );
//		pion_2_ypts[ptBin]  = pion_2_ypts[ptBin]  / pionv2_fit->Eval( pion_2_xpts[ptBin] );
//		pion_1_Eypts[ptBin] = pion_1_Eypts[ptBin] / pionv2_fit->Eval( pion_2_xpts[ptBin] );
//		pion_2_Eypts[ptBin] = pion_2_Eypts[ptBin] / pionv2_fit->Eval( pion_2_xpts[ptBin] );
//	}
//
//	for(int ptBin = 0; ptBin < nPtBins[2]; ptBin++)
//	{
//		kaon_1_ypts[ptBin]  = kaon_1_ypts[ptBin]  / kaonv2_fit->Eval( kaon_2_xpts[ptBin] );
//		kaon_2_ypts[ptBin]  = kaon_2_ypts[ptBin]  / kaonv2_fit->Eval( kaon_2_xpts[ptBin] );
//		kaon_1_Eypts[ptBin] = kaon_1_Eypts[ptBin] / kaonv2_fit->Eval( kaon_2_xpts[ptBin] );
//		kaon_2_Eypts[ptBin] = kaon_2_Eypts[ptBin] / kaonv2_fit->Eval( kaon_2_xpts[ptBin] );
//	}
//
//	for(int ptBin = 0; ptBin < nPtBins[3]; ptBin++)
//	{
//		prot_1_ypts[ptBin]  = prot_1_ypts[ptBin]  / protv2_fit->Eval( prot_2_xpts[ptBin] );
//		prot_2_ypts[ptBin]  = prot_2_ypts[ptBin]  / protv2_fit->Eval( prot_2_xpts[ptBin] );
//		prot_1_Eypts[ptBin] = prot_1_Eypts[ptBin] / protv2_fit->Eval( prot_2_xpts[ptBin] );
//		prot_2_Eypts[ptBin] = prot_2_Eypts[ptBin] / protv2_fit->Eval( prot_2_xpts[ptBin] );
//	}


	for(int ptBin = 0; ptBin < nPtBins[0]; ptBin++)
	{
			  
			  std::cout << "cpar1_y: " << cpar_1_ypts[ptBin] << std::endl;
			  std::cout << "cpar2_y: " << cpar_2_ypts[ptBin] << std::endl;
	}

	double *cpar_diff = new double[nPtBins[0]];
	double *pion_diff = new double[nPtBins[1]];
	double *kaon_diff = new double[nPtBins[2]];
	double *prot_diff = new double[nPtBins[3]];

	double *cpar_reldiff = new double[nPtBins[0]];
	double *pion_reldiff = new double[nPtBins[1]];
	double *kaon_reldiff = new double[nPtBins[2]];
	double *prot_reldiff = new double[nPtBins[3]];

	double *cpar_reldiff_E = new double[nPtBins[0]];
	double *pion_reldiff_E = new double[nPtBins[1]];
	double *kaon_reldiff_E = new double[nPtBins[2]];
	double *prot_reldiff_E = new double[nPtBins[3]];

	double *cpar_1_fitnorm = new double[nPtBins[0]];
	double *pion_1_fitnorm = new double[nPtBins[1]];
	double *kaon_1_fitnorm = new double[nPtBins[2]];
	double *prot_1_fitnorm = new double[nPtBins[3]];

	double *cpar_2_fitnorm = new double[nPtBins[0]];
	double *pion_2_fitnorm = new double[nPtBins[1]];
	double *kaon_2_fitnorm = new double[nPtBins[2]];
	double *prot_2_fitnorm = new double[nPtBins[3]];

	for(int ptBin = 0; ptBin < nPtBins[0]; ptBin++)
	{
		cpar_diff[ptBin] = cpar_2_ypts[ptBin] - cpar_1_ypts[ptBin];
		cpar_reldiff[ptBin]   =   cpar_diff[ptBin] / cpar_1_ypts[ptBin];

		cpar_1_relE[ptBin] = cpar_1_Eypts[ptBin] / cpar_1_ypts[ptBin];
		cpar_2_relE[ptBin] = cpar_2_Eypts[ptBin] / cpar_2_ypts[ptBin];

		cpar_reldiff_E[ptBin] = (cpar_2_ypts[ptBin] / cpar_1_ypts[ptBin]) * sqrt( cpar_1_relE[ptBin] * cpar_1_relE[ptBin] + cpar_2_relE[ptBin] * cpar_2_relE[ptBin] );

		cpar_1_fitnorm[ptBin] = cpar_1_ypts[ptBin]  / cparv2_fit->Eval( cpar_2_xpts[ptBin] );
		cpar_2_fitnorm[ptBin] = cpar_2_ypts[ptBin]  / cparv2_fit->Eval( cpar_2_xpts[ptBin] );

		cpar_2_ypts[ptBin]  = cpar_2_ypts[ptBin]  / cpar_1_ypts[ptBin];
		cpar_1_Eypts[ptBin] = cpar_1_Eypts[ptBin] / cpar_1_ypts[ptBin];
		cpar_2_Eypts[ptBin] = cpar_2_Eypts[ptBin] / cpar_1_ypts[ptBin];
		cpar_1_ypts[ptBin]  = cpar_1_ypts[ptBin]  / cpar_1_ypts[ptBin];

	}

	for(int ptBin = 0; ptBin < nPtBins[1]; ptBin++)
	{
		pion_diff[ptBin] = pion_2_ypts[ptBin] - pion_1_ypts[ptBin];
		pion_reldiff[ptBin] =   pion_diff[ptBin] / pion_1_ypts[ptBin];

		cpar_1_relE[ptBin] = cpar_1_Eypts[ptBin] / cpar_1_ypts[ptBin];
		pion_2_relE[ptBin] = pion_2_Eypts[ptBin] / pion_2_ypts[ptBin];

		pion_reldiff_E[ptBin] = (pion_2_ypts[ptBin] / pion_1_ypts[ptBin]) * sqrt( pion_1_relE[ptBin] * pion_1_relE[ptBin] + pion_2_relE[ptBin] * pion_2_relE[ptBin] );

		pion_1_fitnorm[ptBin] = pion_1_ypts[ptBin]  / pionv2_fit->Eval( pion_2_xpts[ptBin] );
		pion_2_fitnorm[ptBin] = pion_2_ypts[ptBin]  / pionv2_fit->Eval( pion_2_xpts[ptBin] );

		pion_2_ypts[ptBin]  = pion_2_ypts[ptBin]  / pion_1_ypts[ptBin];
		pion_1_Eypts[ptBin] = pion_1_Eypts[ptBin] / pion_1_ypts[ptBin];
		pion_2_Eypts[ptBin] = pion_2_Eypts[ptBin] / pion_1_ypts[ptBin];
		pion_1_ypts[ptBin]  = pion_1_ypts[ptBin]  / pion_1_ypts[ptBin];

		std::cerr << "pion_1_relE: "    << pion_1_relE[ptBin] << std::endl;
		std::cerr << "pion_1_ypts: "    << pion_1_ypts[ptBin] << std::endl;
		std::cerr << "pion_1_relE: "    << pion_1_relE[ptBin] << std::endl;
		std::cerr << "pion_reldiffE: "  << pion_reldiff_E[ptBin] << std::endl;
	}

	for(int ptBin = 0; ptBin < nPtBins[2]; ptBin++)
	{
		kaon_diff[ptBin] = kaon_2_ypts[ptBin] - kaon_1_ypts[ptBin];
		kaon_reldiff[ptBin] =   kaon_diff[ptBin] / kaon_1_ypts[ptBin];

		kaon_1_relE[ptBin] = kaon_1_Eypts[ptBin] / kaon_1_ypts[ptBin];
		kaon_2_relE[ptBin] = kaon_2_Eypts[ptBin] / kaon_2_ypts[ptBin];

		kaon_reldiff_E[ptBin] = (kaon_2_ypts[ptBin] / kaon_1_ypts[ptBin]) * sqrt( kaon_1_relE[ptBin] * kaon_1_relE[ptBin] + kaon_2_relE[ptBin] * kaon_2_relE[ptBin] );

		kaon_1_fitnorm[ptBin] = kaon_1_ypts[ptBin]  / kaonv2_fit->Eval( kaon_2_xpts[ptBin] );
		kaon_2_fitnorm[ptBin] = kaon_2_ypts[ptBin]  / kaonv2_fit->Eval( kaon_2_xpts[ptBin] );

		kaon_2_ypts[ptBin]  = kaon_2_ypts[ptBin]  / kaon_1_ypts[ptBin];
		kaon_1_Eypts[ptBin] = kaon_1_Eypts[ptBin] / kaon_1_ypts[ptBin];
		kaon_2_Eypts[ptBin] = kaon_2_Eypts[ptBin] / kaon_1_ypts[ptBin];
		kaon_1_ypts[ptBin]  = kaon_1_ypts[ptBin]  / kaon_1_ypts[ptBin];
	}

	for(int ptBin = 0; ptBin < nPtBins[3]; ptBin++)
	{
		prot_diff[ptBin] = prot_2_ypts[ptBin] - prot_1_ypts[ptBin];
		prot_reldiff[ptBin] =   prot_diff[ptBin] / prot_1_ypts[ptBin];

		prot_1_relE[ptBin] = prot_1_Eypts[ptBin] / prot_1_ypts[ptBin];
		prot_2_relE[ptBin] = prot_2_Eypts[ptBin] / prot_2_ypts[ptBin];

		prot_reldiff_E[ptBin] = (prot_2_ypts[ptBin] / prot_1_ypts[ptBin]) * sqrt( prot_1_relE[ptBin] * prot_1_relE[ptBin] + prot_2_relE[ptBin] * prot_2_relE[ptBin] );

		prot_1_fitnorm[ptBin] = prot_1_ypts[ptBin]  / protv2_fit->Eval( prot_2_xpts[ptBin] );
		prot_2_fitnorm[ptBin] = prot_2_ypts[ptBin]  / protv2_fit->Eval( prot_2_xpts[ptBin] );

		// Debuggg 
		std::cerr << "prot_1_relE: "    << prot_1_relE[ptBin] << std::endl;
		std::cerr << "prot_1_ypts: "    << prot_1_ypts[ptBin] << std::endl;
		std::cerr << "prot_1_relE: "    << prot_1_relE[ptBin] << std::endl;
		std::cerr << "prot_reldiffE: "  << prot_reldiff_E[ptBin] << std::endl;

		prot_2_ypts[ptBin]  = prot_2_ypts[ptBin]  / prot_1_ypts[ptBin];
		prot_1_Eypts[ptBin] = prot_1_Eypts[ptBin] / prot_1_ypts[ptBin];
		prot_2_Eypts[ptBin] = prot_2_Eypts[ptBin] / prot_1_ypts[ptBin];
		prot_1_ypts[ptBin]  = prot_1_ypts[ptBin]  / prot_1_ypts[ptBin];
	}


	TCanvas canvas_rv2_vs_pT ("rv2 vs pT", "v_{2} values as a function of p_{T}; p_{T} [GeV/c];v_{2}", 1024, 768);

	canvas_rv2_vs_pT.SetLeftMargin (0.10);
   canvas_rv2_vs_pT.SetBottomMargin(0.12);
   canvas_rv2_vs_pT.SetRightMargin(0.05);
   canvas_rv2_vs_pT.SetTopMargin  (0.05);
	
	double rv2vspt_ptmin = 0.0;
	double rv2vspt_ptmax = 2.5;
	double rv2vspt_rmin  = 0.7;
	double rv2vspt_rmax  = 1.3;

	TGraphErrors cpar_1_norm = cpar_v2(nPtBins[0], cpar_2_xpts, cpar_1_fitnorm, 0, cpar_1_Eypts );
	TGraphErrors pion_1_norm = pion_v2(nPtBins[1], pion_2_xpts, pion_1_fitnorm, 0, pion_1_Eypts );
	TGraphErrors kaon_1_norm = kaon_v2(nPtBins[2], kaon_2_xpts, kaon_1_fitnorm, 0, kaon_1_Eypts );
	TGraphErrors prot_1_norm = prot_v2(nPtBins[3], prot_2_xpts, prot_1_fitnorm, 0, prot_1_Eypts );
	TGraphErrors cpar_2_norm = cpar_v2(nPtBins[0], cpar_2_xpts, cpar_2_fitnorm, 0, cpar_2_Eypts );
	TGraphErrors pion_2_norm = pion_v2(nPtBins[1], pion_2_xpts, pion_2_fitnorm, 0, pion_2_Eypts );
	TGraphErrors kaon_2_norm = kaon_v2(nPtBins[2], kaon_2_xpts, kaon_2_fitnorm, 0, kaon_2_Eypts );
	TGraphErrors prot_2_norm = prot_v2(nPtBins[3], prot_2_xpts, prot_2_fitnorm, 0, prot_2_Eypts );

	// Debuggg 
	std::cerr << "cpar_2_Eypts: " << cpar_2_Eypts[0] << std::endl;

	cpar_2_norm.SetTitle("");
   cpar_2_norm.GetXaxis()->SetLimits(rv2vspt_ptmin, rv2vspt_ptmax);
	cpar_2_norm.GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
	cpar_2_norm.GetXaxis()->CenterTitle(1);
	cpar_2_norm.GetXaxis()->SetTitleOffset(1.2);
	cpar_2_norm.GetXaxis()->SetTitleSize(figuretextsize);
	cpar_2_norm.GetYaxis()->SetRangeUser(rv2vspt_rmin,rv2vspt_rmax);
	cpar_2_norm.GetYaxis()->SetLimits(rv2vspt_rmin,rv2vspt_rmax);
	cpar_2_norm.GetYaxis()->SetTitle("v_{2} ratio ");
	cpar_2_norm.GetYaxis()->CenterTitle(1);
	cpar_2_norm.GetYaxis()->SetTitleOffset(1.2);
	cpar_2_norm.GetYaxis()->SetTitleSize(figuretextsize);

	label = Form("ratiov2vspt_MC_comprarison_nTrk_%03d-%03d", mult1, mult2); 
	pngfigure = label+".png";
	pdffigure = label+".pdf";

	cpar_2_norm.SetMarkerStyle(25);
	pion_2_norm.SetMarkerStyle(26);
	kaon_2_norm.SetMarkerStyle(30);
	prot_2_norm.SetMarkerStyle(24);
	cpar_1_norm.SetMarkerStyle(21);
	pion_1_norm.SetMarkerStyle(22);
	kaon_1_norm.SetMarkerStyle(29);
	prot_1_norm.SetMarkerStyle(20);

	tCMSsystemlabel.Draw();

	cpar_2_norm.Draw("AP");
	pion_2_norm.Draw("P");
	kaon_2_norm.Draw("P");
	prot_2_norm.Draw("P");
	cpar_1_norm.Draw("P");
	pion_1_norm.Draw("P");
	kaon_1_norm.Draw("P");
	prot_1_norm.Draw("P");

   double up = 1.05;
   double lo = 0.95;

	TLine *line_up  = new TLine(rv2vspt_ptmin,up,rv2vspt_ptmax,up); 
	TLine *line_mi  = new TLine(rv2vspt_ptmin,1.00,rv2vspt_ptmax,1.00); 
	TLine *line_lo  = new TLine(rv2vspt_ptmin,lo,rv2vspt_ptmax,lo); 

	TF1* funcr_1 = new TF1("func_1", "1", 0, 2);
	TF1* funcr_2 = new TF1("func_2", "1", 0, 2);

	funcr_1->SetLineColor(kBlack);
	funcr_2->SetLineColor(kBlack);
	funcr_1->SetMarkerStyle(21);
	funcr_2->SetMarkerStyle(25);

	double grr_legend_x1=0.6;
	double grr_legend_y1=0.25;
	double grr_legend_x2=grr_legend_x1+.10;
	double grr_legend_y2=grr_legend_y1+.10;

	TLegend grr_legend (grr_legend_x1, grr_legend_y1, grr_legend_x2, grr_legend_y2);
	grr_legend.SetFillStyle(0);
	grr_legend.SetBorderSize(0);
	grr_legend.AddEntry(funcr_1,Form("%s", label1.c_str()) , "P");
   grr_legend.AddEntry(funcr_2,Form("%s", label2.c_str()) , "P");
	grr_legend.SetTextSize(figuretextsize);
	grr_legend.Draw("SAME");

	line_up->SetLineStyle(7);
	line_mi->SetLineStyle(7);
	line_lo->SetLineStyle(7);

	line_up->Draw("SAME");
	line_mi->Draw("SAME");
	line_lo->Draw("SAME");

	canvas_rv2_vs_pT.SaveAs( pngfigure.c_str() );
	canvas_rv2_vs_pT.SaveAs( pdffigure.c_str() );

	// Anna

	TCanvas canvas_reldiff ("rv2 vs pT", "v_{2} values as a function of p_{T}; p_{T} [GeV/c];v_{2}", 1024, 768);

	canvas_reldiff.SetLeftMargin (0.10);
   canvas_reldiff.SetBottomMargin(0.12);
   canvas_reldiff.SetRightMargin(0.05);
   canvas_reldiff.SetTopMargin  (0.05);
	
	double reldiffvspt_ptmin = 0.0;
	double reldiffvspt_ptmax = 2.5;
	double reldiffvspt_rmin  = -0.5;
	double reldiffvspt_rmax  =  0.5;

	TGraphErrors cpar_reldiffgraph = cpar_v2( nPtBins[0], cpar_2_xpts, cpar_reldiff, 0, cpar_reldiff_E );
	TGraphErrors pion_reldiffgraph = pion_v2( nPtBins[1], pion_2_xpts, pion_reldiff, 0, pion_reldiff_E );
	TGraphErrors kaon_reldiffgraph = kaon_v2( nPtBins[2], kaon_2_xpts, kaon_reldiff, 0, kaon_reldiff_E );
	TGraphErrors prot_reldiffgraph = prot_v2( nPtBins[3], prot_2_xpts, prot_reldiff, 0, prot_reldiff_E );

	cpar_reldiffgraph.SetTitle("");
   cpar_reldiffgraph.GetXaxis()->SetLimits(reldiffvspt_ptmin, reldiffvspt_ptmax);
	cpar_reldiffgraph.GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
	cpar_reldiffgraph.GetXaxis()->CenterTitle(1);
	cpar_reldiffgraph.GetXaxis()->SetTitleOffset(1.2);
	cpar_reldiffgraph.GetXaxis()->SetTitleSize(figuretextsize);
	cpar_reldiffgraph.GetYaxis()->SetRangeUser(reldiffvspt_rmin,reldiffvspt_rmax);
	cpar_reldiffgraph.GetYaxis()->SetLimits(reldiffvspt_rmin,reldiffvspt_rmax);
	cpar_reldiffgraph.GetYaxis()->SetTitle("rel. diff");
	cpar_reldiffgraph.GetYaxis()->CenterTitle(1);
	cpar_reldiffgraph.GetYaxis()->SetTitleOffset(1.2);
	cpar_reldiffgraph.GetYaxis()->SetTitleSize(figuretextsize);

	std::string label_reldiff = Form("reldiff_MC_comprarison_nTrk_%03d-%03d", mult1, mult2); 
	std::string pngfigure_reldiff = label_reldiff+".png";
	std::string pdffigure_reldiff = label_reldiff+".pdf";

	cpar_reldiffgraph.SetMarkerStyle(25);
	pion_reldiffgraph.SetMarkerStyle(26);
	kaon_reldiffgraph.SetMarkerStyle(30);
	prot_reldiffgraph.SetMarkerStyle(24);

	cpar_reldiffgraph.Draw("AP");
	pion_reldiffgraph.Draw("P");
	kaon_reldiffgraph.Draw("P");
	prot_reldiffgraph.Draw("P");

   double reldiff_min = -0.05;
   double reldiff_max =  0.05;

	TLine *line_upp  = new TLine(reldiffvspt_ptmin,reldiff_min,reldiffvspt_ptmax,reldiff_min); 
	TLine *line_mid  = new TLine(reldiffvspt_ptmin,0,reldiffvspt_ptmax,0);
	TLine *line_low  = new TLine(reldiffvspt_ptmin,reldiff_max,reldiffvspt_ptmax,reldiff_max); 

	line_upp->SetLineStyle(7);
	line_mid->SetLineStyle(7);
	line_low->SetLineStyle(7);

	line_upp->Draw("SAME");
	line_upp->Draw("SAME");
	line_low->Draw("SAME");

	tCMSsystemlabel.Draw();

	canvas_reldiff.SaveAs( pngfigure_reldiff.c_str() );
	canvas_reldiff.SaveAs( pdffigure_reldiff.c_str() );

}
