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
#include "GraphStyles.h"

int main( int argc, const char *argv[] )
{ 

  if(argc != 3)
  {
    std::cerr << "Usage: MC_CrosscheckFig" << std::endl;
	 exit(1);
  }

 TString inpFilename_RECO = argv[1];
 TString inpFilename_GENE = argv[2];

 // Binning
 int nCorrTyp			  = nCorrTyp_; 
 int *nPtBins = new int[nCorrTyp_];

 for(int TypBin = 0; TypBin < nCorrTyp; TypBin++)
 { nPtBins[TypBin] = nPtBins_[TypBin]; }

 int nMultiplicityBins_Ana = nMultiplicityBins_Ana_HDR;
 int nMultiplicityBins_EvM = nMultiplicityBins_EvM_HDR;
 int nZvtxBins 		      = nZvtxBins_; 


 ///////////////////////////////////////
 //                                   //
 // ****** Opening input files ****** //
 //                                   //
 ///////////////////////////////////////
 
 TFile *fr = new TFile(inpFilename_RECO, "READ");
 if ( fr->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename_RECO << std::endl; exit(-1);}

 TFile *fg = new TFile(inpFilename_GENE, "READ");
 if ( fg->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename_GENE << std::endl; exit(-1);}

 //////////////////////////////
 //                          //
 // ***** Initializing ***** //
 //                          //
 //////////////////////////////

 int mult1 = 0;
 int mult2 = 120;

 TGraphErrors *cpar_RECO = (TGraphErrors*)fr->Get( Form("cpar_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *pion_RECO = (TGraphErrors*)fr->Get( Form("pion_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *kaon_RECO = (TGraphErrors*)fr->Get( Form("kaon_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *prot_RECO = (TGraphErrors*)fr->Get( Form("prot_Ntrk_%03d-%03d", mult1, mult2 ) );

 TGraphErrors *cpar_GENE = (TGraphErrors*)fg->Get( Form("cpar_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *pion_GENE = (TGraphErrors*)fg->Get( Form("pion_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *kaon_GENE = (TGraphErrors*)fg->Get( Form("kaon_Ntrk_%03d-%03d", mult1, mult2 ) );
 TGraphErrors *prot_GENE = (TGraphErrors*)fg->Get( Form("prot_Ntrk_%03d-%03d", mult1, mult2 ) );

 //////////////////////
 //                  //
 // ***** Plot ***** //
 //                  //
 //////////////////////
 
	cpar_GENE->SetMarkerStyle(21);
	pion_GENE->SetMarkerStyle(22);
	kaon_GENE->SetMarkerStyle(29);
	prot_GENE->SetMarkerStyle(20);
	
	cpar_RECO->SetMarkerStyle(25);
	pion_RECO->SetMarkerStyle(26);
	kaon_RECO->SetMarkerStyle(30);
	prot_RECO->SetMarkerStyle(4);

	cpar_GENE->GetFunction("cparv2_fit")->SetBit(TF1::kNotDraw);
	pion_GENE->GetFunction("pionv2_fit")->SetBit(TF1::kNotDraw);
	kaon_GENE->GetFunction("kaonv2_fit")->SetBit(TF1::kNotDraw);
	prot_GENE->GetFunction("protv2_fit")->SetBit(TF1::kNotDraw);
	cpar_RECO->GetFunction("cparv2_fit")->SetBit(TF1::kNotDraw);
	pion_RECO->GetFunction("pionv2_fit")->SetBit(TF1::kNotDraw);
	kaon_RECO->GetFunction("kaonv2_fit")->SetBit(TF1::kNotDraw);
	prot_RECO->GetFunction("protv2_fit")->SetBit(TF1::kNotDraw);

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

	cpar_GENE->SetTitle("");
   cpar_GENE->GetXaxis()->SetLimits(v2vspt_ptmin,v2vspt_ptmax);
	cpar_GENE->GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
	cpar_GENE->GetXaxis()->CenterTitle(1);
	cpar_GENE->GetXaxis()->SetTitleOffset(1.2);
	cpar_GENE->GetXaxis()->SetTitleSize(figuretextsize);
	cpar_GENE->GetYaxis()->SetRangeUser(v2vspt_v2min,v2vspt_v2max);
	cpar_GENE->GetYaxis()->SetTitle("v_{2}");
	cpar_GENE->GetYaxis()->CenterTitle(1);
	cpar_GENE->GetYaxis()->SetTitleOffset(1.2);
	cpar_GENE->GetYaxis()->SetTitleSize(figuretextsize);


	TF1 *cparv2_fit = new TF1 ("cparv2_fit", "[0]+[1]*x+[2]*x*x",0.2,3);
	TF1 *pionv2_fit = new TF1 ("pionv2_fit", "[0]+[1]*x+[2]*x*x",0.2,0.9);
	TF1 *kaonv2_fit = new TF1 ("kaonv2_fit", "[0]+[1]*x+[2]*x*x",0.2,0.9);
	TF1 *protv2_fit = new TF1 ("protv2_fit", "[0]+[1]*x+[2]*x*x",0.2,1.5);

	cparv2_fit->SetLineColor(kBlack);
	pionv2_fit->SetLineColor(kRed);
	kaonv2_fit->SetLineColor(kGreen);
	protv2_fit->SetLineColor(kBlue);

	cparv2_fit->SetLineStyle(9);
	pionv2_fit->SetLineStyle(9);
	kaonv2_fit->SetLineStyle(9);
	protv2_fit->SetLineStyle(9);

	cpar_GENE->Fit("cparv2_fit", "R");
	pion_GENE->Fit("pionv2_fit", "R");
	kaon_GENE->Fit("kaonv2_fit", "R");
	prot_GENE->Fit("protv2_fit", "R");

	// Draw
	cpar_GENE->Draw("AP");
	pion_GENE->Draw("P");
	kaon_GENE->Draw("P");
	prot_GENE->Draw("P");

	cpar_RECO->Draw("P");
	pion_RECO->Draw("P");
	kaon_RECO->Draw("P");
	prot_RECO->Draw("P");

	// Legends, texts
	double legend_x1=.14;
	double legend_y1=0.56;
	double legend_x2=legend_x1+.20;
	double legend_y2=legend_y1+.20;

	double gr_legend_x1=0.8;
	double gr_legend_y1=0.4;
	double gr_legend_x2=gr_legend_x1+.10;
	double gr_legend_y2=gr_legend_y1+.10;

	double CMSsystemlabelposx = 0.14;
	double CMSsystemlabelposy = 0.84;
	double multlabelposx = 0.62;
	double multlabelposy = 0.24;

	TLegend v2vsptlegend (legend_x1, legend_y1, legend_x2, legend_y2);
	v2vsptlegend.SetFillStyle(0);
	v2vsptlegend.SetBorderSize(0);
	v2vsptlegend.AddEntry(cpar_GENE,"charged", "P");
	v2vsptlegend.AddEntry(pion_GENE,"#pi", "P");
	v2vsptlegend.AddEntry(kaon_GENE,"K", "P");
	v2vsptlegend.AddEntry(prot_GENE,"p", "P");
	v2vsptlegend.SetTextSize(figuretextsize);
	v2vsptlegend.Draw("SAME");

	TF1* func_GENE = new TF1("func_GENE", "1", 0, 2);
	TF1* func_RECO = new TF1("func_RECO", "1", 0, 2);

	func_GENE->SetLineColor(kBlack);
	func_RECO->SetLineColor(kBlack);
	func_GENE->SetMarkerStyle(21);
	func_RECO->SetMarkerStyle(25);

	TLegend gr_legend (gr_legend_x1, gr_legend_y1, gr_legend_x2, gr_legend_y2);
	gr_legend.SetFillStyle(0);
	gr_legend.SetBorderSize(0);
	gr_legend.AddEntry(func_RECO,"reco", "P");
	gr_legend.AddEntry(func_GENE,"gen", "P");
	gr_legend.SetTextSize(figuretextsize);
	gr_legend.Draw("SAME");

	std::string CMSsystemlabel = Form("#splitline{CMS (work in progress)}{ pPb EPOS MC }", mult1, mult2);
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

	double *cpar_GENE_Eypts = cpar_GENE->GetEY();
	double *pion_GENE_Eypts = pion_GENE->GetEY();
	double *kaon_GENE_Eypts = kaon_GENE->GetEY();
	double *prot_GENE_Eypts = prot_GENE->GetEY();

	double *cpar_RECO_Eypts = cpar_RECO->GetEY();
	double *pion_RECO_Eypts = pion_RECO->GetEY();
	double *kaon_RECO_Eypts = kaon_RECO->GetEY();
	double *prot_RECO_Eypts = prot_RECO->GetEY();

	double *cpar_GENE_ypts = cpar_GENE->GetY();
	double *pion_GENE_ypts = pion_GENE->GetY();
	double *kaon_GENE_ypts = kaon_GENE->GetY();
	double *prot_GENE_ypts = prot_GENE->GetY();

	double *cpar_RECO_ypts = cpar_RECO->GetY();
	double *pion_RECO_ypts = pion_RECO->GetY();
	double *kaon_RECO_ypts = kaon_RECO->GetY();
	double *prot_RECO_ypts = prot_RECO->GetY();

	double *cpar_RECO_xpts = cpar_RECO->GetX();
	double *pion_RECO_xpts = pion_RECO->GetX();
	double *kaon_RECO_xpts = kaon_RECO->GetX();
	double *prot_RECO_xpts = prot_RECO->GetX();

	for(int ptBin = 0; ptBin < nPtBins[0]; ptBin++)
	{
		cpar_GENE_ypts[ptBin]  = cpar_GENE_ypts[ptBin]  / cparv2_fit->Eval( cpar_RECO_xpts[ptBin] );
		cpar_RECO_ypts[ptBin]  = cpar_RECO_ypts[ptBin]  / cparv2_fit->Eval( cpar_RECO_xpts[ptBin] );
		cpar_GENE_Eypts[ptBin] = cpar_GENE_Eypts[ptBin] / cparv2_fit->Eval( cpar_RECO_xpts[ptBin] );
		cpar_RECO_Eypts[ptBin] = cpar_RECO_Eypts[ptBin] / cparv2_fit->Eval( cpar_RECO_xpts[ptBin] );
	}

	for(int ptBin = 0; ptBin < nPtBins[1]; ptBin++)
	{
		pion_GENE_ypts[ptBin]  = pion_GENE_ypts[ptBin]  / pionv2_fit->Eval( pion_RECO_xpts[ptBin] );
		pion_RECO_ypts[ptBin]  = pion_RECO_ypts[ptBin]  / pionv2_fit->Eval( pion_RECO_xpts[ptBin] );
		pion_GENE_Eypts[ptBin] = pion_GENE_Eypts[ptBin] / pionv2_fit->Eval( pion_RECO_xpts[ptBin] );
		pion_RECO_Eypts[ptBin] = pion_RECO_Eypts[ptBin] / pionv2_fit->Eval( pion_RECO_xpts[ptBin] );
	}

	for(int ptBin = 0; ptBin < nPtBins[2]; ptBin++)
	{
		kaon_GENE_ypts[ptBin]  = kaon_GENE_ypts[ptBin]  / kaonv2_fit->Eval( kaon_RECO_xpts[ptBin] );
		kaon_RECO_ypts[ptBin]  = kaon_RECO_ypts[ptBin]  / kaonv2_fit->Eval( kaon_RECO_xpts[ptBin] );
		kaon_GENE_Eypts[ptBin] = kaon_GENE_Eypts[ptBin] / kaonv2_fit->Eval( kaon_RECO_xpts[ptBin] );
		kaon_RECO_Eypts[ptBin] = kaon_RECO_Eypts[ptBin] / kaonv2_fit->Eval( kaon_RECO_xpts[ptBin] );
	}

	for(int ptBin = 0; ptBin < nPtBins[3]; ptBin++)
	{
		prot_GENE_ypts[ptBin]  = prot_GENE_ypts[ptBin]  / protv2_fit->Eval( prot_RECO_xpts[ptBin] );
		prot_RECO_ypts[ptBin]  = prot_RECO_ypts[ptBin]  / protv2_fit->Eval( prot_RECO_xpts[ptBin] );
		prot_GENE_Eypts[ptBin] = prot_GENE_Eypts[ptBin] / protv2_fit->Eval( prot_RECO_xpts[ptBin] );
		prot_RECO_Eypts[ptBin] = prot_RECO_Eypts[ptBin] / protv2_fit->Eval( prot_RECO_xpts[ptBin] );
	}


	TCanvas canvas_rv2_vs_pT ("rv2 vs pT", "v_{2} values as a function of p_{T}; p_{T} [GeV/c];v_{2}", 1024, 768);

	canvas_rv2_vs_pT.SetLeftMargin(0.10);
   canvas_rv2_vs_pT.SetBottomMargin(0.12);
   canvas_rv2_vs_pT.SetRightMargin(0.05);
   canvas_rv2_vs_pT.SetTopMargin(0.05);
	
	double rv2vspt_ptmin = 0.0;
	double rv2vspt_ptmax = 2.5;
	double rv2vspt_rmin  = 1.2;
	double rv2vspt_rmax  = 0.8;

	TGraphErrors cpar_GENE_norm = cpar_v2(nPtBins[0], cpar_RECO_xpts, cpar_GENE_ypts, 0, cpar_GENE_Eypts );
	TGraphErrors pion_GENE_norm = pion_v2(nPtBins[1], pion_RECO_xpts, pion_GENE_ypts, 0, pion_GENE_Eypts );
	TGraphErrors kaon_GENE_norm = kaon_v2(nPtBins[2], kaon_RECO_xpts, kaon_GENE_ypts, 0, kaon_GENE_Eypts );
	TGraphErrors prot_GENE_norm = prot_v2(nPtBins[3], prot_RECO_xpts, prot_GENE_ypts, 0, prot_GENE_Eypts );
	TGraphErrors cpar_RECO_norm = cpar_v2(nPtBins[0], cpar_RECO_xpts, cpar_RECO_ypts, 0, cpar_RECO_Eypts );
	TGraphErrors pion_RECO_norm = pion_v2(nPtBins[1], pion_RECO_xpts, pion_RECO_ypts, 0, pion_RECO_Eypts );
	TGraphErrors kaon_RECO_norm = kaon_v2(nPtBins[2], kaon_RECO_xpts, kaon_RECO_ypts, 0, kaon_RECO_Eypts );
	TGraphErrors prot_RECO_norm = prot_v2(nPtBins[3], prot_RECO_xpts, prot_RECO_ypts, 0, prot_RECO_Eypts );

	// Debuggg 
	std::cerr << "cpar_RECO_Eypts: " << cpar_RECO_Eypts[0] << std::endl;

	cpar_RECO_norm.SetTitle("");
   cpar_RECO_norm.GetXaxis()->SetLimits(rv2vspt_ptmin,rv2vspt_ptmax);
	cpar_RECO_norm.GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
	cpar_RECO_norm.GetXaxis()->CenterTitle(1);
	cpar_RECO_norm.GetXaxis()->SetTitleOffset(1.2);
	cpar_RECO_norm.GetXaxis()->SetTitleSize(figuretextsize);
	cpar_RECO_norm.GetYaxis()->SetRangeUser(rv2vspt_rmin,rv2vspt_rmax);
	cpar_RECO_norm.GetYaxis()->SetTitle("v_{2} ratio reco/gene");
	cpar_RECO_norm.GetYaxis()->CenterTitle(1);
	cpar_RECO_norm.GetYaxis()->SetTitleOffset(1.2);
	cpar_RECO_norm.GetYaxis()->SetTitleSize(figuretextsize);

	label = Form("ratiov2vspt_MC_comprarison_nTrk_%03d-%03d", mult1, mult2); 
	pngfigure = label+".png";
	pdffigure = label+".pdf";

	cpar_RECO_norm.SetMarkerStyle(25);
	pion_RECO_norm.SetMarkerStyle(26);
	kaon_RECO_norm.SetMarkerStyle(30);
	prot_RECO_norm.SetMarkerStyle(24);
	cpar_GENE_norm.SetMarkerStyle(21);
	pion_GENE_norm.SetMarkerStyle(22);
	kaon_GENE_norm.SetMarkerStyle(29);
	prot_GENE_norm.SetMarkerStyle(20);

	cpar_RECO_norm.Draw("AP");
	pion_RECO_norm.Draw("P");
	kaon_RECO_norm.Draw("P");
	prot_RECO_norm.Draw("P");
	cpar_GENE_norm.Draw("P");
	pion_GENE_norm.Draw("P");
	kaon_GENE_norm.Draw("P");
	prot_GENE_norm.Draw("P");

	TLine *line_12  = new TLine(rv2vspt_ptmin,1.05,rv2vspt_ptmax,1.05); 
	TLine *line_10  = new TLine(rv2vspt_ptmin,1.0,rv2vspt_ptmax,1.0); 
	TLine *line_08  = new TLine(rv2vspt_ptmin,0.95,rv2vspt_ptmax,0.95); 

	line_12->SetLineStyle(7);
	line_10->SetLineStyle(7);
	line_08->SetLineStyle(7);

	line_12->Draw("SAME");
	line_10->Draw("SAME");
	line_08->Draw("SAME");

	canvas_rv2_vs_pT.SaveAs( pngfigure.c_str() );
	canvas_rv2_vs_pT.SaveAs( pdffigure.c_str() );


}
