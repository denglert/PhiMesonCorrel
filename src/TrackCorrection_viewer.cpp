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
#include "TLegend.h"

double figuretextsize = 0.043; 

double TH1Dlabelposx = 0.63;
double TH1Dlabelposy = 0.20;

const double ptGlobalMin = 0.1;

//const int nphi = 70;
const int nphi = 20;
const double phiMin = -TMath::Pi();
const double phiMax =  TMath::Pi();


const int nCorrTyp = 5;

const int npt[nCorrTyp]      = { 15,     8,    7,   14,     6 };
const double ptMin[nCorrTyp] = { 0.3, 0.20, 0.20, 0.20,  1.00 };
const double ptMax[nCorrTyp] = { 3.0, 1.00, 0.90, 1.60,  1.60 };

const int neta[nCorrTyp]      = {  12 ,    4,    4,    4,    4 };
const double etaMin[nCorrTyp] = { -2.4, -0.8, -0.8, -0.8, -0.8 };
const double etaMax[nCorrTyp] = {  2.4,  0.8,  0.8,  0.8,  0.8 };

void plottable( TH3D **table, const char figbasename[] )
{
	gStyle->SetOptStat(0);
	
	double ptbw[nCorrTyp];
	for(int i = 0; i < nCorrTyp; i++)
	{ ptbw[i] = (ptMax[i]-ptMin[i])/npt[i]; }
	
	for (int pid = 0; pid < nCorrTyp; pid++)
	for (int ptBin = 1; ptBin < (npt[pid]+1); ptBin++)
	{
		TCanvas canvas_table ("trkTable", ";Eta;Phi", 900, 600);
	  	canvas_table.Divide(2,1);
	
	  	table[pid]->GetXaxis()->SetRange(ptBin, ptBin);
	  	TH2D *histo = (TH2D*)table[pid]->Project3D("zy");

		histo->SetTitle("");
		histo->GetXaxis()->SetTitleOffset(1.6);
		histo->GetYaxis()->SetTitleOffset(1.6);
		histo->GetZaxis()->SetTitleOffset(1.6);
	
	  	canvas_table.cd(1);

		if ( pid != 0)
		{ histo->GetXaxis()->SetRangeUser(-0.8,0.75); }

	  	histo->Draw("SURF1");
	  	canvas_table.cd(2);
	   histo->Draw("COLZ");

		double pt1 = ptMin[pid] + ptbw[pid] * (ptBin-1);
		double pt2 = ptMin[pid] + ptbw[pid] * (ptBin );
	
		std::string tableFigBase = Form("%s_typ_%d_pt_%.3f-%.3f", figbasename, pid, pt1, pt2);
		
		std::string tableFigPNG = tableFigBase+".png";
		std::string tableFigPDF = tableFigBase+".pdf";
		
		canvas_table.SaveAs(tableFigPNG.c_str() );
		canvas_table.SaveAs(tableFigPDF.c_str() );

//	  	canvas_table.cd(1);
//	  	histo->GetXaxis()->SetRangeUser(-1.0,1.0);
//	  	histo->Draw("SURF1");
//	  	canvas_table.cd(2);
//	   histo->Draw("COLZ");
//
//		tableFigBase = Form("etacut_1_%s_typ_%d_pt_%.3f-%.3f", figbasename, pid, pt1, pt2);
// 	
//		tableFigPNG = tableFigBase+".png";
//		tableFigPDF = tableFigBase+".pdf";
//		
//		canvas_table.SaveAs(tableFigPNG.c_str() );
//		canvas_table.SaveAs(tableFigPDF.c_str() );

//	  	canvas_table.cd(1);
//
//		{ histo->GetXaxis()->SetRangeUser(-0.8,0.75); }
//
//	  	histo->Draw("SURF1");
//	  	canvas_table.cd(2);
//	   histo->Draw("COLZ");
//
//		tableFigBase = Form("%s_typ_%d_pt_%.3f-%.3f_etacut_0.8", figbasename, pid, pt1, pt2);
//		
//		tableFigPNG = tableFigBase+".png";
//		tableFigPDF = tableFigBase+".pdf";
//		
//		canvas_table.SaveAs(tableFigPNG.c_str() );
//		canvas_table.SaveAs(tableFigPDF.c_str() );

//	  	canvas_table.cd(1);
//	  	histo->GetXaxis()->SetRangeUser(-1.5,1.5);
//	  	histo->Draw("SURF1");
//	  	canvas_table.cd(2);
//	   histo->Draw("COLZ");
//
//		tableFigBase = Form("etacut_1.5_%s_typ_%d_pt_%.3f-%.3f", figbasename, pid, pt1, pt2);
//		
//		tableFigPNG = tableFigBase+".png";
//		tableFigPDF = tableFigBase+".pdf";
//		
//		canvas_table.SaveAs(tableFigPNG.c_str() );
//		canvas_table.SaveAs(tableFigPDF.c_str() );

	}
}

void plotTH2D( TH2D *histo, const char figbasename[], const char label[], std::string option )
{
	TCanvas canvas_table ("trkcorrfigure", "lol", 800, 600);

  	canvas_table.Divide(2,1);


   canvas_table.SetRightMargin(0.05);

	histo->SetTitle(label);
	histo->GetXaxis()->SetTitleOffset(1.6);
	histo->GetYaxis()->SetTitleOffset(1.6);
	histo->GetZaxis()->SetTitleOffset(1.6);

	double zmin, zmax;

	if ( option == "0.0-1.0" )
	{
		std::cout << "0.0-1.0" << std::endl;
		zmin = 0.0;
		zmax = 1.0;
		histo->GetZaxis()->SetLimits(zmin, zmax);
		histo->GetZaxis()->SetRangeUser(zmin, zmax);
	}
	else if ( option == "1.0-1.1" )
	{
		std::cout << "1.0-1.1" << std::endl;
		zmin = 1.0;
		zmax = 1.1;
		histo->GetZaxis()->SetLimits(zmin, zmax);
		histo->GetZaxis()->SetRangeUser(zmin, zmax);
	}
	else if ( option == "0.0-0.1" )
	{
		std::cout << "0.0-0.1" << std::endl;
		zmin = 0.0;
		zmax = 0.1;
		histo->GetZaxis()->SetLimits(zmin, zmax);
		histo->GetZaxis()->SetRangeUser(zmin, zmax);
	}
	else if ( option == "0.9-10.0")
	{
		std::cout << "0.9-10.0" << std::endl;
		zmin = 0.9;
		zmax = 10.0;
		histo->GetZaxis()->SetLimits(zmin, zmax);
		histo->GetZaxis()->SetRangeUser(zmin, zmax);
	}
	else if ( option == "auto" ) 
	{
		std::cout << "auto" << std::endl;
	};

  	canvas_table.cd(1);

  	histo->Draw("SURF1");
  	canvas_table.cd(2);
   histo->Draw("COLZ");
	
	std::string tableFigPNG = Form("%s.png", figbasename);
	std::string tableFigPDF = Form("%s.pdf", figbasename);

	canvas_table.SaveAs(tableFigPNG.c_str() );
	canvas_table.SaveAs(tableFigPDF.c_str() );
};


void plotTH1D( TH1D *histo, const char figbasename[], const char xaxisname[], const char yaxisname[], const char label[])
{
	gStyle->SetOptStat(0);
	
	TCanvas canvas_table ("trkcorrection", Form("%s;%s", yaxisname, xaxisname), 800, 600);
	
	histo->GetXaxis()->SetTitle(xaxisname);
	histo->GetYaxis()->SetTitle(yaxisname);

	histo->SetTitle("");

	histo->GetXaxis()->SetTitleOffset(1.6);
	histo->GetYaxis()->SetTitleOffset(1.6);
	histo->GetZaxis()->SetTitleOffset(1.6);


	TLatex tlabel( TH1Dlabelposx, TH1Dlabelposy, label); 
	tlabel.SetTextSize(figuretextsize);
	tlabel.SetNDC(kTRUE);
	
	//if ( pid != 0)
	//{ histo->GetXaxis()->SetRangeUser(-0.8,0.75); }
	
	histo->Draw("");
	tlabel.Draw("SAME");
	
	std::string tableFigPNG = Form("%s.png", figbasename);
	std::string tableFigPDF = Form("%s.pdf", figbasename);
	
	canvas_table.SaveAs(tableFigPNG.c_str() );
	canvas_table.SaveAs(tableFigPDF.c_str() );
}

int main( int argc, const char *argv[] )
{ 

 if(argc != 2)
 {
   std::cerr << "Usage: TrackCorrection_viewer <trkCorr.root file to be displayed>" << std::endl;
   exit(1);
 }

 TString inpFilename     = argv[1];

 // Open file
 TFile *f = TFile::Open(inpFilename);
 if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << "File successfully opened." << std::endl;}

 TH1D *zvtxDistrDATA = (TH1D*)f->Get("zvtxDistrDATA");
 TH1D *zvtxDistrMC   = (TH1D*)f->Get("zvtxDistrMC");
 TH1D *ratiozvtx     = (TH1D*)f->Get("zvtxratio");
;

 TH3D **hgen3D         = Read_TH3D_1Darray( f, "hgen3D typ",  nCorrTyp);
 TH3D **heff3D         = Read_TH3D_1Darray( f, "heff3D part",  nCorrTyp);
 TH3D **trkCorr3D      = Read_TH3D_1Darray( f, "hcorr3D typ",  nCorrTyp);
 TH3D **hmatched3D     = Read_TH3D_1Darray( f, "hmatched3D typ",  nCorrTyp);

 TH1D **hmultrec1D = Read_TH1D_1Darray( f, "hmultrec1D typ", nCorrTyp);
 TH1D **hmatched1D = Read_TH1D_1Darray( f, "hmatched1D typ", nCorrTyp);

 TH2D **hmultrec2D = Read_TH2D_1Darray( f, "hmultrec2D typ", nCorrTyp);
 TH2D **hmatched2D = Read_TH2D_1Darray( f, "hmatched2D typ", nCorrTyp);

 TH2D **hsecondary2D   = Read_TH2D_1Darray( f, "hsecondary2D typ",   nCorrTyp);
 TH2D **hmultrecrate2D = Read_TH2D_1Darray( f, "hmultrecrate2D part", nCorrTyp);
 TH2D **hfake2D        = Read_TH2D_1Darray( f, "hfake2D typ",        nCorrTyp);

 TH1D **hsecondary1D   = Read_TH1D_1Darray( f, "hsecondary1D typ",   nCorrTyp);
 TH1D **hmultrecrate1D = Read_TH1D_1Darray( f, "hmultrecrate1D part", nCorrTyp);
 TH1D **hfake1D        = Read_TH1D_1Darray( f, "hfake1D typ",        nCorrTyp);

 // Initilaizing
 double ptbw[nCorrTyp];
 for(int i = 0; i < nCorrTyp; i++)
 { ptbw[i] = (ptMax[i]-ptMin[i])/npt[i]; 
	std::cout << Form("ptbw: %.2f \n", ptbw[i]);
 }

  double etabw[nCorrTyp]; 
  for(int i = 0; i < nCorrTyp; i++)
  { etabw[i] = (etaMax[i]-etaMin[i])/neta[i];}

  const double phiMin = -TMath::Pi();
  const double phiMax =  TMath::Pi();
  const int nphi = 35;
  const double phibw = (phiMax-phiMin)/nphi;

//  for(int i = 0; i < 1; i++)
//  for(int x = 1; x < npt[i]  +1; x++)
//  for(int y = 1; y < neta[i] +1; y++)
//  for(int z = 1; z < nphi +1; z++)
//  {
//	  double pt  = ptMin[i]+x*ptbw[i]; 
//	  double eta = etaMin[i]+y*etabw[i]; 
//	  double phi = phiMin+z*phibw; 
//
//	  std::cout << Form("pt: %.2f eta: %.2f phi: %.2f", pt, eta, phi) << std::endl;
//	  std::cout << "efff: " << heff3D[i]->GetBinContent(x,y,z) << std::endl;
//	  std::cout << "fake: " << hfake2D[i]->GetBinContent(x,y) << std::endl;
//	  std::cout << "seco: " << hsecondary2D[i]->GetBinContent(x,y) << std::endl;
//	  std::cout << "mult: " << hmultrecrate2D[i]->GetBinContent(x,y) << std::endl;
//
//   double value = (1.0-hfake2D[i]->GetBinContent(x,y))*(1.0-hsecondary2D[i]->GetBinContent(x,y)) / ((heff3D[i]->GetBinContent(x,y,z)) * (hmultrecrate2D[i]->GetBinContent(x,y)));
//
//	std::cout << "trkCorr: " << value << std::endl;
////	  hcorr3D[i]->SetBinContent(x,y,z,value);
// //	  log.wr(Form("%d %.3f %.2f %.2f : %.4f", i, pt, eta, phi, value));
//  }

	gStyle->SetOptStat(0);
	TCanvas canvas_zvtx ("zvtxdistr", ";Eta;Phi", 900, 600);

	canvas_zvtx.SetBottomMargin(0.50);

   zvtxDistrMC->SetMarkerSize(0.9);
   zvtxDistrMC->SetMarkerStyle(8);
   zvtxDistrDATA->SetMarkerSize(0.9);
   zvtxDistrDATA->SetMarkerStyle(8);

   TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   pad1->SetBottomMargin(0);
   pad1->Draw("");
   pad1->cd();
   zvtxDistrDATA->GetYaxis()->SetTitle("dN/d z_{vtx}");
   zvtxDistrDATA->SetMarkerColor(kBlue);
   zvtxDistrMC->SetMarkerColor(kRed);
   zvtxDistrDATA->Draw("P");
   zvtxDistrMC->Draw("SAMEP");

	TLegend zvtxdistr_legend (0.75, 0.6, 0.85, 0.7);
	zvtxdistr_legend.SetFillStyle(1);
	zvtxdistr_legend.SetBorderSize(1);
	zvtxdistr_legend.AddEntry(zvtxDistrDATA,"DATA", "P");
	zvtxdistr_legend.AddEntry(zvtxDistrMC,"MC", "P");
	zvtxdistr_legend.Draw("SAME");
   canvas_zvtx.cd();
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->Draw();
   pad2->cd();
	ratiozvtx->SetMarkerColor(kBlack);
	ratiozvtx->SetLineColor(kBlack);
	ratiozvtx->GetXaxis()->SetTitle("z_{vtx} [cm]");
	ratiozvtx->GetYaxis()->SetTitle("ratio");
	ratiozvtx->GetXaxis()->SetTitleOffset(1.0);
	ratiozvtx->GetYaxis()->SetTitleOffset(0.4);
	ratiozvtx->GetYaxis()->SetLabelSize(0.06);
	ratiozvtx->GetXaxis()->SetLabelSize(0.09);
	ratiozvtx->GetXaxis()->SetTitleSize(0.082);
	ratiozvtx->GetYaxis()->SetTitleSize(0.082);
   ratiozvtx->SetMarkerStyle(8);
	ratiozvtx->Draw("P");



	std::string outPDF = "zvtxdistr.pdf";
	std::string outPNG = "zvtxdistr.png";

	canvas_zvtx.SaveAs(outPDF.c_str());
	canvas_zvtx.SaveAs(outPNG.c_str());


 for(int i = 0; i < nCorrTyp; i++)
 {

 std::string type = particletype(i);

 std::string filebasename_fakerate_ptint_etadep = Form("fakerate_%d_pT_%.2f_%.2f", i, ptMin[i], ptMax[i]);
 std::string filebasename_secondar_ptint_etadep = Form("secondar_%d_pT_%.2f_%.2f", i, ptMin[i], ptMax[i]);
 std::string filebasename_multrate_ptint_etadep = Form("multrate_%d_pT_%.2f_%.2f", i, ptMin[i], ptMax[i]);
 std::string filebasename_multreco_ptint_etadep = Form("multreco_%d_pT_%.2f_%.2f", i, ptMin[i], ptMax[i]);
 std::string filebasename_matchedt_ptint_etadep = Form("matchedt_%d_pT_%.2f_%.2f", i, ptMin[i], ptMax[i]);

 std::string filebasename_fake_ptetadep = Form("fake_%d_pteta", i);
 std::string filebasename_seco_ptetadep = Form("seco_%d_pteta", i);
 std::string filebasename_mult_ptetadep = Form("mult_%d_pteta", i);

 std::string label_ptint_etadep = Form("#splitline{%s}{p_{T} = [%.2f - %.2f] GeV/c}", type.c_str(), ptMin[i], ptMax[i]);

 std::string label_fake_ptetadep = Form("%s fake", type.c_str());
 std::string label_seco_ptetadep = Form("%s seco", type.c_str());
 std::string label_mult_ptetadep = Form("%s mult", type.c_str());

 plotTH1D(hfake1D[i]       ,filebasename_fakerate_ptint_etadep.c_str(),"#eta", "fakerate"      , label_ptint_etadep.c_str());
 plotTH1D(hsecondary1D[i]  ,filebasename_secondar_ptint_etadep.c_str(),"#eta", "secondary rate", label_ptint_etadep.c_str());
 plotTH1D(hmultrecrate1D[i],filebasename_multrate_ptint_etadep.c_str(),"#eta", "multreco rate" , label_ptint_etadep.c_str());
 plotTH1D(hmultrec1D[i],    filebasename_multreco_ptint_etadep.c_str(),"#eta", "multreco"      , label_ptint_etadep.c_str());
 plotTH1D(hmatched1D[i],    filebasename_matchedt_ptint_etadep.c_str(),"#eta", "matched"       , label_ptint_etadep.c_str());

 plotTH2D(hfake2D[i]       ,filebasename_fake_ptetadep.c_str(),Form("%s fake;p_{T};#eta", type.c_str()), "0.0-0.1"  );
 plotTH2D(hsecondary2D[i]  ,filebasename_seco_ptetadep.c_str(),Form("%s seco;p_{T};#eta", type.c_str()), "0.0-0.1" );
 plotTH2D(hmultrecrate2D[i],filebasename_mult_ptetadep.c_str(),Form("%s mult;p_{T};#eta", type.c_str()), "1.0-1.1" );

 }

 for (int i = 0; i < nCorrTyp; i++)
 for (int ptBin = 1; ptBin < (npt[i]+1); ptBin++)
 {

 double pt1 = (ptMin[i] + (ptbw[i] * (ptBin-1)));
 double pt2 = (ptMin[i] + (ptbw[i] * (ptBin  )));
 
 std::cout << Form("i = %d, ptbw[i] = %.2f, ptMin[i] = %.2f, ptBin = %d", i, ptbw[i], ptMin[i], ptBin) << std::endl;
 std::cout << Form("pt = [ %.2f - %.2f ]", pt1, pt2) << std::endl;
 
 trkCorr3D[i]  ->GetXaxis()->SetRange(ptBin, ptBin);
 heff3D[i]     ->GetXaxis()->SetRange(ptBin, ptBin);
 hgen3D[i]     ->GetXaxis()->SetRange(ptBin, ptBin);
 hmatched3D[i] ->GetXaxis()->SetRange(ptBin, ptBin);

 TH2D *trkCorr  = (TH2D*)trkCorr3D  [i]->Project3D("zy");
 TH2D *heff     = (TH2D*)heff3D     [i]->Project3D("zy");
 TH2D *hgen     = (TH2D*)hgen3D     [i]->Project3D("zy");
 TH2D *hmatched = (TH2D*)hmatched3D [i]->Project3D("zy");

 std::string type = particletype(i);

 std::string filebasename_trkCorr = Form("trkCorr_%d_pT_%.2f_%.2f", i, pt1, pt2);
 std::string filebasename_eff     = Form("efficiency_%d_pT_%.2f_%.2f", i, pt1, pt2);
 std::string filebasename_gen     = Form("gen_%d_pT_%.2f_%.2f", i, pt1, pt2);
 std::string filebasename_matched = Form("matchedt_%d_pT_%.2f_%.2f", i, pt1, pt2);

 std::string label_trkCorr        = Form("#splitline{%s, trkCorrections}{p_{T} [ %.2f - %.2f]}", type.c_str(), pt1, pt2);
 std::string label_eff            = Form("#splitline{%s, efficiency}{p_{T} [ %.2f - %.2f]}", type.c_str(), pt1, pt2);
 std::string label_gen            = Form("#splitline{%s, gen}{p_{T} [ %.2f - %.2f]}", type.c_str(), pt1, pt2);
 std::string label_matched        = Form("#splitline{%s, matched}{p_{T} [ %.2f - %.2f]}", type.c_str(), pt1, pt2);

 	plotTH2D(trkCorr,filebasename_trkCorr.c_str(),Form("%s;#eta;#phi", label_trkCorr.c_str()) , "0.9-10.0" );
 	plotTH2D(heff   ,filebasename_eff.c_str()    ,Form("%s;#eta;#phi", label_eff.c_str())     , "0.0-1.0" );
 	plotTH2D(hgen   ,filebasename_gen.c_str()    ,Form("%s;#eta;#phi", label_gen.c_str())     , "auto" );
 	plotTH2D(hmatched   ,filebasename_matched.c_str(),Form("%s;#eta;#phi", label_matched.c_str()), "auto"     );
 }

}
