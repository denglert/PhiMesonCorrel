#include "GraphTools.h"

void plotTH2D (TH2D *histo, const char titlelabels[], const char figbasename[], const char drawmode[] )
{

	gStyle->SetOptStat(0);

	TCanvas canvas ("canvas","", canvas_res_x, canvas_res_y);

	canvas.SetLeftMargin(   canvas_margin_left   );
   canvas.SetBottomMargin( canvas_margin_bottom );
   canvas.SetRightMargin(  canvas_margin_right  );
   canvas.SetTopMargin(    canvas_margin_top    );

	histo->SetTitle( titlelabels );

	histo->GetXaxis()->SetTitleSize( xtitle_size );
	histo->GetYaxis()->SetTitleSize( ytitle_size );

	histo->GetXaxis()->SetTitleOffset( xtitle_offset );
	histo->GetYaxis()->SetTitleOffset( ytitle_offset );

	std::string outPNG = Form("%s.png", figbasename);
	std::string outPDF = Form("%s.pdf", figbasename);

	if ( ! (strcmp( drawmode, "-1") == 0 ) )
	{

	histo->Draw( drawmode );
 	canvas.SaveAs( outPNG.c_str() );
	canvas.SaveAs( outPDF.c_str() ); 

	}

}


TCanvas *GetCanvas ()
{
	TCanvas *canvas = new TCanvas ("canvas", "canvas", canvas_res_x, canvas_res_y);

	canvas->SetLeftMargin  (canvas_margin_left  );
   canvas->SetBottomMargin(canvas_margin_bottom);
   canvas->SetRightMargin (canvas_margin_right );
   canvas->SetTopMargin   (canvas_margin_top   );

	return canvas;
}


///////////////////////////////
// cpar_v2
TGraphErrors *phim_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors *graph = new TGraphErrors(npoint, pt, v2, pt_Error, v2_Error);
	graph->SetTitle("");
	graph->SetLineColor(kMagenta);
	graph->SetMarkerColor(kMagenta);
	// full circle
	graph->SetMarkerStyle(21);
	graph->SetMarkerSize(1);
   graph->GetXaxis()->SetLimits(v2vspt_ptmin,v2vspt_ptmax);
	graph->GetXaxis()->SetTitle("p_{T} [GeV/c]");                              			  
	graph->GetXaxis()->CenterTitle(1);
	graph->GetXaxis()->SetTitleOffset(1.2);
	graph->GetXaxis()->SetTitleSize(figuretextsize);
	graph->GetYaxis()->SetRangeUser(v2vspt_v2min,v2vspt_v2max);
	graph->GetYaxis()->SetTitle("v_{2}");
	graph->GetYaxis()->CenterTitle(1);
	graph->GetYaxis()->SetTitleOffset(1.2);
	graph->GetYaxis()->SetTitleSize(figuretextsize);
	return graph;
};

TGraphErrors *phim_v2_syst (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors *graph = new TGraphErrors (npoint, pt, v2, pt_Error, v2_Error);
	graph->SetTitle("");
	graph->SetFillStyle(0);
	graph->SetMarkerSize(0);
	graph->SetLineColor(kMagenta);
	graph->SetLineWidth(2);
	return graph;
};


///////////////////////////////
// cpar_v2
TGraphErrors cpar_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kBlack);
	graph.SetMarkerColor(kBlack);
	// full circle
	graph.SetMarkerStyle(21);
	graph.SetMarkerSize(1);
	return graph;
};

TGraphErrors cpar_v2_syst (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetFillStyle(0);
	graph.SetMarkerSize(0);
	graph.SetLineColor(16);
	graph.SetLineWidth(2);
	return graph;
};

//////////////////////////////////
// pion_v2
TGraphErrors pion_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{
	// pion
	TGraphErrors graph(npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kRed);
	graph.SetMarkerColor(kRed);
	// full square
	graph.SetMarkerStyle(22);
	graph.SetMarkerSize(1.5);
	return graph;
};

TGraphErrors pion_v2_syst (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetFillStyle(0);
	graph.SetMarkerSize(0);
	graph.SetLineColor(46);
	graph.SetLineWidth(2);
	return graph;
};

//////////////////////////////////
// kaon_v2
TGraphErrors kaon_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{
	// kaon
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kGreen);
	graph.SetMarkerColor(kGreen);
	graph.SetMarkerStyle(29);
	graph.SetMarkerSize(2);
	return graph;
};

TGraphErrors kaon_v2_syst (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetFillStyle(0);
	graph.SetMarkerSize(0);
	graph.SetLineColor(8);
	graph.SetLineWidth(2);
	return graph;
};

//////////////////////////////////
// prot_v2
TGraphErrors prot_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{
	// prot
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kBlue);
	graph.SetMarkerColor(kBlue);
	graph.SetMarkerStyle(20);
	graph.SetMarkerSize(1);
	return graph;
};

TGraphErrors prot_v2_syst (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetFillStyle(0);
	graph.SetMarkerSize(0);
	graph.SetLineColor(38);
	graph.SetLineWidth(2);
	return graph;
};


///////////////////////////////
// cpar_self_v2
TGraphErrors cpar_self_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{ 
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kBlack);
	graph.SetMarkerColor(kBlack);
	graph.SetMarkerStyle(24);
	return graph;
};

//////////////////////////////////
// pion_self_v2
TGraphErrors pion_self_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{
	// pion
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kRed);
	graph.SetMarkerColor(kRed);
	graph.SetMarkerStyle(25);
	return graph;
};

//////////////////////////////////
// kaon_self_v2
TGraphErrors kaon_self_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{
	// kaon
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kGreen);
	graph.SetMarkerColor(kGreen);
	graph.SetMarkerStyle(25);
	return graph;
};

//////////////////////////////////
// prot_self_v2
TGraphErrors prot_self_v2 (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error )
{
	// prot
	TGraphErrors graph (npoint, pt, v2, pt_Error, v2_Error);
	graph.SetLineColor(kBlue);
	graph.SetMarkerColor(kBlue);
	graph.SetMarkerStyle(25);
	return graph;
}

///////////////////////
// label_CMS_pPb
TLatex label_CMS_pPb(double posx, double posy, double figuretextsize)
{
	std::string CMSsystemlabel = "#splitline{CMS (work in progress) pPb}{#sqrt{s_{NN}} = 5.02 TeV L_{int} = 35 nb^{-1}}";

	TLatex tCMSsystemlabel( posx, posy, CMSsystemlabel.c_str()); 
	tCMSsystemlabel.SetTextSize(figuretextsize);
	tCMSsystemlabel.SetNDC(kTRUE);

	return tCMSsystemlabel;
}

/////////////////
// label_multBin
TLatex label_multBin(double posx, double posy, double figuretextsize, int multBin)
{

	int mult1 = multiplicity_Ana(multBin, 0, nMultiplicityBins_Ana_HDR);
	int mult2 = multiplicity_Ana(multBin, 1, nMultiplicityBins_Ana_HDR);

	std::string binlabel = Form("%d #leq N_{trk}^{offline} #leq %d", mult1, mult2);

	TLatex tlabel( posx, posy, binlabel.c_str()); 
	tlabel.SetTextSize(figuretextsize);
	tlabel.SetNDC(kTRUE);

	return tlabel;
}

TLatex label_multBin(double posx, double posy, double figuretextsize, int minNtrk, int maxNtrk, const char pre[])
{

	std::string binlabel = Form("%s %d #leq N_{trk}^{offline} #leq %d", pre, minNtrk, maxNtrk);

	TLatex tlabel( posx, posy, binlabel.c_str()); 
	tlabel.SetTextSize(figuretextsize);
	tlabel.SetNDC(kTRUE);

	return tlabel;
}

///////////////////////
// label_Ntrk_pt
TLatex label_Ntrk_pt(double posx, double posy, double figuretextsize, int TypBin, int multBin, int ptBin)
{

	double pt1 = pt(TypBin, ptBin, 0);
	double pt2 = pt(TypBin, ptBin, 1);
	int mult1 = multiplicity_Ana(multBin, 0, nMultiplicityBins_Ana_HDR);
	int mult2 = multiplicity_Ana(multBin, 1, nMultiplicityBins_Ana_HDR);

	std::string binlabel  = Form("#splitline{#splitline{  %d #leq N_{trk}^{offline} #leq %d}{ %.1f < p_{T}^{trig}  < %.1f GeV/c}}{ %.1f < p_{T}^{assoc} < %.1f GeV/c}", mult1, mult2, pt1, pt2, ptref1, ptref2);

	TLatex tlabel( posx, posy, binlabel.c_str()); 
	tlabel.SetTextSize(figuretextsize);
	tlabel.SetNDC(kTRUE);

	return tlabel;
}

///////////////////////
// label_CorrTyp
TLatex label_CorrTyp(double posx, double posy, double figuretextsize, int TypBin)
{

	std::string rightlabel = Form("%s - charged", particletype(TypBin).c_str());

	std::string binlabel = Form("%s - charged", particletype(TypBin).c_str());

	TLatex tlabel( posx, posy, binlabel.c_str()); 
	tlabel.SetTextSize(figuretextsize);
	tlabel.SetNDC(kTRUE);

	return tlabel;
}


#
