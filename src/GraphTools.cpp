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


