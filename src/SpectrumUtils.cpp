#include "SpectrumUtils.h"

// Spectrum parameters
// 0.5 MeV binwidth
const int nspectrumBins = 140;
const double spectrumMin = 0.985;
const double spectrumMax = 1.05;
const double spectrumbinWidth = (spectrumMax-spectrumMin)/double(nspectrumBins);
 
// Lower sideband window
const double sidebandLowMassMin	= massbins[0][0];
const double sidebandLowMassMax	= massbins[0][1];

// Candidate window
const double candidateMassMin		= massbins[1][0];
const double candidateMassMax		= massbins[1][1];

// Higher sideband window
const double sidebandHighMassMin = massbins[2][0];
const double sidebandHighMassMax = massbins[2][1];

////////////////////////////////////////////
// Spectrum and yield analysis parameters //
////////////////////////////////////////////
	
const double fitmin = 1.000;
const double fitmax = 1.035;
 
const double uppersidebandmin = 1.030;
const double uppersidebandmax = 1.050;

const double sidebandmin = 1.000;
const double sidebandmax = 1.035;
 
const double sidebandfitmin = 0.998;
const double sidebandfitmax = 1.050;
 
const double mass_init  = 1.020;
const double gamma_init = 0.006;
 
// Rejection region
const	double rejectmin = 1.01;
const	double rejectmax = 1.032;	

const bool reject = true;

/////////////////////////////////////////////////////////////
// C type functions used in the spectrum fit

double poly2nd(double *x, double *par)
{

   if (reject && x[0] > rejectmin && x[0] < rejectmax)
	{
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

//////////////////////////

Double_t bkgr(Double_t *x, Double_t *par)
{ return par[0] + par[1]*x[0] + par[2]*x[0]*x[0]; };

//////////////////////////

Double_t lorentzianPeak(Double_t *x, Double_t *par)
{ return (0.5*spectrumbinWidth*par[0]*par[1]/TMath::Pi()) /  TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])  + .25*par[1]*par[1]);	};

//////////////////////////

Double_t fitFunction(Double_t *x, Double_t *par)
{ return bkgr(x,par) + lorentzianPeak(x,&par[3]); };

/////////////////////////////////////////////////////////////
// saveSpectrumIMG
// - exports spectrum as .png and .pdf
void saveSpectrumIMG(TH1D* spectrum, std::string figurebasename)
{
	spectrum->GetXaxis()->SetTitle("m [GeV/c^2]");
	spectrum->GetXaxis()->SetTitleOffset(1.3);
	spectrum->GetYaxis()->SetTitle("Entries"); 
	spectrum->GetYaxis()->SetTitleOffset(1.3);
   spectrum->SetMinimum(0);


   TCanvas *c = new TCanvas("c","c", 1024, 768);
   spectrum->Draw();

	std::string pngfigure = figurebasename+".png";
	std::string pdffigure = figurebasename+".pdf";
   c->SaveAs(pngfigure.c_str());
   c->SaveAs(pdffigure.c_str());
}

/////////////////////////////////////////////////////////////
// fitspectrum
// - used to fit the spectrum

void fitSpectrum(TH1D* spectrum, spectruminfo& spectrinfo, std::string figurebasename, std::string fitlog, std::string label, bool debug)
{

	// need to call AnalysisBinning.h first
	
	// Reroute fit results from stdout to file
	freopen ( fitlog.c_str(), "a", stdout);

	// Calculating number of entries and their uncertainty in signal region
	
	int sidebandLowMassMinBin 		= spectrum->GetXaxis()->FindBin(sidebandLowMassMin);	
	int sidebandLowMassMaxBin 		= spectrum->GetXaxis()->FindBin(sidebandLowMassMax);	

	int candidateMassMinBin 		= spectrum->GetXaxis()->FindBin(candidateMassMin);	
	int candidateMassMaxBin 		= spectrum->GetXaxis()->FindBin(candidateMassMax);	

	int sidebandHighMassMinBin 	= spectrum->GetXaxis()->FindBin(sidebandHighMassMin);	
	int sidebandHighMassMaxBin 	= spectrum->GetXaxis()->FindBin(sidebandHighMassMax);	
	if(sidebandHighMassMaxBin == (nspectrumBins+1)) { sidebandHighMassMaxBin=nspectrumBins;};

	// Integrals
	spectrinfo.sidebandLowEntries  = spectrum->Integral(  sidebandLowMassMinBin,  sidebandLowMassMaxBin);
	spectrinfo.candidateEntries 	 = spectrum->Integral(    candidateMassMinBin,    candidateMassMaxBin);
	spectrinfo.sidebandHighEntries = spectrum->Integral( sidebandHighMassMinBin, sidebandHighMassMaxBin);

	spectrinfo.sidebandLowEntries_error  = sqrt(spectrinfo.sidebandLowEntries); 
	spectrinfo.candidateEntries_error 	 = sqrt(spectrinfo.candidateEntries); 
	spectrinfo.sidebandHighEntries_error = sqrt(spectrinfo.sidebandHighEntries);

	spectrinfo.backgrViaAriEntries		 = ( spectrinfo.sidebandLowEntries + spectrinfo.sidebandHighEntries ) / 2;
	spectrinfo.backgrViaAriEntries_error = sqrt( (  spectrinfo.sidebandLowEntries *  spectrinfo.sidebandLowEntries ) / 4 + 
						 										( spectrinfo.sidebandHighEntries * spectrinfo.sidebandHighEntries ) / 4 );

	spectrinfo.signalViaAriEntries 		 = spectrinfo.candidateEntries - spectrinfo.backgrViaAriEntries;
	spectrinfo.signalViaAriEntries_error = sqrt( (    spectrinfo.candidateEntries_error *    spectrinfo.candidateEntries_error ) +
								 							   ( spectrinfo.backgrViaAriEntries_error * spectrinfo.backgrViaAriEntries_error ) );

	spectrinfo.sigtobkgrViaAri = spectrinfo.signalViaAriEntries / spectrinfo.backgrViaAriEntries;
	
	/////////////////////////////
	// ---  Initialization --- //
	/////////////////////////////
	
	double yield_init = spectrinfo.signalViaAriEntries;

	// Sideband Function
	TF1 sidebandFcn ("sidebandFcn", poly2nd, sidebandfitmin, sidebandfitmax, 3);
	sidebandFcn.SetLineColor(kCyan);
	sidebandFcn.SetLineWidth(3);

	TF1 probeFcn1 ("probeFcn1", "[0]+[1]*x", 1.00, 1.01);
	probeFcn1.SetLineColor(kGreen);

	TF1 probeFcn2 ("probeFcn2", "[0]+[1]*x", 1.03, 1.035);
	probeFcn2.SetLineColor(kGreen);

	// Background function
	TF1 backgroundFcn ("backgroundFcn", bkgr, sidebandmin, sidebandmax, 3);
	backgroundFcn.SetLineColor(kMagenta);
	backgroundFcn.SetLineWidth(2);
	
	// Signal function	
	TF1 signalFcn ("signalFcn", lorentzianPeak, fitmin, fitmax, 3);
	signalFcn.SetLineColor(kBlack);
	signalFcn.SetLineWidth(4);

	// Complete function
	TF1 fitFcn    ("fitFcn", fitFunction, fitmin, fitmax, 6);
	fitFcn.SetNpx(500);
	fitFcn.SetLineWidth(4);
	fitFcn.SetLineColor(kBlue);

	//////////////////////
	// *** FIT AREA *** //
	//////////////////////
	
//	// - Step 1: Fit the upper sideband region
//	double x1 = sidebandHighMassMin;
//	double x2 = sidebandHighMassMax;
//	double y1 = spectrum->GetBinContent( sidebandHighMassMinBin );
//	double y2 = spectrum->GetBinContent( sidebandHighMassMaxBin );
//	double slope = (y2-y1)/(x2-x1);
//	uppersidebandFcn.FixParameter(2, 0);
//	uppersidebandFcn.FixParameter(1, slope);
//	uppersidebandFcn.FixParameter(0, -78000);
//	spectrum->Fit( &uppersidebandFcn, "0");
//
//	// Debuggg 
//	std::cerr << "uppersidebandFcn.GetParameter(0): " << uppersidebandFcn.GetParameter(0) << std::endl;
//	std::cerr << "uppersidebandFcn.GetParameter(1): " << uppersidebandFcn.GetParameter(1) << std::endl;
//	std::cerr << "uppersidebandFcn.GetParameter(2): " << uppersidebandFcn.GetParameter(2) << std::endl;

	// - Step 1: Fit the sideband region
	TFitResultPtr r = spectrum->Fit( &sidebandFcn, "0RS");

	TMatrixDSym cov = r->GetCovarianceMatrix();
	const double *covarr = cov.GetMatrixArray();
	const double *params = r->GetParams();

	spectrum->Fit( &probeFcn1, "0R");
	spectrum->Fit( &probeFcn2, "0R");

	// - Step 2: Pass the parameters of the background interpolated from the sideband
	// 			 and the initial estimations for the signal
	fitFcn.SetParameter(0, sidebandFcn.GetParameter(0) );
	fitFcn.SetParameter(1, sidebandFcn.GetParameter(1) );
	fitFcn.SetParameter(2, sidebandFcn.GetParameter(2) );
	fitFcn.SetParameter(3, yield_init);
	fitFcn.SetParameter(4, gamma_init);
	fitFcn.SetParameter(5, mass_init);

	fitFcn.FixParameter(0, sidebandFcn.GetParameter(0));
	fitFcn.FixParameter(1, sidebandFcn.GetParameter(1));
	fitFcn.FixParameter(2, sidebandFcn.GetParameter(2));

	// - Step 3: Fit with free signal parameters
	spectrum->Fit( &fitFcn, "R");

//	fitFcn.ReleaseParameter(0);
//	fitFcn.ReleaseParameter(1);
//	fitFcn.ReleaseParameter(2);
//	
//	// - Step 4: Fit with all free parameters
//	spectrum->Fit( &fitFcn, "R");

	// Passing the fit parameters to the signal and background
	Double_t par[6];
	fitFcn.GetParameters(par);

	backgroundFcn.SetParameter(0, sidebandFcn.GetParameter(0));
	backgroundFcn.SetParameter(1, sidebandFcn.GetParameter(1));
	backgroundFcn.SetParameter(2, sidebandFcn.GetParameter(2));

	backgroundFcn.SetParameter(0, par[0]);
	backgroundFcn.SetParameter(1, par[1]);
	backgroundFcn.SetParameter(2, par[2]);

	signalFcn.SetParameter(0, par[3]);
	signalFcn.SetParameter(1, par[4]);
	signalFcn.SetParameter(2, par[5]);

	spectrinfo.mass		   = fitFcn.GetParameter(5);
	spectrinfo.mass_error   = fitFcn.GetParError(5);
	
	spectrinfo.gamma		   = fitFcn.GetParameter(4);
	spectrinfo.gamma_error  = fitFcn.GetParError(4);

	spectrinfo.yield		   = fitFcn.GetParameter(3);
	spectrinfo.yield_error  = fitFcn.GetParError(3);
	
	// Debuggg 
	std::cerr << "sidebandHighMassMax = " << sidebandHighMassMax << std::endl;
	std::cerr << "sidebandHighMassMaxBin = " << sidebandHighMassMaxBin << std::endl;
	std::cerr << "spectrinfo.candidateEntries = " << spectrinfo.candidateEntries << std::endl;
	std::cerr << "spectrinfo.sidebandLowEntries = " << spectrinfo.sidebandLowEntries << std::endl;
	std::cerr << "spectrinfo.sidebandHighEntries = " << spectrinfo.sidebandHighEntries << std::endl;
	std::cerr << "spectrinfo.backgrViaAriEntries = " << spectrinfo.backgrViaAriEntries << std::endl;

	// Fit integrals
	spectrinfo.backgrViaFitEntries		 = backgroundFcn.Integral(massbins[1][0], massbins[1][1])/spectrumbinWidth ;
	spectrinfo.backgrViaFitEntries_error = backgroundFcn.IntegralError(massbins[1][0], massbins[1][1], params, covarr )/spectrumbinWidth;


//	spectrinfo.backgrViaFitEntries_error = backgroundFcn.IntegralError(massbins[1][0], massbins[1][1]);
	spectrinfo.signalViaFitEntries 		 = signalFcn.Integral(massbins[1][0], massbins[1][1] )/spectrumbinWidth;
	spectrinfo.signalViaFitEntries_error = spectrinfo.backgrViaFitEntries_error;

	double sig_rel_unc = spectrinfo.signalViaFitEntries_error/spectrinfo.signalViaFitEntries;
	double bkg_rel_unc = spectrinfo.backgrViaFitEntries_error/spectrinfo.backgrViaFitEntries;

	spectrinfo.sigtobkgrViaFit       = spectrinfo.signalViaFitEntries / spectrinfo.backgrViaFitEntries;
	spectrinfo.sigtobkgrViaFit_error = spectrinfo.sigtobkgrViaFit * sqrt( sig_rel_unc*sig_rel_unc + bkg_rel_unc*bkg_rel_unc );

	// --- Create canvas --- //
	TCanvas canvas_spectrum("canvas_spectrum","",800,600);

	canvas_spectrum.SetLeftMargin(   canvas_margin_left   );
   canvas_spectrum.SetBottomMargin( canvas_margin_bottom );
   canvas_spectrum.SetRightMargin(  canvas_margin_right  );
   canvas_spectrum.SetTopMargin(    canvas_margin_top    );

	// No stat box plz. thx.
	gStyle->SetOptStat(0);

	// ---- Spectrum cosmetics ---- //
	spectrum->SetStats(false);
	spectrum->GetXaxis()->SetTitle("m [GeV/c^{2}]");
	spectrum->GetXaxis()->SetTitleOffset(2.5);
	spectrum->GetYaxis()->SetTitle("Entries"); 
	spectrum->GetYaxis()->SetTitleOffset(1.5);
   spectrum->SetMinimum(0);

	// --- Vertical lines indicating the region of integral --- //
	double candidateMassMinBinHeight = spectrum->GetBinContent(candidateMassMinBin);

	// Candidate region
	TLine candidateMassMinline(candidateMassMin,0,candidateMassMin, 0.2*candidateMassMinBinHeight);
	TLine candidateMassMaxline(candidateMassMax,0,candidateMassMax, 0.2*candidateMassMinBinHeight);

	// Sideband regions
	TLine sidebandLowMassMinline(sidebandLowMassMin,0,sidebandLowMassMin, 0.2*candidateMassMinBinHeight);
	TLine sidebandLowMassMaxline(sidebandLowMassMax,0,sidebandLowMassMax, 0.2*candidateMassMinBinHeight);

	TLine sidebandHighMassMinline(sidebandHighMassMin,0,sidebandHighMassMin, 0.2*candidateMassMinBinHeight);
	TLine sidebandHighMassMaxline(sidebandHighMassMax-0.0001,0,sidebandHighMassMax-0.0001, 0.2*candidateMassMinBinHeight);

	double fitinfo_upperleftcornerpositionx = 0.15;
	double fitinfo_upperleftcornerpositiony = 0.88;
	double shift = 0.04;

	// Write FIT info to canvas 
	TLatex tlabel( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony, label.c_str()); 
	tlabel.SetTextSize(0.030);
   tlabel.SetNDC(kTRUE);
	TLatex tmass( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony-2*shift,      Form("m = %.4f #pm %.4f GeV/c^{2}", spectrinfo.mass, spectrinfo.mass_error) ); 
	tmass.SetTextSize(0.030);
   tmass.SetNDC(kTRUE);
	TLatex tgamma( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony-3*shift,    Form("#Gamma   = %.2f #pm %.2f GeV/c^{2}", (spectrinfo.gamma*1000), (1000*spectrinfo.gamma_error)) ); 
	tgamma.SetTextSize(0.030);
   tgamma.SetNDC(kTRUE);

	TLatex tyield( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony-4*shift,     Form("Yield = %.0f #pm %.2f ", spectrinfo.signalViaFitEntries, spectrinfo.signalViaFitEntries_error) ); 
	tyield.SetTextSize(0.030);
   tyield.SetNDC(kTRUE);

	TLatex tsigtobkgrViaAri( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony-5*shift,Form("S/B_{sb} = %.2f", spectrinfo.sigtobkgrViaAri) ); 
	tsigtobkgrViaAri.SetTextSize(0.030);
   tsigtobkgrViaAri.SetNDC(kTRUE);

	TLatex tsigtobkgrViaFit( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony-6*shift,Form("S/B_{fit} = %.2f", spectrinfo.sigtobkgrViaFit) ); 
	tsigtobkgrViaFit.SetTextSize(0.030);
   tsigtobkgrViaFit.SetNDC(kTRUE);


	TLatex tcandentries( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony-7*shift, Form("cand. entries  = %f", spectrinfo.candidateEntries) ); 
	tcandentries.SetTextSize(0.030);
   tcandentries.SetNDC(kTRUE);

	TLatex tprobeFcn1 (0.6, 0.3, Form("#splitline{probeFcn1: a+b*x}{a = %.2f, b = %.2f}", probeFcn1.GetParameter(0), probeFcn1.GetParameter(1)) );
	tprobeFcn1.SetTextSize(0.03);
	tprobeFcn1.SetNDC(kTRUE);

	TLatex tfullyield( fitinfo_upperleftcornerpositionx,fitinfo_upperleftcornerpositiony-8*shift,     Form("Full yield = %.0f #pm %.0f ", spectrinfo.yield, spectrinfo.yield_error) ); 
	tfullyield.SetTextSize(0.030);
   tfullyield.SetNDC(kTRUE);

	spectrum->Draw();

	candidateMassMinline.Draw("SAME");
	candidateMassMaxline.Draw("SAME");
	sidebandLowMassMinline.Draw("SAME");
	sidebandLowMassMaxline.Draw("SAME");
	sidebandHighMassMinline.Draw("SAME");
	sidebandHighMassMaxline.Draw("SAME");
	tlabel.Draw("SAME");
	tmass.Draw("SAME");
	tgamma.Draw("SAME");
	tyield.Draw("SAME");
	if (debug)
	{
		sidebandFcn.Draw("SAME");
		probeFcn1.Draw("SAME");
		probeFcn2.Draw("SAME");
		tcandentries.Draw("SAME");
		tsigtobkgrViaAri.Draw("SAME");
		tsigtobkgrViaFit.Draw("SAME");
		tprobeFcn1.Draw("SAME");
	}

	backgroundFcn.Draw("SAME");
	signalFcn.Draw("SAME");

	fclose (stdout);
	std::string pngfigure;
	std::string pdffigure;

	pngfigure = figurebasename+".png";
	pdffigure = figurebasename+".pdf";

	if (debug)
	{
		pngfigure = figurebasename+"_debug.png";
		pdffigure = figurebasename+"_debug.pdf";
	}

	canvas_spectrum.SaveAs( pngfigure.c_str() );
	canvas_spectrum.SaveAs( pdffigure.c_str() );

}
