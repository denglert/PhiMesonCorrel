#ifndef SPECTRUMUTILS_H
#define SPECTRUMUTILS_H
#include <TH1D.h>
#include "AnalysisBinning.h"
#include <TLatex.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMatrixDSym.h>
#include <TMath.h>
#include "AnalysisFW.h"
#include "GraphTools.h"

// Spectrum parameters
// 0.5 MeV binwidth
extern const int nspectrumBins;
extern const double spectrumMin;
extern const double spectrumMax;
extern const double spectrumbinWidth;

// Lower sideband window
extern const double sidebandLowMassMin;
extern const double sidebandLowMassMax;

// Candidate window
extern const double candidateMassMin;
extern const double candidateMassMax;

// Higher sideband window
extern const double sidebandHighMassMin;
extern const double sidebandHighMassMax;

////////////////////////////////////////////
// Spectrum and yield analysis parameters //
////////////////////////////////////////////
	
extern const double fitmin;
extern const double fitmax;
 
extern const double uppersidebandmin;
extern const double uppersidebandmax;
 
extern const double sidebandmin;
extern const double sidebandmax;
 
extern const double mass_init ;
extern const double gamma_init;
 
// Rejection region
extern const	double rejectmin;
extern const	double rejectmax;	

extern const bool reject;

///////////////////////////////////////
// Spectrum INFO extracted from data //
///////////////////////////////////////

struct spectruminfo
{
   double mass;
	double mass_error;

	double gamma;
	double gamma_error;
	
	double yield;
	double yield_error;

	double sidebandLowEntries;
	double sidebandLowEntries_error;

	double sidebandHighEntries;
	double sidebandHighEntries_error;

	double candidateEntries;
	double candidateEntries_error;

	double signalViaAriEntries;
	double signalViaFitEntries;

	double signalViaAriEntries_error;
	double signalViaFitEntries_error;

	double backgrViaAriEntries;
	double backgrViaFitEntries;

	double backgrViaFitEntries_error;
	double backgrViaAriEntries_error;

	double sigtobkgrViaAri;
	double sigtobkgrViaFit;

	double sigtobkgrViaAri_error;
	double sigtobkgrViaFit_error;
};


/////////////////////////////////
// SpectrumUtils Functions
double   poly2nd         (double *x, double *par);
Double_t bkgr            (Double_t *x, Double_t *par);
Double_t lorentzianPeak  (Double_t *x, Double_t *par);
Double_t fitFunction     (Double_t *x, Double_t *par);
void     saveSpectrumIMG (TH1D* spectrum, std::string figurebasename);
void     fitSpectrum     (TH1D* spectrum, spectruminfo& spectrinfo, std::string figurebasename, std::string fitlog, std::string label, bool debug);

#endif
