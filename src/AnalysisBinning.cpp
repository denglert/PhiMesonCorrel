#include "AnalysisBinning.h"
#include <TMath.h>
#include <string>

////////////////////////////////////
// ------ Analysis binning ------ //
////////////////////////////////////

// Number of correlation types:
// 0 - chadron  - chadron
// 1 - 	pion   - chadron
// 2 - 	kaon   - chadron
// 3 - 	prot   - chadron
// 4 - non-prot - chadron

const int nCorrTyp_ = 3;

////////////////////////
// *** Pt binning *** //
////////////////////////

const int nPtBinsMax_      	=   16; 
const int nPtBins_[nCorrTyp_] = { 6, 6, 6 };
////const int nPtBins_[nCorrTyp_] = { 6, 3, 4, 6 };

const float trigptbins[nCorrTyp_][nPtBinsMax_][2] = 
{
	 {{  0.30,  0.50 },
	  {  0.50,  0.70 },
	  {  0.70,  0.90 },
	  {  0.90,  1.10 },
	  {  1.10,  1.30 },
	  {  1.30,  1.50 }},
	 {{  0.30,  0.50 },
	  {  0.50,  0.70 },
	  {  0.70,  0.90 },
	  {  0.90,  1.10 },
	  {  1.10,  1.30 },
	  {  1.30,  1.50 }},
	 {{  0.30,  0.50 },
	  {  0.50,  0.70 },
	  {  0.70,  0.90 },
	  {  0.90,  1.10 },
	  {  1.10,  1.30 },
	  {  1.30,  1.50 }}
};



// MC
//const int nPtBinsMax_      	= 12; 
//const int nPtBins_[nCorrTyp_] = { 6, 3, 3, 6 };
//
//const float trigptbins[nCorrTyp_][nPtBinsMax_][2] = 
//{
//	 {{  0.30,  0.50 },
//	  {  0.50,  1.00 },
//	  {  1.00,  1.50 },
//	  {  1.50,  2.00 },
//	  {  2.00,  2.50 },
//	  {  2.50,  3.00 },
//	  { -1.00, -1.00 }},
//	 {{  0.20,  0.40 },
//	  {  0.40,  0.60 },
//	  {  0.60,  0.80 },
//	  {  0.80,  1.00 },
//	  { -1.00, -1.00 },
//	  { -1.00, -1.00 },
//	  { -1.00, -1.00 }},
//	 {{  0.20,  0.40 },
//	  {  0.40,  0.60 },
//	  {  0.60,  0.80 },
//	  {  0.80,  0.90 },
//	  { -1.00, -1.00 },
//	  { -1.00, -1.00 },
//	  { -1.00, -1.00 }},
//	 {{  0.20,  0.40 },
//	  {  0.40,  0.60 },
//	  {  0.60,  0.80 },
//	  {  0.80,  1.00 },
//	  {  1.00,  1.20 },
//	  {  1.20,  1.40 },
//	  {  1.40,  1.60 }}
//};

const float ptref1 = 0.3;
const float ptref2 = 3.0;

const float assoptmin = 0.3;
const float assoptmax = 3.0;

//////////////////////////////////
// *** Multiplicity binning *** //
//////////////////////////////////
// Analysis multiplicity binning

////// MinBias and HighMultiplicity
//  const int nMultiplicityBins_Ana_HDR = 9;
//  const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
//  {
//  	  {   0,  30 },
//  	  {  30,  50 },
//  	  {  50,  80 },
//  	  {  80, 100 },
//  	  { 100, 120 },
//  	  { 120, 150 },
//  	  { 150, 185 },
//  	  { 185, 220 },
//  	  { 220, 260 }
// ////	  { 260, 300 }
// ////	  { 300, 350 }
//};

// MINBIAS
//const int nMultiplicityBins_Ana_HDR = 5;
//const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
//{
//	  {   0,  30 },
//	  {  30,  50 },
//	  {  50,  80 },
//	  {  80, 100 },
//	  { 100, 120 }
//};

// New MinBias
// const int nMultiplicityBins_Ana_HDR = 3;
// const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
// {
// 	  {   0,  50 },
// 	  {  50,  80 },
// 	  {  80, 120 }
// };

// MC comparison
const int nMultiplicityBins_Ana_HDR = 1;
const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
{
	  {   0,  120 }
};


// Single HighMult Bin
//const int nMultiplicityBins_Ana_HDR = 1;
//const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
//{
//	  {  120,  260 }
//};


//  HighMult Bin
//const int nMultiplicityBins_Ana_HDR = 1;
//const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
//{
//	  {  120,  260 }
//};


// New MinBias
//const int nMultiplicityBins_Ana_HDR = 2;
//const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
//{
//	  {   0,  60 },
//	  {  60, 120 }
//};

// HIGH-MULTIPLICITY
//const int nMultiplicityBins_Ana_HDR = 2;
//const int multiplicitybins_Ana[nMultiplicityBins_Ana_HDR][2] = 
//{
//	  { 120, 180 },
//	  { 180, 260 },
////	  { 260, 300 }
//////////	  { 300, 350 }
//};

//// Event Mix multiplicity binning
// HIGH-MULTIPLICITY
//const int nMultiplicityBins_EvM_HDR = 4;
//const int multiplicitybins_EvM[nMultiplicityBins_EvM_HDR][2] = 
//{
//	  { 120, 150 },
//	  { 150, 185 },
//	  { 185, 220 },
//	  { 220, 260 },
////	  { 260, 300 }
//////	  { 300, 350 }
//};

//MINBIAS
const int nMultiplicityBins_EvM_HDR = 8;
const int multiplicitybins_EvM[nMultiplicityBins_EvM_HDR][2] = 
{
	  {   0,  20 },
	  {  20,  30 },
	  {  30,  40 },
	  {  40,  50 },
	  {  50,  60 },
	  {  60,  80 },
	  {  80, 100 },
	  { 100, 120 }
};


//////////////////////////
// *** Zvtx binning *** //
//////////////////////////
// Event Mix zvtx binning
const int nZvtxBins_ = 22;
//const int nZvtxBins_ = 2;
const float zvtxbins[nZvtxBins_][2] = 
{
//	  { -15.0, -13.0 },
	  { -13.0, -11.0 },
	  { -11.0,  -9.0 },
	  {  -9.0,  -8.0 },
	  {  -8.0,  -7.0 },
	  {  -7.0,  -6.0 },
	  {  -6.0,  -5.0 },
	  {  -5.0,  -4.0 },
	  {  -4.0,  -3.0 },
	  {  -3.0,  -2.0 },
	  {  -2.0,  -1.0 },
	  {  -1.0,   0.0 },
	  {   0.0,   1.0 },
	  {   1.0,   2.0 },
	  {   2.0,   3.0 },
	  {   3.0,   4.0 },
	  {   4.0,   5.0 },
	  {   5.0,   6.0 },
	  {   6.0,   7.0 },
	  {   7.0,   8.0 },
	  {   8.0,   9.0 },
	  {   9.0,  11.0 },
	  {  11.0,  13.0 },
//	  {  13.0,  15.0 },
};

//////////////////////////
// *** Mass binning *** //
//////////////////////////
const int nMassBins_ = 3;

const double massbins[nMassBins_][2] = 
{
	  {  0.994, 1.004 },
	  {  1.014, 1.024 },
	  {  1.039, 1.049 },
};

//////////////////////////
// *** TH2D binning *** //
//////////////////////////
const double dEtaMin = -4.8;
const double dEtaMax =  4.8;

const double dEtaMin_plot = -3.0;
const double dEtaMax_plot =  2.9;

const double dPhiMin = - TMath::Pi()/2;
const double dPhiMax = 3*TMath::Pi()/2;

const int ndEtaBins = 96;
const int ndPhiBins = 96;

const double dEta_binWidth = (dEtaMax-dEtaMin)/ndEtaBins;
const double dPhi_binWidth = (dPhiMax-dPhiMin)/ndPhiBins;

// ALICE bin
// const int negdEtaCut1Bin_ = 33;
// //const int negdEtaCut1Bin_ = 19;
// const int negdEtaCut2Bin_ = 40;
// const int posdEtaCut1Bin_ = 57;
// //const int posdEtaCut2Bin_ = 78;
// const int posdEtaCut2Bin_ = 64;

// CMS default bins [-4.0 -2.0] & [2.0 - 4.0]
//const int negdEtaCut1Bin_ = 9;
//const int negdEtaCut2Bin_ = 28;
//const int posdEtaCut1Bin_ = 69;
//const int posdEtaCut2Bin_ = 88;

const int negdEtaCut1Bin_ = 12;
const int negdEtaCut2Bin_ = 28;
const int posdEtaCut1Bin_ = 69;
const int posdEtaCut2Bin_ = 85;

// Constrained bins
//const int negdEtaCut1Bin_ = 19;
//const int posdEtaCut2Bin_ = 78;

const double negdEtaCut1 = (dEtaMin                ) + (dEta_binWidth*(negdEtaCut1Bin_-1));
const double negdEtaCut2 = (dEtaMin + dEta_binWidth) + (dEta_binWidth*(negdEtaCut2Bin_-1));
const double posdEtaCut1 = (dEtaMin                ) + (dEta_binWidth*(posdEtaCut1Bin_-1));
const double posdEtaCut2 = (dEtaMin + dEta_binWidth) + (dEta_binWidth*(posdEtaCut2Bin_-1));

////////////////////////////////////////
// ***  Binning Utility functions *** //
////////////////////////////////////////

/////////////////////////////////
// --- Pt binning
// Updated binning for pt with dependence on the correlation typ
double pt(int TypBin, short int bin, int a)
{
	return trigptbins[TypBin][bin][a];
}

short int ptbin(int TypBin, float pt)
{
	int ptbin;

	// Default value, changes if a match is found in the for loop
	ptbin = -1;

	// if TypBin = 99 (unknown)
	if (TypBin == 99)
	{
		for(int ptBin = 0; ptBin < nPtBins_[0]; ptBin++)
		{ if ( (trigptbins[0][ptBin][0] <= pt) && (pt <= trigptbins[0][ptBin][1] ) ) {return ptBin;}; }
	}
	// if TypBin = 1,2,3
	else
	{
		for(int ptBin = 0; ptBin < nPtBins_[TypBin]; ptBin++)
		{ if ( (trigptbins[TypBin][ptBin][0] <= pt) && (pt <= trigptbins[TypBin][ptBin][1] ) ) {return ptBin;}; }
	}

	return ptbin;
}

/////////////////////////////////////////
// --- Multiplicity binning
// Analysis multiplicity binning
int multiplicity_Ana(int bin, int a, int nbins = 1)
{
	if ( (nbins == 1) && (a == 0) ) {return multiplicitybins_Ana[0][0];}
	if ( (nbins == 1) && (a == 1) ) {return multiplicitybins_Ana[nMultiplicityBins_Ana_HDR - 1][1];}
	else
	return multiplicitybins_Ana[bin][a];
}


int multiplicitybin_Ana(int nTrk, int bins = 1)
{
	int multiplicitybin;

	// Default value, changes if a match is found in the for loop
	multiplicitybin = -1;

	for(int multiplicityBin = 0; multiplicityBin < nMultiplicityBins_Ana_HDR; multiplicityBin++)
	{
   	if ( (multiplicitybins_Ana[multiplicityBin][0] <= nTrk) && (nTrk < multiplicitybins_Ana[multiplicityBin][1] ) ) {return multiplicityBin;};
	}

	// Single binning
	if ( (multiplicitybin != -1) && (bins == 1) ) multiplicitybin = 0;
	return multiplicitybin;
}

// Event mix multiplicity binning
int multiplicity_EvM(int bin, int a)
{
	return multiplicitybins_EvM[bin][a];
}


int multiplicitybin_EvM(int nTrk)
{
	int multiplicitybin;

	// Default value, changes if a match is found in the for loop
	multiplicitybin = -1;

	for(int multiplicityBin = 0; multiplicityBin < nMultiplicityBins_EvM_HDR; multiplicityBin++)
	{
   	if ( (multiplicitybins_EvM[multiplicityBin][0] <= nTrk) && (nTrk < multiplicitybins_EvM[multiplicityBin][1] ) ) {return multiplicityBin;};
	}

	return multiplicitybin;
}

///////////////////////////////////////////
// --- zvtx  binning
double zvtx(int bin, int a, int nbins)
{
	if ( (nbins == 1) && (a == 0) ) {return zvtxbins[0][0];}
	if ( (nbins == 1) && (a == 1) ) {return zvtxbins[nZvtxBins_ - 1][1];}
	else
	return zvtxbins[bin][a];
}

int zvtxbin(double zvtx, int bins = 1) 
{
	int zvtxbin;

	// Default value, changes if a match is found in the for loop
	zvtxbin = -1;

	for(int zvtxBin = 0; zvtxBin < nZvtxBins_; zvtxBin++)
	{
   	if ( (zvtxbins[zvtxBin][0] <= zvtx) && (zvtx <= zvtxbins[zvtxBin][1] ) ) {return zvtxBin;};
	}

	// Default to single binning if bins is not specified and the input zvtx is inside our zvtx range
	if ( (zvtxbin != -1 ) && (bins == 1) ) zvtxbin = 0;
	return zvtxbin;
}

short int massbin(double mass)
{
	for(int massBin = 0; massBin < nMassBins_; massBin++)
	{
   	if ( (massbins[massBin][0] <= mass) && (mass <= massbins[massBin][1] ) ) { return massBin; };
	}

	return -1;
}

std::string particletype (int ID)
{

	std::string out;

	switch (ID)
   { 
		case 0: out = "#Phi_{low}"; return out;
		case 1: out = "#Phi_{cand}"; return out;
		case 2: out = "#Phi_{high}"; return out;
   }
}

std::string particletypelabel (int ID)
{

	std::string out;

	switch (ID)
   { 
		case 1: out = "philow"; return out;
		case 2: out = "phican"; return out;
		case 3: out = "phihig"; return out;
   }
}
