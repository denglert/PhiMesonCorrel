#include <string>

//////////////////////////
// Correlation types bins
extern const int nCorrTyp_;
extern const int nPtBinsMax_;

///////////////////////////////////////////////
// pT classes
extern const int nPtBins_[];

extern const float trigptbins[][16][2];

extern const float ptref1;
extern const float ptref2;


extern const float assoptmin;
extern const float assoptmax;

extern const double massbins[][2];

///////////////////////////////////////////////
// Multiplicity bins
extern const int nMultiplicityBins_Ana_HDR;
extern const int nMultiplicityBins_EvM_HDR;
extern const int multiplicitybins[][2];

////////////////////////////////////////
// zvtx bins
extern const int nZvtxBins_ ;
extern const float zvtxbins[][2];

//////////////////////
// TH2D binning
extern const double dEtaMin;
extern const double dEtaMax;

extern const double dEtaMin_plot;
extern const double dEtaMax_plot;
 
extern const double dPhiMin;
extern const double dPhiMax;
 
extern const int ndEtaBins;
extern const int ndPhiBins;
 
extern const double dEta_binWidth;
extern const double dPhi_binWidth;
 
extern const int negdEtaCut1Bin_;
extern const int negdEtaCut2Bin_;
extern const int posdEtaCut1Bin_;
extern const int posdEtaCut2Bin_;
 
extern const double negdEtaCut1;
extern const double negdEtaCut2;
extern const double posdEtaCut1;
extern const double posdEtaCut2;

/////////////////////////////////////
// Binning Utility functions
extern double pt(int TypBin, short int bin, int a);
extern short int ptbin(int TypBin, float pt);
extern int multiplicity_Ana(int bin, int a, int nbins);
extern int multiplicitybin_Ana(int nTrk, int bins);
extern int multiplicity_EvM(int bin, int a);
extern int multiplicitybin_EvM(int nTrk);
extern double zvtx(int bin, int a, int nbins);
extern int zvtxbin(double zvtx, int bins);
extern short int massbin(double mass);


std::string particletype (int ID);
std::string particletypelabel (int ID);
