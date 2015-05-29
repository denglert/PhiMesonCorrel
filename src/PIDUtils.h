#ifndef PIDUTILS_H
#define PIDUTILS_H
#include <cmath>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <fstream>
#include <string>
#include "SetupCustomTrackTree.h"

const int nPIDBins = 4;

const int npBins 		= 300;
const double pMin 	= 0;
const double pMax 	= 3;
const int ndEdxBins 	= 300;
const double dEdxMin = -2;
const double dEdxMax = 30;


const int npBinslog    = 100;
const int ndEdxBinslog = 100;

const double pminlog    = 0.15;
const double pmaxlog    = 2.0;
const double dEdxminlog = 0.2;
const double dEdxmaxlog = 50;
 
const double pminlin    = 0.1;
const double pmaxlin    = 2.0;
const double dEdxminlin = 0.1;
const double dEdxmaxlin = 20;

class dEdxMaps {
public :
   dEdxMaps(const char tag[]);
   ~dEdxMaps();

	// Variables
	TH2D *dEdxvspPID[4];
	TH2D *dEdxvspAll;

	// Functions
	void Fill(int PID, double p, double dedx);
	void PlotFigs(const char tag[]);
};

class PIDUtil {
public :
	
   PIDUtil( );
   ~PIDUtil();

	// Functions
	
	float BBcurve(float *x, const float *par);
	bool  isPion(float p, float dEdx);
	bool  isKaon(float p, float dEdx);
	bool  isProt(float p, float dEdx);

	// Variables
	float BB_Pion_low_par[3];
	float BB_Pion_hig_par[3];
	float BB_Pion_mindEdxcut;
	float BB_Pion_minpcut; 
	float BB_Pion_maxpcut;
	                        
	float BB_Kaon_low_par[3];
	float BB_Kaon_hig_par[3];
	float BB_Kaon_mindEdxcut;
	float BB_Kaon_maxpcut;
	                        
	float BB_Prot_low_par[3];
	float BB_Prot_hig_par[3];
	float BB_Prot_mindEdxcut;
	float BB_Prot_maxpcut;

	float BB_NonProt_mindEdxcut;
	float BB_NonProt_maxpcut;

	std::string configfile;

	int unIDcode;
	double unIDcode_cm;
	double etaMax;
	
	void ReadInConfig( std::string PIDconfigfile_str );

	int GetID      (const Tracks &tTracks, int iTrk);
	bool IsKaon    (const Tracks &tTracks, int iTrk);
	int GetID      (const Tracks_c &tTracks, int iTrk);
	int GetIDmTrk_trkCorr(const Tracks_c &tTracks, int iTrk);
	int GetIDgenPart_trkCorr(const Tracks_c &tTracks, int iPart);
	int GetIDgenPart_trkCorr(const Particles &tTracks, int iPart);
	double GetID_cm(const Tracks_c &tTracks, int iTrk);
};

// Bethe-Bloch Curve function
double poly2nd(double *x, double *par);

extern const float BB_Pion_low_par[3];
extern const float BB_Pion_maxpcut;
extern const float BB_Pion_hig_par[3];
extern const float BB_Kaon_low_par[3];
extern const float BB_Kaon_hig_par[3];
extern const float BB_Kaon_mindEdxcut;
extern const float BB_Kaon_maxpcut;

extern const float BB_Prot_low_par[3];
extern const float BB_Prot_hig_par[3];
extern const float BB_Prot_mindEdxcut;
extern const float BB_Prot_maxpcut;

extern const double delta;

// PIDUtils parameters
extern const int nPIDBins;
extern const int npBins;
extern const double pMin;
extern const double pMax;
extern const int ndEdxBins;
extern const double dEdxMin;
extern const double dEdxMax;

extern const float asymmetry;
extern const float mindEdx;

/////////////////////////////
// PIDUtils functions

float BBcurve1c(float *x, const float *par);
double BBcurve1 (double *x, double *par);
int GetPID(float p, float dEdx, float eta);

int McPID2AnaPID ( const Particles &tTracks, int iPart);
int McPID2AnaPID ( const Tracks_c  &tTracks, int iPart);
double McPID2AnaPID_cm(const Tracks_c &tTracks, int iPart);

bool isPion(float p, float dEdx);
bool isKaon(float p, float dEdx);
bool isProt(float p, float dEdx);
void makedEdxvspFiglinlin(TH2D* dEdxvsP, std::string PIDconfig, std::string figurename);
void makedEdxvspFigloglog(TH2D* dEdxvsP, std::string PIDconfig, std::string figurename);

#endif
