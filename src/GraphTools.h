#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <string.h>
#include "TGraphErrors.h"
#include "TLatex.h"
#include "AnalysisFW.h"
#include "AnalysisBinning.h"

// General Graph variables
const int canvas_res_x = 800;
const int canvas_res_y = 600;

const double canvas_margin_left   = 0.1;
const double canvas_margin_right  = 0.15;
const double canvas_margin_bottom = 0.12;
const double canvas_margin_top    = 0.05;

const double xtitle_offset = 1.3;
const double ytitle_offset = 1.1;

const double xtitle_size = 0.043;
const double ytitle_size = 0.043;

const double figuretextsize = 0.043;
const double figuretextsize_correl1D = 0.035;
const int fillstyle = 3009;

const int Colors_Mult[4] = {kBlack, kBlue, kRed, kMagenta};

//////////////
// sig2bkgr //
const double sig2bkgr_ptmin = 0.0;
const double sig2bkgr_ptmax = 2.5;
const double sig2bkgr_min = 0;
const double sig2bkgr_max = 2;

const double sig2bkgr_multlabel_x1 = 0.58;
const double sig2bkgr_multlabel_y1 = 0.36;
const double sig2bkgr_multlabel_shift = 0.06;

///////////
// dmass //
const double dmass_ptmin = 0.0;
const double dmass_ptmax = 2.5;
const double dmass_min = -5;
const double dmass_max =  5;

const double dmass_multlabel_x1 = 0.58;
const double dmass_multlabel_y1 = 0.34;
const double dmass_multlabel_shift = 0.06;

const double dmass_CMSsystemlabel_x1 = 0.14;
const double dmass_CMSsystemlabel_y1 = 0.84;

//////////////
// v2 vs pt //
const double v2vspt_ptmin = 0.0;
const double v2vspt_ptmax = 2.5;
const double v2vspt_v2min = 0.0;
const double v2vspt_v2max = 0.16;

const double v2vspt_legend_x1 = 0.12;
const double v2vspt_legend_y1 = 0.56;
const double v2vspt_legend_xw = 0.20;
const double v2vspt_legend_yw = 0.18;

const double v2vspt_CMSsystemlabel_x1 = 0.14;
const double v2vspt_CMSsystemlabel_y1 = 0.84;

const double v2vspt_multlabel_x1 = 0.58;
const double v2vspt_multlabel_y1 = 0.24;

/////////////////////
// Graph Functions //

TGraphErrors *Get_TGraph_sig2bkgr (int npoint, double *pt, double *sig2bkgr, double *pt_Error, double *sig2bkgr_Error );
TGraphErrors *Get_TGraph_dmass (int npoint, double *pt, double *dmass, double *pt_Error, double *dmass_Error );

TGraphErrors *phim_v2      (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error );
TGraphErrors *phim_v2_syst (int npoint, double *pt, double *v2, double *pt_Error, double *v2_Error );

TCanvas *GetCanvas ();

void plotTH2D (TH2D *histo, const char titlelabels[], const char figbasename[], const char drawmode[] );
TLatex label_CorrTyp(double posx, double posy, double figuretextsize, int TypBin);
TLatex *label_multBinPtr(double posx, double posy, double figuretextsize, int multBin);
TLatex label_multBin(double posx, double posy, double figuretextsize, int multBin);
TLatex label_multBin(double posx, double posy, double figuretextsize, int minNtrk, int maxNtrk, const char pre[]);
TLatex label_Ntrk_pt(double posx, double posy, double figuretextsize, int TypBin, int multBin, int ptBin);
TLatex label_CMS_pPb(double posx, double posy, double figuretextsize);

#endif
