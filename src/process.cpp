#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <iostream>
#include "AnalysisFW.h"
#include "AnalysisBinning.h"
#include "PIDUtils.h"
#include "CorrelationUtils.h"


int main( int argc, const char *argv[] )
{ 


 if(argc != 4)
 { std::cerr << "Usage: process <.root file to be processed> <tag> <cont matrix>" << std::endl; }

 std::string inpFilename 		   = argv[1];
 std::string tag		            = argv[2];
 std::string contmatrix_filename = argv[3];

 // Binning
 int nCorrTyp 	    = nCorrTyp_;
 int *nPtBins      = new int[nCorrTyp_];

 for(int TypBin = 0; TypBin < nCorrTyp; TypBin++)
 { nPtBins[TypBin] = nPtBins_[TypBin]; }
 int nMultiplicityBins_Ana = nMultiplicityBins_Ana_HDR; 
 int nMultiplicityBins_EvM = nMultiplicityBins_EvM_HDR; 
 int nZvtxBins 	 = nZvtxBins_; 

 prepareDIR(tag);

 //////////////////////////////////////
 // ****** Opening input file ****** //
 //////////////////////////////////////
 
// TFile *f = new TFile(inpFilename.c_str(), "READ");
// if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename.c_str() << std::endl; }
// else{ std::cout << "TFile seems to be loaded." << std::endl; };
  
 //////////////////////////////
 // ***** Initializing ***** //
 //////////////////////////////
 
 TH1::SetDefaultSumw2( );
 TH2::SetDefaultSumw2( );

 CorrelationFramework CFW(nCorrTyp, nPtBins, nMultiplicityBins_Ana, nMultiplicityBins_EvM);
 CFW.tag = tag;
 CFW.contmatrix_filename   = contmatrix_filename;
 CFW.preprocessed_filename = inpFilename;
 CFW.DoSelfCorrelation = false;
 CFW.SetupForProcess();
 // Old Read-In, now integrated with CFW.SetupForProcess();
 // CFW.ReadIn_CorrelationFuncs( f );
 CFW.DeContaminate();
 CFW.Set_dEtacut();

 ///////////////////////////
 // ***** ANALYISIS ***** //
 ///////////////////////////
 
 CFW.doAnalysis();

 CFW.ReBin();

// std::string signal = "correl2D_signal_pion_pT_0.40-0.60_nTrk_000-300.png";
// std::string backgr = "correl2D_backgr_pion_pT_0.40-0.60_nTrk_000-300.png";
// std::string functi = "correl2D_functi_pion_pT_0.40-0.60_nTrk_000-300.png";

// std::string signal_title = "sign";
// std::string backgr_title = "bkgr";
// std::string functi_title = "func";

// plot2DCorrelation_custom( CFW.correl2D_signal[1][1][0], signal,  signal_title);
// plot2DCorrelation_custom( CFW.correl2D_functi[1][1][0], functi,  functi_title);
// plot2DCorrelation_custom( CFW.correl2D_backgr[1][1][0], backgr,  backgr_title);

 CFW.Calcvns();
 CFW.display_v2s();


 // Make figures
 for(int multBin=0; multBin < nMultiplicityBins_Ana; multBin++)
 {
	CFW.makeFigv2vspT_allparticles(multBin, tag);
	CFW.makeFigv2vspT_allparticles_ALICE_comparison(multBin, tag);
	CFW.makeFigv3vspT_allparticles(multBin, tag);
//	CFW.makeFigv2vspT_allparticles_with_selfcorrelation(multBin, tag);

 }

 CFW.makeFigv2vspT_HIN13002(tag);
 CFW.makeFigv2vsnTrk_cpar_ref(tag);

 CFW.makeFigCorrel2D( tag );
 CFW.ReBin();
 CFW.makeFigCorrel1D( tag );


////////////////////////////////
 // ***** Setting output ***** //
 ////////////////////////////////

 CFW.Corr_Results->Close();



}
