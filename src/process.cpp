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


 if(argc != 5)
 { std::cerr << "Usage: process <.root file to be processed> <tag> <cont matrix>" << std::endl; }

 std::string inpFilename 		   = argv[1];
 std::string tag		            = argv[2];
 std::string contmatrix_filename = argv[3];
 std::string dataset_name 			= argv[4];

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
 CFW.dataset_name = dataset_name;
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
//

 CFW.fitSpectra();

 CFW.Calcvns();
 CFW.SetupTGraphs();

 CFW.RemovePoint(0);

 CFW.makeFigsig2bkgr();
 CFW.makeFigdMass( );

 // Make figures
 for(int multBin=0; multBin < nMultiplicityBins_Ana; multBin++)
 {
 	CFW.makeFigPhiv2vspT_lowavghig( multBin );
 	CFW.makeFigPhiv2vspT_control( multBin );
 }


 CFW.display_v2s();

 CFW.makeFigCorrel2D( tag );
 CFW.ReBin();
 CFW.makeFigCorrel1D( tag );


////////////////////////////////
 // ***** Setting output ***** //
 ////////////////////////////////

 CFW.Save();

 // Combined results
 std::string pkp_res = "/afs/cern.ch/work/d/denglert/public/projects/PKPCorrelation_SLC6/CMSSW_5_3_20/src/denglert/PKPCorrelationAna/results/Fiji_MinBias_HighMult_Full_trkCorr_yes_ContMatrix_yes/dump.root";
// std::string phi_res = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/Control_HighMult_Full_singlebin_trkCorr_no_temp_fit/dump.root";
// std::string phi_res = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/PIDscan_MinBias_dEdxminsweep_full_2nd_config_0_trkCorr_no_fit/dump.root";
 std::string phi_res = Form("/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/results/%s/dump.root", tag.c_str());
 
 int mult_pkp1;
 int mult_pkp2;
 int mult_pkp1_l;
 int mult_pkp2_l;
 int mult_phi1;
 int mult_phi2;

 mult_pkp1 = 100;
 mult_pkp2 = 120;
 mult_pkp1_l = mult_pkp1;
 mult_pkp2_l = mult_pkp2;
 mult_phi1 = 0;
 mult_phi2 = 120;
 CFW.makeFigCombinedResults(pkp_res.c_str(), phi_res.c_str(), mult_pkp1, mult_pkp2, mult_pkp1_l, mult_pkp2_l, mult_phi1, mult_phi2, Form("_pkp_%03d-%03d_phi_%03d-%03d", mult_pkp1, mult_pkp2, mult_phi1, mult_phi2) );
 mult_pkp1 = 120;
 mult_pkp2 = 150;
 mult_pkp1_l = mult_pkp1;
 mult_pkp2_l = mult_pkp2;
 mult_phi1 = 120;
 mult_phi2 = 180;
 CFW.makeFigCombinedResults(pkp_res.c_str(), phi_res.c_str(), mult_pkp1, mult_pkp2, mult_pkp1_l, mult_pkp2_l, mult_phi1, mult_phi2, Form("_pkp_%03d-%03d_phi_%03d-%03d", mult_pkp1, mult_pkp2, mult_phi1, mult_phi2) );
 mult_pkp1 = 185;
 mult_pkp2 = 220;
 mult_pkp1_l = 180;
 mult_pkp2_l = 220;
 mult_phi1 = 180;
 mult_phi2 = 260;
 CFW.makeFigCombinedResults(pkp_res.c_str(), phi_res.c_str(), mult_pkp1, mult_pkp2, mult_pkp1_l, mult_pkp2_l, mult_phi1, mult_phi2, Form("_pkp_%03d-%03d_phi_%03d-%03d", mult_pkp1, mult_pkp2, mult_phi1, mult_phi2) );

}
