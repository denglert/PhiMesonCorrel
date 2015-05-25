#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <iostream>
#include <TLorentzVector.h>
#include "CorrelationUtils.h"
#include "SetupCustomTrackTree.h"
#include "AnalysisFW.h"
#include "PIDUtils.h"


int main(int argc, const char *argv[])
{ 

  if(argc != 3)
  {
    std::cerr << "Usage: PIDSetup <.root file to be processed>" << std::endl;
	 exit(1);
  }

 TString inpFilename   = argv[1];
// std::string PIDconfig = argv[2];

 const int nFiles = 5;
 std::string PIDconfig[nFiles]; 
 PIDconfig[0] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/PIDUtils/config/dEdxmin_sweep/config_0";
 PIDconfig[1] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/PIDUtils/config/dEdxmin_sweep/config_1";
 PIDconfig[2] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/PIDUtils/config/dEdxmin_sweep/config_2";
 PIDconfig[3] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/PIDUtils/config/dEdxmin_sweep/config_3";
 PIDconfig[4] = "/afs/cern.ch/work/d/denglert/public/projects/PhiMesonCorrel_Nemo/CMSSW_5_3_20/src/denglert/PhiMesonCorrel/PIDUtils/config/dEdxmin_sweep/config_4";

 ///////////////////////////////////////////
 // Opening input file, setting up TTrees //
 ///////////////////////////////////////////
 // Open file
 TFile *f = TFile::Open(inpFilename);
 if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << "Successfully opened file " << inpFilename << std::endl ;}

 for(int i=0; i < nFiles; i++)
 {

 // dEdxvsP
 TH2D *dEdxvsP_lin = (TH2D*)f->Get("dEdxVsP lin-lin");
 TH2D *dEdxvsP_log = (TH2D*)f->Get("dEdxVsP log-log");
 
 // Save dEdx value map
 std::string figurename_lin = Form("dEdx_lin_%d.png", i);
 std::string figurename_log = Form("dEdx_log_%d.png", i);
 
 makedEdxvspFiglinlin(dEdxvsP_lin, PIDconfig[i], figurename_lin);
 makedEdxvspFigloglog(dEdxvsP_log, PIDconfig[i], figurename_log);

 }

}
