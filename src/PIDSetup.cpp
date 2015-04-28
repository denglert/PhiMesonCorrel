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

  if(argc != 2)
  {
    std::cerr << "Usage: PIDSetup <.root file to be processed>" << std::endl;
	 exit(1);
  }

 TString inpFilename   = 		 argv[1];

 ///////////////////////////////////////////
 // Opening input file, setting up TTrees //
 ///////////////////////////////////////////
 
 // Open file
 TFile *f = TFile::Open(inpFilename);
 if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << "Successfully opened file " << inpFilename << std::endl ;}

 // dEdxvsP
 TH2D *dEdxvsP = (TH2D*)f->Get("dEdxVsP");
 
 // Save dEdx value map
 std::string figurename = "dEdxtest.png";

 makedEdxvspFigloglog(dEdxvsP, figurename);

}
