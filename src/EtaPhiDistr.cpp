#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TNtuple.h>
#include <iostream>
#include <deque>
#include "AnalysisFW.h"
#include "AnalysisBinning.h"
#include "PIDUtils.h"
#include "SetupCustomTrackTree.h"

const double latextxtsize_std = 0.043;


TLatex MakeLabel ( double posx, double posy, double textsize, const char text[] )
{
 TLatex tlabel( posx, posy, text); 
 tlabel.SetTextSize(textsize);
 tlabel.SetNDC(kTRUE);
 return tlabel;
}



void plotEtaDistr ( TH1D *EtaDistr, const char figurebasename[], const char text[] )
{
 TCanvas canvas_etadistr ("EtaDistr", ";Eta;Entries", 800, 600);
 double posx = 0.15;
 double posy = 0.20;
 TLatex label = MakeLabel(posx, posy, latextxtsize_std, text);
 EtaDistr->SetTitleOffset(1.2,"X");
 EtaDistr->SetTitleOffset(1.6,"Y");
 EtaDistr->GetXaxis()->SetTitleSize(0.052);
 EtaDistr->GetYaxis()->SetTitleSize(0.052);
 EtaDistr->Draw();
 label.Draw();
 canvas_etadistr.SaveAs(Form("%s.png", figurebasename));
 canvas_etadistr.SaveAs(Form("%s.pdf", figurebasename));

}

int main( int argc, const char *argv[] )
{ 

  if(argc != 7)
  {
    std::cerr << "Usage: process <.root file to be preprocessed> <dotrkCorr> <trkCorrFilename> <tag> <nEvents>" << std::endl;
	 exit(1);
  }

 TString inpFilename     = argv[1];
 TString dotrkCorr_str 	 = argv[2];
 TString trkCorrFilename = argv[3];
 std::string PIDconfig   = argv[4];
 std::string tag		    = argv[5];
 int nEvMax 	  		    = atoi( argv[6] );

 // Binning
 int nCorrTyp			  = nCorrTyp_; 
 int *nPtBins = new int[nCorrTyp_];

 for(int TypBin = 0; TypBin < nCorrTyp; TypBin++)
 { nPtBins[TypBin] = nPtBins_[TypBin]; }

 int nMultiplicityBins_Ana = nMultiplicityBins_Ana_HDR;
 int nMultiplicityBins_EvM = nMultiplicityBins_EvM_HDR;
 int nZvtxBins 		      = nZvtxBins_; 

 PIDUtil *pidutil = new PIDUtil();
 pidutil->ReadInConfig( PIDconfig );

 bool dotrkCorr;
      if( dotrkCorr_str == "yes") { dotrkCorr = true;}
 else if( dotrkCorr_str == "no") { dotrkCorr = false;}
 else {std::cerr << "dotrkCorr not defined." << std::endl; exit(-1);}

 //////////////////////////////////////
 //                                  //
 // ****** Opening input file ****** //
 //                                  //
 //////////////////////////////////////
 
 // Open trkCorr file
 TFile *f_trkCorr = NULL;
 f_trkCorr = new TFile(trkCorrFilename, "READ");
 if ( (f_trkCorr->IsZombie() ) || f_trkCorr == NULL ) {std::cerr << "Error opening file: " << trkCorrFilename << std::endl; exit(-1);}
 else {std::cout << "trkCorr File successfully opened." << std::endl;}

 TH3D **trkEff = Read_TH3D_1Darray(f_trkCorr, "hcorr3D typ", nCorrTyp);


 for (int i = 0; i < nCorrTyp; i++)
 {
  trkEff[i]->SetDirectory(0);
 }

 f_trkCorr->Close();

 // Open file
 TFile *f = TFile::Open(inpFilename);
 // TFile *f = new TFile(inpFilename, "READ");
 if ( f->IsZombie() ) {std::cerr << "Error opening file: " << inpFilename << std::endl; exit(-1);}
 else {std::cout << "File successfully opened." << std::endl;}

 // trackTree
 // DATA - pptracks
 // MC   - ppTrack
 TTree *trackTree = (TTree*)f->Get("pptracks/trackTree");
 //TTree *trackTree = (TTree*)f->Get("ppTrack/trackTree");
 Tracks_c tTracks;
 bool isMC = true;
 setupTrackTree_c(trackTree, tTracks, isMC);

 // hiEvtAnalyzer
 TTree *EvtAna= (TTree*)f->Get("hiEvtAnalyzer/HiTree");
 int hiNtracks; EvtAna->SetBranchAddress("hiNtracks", &hiNtracks);
 float vz; EvtAna->SetBranchAddress("vz", &vz);

 // skimanalysis for event selection
 TTree *SkimAna= (TTree*)f->Get("skimanalysis/HltTree");
 // DATA  -pPAcollisionEventSelection
 // MC   - pPAcollisionEventSelectionPA
// int pPAcollisionEventSelection; SkimAna->SetBranchAddress("pPAcollisionEventSelectionPA", &pPAcollisionEventSelection);
 int pPAcollisionEventSelection; SkimAna->SetBranchAddress("pPAcollisionEventSelection", &pPAcollisionEventSelection);
 int pileUpBit;                  SkimAna->SetBranchAddress("pVertexFilterCutGplus", &pileUpBit);

 ////////////////////////////////
 //                            //
 // ***** Setting output ***** //
 //                            //
 ////////////////////////////////
 
 TFile *output = new TFile(Form("./EtaPhi_distr_%s.root", tag.c_str() ),"RECREATE");
 output->cd();

 // Debuggg 
 std::cerr << "chicks " << std::endl;

 //////////////////////////////
 //                          //
 // ***** Initializing ***** //
 //                          //
 //////////////////////////////
 
 TH1::SetDefaultSumw2( );
 TH2::SetDefaultSumw2( );

 AnalysisFW *AnaFW = new AnalysisFW();
 AnaFW->Setup();

 // dEdxvsp map
 TH2D *dEdxvsplinlintyp[nCorrTyp];
 TH2D *dEdxvsploglogtyp[nCorrTyp];

 double l10 = TMath::Log(10);

 int    npBinslog    = 100;
 int    ndEdxBinslog = 100;
 double pBins[npBinslog+1];
 double dEdxBins[ndEdxBinslog+1];

 double pminlog    = 0.1;
 double pmaxlog    = 2.0;
 double dEdxminlog = 0.2;
 double dEdxmaxlog = 50;

 double pminlin    = 0.1;
 double pmaxlin    = 2.0;
 double dEdxminlin = 0.1;
 double dEdxmaxlin = 20;

 double dplog   = ( TMath::Log(pmaxlog)   -TMath::Log(pminlog)    )/npBinslog/l10;
 double ddEdxlog = ( TMath::Log(dEdxmaxlog)-TMath::Log(dEdxminlog) )/ndEdxBinslog/l10;

 for (int i=0; i<=npBinslog; i++)
 { pBins[i] = TMath::Exp(l10*(i*dplog + TMath::Log(pminlog)/l10)); }

 for (int i=0; i<=ndEdxBinslog; i++)
 { dEdxBins[i] = TMath::Exp(l10*(i*ddEdxlog+ TMath::Log(dEdxminlog)/l10)); }

 for (int i = 0; i < nCorrTyp ; i++)
 {
	 dEdxvsplinlintyp[i] = new TH2D (Form("dEdxVsP lin-lin %s", particletype(i).c_str()) ,";p(GeV/c);dE/dx [MeV/cm]", npBins, pminlin, pmaxlin, ndEdxBins, dEdxminlin, dEdxmaxlin);
	 dEdxvsploglogtyp[i] = new TH2D (Form("dEdxVsP log-log %s", particletype(i).c_str()) ,";p(GeV/c);dE/dx [MeV/cm]", npBinslog, pBins, ndEdxBinslog, dEdxBins);
 }

 TH2D* dEdxvsplinlinall = new TH2D ("dEdxVsP lin-lin " ,";p(GeV/c);dE/dx [MeV/cm]", npBins, pminlin, pmaxlin, ndEdxBins, dEdxminlin, dEdxmaxlin);
 TH2D* dEdxvsploglogall = new TH2D ("dEdxVsP log-log " ,";p(GeV/c);dE/dx [MeV/cm]", npBinslog, pBins, ndEdxBinslog, dEdxBins);

 // Eta and (Eta,Phi) distributions

 int nParticles = nCorrTyp;

 int nEtaBins  = 100;
 int nPhiBins  = 100;
 double EtaMin = -2.4;
 double EtaMax =  2.4;
 double PhiMin = -TMath::Pi();
 double PhiMax =  TMath::Pi();

 TH2D *EtaPhiDistr_cpar = new TH2D (Form("EtaPhiDistribution id = %d, pt [%.2f - %.2f]", 0, ptref1, ptref2),";#eta;#Phi [rad]", nEtaBins, EtaMin, EtaMax, nPhiBins, PhiMin, PhiMax);
 TH1D *EtaDistr_cpar    = new TH1D (Form("EtaDistribution id = %d, pt [%.2f - %.2f]", 0, ptref1, ptref2),";#eta;Entries", nEtaBins, EtaMin, EtaMax);


 TH1D ***EtaDistr;
 TH2D ***EtaPhiDistr;

 EtaPhiDistr = new TH2D **[nParticles];
 EtaDistr    = new TH1D **[nParticles];

 for (int pid = 0; pid < nParticles; pid++)
 {
	EtaPhiDistr[pid] = new TH2D*[nPtBins[pid]];
	EtaDistr   [pid] = new TH1D*[nPtBins[pid]];

 	for (int ptBin = 0; ptBin < nPtBins[pid]; ptBin++)
	{ EtaPhiDistr[pid][ptBin]  = new TH2D (Form("EtaPhiDistribution id = %d, pt [%.2f - %.2f]", pid, pt(pid, ptBin, 0), pt(pid, ptBin, 1)),";#eta;#Phi [rad]", nEtaBins, EtaMin, EtaMax, nPhiBins, PhiMin, PhiMax); 
  	     EtaDistr[pid][ptBin]  = new TH1D (Form("EtaDistribution id = %d, pt [%.2f - %.2f]", pid, pt(pid, ptBin, 0), pt(pid, ptBin, 1)),";#eta;Entries", nEtaBins, EtaMin, EtaMax); }
 }

 // pT distribution
 const int nPtDistrBins = 15;
 const double pTMin = 0.1;
 const double pTMax = 1.6;

 TH1D *pTDistr[nParticles];
 for (int pid = 0; pid < nParticles; pid++)
 {
 	pTDistr[pid]  = new TH1D (Form("pTDistribution id = %d", pid),";pT;Entries", nPtDistrBins, pTMin, pTMax);
 }


 // vertez z distribution
 int nZvtxDistrBins = 30;
 double zVtxMin = -15.;
 double zVtxMax =  15.;
 TH1D *zVtxDistr = new TH1D ("zv",";vz;nEvents", nZvtxDistrBins, zVtxMin, zVtxMax);

 // Init finish
 std::cout << "Initialization completed." << std::endl;

 ///////////////////////////
 //                       //
 // ***** ANALYISIS ***** //
 //                       //
 ///////////////////////////
 
 if (nEvMax == -1) {nEvMax = trackTree->GetEntries();}
 for (int iEvA = 0; iEvA < nEvMax; iEvA++)
 {
	
	// Event counter info
	if ( (iEvA % 1000) == 0 )
	{ std::cout << "Event: " << iEvA << std::endl; }

	// Get current event info
	EvtAna->GetEntry(iEvA);
	SkimAna->GetEntry(iEvA);

	AnaFW->CountProcessedEvent();

	// Event Selection
	if ( !EventSelection( pPAcollisionEventSelection, pileUpBit) ) continue;
	if ( zvtxbin(vz, nZvtxBins) == -1 ) continue;

	zVtxDistr->Fill(vz);

	AnaFW->CountPassedEvent();
	AnaFW->FillnTrk(hiNtracks);

	// Load in tracks
	trackTree->GetEntry(iEvA);

	int nTrk = tTracks.nTrk;

	for (int iTrk = 0; iTrk < nTrk; iTrk++)
	{

  		float p = tTracks.trkPt[iTrk] * cosh(tTracks.trkEta[iTrk]);

		// *** Track selection *** //
		if ( !TrackSelection(tTracks, iTrk ) ) continue;

		int PID = pidutil->GetID(tTracks, iTrk);

		int ptBin_ID = ptbin( PID , tTracks.trkPt[iTrk]);
		int ptBin_CH = ptbin(   0 , tTracks.trkPt[iTrk]);

		bool isOutsideIdentifHadronPtRange = ( ptBin_ID == -1);
		bool isOutsideChargedHadronPtRange = ( ptBin_CH == -1);

		double Eta = tTracks.trkEta[iTrk];
		double Phi = tTracks.trkPhi[iTrk];
		double pt = tTracks.trkPt[iTrk];
		double w0 = trackWeight(trkEff,   0, pt, Eta, Phi, dotrkCorr);
		double w  = trackWeight(trkEff, PID, pt, Eta, Phi, dotrkCorr);

		if ( (ptref1 < tTracks.trkPt[iTrk]) && (tTracks.trkPt[iTrk] < ptref2) )
		{
			EtaPhiDistr_cpar->Fill(Eta,Phi,w0);
			EtaDistr_cpar->Fill(Eta,w0);
		}

		// chadron
		if( !isOutsideChargedHadronPtRange )
		{ 
			EtaPhiDistr[0][ptBin_CH]->Fill(Eta,Phi,w0);
		   EtaDistr[0][ptBin_CH]->Fill(Eta,w0);
			pTDistr[0]->Fill(pt, w0);
			dEdxvsplinlintyp[0]->Fill( p, tTracks.dedx[iTrk] );
			dEdxvsploglogtyp[0]->Fill( p, tTracks.dedx[iTrk] );
		}

		// pid particle
		if( (PID != 99) && !isOutsideIdentifHadronPtRange )
		{
			EtaPhiDistr[PID][ptBin_ID]->Fill(Eta,Phi, w);
		  	EtaDistr[PID][ptBin_ID]->Fill(Eta, w);
			pTDistr[PID]->Fill(pt, w);
			dEdxvsplinlintyp[PID]->Fill( p, tTracks.dedx[iTrk] );
			dEdxvsploglogtyp[PID]->Fill( p, tTracks.dedx[iTrk] );
		}

		dEdxvsploglogall->Fill( p, tTracks.dedx[iTrk] );
		dEdxvsplinlinall->Fill( p, tTracks.dedx[iTrk] );
		// *** Track selection *** //


	}

 }

 // Normalization
 for (int pid = 0; pid < nParticles; pid++)
 for (int ptBin = 0; ptBin < nPtBins[pid]; ptBin++)
 {
	EtaPhiDistr[pid][ptBin]->Scale(1./nEvMax);
	   EtaDistr[pid][ptBin]->Scale(1./nEvMax);
 }


 ////////////////////////////
 //                        //
 // **** Plot Figures **** //
 //                        //
 ////////////////////////////

 gStyle->SetOptStat(0);

 /////////////////////////
 // (Eta,Phi) distribution
 for (int pid = 0; pid < nParticles; pid++)
 for (int ptBin = 0; ptBin < nPtBins[pid]; ptBin++)
 {
	TCanvas canvas_EtaPhiDistr ("EtaPhiDistr", ";Eta;Phi", 800, 600);
 	EtaPhiDistr[pid][ptBin]->Draw("COLZ");

	canvas_EtaPhiDistr.SetRightMargin(0.15);

 	std::string EtaPhiDistrFigBase = Form("EtaPhiDistr_typ_%d_pt_%.2f-%.2f", pid, pt(pid, ptBin, 0), pid, pt(pid, ptBin, 1));

 	std::string EtaPhiDistrFigPNG = EtaPhiDistrFigBase+".png";
 	std::string EtaPhiDistrFigPDF = EtaPhiDistrFigBase+".pdf";

 	canvas_EtaPhiDistr.SaveAs(EtaPhiDistrFigPNG.c_str() );
 	canvas_EtaPhiDistr.SaveAs(EtaPhiDistrFigPDF.c_str() );
 }

 TCanvas canvas_EtaPhiDistr_cpar ("EtaPhiDistr", ";Eta;Phi", 800, 600);
 canvas_EtaPhiDistr_cpar.SetRightMargin(0.15);

 EtaPhiDistr_cpar->Draw("COLZ");
 EtaPhiDistr_cpar->Draw("COLZ");

 canvas_EtaPhiDistr_cpar.SaveAs(Form("EtaPhiDistr_typ_0_pt_%.2f-%.2f.png", ptref1, ptref2));
 canvas_EtaPhiDistr_cpar.SaveAs(Form("EtaPhiDistr_typ_0_pt_%.2f-%.2f.pdf", ptref1, ptref2));

 ////////////////////
 // Eta distribution
 for (int pid = 0; pid < nParticles; pid++)
 for (int ptBin = 0; ptBin < nPtBins[pid]; ptBin++)
 {
	TCanvas canvas_EtaDistr ("EtaDistr", ";Eta;Phi", 800, 600);
 	EtaDistr[pid][ptBin]->Draw("COLZ");

 	std::string EtaDistrFigBase = Form("EtaDistr_typ_%d_pt_%.2f-%.2f", pid, pt(pid, ptBin, 0), pid, pt(pid, ptBin, 1));

 	std::string EtaDistrFigPNG = EtaDistrFigBase+".png";
 	std::string EtaDistrFigPDF = EtaDistrFigBase+".pdf";

 	canvas_EtaDistr.SaveAs(EtaDistrFigPNG.c_str() );
 	canvas_EtaDistr.SaveAs(EtaDistrFigPDF.c_str() );
 }

 plotEtaDistr( EtaDistr_cpar, Form("EtaDistr_typ_0_pt_%.2f-%.2f", ptref1, ptref2), Form("#splitline{CMS Preliminary pPb, #sqrt{s_{NN}} = 5.02 TeV}{charged %.2f < p_{T} < %.2f GeV/c}", ptref1, ptref2) );

// TCanvas canvas_EtaDistr_cpar ("EtaDistr", ";Eta;Entries", 800, 600);
// EtaDistr_cpar->Draw("COLZ");
// canvas_EtaPhiDistr_cpar.SetTitleOffset(1.2,"X")
// canvas_EtaPhiDistr_cpar.SetTitleOffset(1.6,"Y")
//
//
// canvas_EtaDistr_cpar.SaveAs(Form("EtaDistr_typ_0_pt_%.2f-%.2f.png", ptref1, ptref2));
// canvas_EtaDistr_cpar.SaveAs(Form("EtaDistr_typ_0_pt_%.2f-%.2f.pdf", ptref1, ptref2));

 /////////////////
 // dEdx vs p plots
 for (int i = 0; i < 4; i++)
 {
 	TCanvas canvas_dEdxvsplin ("dEdx", ";p [GeV/c];dE/dx", 800, 600);

	//canvas_dEdxvsplin.GetXaxis()->SetRangeUser(0.15,5);
	//canvas_dEdxvsplin.GetYaxis()->SetRangeUser(1,20);

	canvas_dEdxvsplin.SetLogz(1);
	dEdxvsplinlintyp[i]->Draw("COLZ");

 	std::string dEdxvsplinlinFigBase = Form("dEdxvspLin_typ_%d", i);

 	std::string dEdxvspFigPNG = dEdxvsplinlinFigBase+".png";
 	std::string dEdxvspFigPDF = dEdxvsplinlinFigBase+".pdf";

 	canvas_dEdxvsplin.SaveAs(dEdxvspFigPNG.c_str() );
 	canvas_dEdxvsplin.SaveAs(dEdxvspFigPDF.c_str() );
 }

 TF1 *ptexpfit[nCorrTyp];

 for (int i = 0; i < nCorrTyp; i++)
 {
 	ptexpfit[i] = new TF1(Form("ptdistr_exp_%d", i), "[0]*exp(- [1]*(x-[2]) )", 0.4, 1.0 );
 }

 ////////////////////
 // pT distribution
 for (int pid = 0; pid < nParticles; pid++)
 {
	TCanvas canvas_pTDistr ("pTDistr", ";Eta;Phi", 800, 600);
	pTDistr[pid]->Fit( ptexpfit[pid] );
 	pTDistr[pid]->Draw("");

 	std::string pTDistrFigBase = Form("pTDistr_typ_%d", pid);

 	std::string pTDistrFigPNG = pTDistrFigBase+".png";
 	std::string pTDistrFigPDF = pTDistrFigBase+".pdf";

 	canvas_pTDistr.SaveAs(pTDistrFigPNG.c_str() );
 	canvas_pTDistr.SaveAs(pTDistrFigPDF.c_str() );
 }

 // zVtxDistr
 gStyle->SetOptStat(1);
 	TCanvas canvas_zVtxDistr ("canvas_zvtx", ";vz;nEvents", 800, 600);
	zVtxDistr->Draw();

 	std::string zVtxDistrFigBase = "zVtxDistr";

 	std::string zVtxDistrFigPNG = zVtxDistrFigBase+".png";
 	std::string zVtxDistrFigPDF = zVtxDistrFigBase+".pdf";

	canvas_zVtxDistr.SaveAs(zVtxDistrFigPNG.c_str());
	canvas_zVtxDistr.SaveAs(zVtxDistrFigPDF.c_str());

 gStyle->SetOptStat(0);
 // dEdx vs. p map
 for (int i = 0; i < 4; i++)
 {
 	TCanvas canvas_dEdxvsplog ("dEdx", ";p [GeV/c];dE/dx", 800, 600);
	canvas_dEdxvsplog.SetLogz(1);
	canvas_dEdxvsplog.SetLogx(1);
	canvas_dEdxvsplog.SetLogy(1);

	dEdxvsploglogtyp[i]->Draw("COLZ");
//	dEdxvsploglogtyp[i]->GetYaxis()->SetUserRange(0.2,50);

	std::string PID = particletype (i);
	std::string label = PID;

	TLatex tlabel ( 0.25,0.25, label.c_str() ); 
	tlabel.SetTextSize(0.062);
	tlabel.SetNDC(kTRUE);
	tlabel.Draw();

 	std::string dEdxvsploglogFigBase = Form("dEdxvspLog_typ_%d", i);

 	std::string dEdxvspFigPNG = dEdxvsploglogFigBase+".png";
 	std::string dEdxvspFigPDF = dEdxvsploglogFigBase+".pdf";

 	canvas_dEdxvsplog.SaveAs(dEdxvspFigPNG.c_str() );
 	canvas_dEdxvsplog.SaveAs(dEdxvspFigPDF.c_str() );
 }

 	std::string dEdxvsplinlinallFigBase = "dEdxvspLin_typ_all";
 	std::string dEdxvsploglogallFigBase = "dEdxvspLog_typ_all";
 	std::string dEdxvsplinFigallPNG = dEdxvsplinlinallFigBase+".png";
 	std::string dEdxvsplinFigallPDF = dEdxvsplinlinallFigBase+".pdf";
 	std::string dEdxvsplogFigallPNG = dEdxvsploglogallFigBase+".png";
 	std::string dEdxvsplogFigallPDF = dEdxvsploglogallFigBase+".pdf";

	makedEdxvspFigloglog( dEdxvsploglogall, PIDconfig, dEdxvsplogFigallPNG);
	makedEdxvspFigloglog( dEdxvsploglogall, PIDconfig ,dEdxvsplogFigallPDF);
	makedEdxvspFiglinlin( dEdxvsplinlinall, PIDconfig, dEdxvsplinFigallPNG);
	makedEdxvspFiglinlin( dEdxvsplinlinall, PIDconfig, dEdxvsplinFigallPDF);

 //////////////////////
 //                  //
 // **** OUTPUT **** //
 //                  //
 //////////////////////
 
 output->Write();
 output->Close();

}
