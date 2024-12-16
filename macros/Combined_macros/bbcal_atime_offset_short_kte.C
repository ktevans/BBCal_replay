/*
  This script Script to determine ADC time offsets for all SH and PS channels w.r.t
  BBTH cluster mean time (bbhodo.clus.tmean). One needs a configfile, similar to /cfg/atimeOff-example.cfg,
  to execute this script. To execute, do:
  ----
  [a-onl@aonl2 macros]$ pwd
  /adaqfs/home/a-onl/sbs/BBCal_replay/macros
  [a-onl@aonl2 macros]$ root -l
  root [0] .x Combined_macros/bbcal_atime_offset.C("Combined_macros/cfg/atimeOff-example.cfg")
  ----
  P. Datta  <pdbforce@jlab.org>  Created  16 Feb 2022
*/

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <TH2F.h>
#include <TGraph.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStopwatch.h>

const Double_t Mp = 0.938272081;  // +/- 6E-9 GeV

const Int_t kNcolsSH = 7;   // SH columns
const Int_t kNrowsSH = 27;  // SH rows
const Int_t kNblksSH = 189; // Total # SH blocks/PMTs
const Int_t kNcolsPS = 2;   // PS columns
const Int_t kNrowsPS = 26;  // PS rows
const Int_t kNblksPS = 52;  // Total # PS blocks/PMTs

string getDate();
void Custm1DHisto(TH1F*);
TString GetOutFileBase(TString);
void ReadOffset(TString, Double_t*);
TH1F* MakeHisto(Int_t, Int_t, char const *, Int_t, Double_t, Double_t);
std::vector<std::string> SplitString(char const delim, std::string const myStr);
void Custm2DRnumHisto(TH2F*, std::vector<std::string> const & lrnum);

TH1F *h_atime_sh[kNrowsSH][kNcolsSH];
TH1F *h_atime_sh_corr[kNrowsSH][kNcolsSH];
TH1F *h_atime_ps[kNrowsPS][kNcolsPS];
TH1F *h_atime_ps_corr[kNrowsPS][kNcolsPS];

vector<Long64_t> goodevents;
Double_t ash_atimeOffs[kNblksSH];
Double_t aps_atimeOffs[kNblksPS];

TCanvas *subCanv[8];
const Int_t kCanvSize = 100;
/*
namespace shgui {
  TGMainFrame *main = 0;
  TGHorizontalFrame *frame1 = 0;
  TGTab *fTab;
  TGLayoutHints *fL3;
  TGCompositeFrame *tf;
  TGTextButton *exitButton;
  TGNumberEntry *entryInput;
  TGLabel *runLabel;

  TRootEmbeddedCanvas *canv[8];

  TGCompositeFrame* AddTabSub(Int_t sub) {
    tf = fTab->AddTab(Form("Tab %d",sub+1));

    TGCompositeFrame *fF5 = new TGCompositeFrame(tf, (12+1)*kCanvSize,(6+1)*kCanvSize , kHorizontalFrame);
    TGLayoutHints *fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
					   kLHintsExpandY, 5, 5, 5, 5);
    TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(Form("shSubCanv%d",sub), fF5, 6*kCanvSize,8*kCanvSize);
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }

  void SetupGUI() {
    if(!main) {
      main = new TGMainFrame(gClient->GetRoot(), 1200, 1100);
      frame1 = new TGHorizontalFrame(main, 150, 20, kFixedWidth);
      runLabel = new TGLabel(frame1,"Run ");
      exitButton = new TGTextButton(frame1, "&Exit",
				    "gApplication->Terminate(0)");
      TGLayoutHints *frame1LH = new TGLayoutHints(kLHintsTop|kLHintsLeft|
						  kLHintsExpandX,2,2,2,2);
      frame1->AddFrame(runLabel,frame1LH);
      frame1->AddFrame(exitButton,frame1LH);
      main->AddFrame(frame1, new TGLayoutHints(kLHintsBottom |
					       kLHintsRight, 2, 2, 5, 1));

      // Create the tab widget
      fTab = new TGTab(main, 300, 300);
      fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

      // Create Tab1 (SH Sub1)
      for(Int_t i = 0; i < 8; i++) {
        tf = AddTabSub(i);
      }
      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
					     kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();

      for(Int_t i = 0; i < 4; i++) {
        subCanv[i] = canv[i]->GetCanvas();
	subCanv[i]->Divide(kNcolsSH,7,0.001,0.001);
      }
      for(Int_t i = 4; i < 8; i++) {
        subCanv[i] = canv[i]->GetCanvas();
	subCanv[i]->Divide(kNcolsPS,7,0.001,0.001);
      }
    }
  }
};
*/
// main function
void bbcal_atime_offset_short (char const * configfilename, bool isdebug = 1) {

  gErrorIgnoreLevel = kError; // Ignore all ROOT warnings

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  TStopwatch *sw2 = new TStopwatch();
  sw->Start(); sw2->Start();

  //gui setup
  //shgui::SetupGUI();

  TChain *C = new TChain("T");
  //creating base for outfile names
  TString cfgfilebase = GetOutFileBase(configfilename);

  TString exp = "unknown";
  Int_t config = -1;       // Experimental configuration
  TString set = "N/A";     // Needed when we have multiple calibration sets within a config
  Int_t ppass = -1;        // Replay pass to get ready for
  Double_t Ebeam = 0.;     // GeV
  Double_t atppos_nom = -40.; // ns Nominal ADC time peak position determined by the latency in FADC config file
  Double_t atppos_old = 0.;   // ns Current BBCAL ADC time peak position
  Double_t atppos_new = 0.;   // ns Desired BBCAL ADC time peak position after calibration
  Double_t hcal_atppos = 0.;  // ns HCAL ADC time peak position

  // Reading configfile
  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endRunlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }
  TCut globalcut = ""; TString gcutstr;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
      gcutstr += currentline;
    }
  }
  std::vector<std::string> gCutList = SplitString('&', gcutstr.Data());
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);
  while( currentline.ReadLine( configfile ) ){
    if( currentline.BeginsWith("#") ) continue;
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if (skey == "exp") exp = ((TObjString*)(*tokens)[1])->GetString();
      if (skey == "config") config = ((TObjString*)(*tokens)[1])->GetString().Atoi();
      if (skey == "set") set = ((TObjString*)(*tokens)[1])->GetString();
      if (skey == "pre_pass") ppass = ((TObjString*)(*tokens)[1])->GetString().Atoi();
      if (skey == "E_beam") Ebeam = ((TObjString*)(*tokens)[1])->GetString().Atof();
      if (skey == "atppos_nom") atppos_nom = ((TObjString*)(*tokens)[1])->GetString().Atof();
      if (skey == "atppos_old") atppos_old = ((TObjString*)(*tokens)[1])->GetString().Atof();
      if (skey == "atppos_new") atppos_new = ((TObjString*)(*tokens)[1])->GetString().Atof();
      if (skey == "hcal_atppos") hcal_atppos = ((TObjString*)(*tokens)[1])->GetString().Atof();
      if( skey == "*****" ){
	break;
      }
    }
    delete tokens;
  }
  atppos_nom = -1.*abs(atppos_nom); // Foolproofing: atppos_nom should always be a positive no.

  // Check for empty rootfiles and set tree branches
  if(C->GetEntries()==0) {cerr << endl << " --- No ROOT file found!! --- " << endl << endl; exit(1);}
  else cout << endl << "Found " << C->GetEntries() << " events. Starting analysis.. " << endl;

  // Setting useful ROOT tree branch addresses
  C->SetBranchStatus("*",0);
  Int_t maxtr=500;
  //shower
  Double_t sh_nclus;     C->SetBranchStatus("bb.sh.nclus",1); C->SetBranchAddress("bb.sh.nclus",&sh_nclus);
  Double_t sh_e;         C->SetBranchStatus("bb.sh.e",1); C->SetBranchAddress("bb.sh.e",&sh_e);
  Double_t sh_rowblk;    C->SetBranchStatus("bb.sh.rowblk",1); C->SetBranchAddress("bb.sh.rowblk",&sh_rowblk);
  Double_t sh_colblk;    C->SetBranchStatus("bb.sh.colblk",1); C->SetBranchAddress("bb.sh.colblk",&sh_colblk);
  Double_t sh_atimeblk;  C->SetBranchStatus("bb.sh.atimeblk",1); C->SetBranchAddress("bb.sh.atimeblk",&sh_atimeblk);
  Double_t sh_nblk;      C->SetBranchStatus("bb.sh.nblk",1); C->SetBranchAddress("bb.sh.nblk",&sh_nblk);
  Double_t sh_clblk_e;   C->SetBranchStatus("bb.sh.clus_blk.e",1); C->SetBranchAddress("bb.sh.clus_blk.e",&sh_clblk_e);
  Double_t sh_idblk;     C->SetBranchStatus("bb.sh.idblk",1); C->SetBranchAddress("bb.sh.idblk",&sh_idblk);
  Double_t sh_clblk_atime[maxtr]; C->SetBranchStatus("bb.sh.clus_blk.atime",1); C->SetBranchAddress("bb.sh.clus_blk.atime",&sh_clblk_atime);
  //preshower
  Double_t ps_nclus;     C->SetBranchStatus("bb.ps.nclus",1); C->SetBranchAddress("bb.ps.nclus",&ps_nclus);
  Double_t ps_e;         C->SetBranchStatus("bb.ps.e",1); C->SetBranchAddress("bb.ps.e",&ps_e);
  Double_t ps_rowblk;    C->SetBranchStatus("bb.ps.rowblk",1); C->SetBranchAddress("bb.ps.rowblk",&ps_rowblk);
  Double_t ps_colblk;    C->SetBranchStatus("bb.ps.colblk",1); C->SetBranchAddress("bb.ps.colblk",&ps_colblk);
  Double_t ps_atimeblk;  C->SetBranchStatus("bb.ps.atimeblk",1); C->SetBranchAddress("bb.ps.atimeblk",&ps_atimeblk);
  Double_t ps_idblk;     C->SetBranchStatus("bb.ps.idblk",1); C->SetBranchAddress("bb.ps.idblk",&ps_idblk);
  Double_t ps_clblk_atime[maxtr]; C->SetBranchStatus("bb.ps.clus_blk.atime",1); C->SetBranchAddress("bb.ps.clus_blk.atime",&ps_clblk_atime);
  //hcal
  Double_t hcal_atimeblk;C->SetBranchStatus("sbs.hcal.atimeblk",1); C->SetBranchAddress("sbs.hcal.atimeblk",&hcal_atimeblk);
  //gem
  Double_t p[maxtr];     C->SetBranchStatus("bb.tr.p",1); C->SetBranchAddress("bb.tr.p",&p);
  Double_t pz[maxtr];    C->SetBranchStatus("bb.tr.pz",1); C->SetBranchAddress("bb.tr.pz",&pz);
  Double_t tg_th[maxtr]; C->SetBranchStatus("bb.tr.tg_th",1); C->SetBranchAddress("bb.tr.tg_th",&tg_th);
  Double_t tg_ph[maxtr]; C->SetBranchStatus("bb.tr.tg_ph",1); C->SetBranchAddress("bb.tr.tg_ph",&tg_ph);
  //hodo
  Double_t hodo_nclus;   C->SetBranchStatus("bb.hodotdc.nclus",1); C->SetBranchAddress("bb.hodotdc.nclus",&hodo_nclus);
  Double_t hodo_tmean[maxtr];   C->SetBranchStatus("bb.hodotdc.clus.tmean",1); C->SetBranchAddress("bb.hodotdc.clus.tmean",&hodo_tmean);
  Double_t hodo_trIndex[maxtr]; C->SetBranchStatus("bb.hodotdc.clus.trackindex",1); C->SetBranchAddress("bb.hodotdc.clus.trackindex",&hodo_trIndex);
  // Event info
  C->SetMakeClass(1);
  C->SetBranchStatus("fEvtHdr.*", 1);
  UInt_t rnum;           C->SetBranchAddress("fEvtHdr.fRun", &rnum);
  UInt_t trigbits;       C->SetBranchAddress("fEvtHdr.fTrigBits", &trigbits);
  ULong64_t gevnum;      C->SetBranchAddress("fEvtHdr.fEvtNum", &gevnum);
  // turning on additional branches for the global cut
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.gem.track.nhits",1);
  //if (exp=="gmn" && ppass<=2 && config>7)
  C->SetBranchStatus("g.trigbits",1);

  // creating atimeOff histograms per BBCal block
  Double_t h_atime_blk_bin = 240, h_atime_blk_min = atppos_nom-60., h_atime_blk_max = atppos_nom+60.;
  Double_t h_atime_blk_corr_bin = 240, h_atime_blk_corr_min = -60., h_atime_blk_corr_max = 60.;
  
  // Let's read in old ADC time offsets for both SH and PS
  //std::cout << std::endl;
  //Double_t old_ash_atimeOffs[kNblksSH];
  //Double_t old_aps_atimeOffs[kNblksPS];
  //for (int i=0; i<kNblksSH; i++) { old_ash_atimeOffs[i] = -1000; }
  //for (int i=0; i<kNblksPS; i++) { old_aps_atimeOffs[i] = -1000; }
  //TString atimeOff_sh, atimeOff_ps;
  char const * exptag = "unknown";
  if (exp=="gmn") exptag = "sbs";
  else if (exp=="gen") exptag = "gen";
  //int cpass = ppass-1;
  if (exp=="gmn" && ppass==2) cpass = 0;
  char const * setno = set.Atoi() < 0 ? "" : ("_set" + set).Data();
  //atimeOff_sh = Form("Output/%s%d%s_prepass%d_atimeOff_sh.txt",exptag,config,setno,cpass);
  //atimeOff_ps = Form("Output/%s%d%s_prepass%d_atimeOff_ps.txt",exptag,config,setno,cpass);
  //ReadOffset(atimeOff_sh, old_ash_atimeOffs);
  //ReadOffset(atimeOff_ps, old_aps_atimeOffs);

  // define output files
  char const * debug = isdebug ? "_test" : "";
  //TString outFile, outPeaks, toffset_ps, toffset_sh;
  //outPeaks = Form("plots/%s%d%s_prepass%d_atimeOff%s.pdf",exptag,config,setno,ppass,debug);
  //outFile = Form("hist/%s%d%s_prepass%d_atimeOff%s.root",exptag,config,setno,ppass,debug);
  //toffset_sh = Form("Output/%s%d%s_prepass%d_atimeOff_sh%s.txt",exptag,config,setno,ppass,debug);
  //toffset_ps = Form("Output/%s%d%s_prepass%d_atimeOff_ps%s.txt",exptag,config,setno,ppass,debug);
  //ofstream toffset_psdata, toffset_shdata;
  //toffset_psdata.open(toffset_ps);
  //toffset_shdata.open(toffset_sh);
  //TFile *fout = new TFile(outFile,"RECREATE");
  //fout->cd();

  // defining important histograms
  //TH1F *h_W = new TH1F("h_W","W distribution",200,0.,5.);
  //TH1F *h_Q2 = new TH1F("h_Q2","Q2 distribution",300,1.,15.);

  // determining histogram ranges --
  Double_t h_atime_bin = 250, h_atime_min = atppos_old-20., h_atime_max = atppos_old+20;
  Double_t h_atime_corr_bin = 250, h_atime_corr_min = atppos_new-20., h_atime_corr_max = atppos_new+20;
  Double_t h_atime_off_bin = 250, h_atime_off_min = -atppos_old-20, h_atime_off_max = -atppos_old+20;
  Double_t h_atime_off_corr_bin = 250, h_atime_off_corr_min = -atppos_new-20, h_atime_off_corr_max = -atppos_new+20;
  // --

  Double_t Nruns = 1000; // Max # runs we anticipate to analyze
  // Double_t coin_ppos = hcal_atppos - atppos_old;
  //Double_t coin_ppos_corr = hcal_atppos - atppos_new;

  // diagnostic plots
  TH2F *h2_tmean_vs_rnum = new TH2F("h2_tmean_vs_rnum","Before offset correction;Run no.;TH ClusTmean (ns)",Nruns,0.5,Nruns+0.5,20,-10,10);
  TH2F *h2_trigbit_vs_rnum = new TH2F("h2_trigbit_vs_rnum","Before offset correction;Run no.;TrigBits Value",Nruns,0.5,Nruns+0.5,5,0,5);
  TH2F *h2_bbsh_time_vs_rnum = new TH2F("h2_bbsh_time_vs_rnum","Before offset correction;Run no.;SH ADCtime (ns)",Nruns,0.5,Nruns+0.5,h_atime_bin,h_atime_min,h_atime_max);
  TH2F *h2_hcal_time_vs_rnum = new TH2F("h2_hcal_time_vs_rnum","Before offset correction;Run no.;HCAL ADCtime (ns)",Nruns,0.5,Nruns+0.5,300,0,300);
  TH2F *h2_coin_time_vs_rnum = new TH2F("h2_coin_time_vs_rnum","Before offset correction;Run no.;HCAL ADCtime - SH ADCtime (ns)",Nruns,0.5,Nruns+0.5,400,50,250);

  ///////////////////////////////////////////
  // 1st Loop over all events to calibrate //
  ///////////////////////////////////////////

  cout << endl;
  Long64_t nevent=0, nevents=C->GetEntries(); UInt_t runnum=0;
  Double_t timekeeper = 0., timeremains = 0.;
  int treenum = 0, currenttreenum = 0, itrrun=0;
  std::vector<std::string> lrnum;    // list of run numbers
  lrnum.reserve(100);

  while( C->GetEntry( nevent++ ) ){
    // Calculating remaining time
    sw2->Stop();
    timekeeper += sw2->RealTime();
    if (nevent % 25000 == 0 && nevent != 0)
      timeremains = timekeeper * (double(nevents) / double(nevent) - 1.);
    sw2->Reset();
    sw2->Continue();

    if(nevent % 100 == 0) cout << nevent << "/" << nevents  << ", " << int(timeremains/60.) << "m \r";;
    cout.flush();
    // ------

    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();

      // track change of runnum
      if (nevent == 1 || rnum != runnum) {
	runnum = rnum; itrrun++;
	lrnum.push_back(to_string(rnum));
      }
    }
    //lrnum.push_back(to_string(rnum));
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;
    if (passedgCut) {

      //calculating physics parameters
      //Double_t P_ang = 57.3*TMath::ACos(pz[0]/p[0]);
      //Double_t Q2 = 4.*Ebeam*p[0]*pow( TMath::Sin(P_ang/57.3/2.),2. );
      //Double_t W2 = Mp*Mp + 2.*Mp*(Ebeam-p[0]) - Q2;
      //Double_t W = 0.;

      //h_Q2->Fill(Q2);
      //if( W2>0. ){
      //W = TMath::Sqrt(W2);
      //h_W->Fill(W);
      //}

      //hodo cut
      //if( hodo_trIndex[0]!=0 ) continue;
     
      h2_tmean_vs_rnum -> Fill(itrrun, hodo_tmean[0]);
      h2_trigbit_vs_rnum -> Fill(itrrun, trigbits);
      h2_bbsh_time_vs_rnum -> Fill(itrrun, sh_clblk_atime[0]);
      h2_hcal_time_vs_rnum -> Fill(itrrun, hcal_atimeblk);
      h2_coin_time_vs_rnum -> Fill(itrrun, hcal_atimeblk-sh_clblk_atime[0]);

      // storing good event numbers
      goodevents.push_back(nevent);
  } //global cut
  } //while
  cout << endl << endl;

  // customizing histos with run # on the x-axis
  Custm2DRnumHisto(h2_tmean_vs_rnum, lrnum); Custm2DRnumHisto(h2_trigbit_vs_rnum, lrnum);
  Custm2DRnumHisto(h2_bbsh_time_vs_rnum, lrnum); Custm2DRnumHisto(h2_hcal_time_vs_rnum, lrnum);
  Custm2DRnumHisto(h2_coin_time_vs_rnum, lrnum);

  /////////////////////////////////
  // Generating diagnostic plots //
  /////////////////////////////////
  
  /**** Canvas 11 (Diagnostics) ****/
  TCanvas *c11 = new TCanvas("c11", "Hodo Diagnostic Plots", 1200, 800);
  c11->Divide(1,2);
  c11->cd(1); //
  gPad->SetGridy();
  gPad->SetLogz();
  h2_tmean_vs_rnum->SetStats(0);
  h2_tmean_vs_rnum->Draw("colz");
  c11->cd(2); //
  gPad->SetGridy();
  h2_trigbit_vs_rnum->SetStats(0);
  h2_trigbit_vs_rnum->Draw("colz");
  c11->SaveAs("DiagnosticPlots_hodo.pdf");
  //**** -- ***//

  /**** Canvas 12 (Diagnostics) ****/
  TCanvas *c12 = new TCanvas("c12", "Diagnostic Plots", 1200, 800);
  c12->Divide(1,2);
  c12->cd(1); //
  gPad->SetGridy();
  gPad->SetLogz();
  h2_bbsh_time_vs_rnum->SetStats(0);
  h2_bbsh_time_vs_rnum->Draw("colz");
  c12->cd(2); //
  gPad->SetGridy();
  gPad->SetLogz();
  h2_hcal_time_vs_rnum->SetStats(0);
  h2_hcal_time_vs_rnum->Draw("colz");
  c12->SaveAs("DiagnosticPlots.pdf");
  //**** -- ***//

  /**** Canvas 13 (Coincidence Diagnostics) ****/
  TCanvas *c13 = new TCanvas("c13", "Coin Plot", 1200, 800);
  c13->Divide(1,1);
  c13->cd(1);
  gPad->SetGridy();
  gPad->SetLogz();
  h2_coin_time_vs_rnum->SetStats(0);
  h2_coin_time_vs_rnum->Draw("colz");
  c13->SaveAs("CoinPlot.pdf");
  //**** -- ***//

  cout << "CPU time = " << sw->CpuTime() << "s. Real time = " << sw->RealTime() << "s. \n\n";
}

// **** ========== Useful functions ========== ****
// Returns today's date
string getDate(){
  // returns today's date
  time_t now = time(0);
  tm ltm = *localtime(&now);
  string yyyy = to_string(1900 + ltm.tm_year);
  string mm = to_string(1 + ltm.tm_mon);
  string dd = to_string(ltm.tm_mday);
  string date = mm + '/' + dd + '/' + yyyy;
  return date;
}

/*
// reads old ADC gain coefficients from TXT files
void ReadOffset(TString atimeOff_rfile, Double_t* atimeOff){
  ifstream atimeOff_data;
  atimeOff_data.open(atimeOff_rfile);
  std::string readline;
  Int_t elemID=0;
  if(atimeOff_data.is_open()){
    std::cout << " Reading ADC gain from : "<< atimeOff_rfile << "\n";
    while(getline(atimeOff_data,readline)){
      istringstream tokenStream(readline);
      std::string token;
      char delimiter = ' ';
      while(getline(tokenStream,token,delimiter)){
  	TString temptoken=token;
  	atimeOff[elemID] = temptoken.Atof();
  	elemID++;
      }
    }
  }else{
    std::cerr << " **!** No file : " << atimeOff_rfile << "\n\n";
    std::exit(1);
  }
  atimeOff_data.close();
}
*/

// Create generic histogram function
TH1F* MakeHisto(Int_t row, Int_t col, char const * suf, Int_t bins, Double_t min=0., Double_t max=50.)
{
  Int_t blkid = row*kNcolsSH+col;
  TH1F *h = new TH1F(TString::Format("h_R%d_C%d_Blk_%d%s",row,col,blkid,suf),"",bins,min,max);
  //h->SetStats(0);
  h->SetLineWidth(1);
  h->SetLineColor(2);
  h->GetYaxis()->SetLabelSize(0.06);
  //h->GetYaxis()->SetLabelOffset(-0.17);
  h->GetYaxis()->SetMaxDigits(3);
  h->GetYaxis()->SetNdivisions(5);
  return h;
}

// Customizes 1D histos
void Custm1DHisto(TH1F* h)
{
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetLabelSize(0.045);
  h->SetLineWidth(1);
  h->SetLineColor(kBlack);
  h->SetMarkerSize(0.9);
  h->SetMarkerColor(kRed);
  h->SetMarkerStyle(kFullCircle);
}

// Customizes 2D histos with run # on the X-axis
void Custm2DRnumHisto(TH2F* h, std::vector<std::string> const & lrnum)
{
  h->GetXaxis()->SetRange(1,lrnum.size());
  for (int i=0; i<lrnum.size(); i++) h->GetXaxis()->SetBinLabel(i+1,lrnum[i].c_str());
  if (lrnum.size()>15) h->LabelsOption("v", "X");
}

// splits a string by a delimiter (doesn't include empty sub-strings)
std::vector<std::string> SplitString(char const delim, std::string const myStr) {
  std::stringstream ss(myStr);
  std::vector<std::string> out;
  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, delim);
    if (!substr.empty()) out.push_back(substr);
  }
  if (out.empty()) std::cerr << "WARNING! No substrings found!\n";
  return out;
}

// returns output file base from configfilename
TString GetOutFileBase(TString configfilename) {
  std::vector<std::string> result;
  result = SplitString('/',configfilename.Data());
  TString temp = result[result.size() - 1];
  return temp.ReplaceAll(".cfg", "");
}
