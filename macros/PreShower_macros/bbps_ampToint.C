/*
This macro computes the signal amplitude to integral ratios for individual PS blocks.
To execute:
----
[a-onl@aonl2 macros]$ pwd
/adaqfs/home/a-onl/sbs/BBCal_replay/macros
[a-onl@aonl2 macros]$ root -l
root [0] .x PreShower_macros/bbps_ampToint.C
----
- P. Datta
*/

#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include "fadc_data.h"

using namespace std;

// Detector parameters
const int kNrows = 26;
const int kNcols = 2;

const double TargetADC = 5.; //pC

TChain *T = 0;
int gCurrentEntry = -1;
TCanvas *subCanv[4];

// Declare necessary histograms
TH1F *hADCamp[kNrows][kNcols];
TH1F *hamptointratio[kNrows][kNcols];

// Declare necessary arrayes 
bool gPulse[kNrows+2][kNcols+2];
double Pars[3];
double ParErrs[3];

// Declare necessary functions
string getDate();
TH1F* MakeHisto( Int_t, Int_t, Int_t, const char*, Double_t, Double_t );
void processEvent( int, bool );
void goodHistoTest( double, int, int );
void makeSummaryPlots( int, string, bool );
void GetTrigtoFADCratio();

// Declare vectors necessary to make diagnostic plots
double blocks[kNrows*kNcols], peakPos[kNrows*kNcols], peakPosErr[kNrows*kNcols];
double RMS[kNrows*kNcols], RMSErr[kNrows*kNcols], NinPeak[kNrows*kNcols], HVCrrFact[kNrows*kNcols];
vector<double> trigTofadcRatio;

const Int_t kCanvSize = 100;
namespace psgui {
  TGMainFrame *main = 0;
  TGHorizontalFrame *frame1 = 0;
  TGTab *fTab;
  TGLayoutHints *fL3;
  TGCompositeFrame *tf;
  TGTextButton *exitButton;
  // TGTextButton *displayEntryButton;
  // TGTextButton *displayNextButton;
  TGNumberEntry *entryInput;
  // TGLabel *ledLabel;

  TRootEmbeddedCanvas *canv[4];

  TGCompositeFrame* AddTabSub(Int_t sub) {
    tf = fTab->AddTab(Form("PS Sub%d",sub+1));

    TGCompositeFrame *fF5 = new TGCompositeFrame(tf, (12+1)*kCanvSize,(6+1)*kCanvSize , kHorizontalFrame);
    TGLayoutHints *fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
        kLHintsExpandY, 5, 5, 5, 5);
    TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(Form("psSubCanv%d",sub), fF5, 6*kCanvSize,8*kCanvSize);
    //TRootEmbeddedCanvas *fEc1 = new TRootEmbeddedCanvas(0, fF5, 600, 600);
    //Int_t wid = fEc1->GetCanvasWindowId();
    //subCanv[sub] = new TCanvas(Form("subCanv%d",sub),10,10,wid);
    //subCanv[sub]->Divide(12,6);
    //fEc1->AdoptCanvas(subCanv[sub]);
    canv[sub] = fEc1;
    fF5->AddFrame(fEc1,fL4);
    tf->AddFrame(fF5,fL4);
    return tf;
  }

  void SetupGUI() {
    if(!main) {
      main = new TGMainFrame(gClient->GetRoot(), 1200, 1100);
      frame1 = new TGHorizontalFrame(main, 150, 20, kFixedWidth);
      //ledLabel = new TGLabel(frame1,"LED Bit:    , Count:      ");
      //displayEntryButton = new TGTextButton(frame1,"&Display Entry","clicked_displayEntryButton()");
      //entryInput = new TGNumberEntry(frame1,0,5,-1,TGNumberFormat::kNESInteger);
      //displayNextButton = new TGTextButton(frame1,"&Next Entry","clicked_displayNextButton()");
      exitButton = new TGTextButton(frame1, "&Exit", 
          "gApplication->Terminate(0)");
      TGLayoutHints *frame1LH = new TGLayoutHints(kLHintsTop|kLHintsLeft|
          kLHintsExpandX,2,2,2,2);
      //frame1->AddFrame(ledLabel,frame1LH);
      //frame1->AddFrame(displayEntryButton,frame1LH);
      //frame1->AddFrame(entryInput,frame1LH);
      //frame1->AddFrame(displayNextButton,frame1LH);
      frame1->AddFrame(exitButton,frame1LH);
      //frame1->Resize(800, displayNextButton->GetDefaultHeight());
      main->AddFrame(frame1, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1));

      // Create the tab widget
      fTab = new TGTab(main, 300, 300);
      fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

      // Create Tab1 (PS Sub1)
      for(Int_t i = 0; i < 4; i++) {
        tf = AddTabSub(i);
      }
      main->AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1));
      main->MapSubwindows();
      main->Resize();   // resize to default size
      main->MapWindow();

      for(Int_t i = 0; i < 4; i++) {
        subCanv[i] = canv[i]->GetCanvas();
	subCanv[i]->Divide(kNcols,7,0.001,0.001);
      }
    }
  }
};

// Main
void bbps_ampToint ( int run = 366, int event = -1 ){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  psgui::SetupGUI();
  gStyle->SetLabelSize(0.08,"XY");
  gStyle->SetTitleFontSize(0.08);

  bool diagPlots = 0;
  bool trigAmp = 0;
  string date = getDate();

  // Take user inputs
  // cout << " Run number? " << endl;
  // cin >> run;
  // cout << " No. of events replayed? [-1 => All] " << endl;
  // cin >> event;

  // Out files
  string OutFile = Form("Output/fit_results/run_%d_ps_ampToint.txt",run);
  ofstream outfile_data;
  outfile_data.open(OutFile);

  int hADCamp_bin = 60, hamptointratio_bin = 18;
  double hADCamp_min = 0., hamptointratio_min = 2.;
  double hADCamp_max = 20., hamptointratio_max = 4.2;

  // Read in data produced by analyzer in root format
  if(!T) { 
    T = new TChain("T");
    //TString dataDIR = gSystem->Getenv("OUT_DIR");
    TString filename = Form("../../Rootfiles/bbshower_%d_%d.root", run, event);
    TString filename_seg = Form("../../Rootfiles/bbshower_%d_%d_*.root", run, event);
    // TString filename_seg = Form("../../Rootfiles/bbshower_%d_%d_seg_*.root", run, event); // use this for runs replayed during GMn time


    T->Add(filename);
    T->Add(filename_seg);
    T->SetBranchStatus("*",0);
    T->SetBranchStatus("bb.ps.*",1);
    T->SetBranchAddress("bb.ps.a_p",fadc_datat::a);
    T->SetBranchAddress("bb.ps.a_amp_p",fadc_datat::amp);
    T->SetBranchAddress("bb.ps.a_time",fadc_datat::tdc);
    T->SetBranchAddress("bb.ps.adcrow",fadc_datat::row);
    T->SetBranchAddress("bb.ps.adccol",fadc_datat::col);
    T->SetBranchStatus("Ndata.bb.ps.adcrow",1);
    T->SetBranchAddress("Ndata.bb.ps.adcrow",&fadc_datat::ndata);
    for(int r = 0; r < kNrows; r++) {
      for(int c = 0; c < kNcols; c++) {
	hADCamp[r][c] = MakeHisto(r, c, hADCamp_bin, "_i", hADCamp_min, hADCamp_max);
	hamptointratio[r][c] = MakeHisto(r, c, hamptointratio_bin, "_r", hamptointratio_min, hamptointratio_max);
      }
    }
  }
  
  cout << endl;
  Long64_t nevents = T->GetEntries();

  // Set default values for pulse check bool arrays. Arrays contain false-valued buffer for possible future diagonal track implementation.
  for(int r = 0; r < kNrows+2; r++) {
    for(int c = 0; c < kNcols+2; c++) {
      if( r==0 || r==kNrows+1 ){
	gPulse[r][c] = true;
      } else {
	gPulse[r][c] = false;
      }
      
      blocks[kNcols*r+c] = 0.;
      peakPos[kNcols*r+c] = 0.;
      peakPosErr[kNcols*r+c] = 0.;
      RMS[kNcols*r+c] = 0.;
      RMSErr[kNcols*r+c] = 0.;
      NinPeak[kNcols*r+c] = 0.;
      HVCrrFact[kNcols*r+c] = 0.;  
    }  
  }

  cout << "Processing " << nevents << " events ...." << endl;

  // Looping through events
  double progress = 0.;
  while(progress<1.0){
    int barwidth = 70;
    for (int nev = 0; nev < nevents; nev++){ 
      processEvent( nev, trigAmp );
      
      cout << "[";
      int pos = barwidth * progress;
      for(int i=0; i<barwidth; ++i){
    	if(i<pos) cout << "=";
    	else if(i==pos) cout << ">";
    	else cout << " ";
      }
      progress = (double)((nev+1.)/nevents);
      cout << "] " << int(progress*100.) << "%\r";
      cout.flush();
    }
  }
  cout << endl << endl;

  // Let's fit the histograms with Gauss (twice)
  TF1 *fgaus = new TF1("fgaus","gaus");
  TF1 *fgaus2 = new TF1("fgaus2","gaus");

  int sub = 0;
  for(int r=0; r<kNrows; r++){
    for(int c=0; c<kNcols; c++){
      string Flag;
      double HVcorrection = 1.0;
      double_t NinPeakPS =0.0;
      for( int i=0; i<4; i++ ) {
      	Pars[i] = 0.;
      	ParErrs[i] = 0.;
      }
      
      if( hADCamp[r][c]->GetEntries() > 0. ){
	// Create fit functions for each module
	fgaus->SetLineColor(2);
	fgaus->SetNpx(1000);

	// first fit
	int maxBin = hADCamp[r][c]->GetMaximumBin();
	double maxBinCenter = hADCamp[r][c]->GetXaxis()->GetBinCenter( maxBin );
	double maxCount = hADCamp[r][c]->GetMaximum();
	double binWidth = hADCamp[r][c]->GetBinWidth(maxBin);
	double stdDev = hADCamp[r][c]->GetStdDev();

	// Reject low energy peak
	if( hADCamp[r][c]->GetBinContent(maxBin-2) < 0.02*hADCamp[r][c]->GetBinContent(maxBin) ){
	  while ( hADCamp[r][c]->GetBinContent(maxBin+1) < hADCamp[r][c]->GetBinContent(maxBin) || hADCamp[r][c]->GetBinContent(maxBin+1) == hADCamp[r][c]->GetBinContent(maxBin) ) 
	    {
		maxBin++;
	    };
	  hADCamp[r][c]->GetXaxis()->SetRange( maxBin+1 , hADCamp[r][c]->GetNbinsX() );
	  maxBin = hADCamp[r][c]->GetMaximumBin();
	  maxBinCenter = hADCamp[r][c]->GetXaxis()->GetBinCenter( maxBin );
	  maxCount = hADCamp[r][c]->GetMaximum();
	  binWidth = hADCamp[r][c]->GetBinWidth(maxBin);
	  stdDev = hADCamp[r][c]->GetStdDev();
	}

	fgaus->SetParameters( maxCount,maxBinCenter,stdDev );
	fgaus->SetRange( hADCamp_min, hADCamp_max );
	hADCamp[r][c]->Fit(fgaus,"NO+RQ");
	fgaus->GetParameters(Pars);

  
	// Second fit with tailored range
	int lowerBinC = 0;
	int upperBinC = 0;
	if( r==0 || r==(kNrows-1)){ // To make the fits better for top and bottom rows
	  lowerBinC = hADCamp_min + (maxBin)*binWidth - (0.8*Pars[2]);
	  upperBinC = hADCamp_min + (maxBin)*binWidth + (1.8*Pars[2]);
	}// else if( r==0&&c==1 ){
	//   lowerBinC = hADCamp_min + (maxBin)*binWidth - (0.5*Pars[2]);
	//   upperBinC = hADCamp_min + (maxBin)*binWidth + (0.5*Pars[2]);
	// }
	else{
	  lowerBinC = hADCamp_min + (maxBin)*binWidth - (2.0*Pars[2]);
	  upperBinC = hADCamp_min + (maxBin)*binWidth + (2.0*Pars[2]);
	}
	fgaus->SetParameters( Pars[0],Pars[1],Pars[2] );
	fgaus->SetRange( lowerBinC, upperBinC );
	
	sub = r/7;
	subCanv[sub]->cd((r%7)*kNcols + c + 1);
	
	hADCamp[r][c]->Fit( fgaus,"NO+RQ" );
	fgaus->GetParameters(Pars);
	for ( int i=0; i<3; i++ ) ParErrs[i] = fgaus->GetParError(i); 

	Flag = "Good"; // States quality of fit
	HVcorrection = pow( (TargetADC/Pars[1]), 0.10); // Correction term for HV
	NinPeakPS = hADCamp[r][c]->Integral( hADCamp[r][c]->FindFixBin(lowerBinC),hADCamp[r][c]->FindFixBin(upperBinC),"" ); // # of good events
	  
	// hADCamp[r][c]->GetXaxis()->SetRange( 2 , hADCamp[r][c]->GetNbinsX() );
	// if(c%2==0){
	//   hADCamp[r][c]->SetTitle(TString::Format("PS-R%d | Amp  ",r+1));
	// }else{
	//   hADCamp[r][c]->SetTitle(TString::Format("PS-L%d | Amp  ",r+1));
	// }
	// hADCamp[r][c]->GetYaxis()->SetLabelSize(0.06);
	// hADCamp[r][c]->SetLineColor(kBlue+1);
	// hADCamp[r][c]->Draw();
	// gPad->Update();

	fgaus2->SetLineColor(2);
	fgaus2->SetNpx(1000);
	fgaus2->SetParameters( maxCount,maxBinCenter,stdDev );
	fgaus2->SetRange( hamptointratio_min, hamptointratio_max );
	hamptointratio[r][c]->Fit(fgaus2,"+RQ");

	// outfile_data << r*kNcols+c << "\t" << fgaus2->GetParameter(1) << endl;
	outfile_data << fgaus2->GetParameter(1) << ' ';
	cout <<  fgaus2->GetParameter(1) << ' ';

	if(c%2==0){
	  hamptointratio[r][c]->SetTitle(TString::Format("PS-R%d | Amp/Int  ",r+1));
	}else{
	  hamptointratio[r][c]->SetTitle(TString::Format("PS-L%d | Amp/Int  ",r+1));
	}
	//hamptointratio[r][c]->GetYaxis()->SetLabelSize(0.06);
        hamptointratio[r][c]->GetXaxis()->SetLabelSize(0.04);
	hamptointratio[r][c]->GetYaxis()->SetLabelSize(0.04);
	hamptointratio[r][c]->SetLineColor(kBlue+1);
	hamptointratio[r][c]->Draw();
	gPad->Update();

	// Let's determine how good is the fit
	if( hADCamp[r][c]->GetEntries() < 20 ){
	  cout << " ** The histogram for module # " << r+1 << "." << c+1 << " was empty!! " << endl;
	  Flag = "No_Data";
	}else if( Pars[2] > 60.0 ){
	  cout << " ** Fit for module # " << r+1 << "." << c+1 << " was too wide!! " << endl;
	  Flag = "Wide";
	  HVcorrection = 1.0;
	}else if( Pars[2] < 0.1 ){
	  cout << " ** Fit for module # " << r+1 << "." << c+1 << " was too narrow!! " << endl;
	  Flag = "Narrow";
	  HVcorrection = 1.0;
	}else if( ParErrs[1] > 100. || ParErrs[2] > 20. ){
	  cout << " ** For module # " << r+1 << "." << c+1 << ", error bar was too high!! " << endl;
	  Flag = "Big_error";
	  HVcorrection = 1.0;
	}

	if ( Flag != "Good" ){
	  for( int i=0; i<4; i++ ) { Pars[i] = 0.; ParErrs[i] = 0.; }
	  NinPeakPS = 0.;
	}

	// Fill in the arrays to create diagnostic plots
	blocks[kNcols*r+c] = kNcols*r+c;
	peakPos[kNcols*r+c] = Pars[1];
	peakPosErr[kNcols*r+c] = ParErrs[1];
	RMS[kNcols*r+c] = Pars[2];
	RMSErr[kNcols*r+c] = ParErrs[2];
	NinPeak[kNcols*r+c] = NinPeakPS;
	HVCrrFact[kNcols*r+c] = HVcorrection;
	
      }
    }
    cout << endl;
    outfile_data << endl;
  }

  TString OutFile_peaks = Form("plots/PS_ampToint_%d.pdf",run);
  subCanv[0]->SaveAs(Form("%s[", OutFile_peaks.Data()));
  for( int canC=0; canC<4; canC++ ) subCanv[canC]->SaveAs(Form("%s", OutFile_peaks.Data()));
  subCanv[3]->SaveAs(Form("%s]", OutFile_peaks.Data()));

  // Making diagnostic plots
  if( diagPlots ){
    cout << "Creating diagnostic plots.." << endl;
    makeSummaryPlots( run, date, trigAmp );
  }

  // Close all the outFiles
  outfile_data.close();

  // Post analysis reporting
  cout << "-------" << endl;
  cout << " Amplitude to Intigral ratios written to : " << OutFile << endl;
  cout << " Peaks saved to : " << OutFile_peaks << endl;
  cout << "-------" << endl;
  
} //main



// Get today's date
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);

  string yyyy = to_string(1900 + ltm.tm_year);
  string mm = to_string(1 + ltm.tm_mon);
  string dd = to_string(ltm.tm_mday);
  string date = mm + '/' + dd + '/' + yyyy;

  return date;
} // getDate

// Create generic histogram function
TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins, const char* suf="", Double_t min=0., Double_t max=50.)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d%s",row,col, suf),
		     TString::Format("%d-%d",row+1,col+1),bins,min,max);
  h->SetStats(0);
  h->SetLineWidth(2);
  h->GetYaxis()->SetLabelSize(0.1);
  //h->GetYaxis()->SetLabelOffset(-0.17);
  h->GetYaxis()->SetNdivisions(5);
  return h;
}

// Process events
void processEvent( int entry, bool trigAmp ){
  // Check event increment and increment
  if(entry == -1) {
    gCurrentEntry++;
  } else {
    gCurrentEntry = entry;
  }
  
  if(gCurrentEntry<0) {
    gCurrentEntry = 0;
  }

  // Get the event from the TTree
  T->GetEntry(gCurrentEntry);

  int r,c,idx,n,sub;
  // Clear old signal histograms, just in case modules are not in the tree
  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      gPulse[r+1][c+1] = false;
    }
  }
  
  
  // Reset signal peak, adc, and tdc arrays
  double adc[kNrows][kNcols];
  double adc_amp[kNrows][kNcols];
  double tdc[kNrows][kNcols];
  for(r  = 0; r < kNrows; r++) {
    for(c  = 0; c < kNcols; c++) {
      adc[r][c] = 0.0;
      adc_amp[r][c] = 0.0;
      tdc[r][c] = 0.0;
    }
  }

  // Process events with module data
  for(int m = 0; m < fadc_datat::ndata; m++) {
    // Define row and column
    r = fadc_datat::row[m]; 
    c = fadc_datat::col[m]; 
    if(r < 0 || c < 0) {
      cerr << "Why is row negative? Or col?" << endl;
      continue;
    }
    
    if(r>= kNrows || c >= kNcols) continue;
    
    // Define index, number of samples, fill adc and tdc arrays, and switch processed marker for error reporting
    adc[r][c] = fadc_datat::a[m];
    adc_amp[r][c] = fadc_datat::amp[m];
    tdc[r][c] = fadc_datat::tdc[m];
    goodHistoTest( tdc[r][c],r,c ); // Setting flag for good events
  }

  if(trigAmp){
    trigTofadcRatio.clear();
    GetTrigtoFADCratio();
  }

  for(r = 0; r < kNrows; r++) {
    for(c = 0; c < kNcols; c++) {
      Int_t m = r*kNcols+c;
      if( c==0 ){ // PS R
	if( r==0 ){ // Bottom row     
	  if( (gPulse[r][c+1]==true)&&(gPulse[r+2][c+1]==true&&gPulse[r+3][c+1]==true)&&(gPulse[r+1][c+2]==false) ){
	    if(trigAmp) {
	      hADCamp[r][c]->Fill( trigTofadcRatio.at(m)*adc[r][c] ); 
	    }else{
	      hADCamp[r][c]->Fill( adc[r][c] );
	    }
	    hamptointratio[r][c]->Fill( adc_amp[r][c]/adc[r][c] ); 	
	  }	
	}else if( r==(kNrows-1) ){ // Top row
	  if( (gPulse[r-1][c+1]==true&&gPulse[r][c+1]==true)&&(gPulse[r+2][c+1]==true)&&(gPulse[r+1][c+2]==false) ){
	    if(trigAmp) {
	      hADCamp[r][c]->Fill( trigTofadcRatio.at(m)*adc[r][c] ); 
	    }else{
	      hADCamp[r][c]->Fill( adc[r][c] );
	    }
	    hamptointratio[r][c]->Fill( adc_amp[r][c]/adc[r][c] ); 	      	
	  }
	}else{ // All the other rows
	  if( (gPulse[r-1][c+1]==true&&gPulse[r][c+1]==true)&&(gPulse[r+2][c+1]==true&&gPulse[r+3][c+1]==true)&&(gPulse[r+1][c+2]==false) ){
	    if(trigAmp) {
	      hADCamp[r][c]->Fill( trigTofadcRatio.at(m)*adc[r][c] ); 
	    }else{
	      hADCamp[r][c]->Fill( adc[r][c] );
	    }
	    hamptointratio[r][c]->Fill( adc_amp[r][c]/adc[r][c] ); 	      	
	  }
	}      
      }else{ // col 1 : PS L 
	if( r==0 ){      
	  if( (gPulse[r][c+1]==true)&&(gPulse[r+2][c+1]==true&&gPulse[r+3][c+1]==true)&&(gPulse[r+1][c]==false) ){
	    if(trigAmp) {
	      hADCamp[r][c]->Fill( trigTofadcRatio.at(m)*adc[r][c] ); 
	    }else{
	      hADCamp[r][c]->Fill( adc[r][c] );
	    }
	    hamptointratio[r][c]->Fill( adc_amp[r][c]/adc[r][c] ); 	
	  }	
	}else if( r==(kNrows-1) ){
	  if( (gPulse[r-1][c+1]==true&&gPulse[r][c+1]==true)&&(gPulse[r+2][c+1]==true)&&(gPulse[r+1][c]==false) ){
	    if(trigAmp) {
	      hADCamp[r][c]->Fill( trigTofadcRatio.at(m)*adc[r][c] ); 
	    }else{
	      hADCamp[r][c]->Fill( adc[r][c] );
	    }
	    hamptointratio[r][c]->Fill( adc_amp[r][c]/adc[r][c] ); 	      	
	  }
	}else{
	  if( (gPulse[r-1][c+1]==true&&gPulse[r][c+1]==true)&&(gPulse[r+2][c+1]==true&&gPulse[r+3][c+1]==true)&&(gPulse[r+1][c]==false) ){
	    if(trigAmp) {
	      hADCamp[r][c]->Fill( trigTofadcRatio.at(m)*adc[r][c] ); 
	    }else{
	      hADCamp[r][c]->Fill( adc[r][c] );
	    }
	    hamptointratio[r][c]->Fill( adc_amp[r][c]/adc[r][c] ); 	      	
	  }
	}      
      }
    }
  }
  
} //processEvent

// Determine the cut using RMS deviation of baseline fluctuation
void goodHistoTest( double tdcVal, int row, int col ){
    
  //Switch if passes FADC time cut
  if(tdcVal!=0) {
    gPulse[row+1][col+1]=true;
  }
  
} //goodHistoTest

// Create diagnostic plots
void makeSummaryPlots( int run, string date, bool trigAmp = 0 ){
    char CName[9], CTitle[100];
    TCanvas *CGr[4];
    TGraph *Gr[4];

    // gROOT->SetBatch();
    CGr[0] = new TCanvas("c1PS","pPosvsBlocks",100,10,700,500);
    CGr[1] = new TCanvas("c2PS","pRMSvsBlocks",100,10,700,500);
    CGr[2] = new TCanvas("c3PS","#EvInPeakvsBlocks",100,10,700,500);
    CGr[3] = new TCanvas("c4PS","HVcorrvsBlocks",100,10,700,500);
 
    int totalBlocks = kNrows*kNcols;
    double xErr[totalBlocks]; // Holds x-error, which is essentially zero for our case
    for( int i=0; i<totalBlocks; i++ ) xErr[i] = 0. ;
	 
    // Gr[0] = new TGraphErrors( totalBlocks, &(blocks[0]), &(peakPos[0]), xErr, &(peakPosErr[0]) );
    // Gr[1] = new TGraphErrors( totalBlocks, &(blocks[0]), &(RMS[0]), xErr, &(RMSErr[0]) );
    // Gr[2] = new TGraph( totalBlocks, &(blocks[0]), &(NinPeak[0]) );
    // Gr[3] = new TGraph( totalBlocks, &(blocks[0]), &(HVCrrFact[0]) );

    Gr[0] = new TGraphErrors( totalBlocks, blocks, peakPos, xErr, peakPosErr );
    Gr[1] = new TGraphErrors( totalBlocks, blocks, RMS, xErr, RMSErr );
    Gr[2] = new TGraph( totalBlocks, blocks, NinPeak );
    Gr[3] = new TGraph( totalBlocks, blocks, HVCrrFact );

    for( int i = 0; i < 4; i++ ){
      CGr[i]->cd();
      gPad->SetGridy();
      Gr[i]->SetLineColor(2);
      Gr[i]->SetLineWidth(2);
      Gr[i]->SetMarkerColor(1);
      Gr[i]->SetMarkerStyle(20);
      Gr[i]->GetXaxis()->SetLabelSize(0.04);
      Gr[i]->GetYaxis()->SetLabelSize(0.04);
      if(i == 0 ){
    	Gr[i]->SetTitle(Form("Run# %d | Peak Position vs. Block No. for PS Blocks | %s",run,date.c_str()) );
    	Gr[i]->GetXaxis()->SetTitle("Block Number");
    	Gr[i]->GetYaxis()->SetTitle("Peak Position (mV)");
    	Gr[i]->GetYaxis()->SetRangeUser(0.,100.);
      } else if (i == 1){
    	Gr[i]->SetTitle( Form("Run# %d | Peak RMS vs. Block No. for PS Blocks | %s",run,date.c_str()) ); 
    	Gr[i]->GetXaxis()->SetTitle("Block Number");
    	Gr[i]->GetYaxis()->SetTitle("Peak RMS (mV)");
    	Gr[i]->GetYaxis()->SetTitleOffset(1.4);
    	Gr[i]->GetYaxis()->SetRangeUser(0.,20.);
      }else if (i == 2){
    	Gr[i]->SetTitle( Form("Run# %d | N of Events in Peak(fitted region) vs Block No. for PS | %s",run,date.c_str()) );
    	Gr[i]->GetXaxis()->SetTitle("Block Number");
    	Gr[i]->GetYaxis()->SetTitle("N of Events in Peak(fitted region)");
    	Gr[i]->GetYaxis()->SetTitleOffset(1.4);
      }else if (i == 3){
    	Gr[i]->SetTitle( Form("Run# %d | HV Correction Factor vs Block No. for PS blocks | %s",run,date.c_str()) );
    	Gr[i]->GetXaxis()->SetTitle("Block Number");
    	Gr[i]->GetYaxis()->SetTitle("HV Correction Factor");
    	Gr[i]->GetYaxis()->SetTitleOffset(1.4);
	Gr[i]->GetYaxis()->SetRangeUser(0.4,1.6);
      }   
      Gr[i]->Draw("AP");
      // CGr[i]->Write();
    }  

    // CGr[0]->SaveAs( Form("plots/BBPS_Cosmic_Cal_%d.pdf[",run) );
    // for( int i=0; i<4; i++ ) CGr[i]->SaveAs( Form("plots/BBPS_Cosmic_Cal_%d.pdf",run) );
    // CGr[3]->SaveAs( Form("plots/BBPS_Cosmic_Cal_%d.pdf]",run) );
}


void GetTrigtoFADCratio(){
  trigTofadcRatio.clear();
  string InFile = "Coefficients/trigtoFADCcoef_PS.txt";
  ifstream infile_data;
  infile_data.open(InFile);
  TString currentline;
  if (infile_data.is_open() ) {
    // cout << " Reading trigger to FADC ratios from " << InFile << endl;
    TString temp;
    while( currentline.ReadLine( infile_data ) ){
      TObjArray *tokens = currentline.Tokenize("\t");
      int ntokens = tokens->GetEntries();
      if( ntokens > 1 ){
	// temp = ( (TObjString*) (*tokens)[0] )->GetString();
	// elemID.push_back( temp.Atof() );
	temp = ( (TObjString*) (*tokens)[1] )->GetString();
	double theRatio = temp.Atof();
	trigTofadcRatio.push_back( theRatio );
      }
      delete tokens;
    }
    infile_data.close();
  } else {
    cerr << " No file : " << InFile << endl;
  }
}
