#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStopwatch.h>

//const Double_t Mp = 0.938272081;  // +/- 6E-9 GeV
//const Int_t kNcolsSH = 7;   // SH columns
//const Int_t kNrowsSH = 27;  // SH rows
//const Int_t kNblksSH = 189; // Total # SH blocks/PMTs
//const Int_t kNcolsPS = 2;   // PS columns
//const Int_t kNrowsPS = 26;  // PS rows
//const Int_t kNblksPS = 52;  // Total # PS blocks/PMTs

TH1F* MakeHisto(Int_t, char const *, Int_t, Double_t, Double_t);

double SHblocks[189] {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,61.0,62.0,63.0,64.0,65.0,66.0,67.0,68.0,69.0,70.0,71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0,101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0,111.0,112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0,120.0,121.0,122.0,123.0,124.0,125.0,126.0,127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0,136.0,137.0,138.0,139.0,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,161.0,162.0,163.0,164.0,165.0,166.0,167.0,168.0,169.0,170.0,171.0,172.0,173.0,174.0,175.0,176.0,177.0,178.0,179.0,180.0,181.0,182.0,183.0,184.0,185.0,186.0,187.0,188.0};
double PSblocks[52] {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0};
double HCALblocks[288] {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,61.0,62.0,63.0,64.0,65.0,66.0,67.0,68.0,69.0,70.0,71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0,101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0,111.0,112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0,120.0,121.0,122.0,123.0,124.0,125.0,126.0,127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0,136.0,137.0,138.0,139.0,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,161.0,162.0,163.0,164.0,165.0,166.0,167.0,168.0,169.0,170.0,171.0,172.0,173.0,174.0,175.0,176.0,177.0,178.0,179.0,180.0,181.0,182.0,183.0,184.0,185.0,186.0,187.0,188.0,189.0,190.0,191.0,192.0,193.0,194.0,195.0,196.0,197.0,198.0,199.0,200.0,201.0,202.0,203.0,204.0,205.0,206.0,207.0,208.0,209.0,210.0,211.0,212.0,213.0,214.0,215.0,216.0,217.0,218.0,219.0,220.0,221.0,222.0,223.0,224.0,225.0,226.0,227.0,228.0,229.0,230.0,231.0,232.0,233.0,234.0,235.0,236.0,237.0,238.0,239.0,240.0,241.0,242.0,243.0,244.0,245.0,246.0,247.0,248.0,249.0,250.0,251.0,252.0,253.0,254.0,255.0,256.0,257.0,258.0,259.0,260.0,261.0,262.0,263.0,264.0,265.0,266.0,267.0,268.0,269.0,270.0,271.0,272.0,273.0,274.0,275.0,276.0,277.0,278.0,279.0,280.0,281.0,282.0,283.0,284.0,285.0,286.0,287.0};

void bbcal_TimingSelfAlignment(char const * filename)
{

  gErrorIgnoreLevel = kError; // Ignore all ROOT warnings
  //gROOT -> SetBatch(kTRUE); // Don't display every plot you draw

  TChain *C = new TChain("Tout");
  C->Add(filename);

  // Check for empty rootfiles and set tree branches
  if(C->GetEntries()==0)
  {
    cerr << endl << " --- No ROOT file found!! --- " << endl << endl;
    exit(1);
  }
  else cout << endl << "Found " << C->GetEntries() << " events. Starting analysis.. " << endl;

  // Define branches

  // Cuts - when a cut==1 it is activated, when it ==0 it is not
  char WCut;         C->SetBranchAddress("WCut",&WCut);
  char coinCut;      C->SetBranchAddress("coinCut",&coinCut);
  char nCut;         C->SetBranchAddress("nCut",&nCut);
  char pCut;         C->SetBranchAddress("pCut",&pCut);

  // Analysis variables
  Int_t runnum;      C->SetBranchAddress("runnum",&runnum);
  double atimePS;    C->SetBranchAddress("atimePS",&atimePS);
  double atimeSH;    C->SetBranchAddress("atimeSH",&atimeSH);
  double idPS;       C->SetBranchAddress("idPS",&idPS);
  double idSH;       C->SetBranchAddress("idSH",&idSH);
  double idblkHCAL;  C->SetBranchAddress("idblkHCAL",&idblkHCAL);
  double hcal_time;  C->SetBranchAddress("hcal_time",&hcal_time);
  double hodo_time;  C->SetBranchAddress("hodo_time",&hodo_time);
  double dx;         C->SetBranchAddress("dx",&dx);
  double dy;         C->SetBranchAddress("dy",&dy);

  // Define histograms ---------------------------------------------------------

  // Preshower Histograms ------------------------------------------------------

  // 2D histos for ADC time vs block id - before and after alignment

  TH2F *h_psTime_before = new TH2F("h_psTime_before", "ADC Time of the Seed Block for the Best Preshower Cluster [ns] vs. Preshower Block ID BEFORE Alignment", 52, 0, 52, 300, -15, 15);
  h_psTime_before -> GetXaxis() -> SetTitle("Preshower Block ID");
  h_psTime_before -> GetYaxis() -> SetTitle("Preshower ADC Time [ns]");

  TH2F *h_psTime_before_scatter = new TH2F("h_psTime_before_scatter", "ADC Time of the Seed Block for the Best Preshower Cluster [ns] vs. Preshower Block ID BEFORE Alignment", 52, 0, 52, 300, -15, 15);
  h_psTime_before_scatter->SetMarkerSize(0.9);
  h_psTime_before_scatter->SetMarkerColor(kRed);
  h_psTime_before_scatter->SetMarkerStyle(kFullCircle);
  h_psTime_before_scatter->SetLineWidth(1);
  h_psTime_before_scatter->SetLineColor(kBlack);

  TH2F *h_psTime_after = new TH2F("h_psTime_after", "ADC Time of the Seed Block for the Best Preshower Cluster [ns] vs. Preshower Block ID AFTER Alignment", 52, 0, 52, 300, -15, 15);
  h_psTime_after -> GetXaxis() -> SetTitle("Preshower Block ID");
  h_psTime_after -> GetYaxis() -> SetTitle("Preshower ADC Time [ns]");

  TH2F *h_psTime_after_scatter = new TH2F("h_psTime_after_scatter", "ADC Time of the Seed Block for the Best Preshower Cluster [ns] vs. Preshower Block ID AFTER Alignment", 52, 0, 52, 300, -15, 15);
  h_psTime_after_scatter->SetMarkerSize(0.9);
  h_psTime_after_scatter->SetMarkerColor(kRed);
  h_psTime_after_scatter->SetMarkerStyle(kFullCircle);

  // 1D histos for ADC time vs block id - before and after alignment

  TH1F *h_ps_before = new TH1F("h_ps_before", "Preshower ADC Time BEFORE Alignment",300,-15,15);
  h_ps_before -> GetXaxis() -> SetTitle("Preshower ADC Time [ns]");
  TH1F *h_ps_after = new TH1F("h_ps_after", "Preshower ADC Time AFTER Alignment",300,-15,15);
  h_ps_after -> GetXaxis() -> SetTitle("Preshower ADC Time [ns]");

  // 1D histogram for the main block that we will align to
  TH1F *h_mainPS = new TH1F("h_mainPS", "ADC Time of Main Preshower Block",300,-15,15);

  // Make histograms for each individual block
  TH1F *h_atime_psBlocks[52];
  for(int i=0; i<=52; i++)
  {
    h_atime_psBlocks[i] = MakeHisto(i,"", 300, -15, 15);
  }

  // Shower Histograms ---------------------------------------------------------

  // 2D histos for ADC time vs block id - before and after alignment

  TH2F *h_shTime_before = new TH2F("h_shTime_before", "ADC Time of the Seed Block for the Best Shower Cluster [ns] vs. Shower Block ID BEFORE Alignment", 189, 0, 189, 300, -15, 15);
  h_shTime_before -> GetXaxis() -> SetTitle("Shower Block ID");
  h_shTime_before -> GetYaxis() -> SetTitle("Shower ADC Time [ns]");

  TH2F *h_shTime_before_scatter = new TH2F("h_shTime_before_scatter", "ADC Time of the Seed Block for the Best Shower Cluster [ns] vs. Shower Block ID BEFORE Alignment", 189, 0, 189, 300, -15, 15);
  h_shTime_before_scatter->SetMarkerSize(0.9);
  h_shTime_before_scatter->SetMarkerColor(kRed);
  h_shTime_before_scatter->SetMarkerStyle(kFullCircle);

  TH2F *h_shTime_after = new TH2F("h_shTime_after", "ADC Time of the Seed Block for the Best Shower Cluster [ns] vs. Shower Block ID AFTER Alignment", 189, 0, 189, 300, -15, 15);
  h_shTime_after -> GetXaxis() -> SetTitle("Shower Block ID");
  h_shTime_after -> GetYaxis() -> SetTitle("Shower ADC Time [ns]");

  TH2F *h_shTime_after_scatter = new TH2F("h_shTime_after_scatter", "ADC Time of the Seed Block for the Best Shower Cluster [ns] vs. Shower Block ID AFTER Alignment", 189, 0, 189, 300, -15, 15);
  h_shTime_after_scatter->SetMarkerSize(0.9);
  h_shTime_after_scatter->SetMarkerColor(kRed);
  h_shTime_after_scatter->SetMarkerStyle(kFullCircle);

  // 1D histos for ADC time vs block id - before and after alignment

  TH1F *h_sh_before = new TH1F("h_sh_before", "Shower ADC Time BEFORE Alignment",300,-15,15);
  h_sh_before -> GetXaxis() -> SetTitle("Shower ADC Time [ns]");
  TH1F *h_sh_after = new TH1F("h_sh_after", "Shower ADC Time AFTER Alignment",300,-15,15);
  h_sh_after -> GetXaxis() -> SetTitle("Shower ADC Time [ns]");

  // 1D histogram for the main block that we will align to
  TH1F *h_mainSH = new TH1F("h_mainSH", "ADC Time of Main Shower Block",300,-15,15);

  // Make histograms for each individual block
  TH1F *h_atime_shBlocks[189];
  for(int i=0; i<=189; i++)
  {
    h_atime_shBlocks[i] = MakeHisto(i,"", 300, -15, 15);
  }

  // HCal Histograms -----------------------------------------------------------

  // 2D histos for ADC time vs block id - before and after alignment

  TH2F *h_hcalTime_before = new TH2F("h_hcalTime_before", "ADC Time of the Seed Block for the Best HCal Cluster [ns] vs. HCal Block ID BEFORE Alignment", 288, 0, 288, 700, 60, 200);
  h_hcalTime_before -> GetXaxis() -> SetTitle("HCal Block ID");
  h_hcalTime_before -> GetYaxis() -> SetTitle("HCal ADC Time [ns]");

  TH2F *h_hcalTime_after = new TH2F("h_hcalTime_after", "ADC Time of the Seed Block for the Best HCal Cluster [ns] vs. HCal Block ID AFTER Alignment", 288, 0, 288, 700, 60, 200);
  h_hcalTime_after -> GetXaxis() -> SetTitle("HCal Block ID");
  h_hcalTime_after -> GetYaxis() -> SetTitle("HCal ADC Time [ns]");

  // 1D histos for ADC time vs block id - before and after alignment

  TH1F *h_hcal_before = new TH1F("h_hcal_before", "HCal ADC Time BEFORE Alignment",700,60,200);
  h_hcal_before -> GetXaxis() -> SetTitle("HCal ADC Time [ns]");
  TH1F *h_hcal_after = new TH1F("h_hcal_after", "HCal ADC Time AFTER Alignment",700,60,200);
  h_hcal_after -> GetXaxis() -> SetTitle("HCal ADC Time [ns]");

  // 1D histogram for the main block that we will align to
  TH1F *h_mainHCal = new TH1F("h_mainHCal", "ADC Time of Main HCal Block",700,60,200);

  // Make histograms for each individual block
  TH1F *h_atime_hcalBlocks[288];
  for(int i=0; i<=288; i++)
  {
    h_atime_hcalBlocks[i] = MakeHisto(i,"", 700, 60, 200);
  }

  // Coincidence Time Histograms -----------------------------------------------

  // 1D histos for ADC time vs block id - before and after alignment

  TH1F *h_coin_before = new TH1F("h_coin_before", "(HCal ADC - Shower ADC) Coincidence Time BEFORE Alignment [ns]",400,110,150);
  TH1F *h_coin_after = new TH1F("h_coin_after", "(HCal ADC - Shower ADC) Coincidence Time AFTER Alignment [ns]",400,110,150);

  // Start Looping Through Events ----------------------------------------------

  int itrrun=0; // Assign an index value to each run number
  int rnum=0; // Track the actual run number

  // Quick loop through to get the number of runs we have
  for (size_t iev = 0; iev < C->GetEntries(); iev++)
  {
    C->GetEntry(iev);

    if (runnum != rnum)
    {
      rnum = runnum;
      itrrun++;
    }
  } // end loop over events #1


  // Use the number of runs to define histograms that plot again run number or run number index
  TH2F *h2_atimePS_vs_rnum_before = new TH2F("h2_atimePS_vs_rnum_before","BEFORE Alignment;Run Number;Preshower ADC Time [ns]",itrrun,0.5,itrrun+0.5,300,-15,15);
  TH2F *h2_atimePS_vs_rnum_after = new TH2F("h2_atimePS_vs_rnum_after","AFTER Alignment;Run Number;Preshower ADC Time [ns]",itrrun,0.5,itrrun+0.5,300,-15,15);

  TH2F *h2_atimeSH_vs_rnum_before = new TH2F("h2_atimeSH_vs_rnum_before","BEFORE Alignment;Run Number;Shower ADC Time [ns]",itrrun,0.5,itrrun+0.5,300,-15,15);
  TH2F *h2_atimeSH_vs_rnum_after = new TH2F("h2_atimeSH_vs_rnum_after","AFTER Alignment;Run Number;Shower ADC Time [ns]",itrrun,0.5,itrrun+0.5,300,-15,15);

  TH2F *h2_atimeHCal_vs_rnum_before = new TH2F("h2_atimeHCal_vs_rnum_before","BEFORE Alignment;Run Number;HCal ADC Time [ns]",itrrun,0.5,itrrun+0.5,700,60,200);
  TH2F *h2_atimeHCal_vs_rnum_after = new TH2F("h2_atimeHCal_vs_rnum_after","AFTER Alignment;Run Number;HCal ADC Time [ns]",itrrun,0.5,itrrun+0.5,700,60,200);

  TH2F *h2_cointime_vs_rnum_before = new TH2F("h2_cointime_vs_rnum_before","BEFORE Alignment;Run Number;HCal ADC Time - Shower ADC Time [ns]",itrrun,0.5,itrrun+0.5,400,110,150);
  TH2F *h2_cointime_vs_rnum_after = new TH2F("h2_cointime_vs_rnum_after","AFTER Alignment;Run Number;HCal ADC Time - Shower ADC Time [ns]",itrrun,0.5,itrrun+0.5,400,110,150);

  // Start another loop where we will be using the run index to plot against, so we reset this variable to zero
  itrrun = 0;

  // Loop through the entries in the TChain
  for (size_t iev = 0; iev < C->GetEntries(); iev++)
  {
    C->GetEntry(iev);

    // Cast the block IDs as integers
    Int_t psBlkIDint = (Int_t) idPS;
    Int_t shBlkIDint = (Int_t) idSH;
    Int_t hcalBlkIDint = (Int_t) idblkHCAL;

    if (runnum != rnum)
    {
      rnum = runnum;
      itrrun++; // move to the next index when the run number changes
    }

    if(nCut==1||pCut==1) // apply spot cuts
    {
      // fill in adc time vs block id
      h_psTime_before->Fill(idPS,atimePS);
      h_shTime_before->Fill(idSH,atimeSH);
      h_hcalTime_before->Fill(idblkHCAL,hcal_time);

      // fill in 1d adc time
      h_ps_before->Fill(atimePS);
      h_sh_before->Fill(atimeSH);
      h_coin_before->Fill(hcal_time-atimeSH);

      // fill in adc time vs run index
      h2_atimePS_vs_rnum_before->Fill(itrrun,atimePS);
      h2_atimeSH_vs_rnum_before->Fill(itrrun,atimeSH);
      h2_atimeHCal_vs_rnum_before->Fill(itrrun,hcal_time);
      h2_cointime_vs_rnum_before->Fill(itrrun,(hcal_time-atimeSH));

      // Fill a histogram for just the main block ADC time for SH and PS and HCal
      if(idPS==25)
      {
        h_mainPS->Fill(atimePS);
      }
      if(idSH==73)
      {
        h_mainSH->Fill(atimePS);
      }
      if(idblkHCAL==150)
      {
        h_mainHCal->Fill(hcal_time);
      }

      // Loop through the blocks and fill their 1D histos
      for (int i = 0; i < 52; i++)
      {
        if(psBlkIDint == i)
        {
  	       h_atime_psBlocks[i]->Fill(atimePS);
        }
      }
      for (int i = 0; i < 189; i++)
      {
        if(shBlkIDint == i)
        {
  	       h_atime_shBlocks[i]->Fill(atimeSH);
        }
      }
      for (int i = 0; i < 288; i++)
      {
        if(hcalBlkIDint == i)
        {
  	       h_atime_hcalBlocks[i]->Fill(hcal_time);
        }
      }

    } // end spot cuts

  } // end loop over entries

  // Define the mean values of the "main" blocks
  double meanSH = h_mainSH->GetMean();
  double meanPS = h_mainPS->GetMean();
  double meanHCal = h_mainHCal->GetMean();
  //cout << meanHCal << endl;

  // Define lists to store all the offsets from the main block means
  double meanOffsetPS[52];
  double meanOffsetSH[189];
  double meanOffsetHCal[288];

  // Define lists to store all the means and standard deviations for all of the blocks
  double meanAllPS[52];
  double meanAllSH[189];
  double meanAllHCal[288];
  double rmsAllPS[52];
  double rmsAllSH[189];
  double rmsAllHCal[288];

  for (int i = 0; i < 52; i++)
  {
    //meanOffsetPS[i] = (h_atime_psBlocks[i] -> GetMean()) - meanPS;
    meanAllPS[i] = h_atime_psBlocks[i] -> GetMean();
    rmsAllPS[i] = h_atime_psBlocks[i] -> GetRMS();
    meanOffsetPS[i] = meanAllPS[i] - meanPS;
    h_psTime_before_scatter -> Fill(PSblocks[i]+0.5,meanAllPS[i]);
    h_psTime_after_scatter -> Fill(PSblocks[i]+0.5,meanAllPS[i]-meanOffsetPS[i]);
    //cout << endl << meanOffsetPS[i] << endl;
  }
  for (int i = 0; i < 189; i++)
  {
    //meanOffsetSH[i] = (h_atime_shBlocks[i] -> GetMean()) - meanSH;
    meanAllSH[i] = h_atime_shBlocks[i] -> GetMean();
    rmsAllSH[i] = h_atime_shBlocks[i] -> GetRMS();
    meanOffsetSH[i] = meanAllSH[i] - meanSH;
    h_shTime_before_scatter -> Fill(SHblocks[i]+0.5,meanAllSH[i]);
    h_shTime_after_scatter -> Fill(SHblocks[i]+0.5,meanAllSH[i]-meanOffsetSH[i]);
    //cout << endl << meanOffsetSH[i] << endl;
  }
  for (int i = 0; i < 288; i++)
  {
    //meanOffsetHCal[i] = (h_atime_hcalBlocks[i] -> GetMean()) - meanHCal;
    meanAllHCal[i] = h_atime_hcalBlocks[i] -> GetMean();
    rmsAllHCal[i] = h_atime_hcalBlocks[i] -> GetRMS();
    meanOffsetHCal[i] = meanAllHCal[i] - meanHCal;
    //cout << meanOffsetHCal[i] << endl;
  }

  // Assign error bars to the scatter plots based on the sigma of the adc time data
  //for (int bin = 1; bin <= h_psTime_before_scatter->GetNbinsX(); bin++)
  //{
    //h_psTime_before_scatter->SetBinError(bin, rmsAllPS[bin-1]);
  //}

  // Start another loop where we will be using the run index to plot against, so we reset this variable to zero
  itrrun = 0;

  // Loop through the entries in the TChain AGAIN
  for (size_t iev = 0; iev < C->GetEntries(); iev++)
  {
    C->GetEntry(iev);

    // Cast the block IDs as integers
    Int_t psBlkIDint2 = (Int_t) idPS;
    Int_t shBlkIDint2 = (Int_t) idSH;
    Int_t hcalBlkIDint2 = (Int_t) idblkHCAL;

    if (runnum != rnum)
    {
      rnum = runnum;
      itrrun++;
    }

    if(nCut==1||pCut==1) // spot cuts
    {
      // Apply offsets based on block ID and fill histos
      for (int i = 0; i < 52; i++)
      {
        if(psBlkIDint2==i)
        {
          h_psTime_after->Fill(idPS,(atimePS-meanOffsetPS[i]));
          h_ps_after->Fill(atimePS-meanOffsetPS[i]);
          h2_atimePS_vs_rnum_after->Fill(itrrun,atimePS-meanOffsetPS[i]);
        }
      }
      for (int i = 0; i < 189; i++)
      {
        if(shBlkIDint2==i)
        {
          h_shTime_after->Fill(idSH,(atimeSH-meanOffsetSH[i]));
          h_sh_after->Fill(atimeSH-meanOffsetSH[i]);
          h_coin_after->Fill((hcal_time-meanOffsetHCal[i])-(atimeSH-meanOffsetSH[i]));
          h2_atimeSH_vs_rnum_after->Fill(itrrun,atimeSH-meanOffsetSH[i]);
          h2_cointime_vs_rnum_after->Fill(itrrun,((hcal_time-meanOffsetHCal[i])-(atimeSH-meanOffsetSH[i])));
        }
      }
      for (int i = 0; i < 288; i++)
      {
        if(hcalBlkIDint2==i)
        {
          h_hcalTime_after->Fill(idblkHCAL,(hcal_time-meanOffsetHCal[i]));
          h_hcal_after->Fill(hcal_time-meanOffsetHCal[i]);
          h2_atimeHCal_vs_rnum_after->Fill(itrrun,hcal_time-meanOffsetHCal[i]);
        }
      }

    } // end spot cuts

  } // end third loop over entries

  // Make Canvases

  // Confirmation plots to make sure we see something in all the blocks and that I've made my loops correctly - not super useful for interpretting results
  TCanvas *timePSblk = new TCanvas("timePSblk","time of PS blocks", 1000, 1000, 1000, 1000);
  timePSblk -> Divide(2,26);
  for(int i=1; i<53; i++)
  {
    timePSblk -> cd(i);
    h_atime_psBlocks[i-1]->Draw();
  }
  timePSblk -> SaveAs("../plots/preshowerTimes.pdf");

  TCanvas *timeSHblk = new TCanvas("timeSHblk","time of SH blocks", 1000, 1000, 1000, 1000);
  timeSHblk -> Divide(7,27);
  for(int i=1; i<190; i++)
  {
    timeSHblk -> cd(i);
    h_atime_shBlocks[i-1]->Draw();
  }
  timeSHblk -> SaveAs("../plots/showerTimes.pdf");

  // Show times vs block IDs before and after alignment
  TCanvas *psADCtime_BlockID = new TCanvas("psADCtime_BlockID","Preshower ADC time vs block IDs", 1000, 1000, 1000, 1000);
  psADCtime_BlockID -> Divide(1,2);
  psADCtime_BlockID -> cd(1);
  h_psTime_before -> Draw("colz");
  h_psTime_before_scatter->Draw("SAME SCAT E1");
  //gr_psTime_before_scatter->Draw("SAME ALP");
  psADCtime_BlockID -> cd(2);
  h_psTime_after -> Draw("colz");
  h_psTime_after_scatter->Draw("SAME SCAT");
  psADCtime_BlockID -> SaveAs("../plots/preshowerAlignmentResults.pdf");

  TCanvas *shADCtime_BlockID = new TCanvas("shADCtime_BlockID","Shower ADC time vs block IDs", 1000, 1000, 1000, 1000);
  shADCtime_BlockID -> Divide(1,2);
  shADCtime_BlockID -> cd(1);
  h_shTime_before -> Draw("colz");
  h_shTime_before_scatter->Draw("SAME SCAT");
  shADCtime_BlockID -> cd(2);
  h_shTime_after -> Draw("colz");
  h_shTime_after_scatter->Draw("SAME SCAT");
  shADCtime_BlockID -> SaveAs("../plots/showerAlignmentResults.pdf");

  TCanvas *hcalADCtime_BlockID = new TCanvas("hcalADCtime_BlockID","HCal ADC time vs block IDs", 1000, 1000, 1000, 1000);
  hcalADCtime_BlockID -> Divide(1,2);
  hcalADCtime_BlockID -> cd(1);
  h_hcalTime_before -> Draw("colz");
  hcalADCtime_BlockID -> cd(2);
  h_hcalTime_after -> Draw("colz");
  hcalADCtime_BlockID -> SaveAs("../plots/hcalAlignmentResults.pdf");

  TCanvas *ADCtimeAlignment = new TCanvas("ADCtimeAlignment","ADC Time Alignment", 1000, 1000, 1000, 1000);
  ADCtimeAlignment -> Divide(2,2);
  ADCtimeAlignment -> cd(1);
  h_ps_before -> Draw();
  ADCtimeAlignment -> cd(2);
  h_ps_after -> Draw();
  ADCtimeAlignment -> cd(3);
  h_sh_before -> Draw();
  ADCtimeAlignment -> cd(4);
  h_sh_after -> Draw();
  ADCtimeAlignment -> SaveAs("../plots/1D_AlignmentResults.pdf");

  TCanvas *T1DcointimeAlignment = new TCanvas("T1DcointimeAlignment","ADC Time Alignment", 1000, 1000, 1000, 1000);
  T1DcointimeAlignment -> Divide(2,1);
  T1DcointimeAlignment -> cd(1);
  h_coin_before -> Draw();
  T1DcointimeAlignment -> cd(2);
  h_coin_after -> Draw();
  T1DcointimeAlignment -> SaveAs("../plots/1D_coin_AlignmentResults.pdf");

  // Show times vs block IDs before and after alignment
  TCanvas *psADCtime_RunNum = new TCanvas("psADCtime_RunNum","Preshower ADC time vs Run Number", 1000, 1000, 1000, 1000);
  psADCtime_RunNum -> Divide(1,2);
  psADCtime_RunNum -> cd(1);
  h2_atimePS_vs_rnum_before -> Draw("colz");
  psADCtime_RunNum -> cd(2);
  h2_atimePS_vs_rnum_after -> Draw("colz");
  psADCtime_RunNum -> SaveAs("../plots/preshowerAlignmentResults_RunNum.pdf");

  TCanvas *shADCtime_RunNum = new TCanvas("shADCtime_RunNum","Shower ADC time vs Run Number", 1000, 1000, 1000, 1000);
  shADCtime_RunNum -> Divide(1,2);
  shADCtime_RunNum -> cd(1);
  h2_atimeSH_vs_rnum_before -> Draw("colz");
  shADCtime_RunNum -> cd(2);
  h2_atimeSH_vs_rnum_after -> Draw("colz");
  shADCtime_RunNum -> SaveAs("../plots/showerAlignmentResults_RunNum.pdf");

  TCanvas *hcalADCtime_RunNum = new TCanvas("hcalADCtime_RunNum","HCal ADC time vs Run Number", 1000, 1000, 1000, 1000);
  hcalADCtime_RunNum -> Divide(1,2);
  hcalADCtime_RunNum -> cd(1);
  h2_atimeHCal_vs_rnum_before -> Draw("colz");
  hcalADCtime_RunNum -> cd(2);
  h2_atimeHCal_vs_rnum_after -> Draw("colz");
  hcalADCtime_RunNum -> SaveAs("../plots/hcalAlignmentResults_RunNum.pdf");

  TCanvas *cointime_RunNum = new TCanvas("cointime_RunNum","Coincidence time vs Run Number", 1000, 1000, 1000, 1000);
  cointime_RunNum -> Divide(1,2);
  cointime_RunNum -> cd(1);
  h2_cointime_vs_rnum_before -> Draw("colz");
  cointime_RunNum -> cd(2);
  h2_cointime_vs_rnum_after -> Draw("colz");
  cointime_RunNum -> SaveAs("../plots/coinAlignmentResults_RunNum.pdf");

} // end main function

TH1F* MakeHisto(Int_t index, char const * suf, Int_t bins=300, Double_t min=-15., Double_t max=15.)
{
  Int_t blkid = index;
  TH1F *h = new TH1F(TString::Format("h_Blk_%d%s",blkid,suf),"",bins,min,max);
  //h->SetStats(0);
  h->SetLineWidth(1);
  h->SetLineColor(2);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetMaxDigits(3);
  h->GetYaxis()->SetNdivisions(5);
  return h;
}
