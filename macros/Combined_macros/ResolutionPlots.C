#include <iostream>
#include <sstream>
#include <fstream>

#include "TCut.h"
#include "TH2F.h"
#include "TMath.h"
#include "TChain.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"

const Int_t kNcolsSH = 7;   // SH columns
const Int_t kNrowsSH = 27;  // SH rows
const Int_t kNblksSH = 189; // Total # SH blocks/PMTs
const Int_t kNcolsPS = 2;   // PS columns
const Int_t kNrowsPS = 26;  // PS rows
const Int_t kNblksPS = 52;  // Total # PS blocks/PMTs

const Double_t Mp = 0.938272081;  // +/- 6E-9 GeV

vector<Long64_t> goodevents;

//string getDate();
//void Custm1DHisto(TH1F*);
TString GetOutFileBase(TString);
//void ReadOffset(TString, Double_t*);
//TH1F* MakeHisto(Int_t, Int_t, char const *, Int_t, Double_t, Double_t);
std::vector<std::string> SplitString(char const delim, std::string const myStr);
//void Custm2DRnumHisto(TH2F*, std::vector<std::string> const & lrnum);

void ResolutionPlots(char const * configfilename)
{

  gErrorIgnoreLevel = kError; // ignore all root warnings
  gStyle -> SetOptFit(1);

  TChain *C = new TChain("T");
  TString cfgfilebase = GetOutFileBase(configfilename);

  TString exp = "unknown";
  Int_t config = -1;
  TString set = "N/A";

  ifstream configfile(configfilename);
  TString currentline;

  while(currentline.ReadLine( configfile ) && !currentline.BeginsWith("endRunlist"))
  {
    if(!currentline.BeginsWith("#"))
    {
      C->Add(currentline);
    }// end if
  }// end while

  TCut globalcut = "";
  TString gcutstr;

  while(currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut"))
  {
    if(!currentline.BeginsWith("#") )
    {
      globalcut += currentline;
      gcutstr += currentline;
    }// end if
  }// end while

  std::vector<std::string> gCutList = SplitString('&', gcutstr.Data());

  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  while( currentline.ReadLine( configfile ) )
  {
    if( currentline.BeginsWith("#") ) continue;

    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();

    if( ntokens>1 )
    {
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if (skey == "exp") exp = ((TObjString*)(*tokens)[1])->GetString();
      if (skey == "config") config = ((TObjString*)(*tokens)[1])->GetString().Atoi();
      if (skey == "set") set = ((TObjString*)(*tokens)[1])->GetString();
      if( skey == "*****" )
      {
	       break;
      }// end if
    }// end if ntokens

    delete tokens;
  }// end while

  // Check for empty rootfiles and set tree branches
  if(C->GetEntries()==0) {cerr << endl << " --- No ROOT file found!! --- " << endl << endl; exit(1);}
  else cout << endl << "Found " << C->GetEntries() << " events. Starting analysis.. " << endl;

  // Setting useful ROOT tree branch addresses
  C->SetBranchStatus("*",0);
  Int_t maxtr=500;

  //shower
  Double_t sh_nclus;              C->SetBranchStatus("bb.sh.nclus",1); C->SetBranchAddress("bb.sh.nclus",&sh_nclus);
  Double_t sh_e;                  C->SetBranchStatus("bb.sh.e",1); C->SetBranchAddress("bb.sh.e",&sh_e);
  Double_t sh_clblk_atime[maxtr]; C->SetBranchStatus("bb.sh.clus_blk.atime",1); C->SetBranchAddress("bb.sh.clus_blk.atime",&sh_clblk_atime);
  Double_t sh_nblk;               C->SetBranchStatus("bb.sh.nblk",1); C->SetBranchAddress("bb.sh.nblk",&sh_nblk);
  Double_t sh_clus_nblk;          C->SetBranchStatus("bb.sh.clus.nblk",1); C->SetBranchAddress("bb.sh.clus.nblk",&sh_clus_nblk);

  //preshower
  Double_t ps_nclus;              C->SetBranchStatus("bb.ps.nclus",1); C->SetBranchAddress("bb.ps.nclus",&ps_nclus);
  Double_t ps_e;                  C->SetBranchStatus("bb.ps.e",1); C->SetBranchAddress("bb.ps.e",&ps_e);
  Double_t ps_clblk_atime[maxtr]; C->SetBranchStatus("bb.ps.clus_blk.atime",1); C->SetBranchAddress("bb.ps.clus_blk.atime",&ps_clblk_atime);
  Double_t ps_nblk;               C->SetBranchStatus("bb.ps.nblk",1); C->SetBranchAddress("bb.ps.nblk",&ps_nblk);
  Double_t ps_clus_nblk;          C->SetBranchStatus("bb.ps.clus.nblk",1); C->SetBranchAddress("bb.ps.clus.nblk",&ps_clus_nblk);

  //global cut branches
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.gem.track.nhits",1);
  Double_t tg_th[maxtr]; C->SetBranchStatus("bb.tr.tg_th",1); C->SetBranchAddress("bb.tr.tg_th",&tg_th);
  Double_t tg_ph[maxtr]; C->SetBranchStatus("bb.tr.tg_ph",1); C->SetBranchAddress("bb.tr.tg_ph",&tg_ph);

  //event info
  C->SetMakeClass(1);
  C->SetBranchStatus("fEvtHdr.*", 1);
  UInt_t rnum;           C->SetBranchAddress("fEvtHdr.fRun", &rnum);
  UInt_t trigbits;       C->SetBranchAddress("fEvtHdr.fTrigBits", &trigbits);
  ULong64_t gevnum;      C->SetBranchAddress("fEvtHdr.fEvtNum", &gevnum);

  //define histograms
  TH1F *h_sh_time_resolution = new TH1F("h_sh_time_resolution", "BBCal Shower Time Resolution", 2000, -10., 10.);
  TH1F *h_ps_time_resolution = new TH1F("h_ps_time_resolution", "BBCal Preshower Time Resolution", 2000, -10., 10.);

  TH1F *h_blocks_in_cluster = new TH1F("h_blocks_in_cluster", "Number of SH and PS blocks in Each Cluster", 20, 0., 20.);

  TH1F *h_num_clusters = new TH1F("h_num_clusters", "Number of SH clusters in Each Event", 6, 0., 6.);

  // Loop over all events //

  Long64_t nevent=0, nevents=C->GetEntries();
  UInt_t runnum=0;

  int treenum = 0;
  int currenttreenum = 0;
  int itrrun=0;
  std::vector<std::string> lrnum;    // list of run numbers
  lrnum.reserve(100);

  while(C->GetEntry(nevent++))
  {
    //apply global cuts using andrew's method
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum)
    {
      treenum = currenttreenum;
      GlobalCut -> UpdateFormulaLeaves();

      //track change of runnum
      if (nevent == 1 || rnum != runnum)
      {
        runnum = rnum;
        itrrun++;
        lrnum.push_back(to_string(rnum));
      }// end if runnum

    }// end if nevent == 1

    bool passedgCut = GlobalCut -> EvalInstance(0) != 0;

    if (passedgCut)
    {
      goodevents.push_back(nevent);

      h_sh_time_resolution -> Fill(sh_clblk_atime[0]-sh_clblk_atime[1]);
      h_ps_time_resolution -> Fill(ps_clblk_atime[0]-ps_clblk_atime[1]);

      h_blocks_in_cluster -> Fill(sh_clus_nblk + ps_clus_nblk);

      h_num_clusters -> Fill(sh_nclus);

    }// end if passed cut

  }// end while through events

  TF1 *sh_fgaus = new TF1("sh_fgaus","gaus");
  TF1 *ps_fgaus = new TF1("ps_fgaus","gaus");

  double sh_Mean = h_sh_time_resolution -> GetMean();
  double sh_stdDev = h_sh_time_resolution -> GetStdDev();
  double sh_maxCount = h_sh_time_resolution -> GetMaximum();
  double sh_maxBinCenter = h_sh_time_resolution -> GetXaxis() -> GetBinCenter(h_sh_time_resolution -> GetMaximumBin());

  double ps_Mean = h_ps_time_resolution -> GetMean();
  double ps_stdDev = h_ps_time_resolution -> GetStdDev();
  double ps_maxCount = h_ps_time_resolution -> GetMaximum();
  double ps_maxBinCenter = h_ps_time_resolution -> GetXaxis() -> GetBinCenter(h_ps_time_resolution -> GetMaximumBin());

  sh_fgaus -> SetParameters(sh_maxCount, sh_maxBinCenter, sh_stdDev);
  sh_fgaus -> SetRange(-5.0, 5.0);

  ps_fgaus -> SetParameters(ps_maxCount, ps_maxBinCenter, ps_stdDev);
  ps_fgaus -> SetRange(-5.0, 5.0);

  h_sh_time_resolution -> Fit(sh_fgaus, "RQ");
  //gPad -> Modified();
  //gPad -> Update();

  for(int i = 0; i < 3; i++)
  {
    h_sh_time_resolution -> Fit(sh_fgaus, "RQ", "", sh_fgaus->GetParameter(1)-1.0*sh_fgaus->GetParameter(2), sh_fgaus->GetParameter(1)+1.0*sh_fgaus->GetParameter(2));
  }

  gPad -> Modified();
  gPad -> Update();
  
  h_ps_time_resolution -> Fit(ps_fgaus, "RQ");
  //gPad -> Modified();
  //gPad -> Update();

  for(int j = 0; j < 3; j++)
  {
    h_ps_time_resolution -> Fit(ps_fgaus, "RQ", "", ps_fgaus->GetParameter(1)-1.0*ps_fgaus->GetParameter(2), ps_fgaus->GetParameter(1)+1.0*ps_fgaus->GetParameter(2));
  }

  gPad -> Modified();
  gPad -> Update();

  TCanvas *cTimeResolution = new TCanvas ("cTimeResolution", "Time Resolution", 1200, 800);
  cTimeResolution -> Divide(1,2);
  cTimeResolution -> cd(1);
  h_sh_time_resolution -> Draw();
  cTimeResolution -> cd(2);
  h_ps_time_resolution -> Draw();
  cTimeResolution -> SaveAs("Output/TimeResolutionPeaks.pdf");

  TCanvas *cBlocksCluster = new TCanvas ("cBlocksCluster", "Blocks in Each Cluster", 1200, 800);
  h_blocks_in_cluster -> Draw();
  cBlocksCluster -> SaveAs("Output/BlocksCluster.pdf");

  TCanvas *cNumClusters = new TCanvas ("cNumClusters", "Number of Clusters in Each Event", 1200, 800);
  h_num_clusters -> Draw();
  cNumClusters -> SaveAs("Output/NumClusters.pdf");

}// end main

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
