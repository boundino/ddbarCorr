#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

namespace xjjc
{
  void progressbar(int index_, int total_, int morespace_=0);
  void progressbar_summary(int total_);
  bool str_contains(std::string str1, std::string str2) { return str1.find(str2)!=std::string::npos; }
}

void merge(std::string outputname, std::string filelist, std::string ntname = "ntDkpi", bool isskim = false, int ntotal = -1)
{
  TTree::SetMaxTreeSize(1LL * 1024 * 1024 * 1024 * 1024);

  std::vector<std::string> files;
  std::ifstream getfile(filelist.c_str());
  int count = 0;
  for(int i=0; (i<ntotal || ntotal<0); i++)
    {
      std::string filename;
      getfile >> filename;
      if(getfile.eof()) { break; }
      files.push_back(filename);
      count++;
    }

  bool isbnt = (ntname!="empty"); bool isnt = xjjc::str_contains(ntname, "nt");

  // TChain* hlt = new TChain("hltanalysis/HltTree");
  // TChain* hltobj = new TChain("hltobject/HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v");
  TChain* skim = new TChain("skimanalysis/HltTree");
  TChain* hi = new TChain("hiEvtAnalyzer/HiTree");
  // TChain* forest = new TChain("HiForest/HiForestInfo");
  TChain *nt = 0;
  if(isbnt) 
    {
      nt = new TChain(Form("Dfinder/%s", ntname.c_str()));
    }
  for(auto& i : files)
    {
      skim->Add(i.c_str());
      hi->Add(i.c_str());
      // forest->Add(i.c_str());
      if(isbnt) 
        {
          nt->Add(i.c_str());
        }
      std::cout<<i<<std::endl;
    }
  std::cout<<"Merged \e[31;1m"<<count<<"\e[0m files"<<std::endl;

  if(ntname == "root")
    {
      nt->SetBranchStatus("*", 0);
      nt->SetBranchStatus("TrackInfo.*", 1);
    }
  else if(isbnt)
    {
      // blind branch
      nt->SetBranchStatus("Dtrk3*", 0);
      nt->SetBranchStatus("Dtrk4*", 0);
      nt->SetBranchStatus("Dtktk*", 0);
      nt->SetBranchStatus("DRes*", 0);
      nt->SetBranchStatus("Dgen*", 0);
      nt->SetBranchStatus("Dgen", 1);
      nt->SetBranchStatus("DsGen", 0);
      nt->SetBranchStatus("D_unfitted_mass", 0);
      nt->SetBranchStatus("D_unfitted_pt", 0);
      nt->SetBranchStatus("Dvtx*", 0);
      nt->SetBranchStatus("Dtrk1theta*", 0);
      nt->SetBranchStatus("Dtrk2theta*", 0);
      nt->SetBranchStatus("Dtrk*MVAVal", 0);
      nt->SetBranchStatus("Dtrk*PhiErr", 0);
      nt->SetBranchStatus("Dtrk*EtaErr", 0);
      nt->SetBranchStatus("Dtrk*Idx", 0);
      nt->SetBranchStatus("Dtrk*Quality", 0);
      nt->SetBranchStatus("Dtrk*Y", 0);
      nt->SetBranchStatus("Dd0*", 0);
      nt->SetBranchStatus("Ddxyz*", 0);
      nt->SetBranchStatus("Dindex*", 0);
      nt->SetBranchStatus("BS*", 0);
      nt->SetBranchStatus("DisSequentialFit", 0);
      nt->SetBranchStatus("Dtrk*Algo", 0);
      nt->SetBranchStatus("Dtrk*originalAlgo", 0);
      // nt->SetBranchStatus("", 0);
      // nt->SetBranchStatus("", 0);
      // nt->SetBranchStatus("", 0);
      // nt->SetBranchStatus("", 0);
      // nt->SetBranchStatus("", 0);
      // nt->SetBranchStatus("", 0);
    }

  TFile* outf = new TFile(outputname.c_str(), "recreate");
  TDirectory* dskim = outf->mkdir("skimanalysis","");
  TDirectory* dhi = outf->mkdir("hiEvtAnalyzer","");
  // TDirectory* dforest = outf->mkdir("HiForest","");
  TDirectory* dDfinder = 0;
  if(isbnt) dDfinder = outf->mkdir("Dfinder","");
      
  dskim->cd();
  TTree* new_skim = skim->CloneTree(0);
  dhi->cd();
  TTree* new_hi = hi->CloneTree(0);
  // dforest->cd();
  // TTree* new_forest = forest->CloneTree(0);
  int Dsize = 0;
  TTree *new_nt = 0;
  if(isbnt)
    {
      dDfinder->cd();
      new_nt = nt->CloneTree(0);
      nt->SetBranchAddress("Dsize", &Dsize);
    }

  Long64_t nentries = hi->GetEntries();
  std::cout<<" -- Event reading"<<std::endl;
  Long64_t j = 0;
  for(Long64_t i=0;i<nentries;i++)
    {
      if(i%10000==0) { xjjc::progressbar(i, nentries); }

      skim->GetEntry(i);
      hi->GetEntry(i);
      // forest->GetEntry(i);
      if(isbnt) 
        {
          nt->GetEntry(i);
        }

      if(isskim)
        {
          if(isbnt && xjjc::str_contains(ntname, "nt") && !Dsize) continue; //
        }

      dskim->cd();
      new_skim->Fill();
      dhi->cd();
      new_hi->Fill();
      // dforest->cd();
      // new_forest->Fill();
      if(isbnt)
        {
          dDfinder->cd();
          new_nt->Fill();
        }
    }
  xjjc::progressbar_summary(nentries);
      
  outf->Write();
  std::cout<<" -- Writing new trees done"<<std::endl;
  outf->Close();
      
}

int main(int argc, char* argv[])
{
  if(argc==6) { merge(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5])); return 0; }
  if(argc==5) { merge(argv[1], argv[2], argv[3], atoi(argv[4])); return 0; }
  std::cout<<__FUNCTION__<<": ./merge.exe [outputname] [filelist] [treename](ntKp, ntphi, ntmix, root, empty) [skim] (optional)[number of events]";
  return 1;
}

void xjjc::progressbar(int index_, int total_, int morespace_/*=0*/)
{
  std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10+morespace_)<<index_<<"\033[0m"<<" / "<<std::setw(10+morespace_)<<total_<<" ] "<<"\033[1;36m"<<(int)(100.*index_/total_)<<"%\033[0m"<<"\r"<<std::flush;
}

void xjjc::progressbar_summary(int total_)
{
  std::cout<<std::endl<<"  Processed "<<"\033[1;31m"<<total_<<"\033[0m event(s)."<<std::endl;
}

// hlt->SetBranchStatus("HLT_HIIslandPhoton*", 0);
// hlt->SetBranchStatus("HLT_HIMinimumBias*", 0);
// hlt->SetBranchStatus("HLT_HIL3Mu3Eta2p5*", 0);
// hlt->SetBranchStatus("HLT_HICastor*", 0);
// hlt->SetBranchStatus("HLT_HIGEDPhoton*", 0);
// hlt->SetBranchStatus("HLT_HIPuAK4CaloJet*", 0);
// hlt->SetBranchStatus("HLT_HIL1NotBptxOR*", 0);
// hlt->SetBranchStatus("HLT_HIL1UnpairedBunch*", 0);
// hlt->SetBranchStatus("HLTriggerFinal*", 0);
// hlt->SetBranchStatus("L1_*", 0);

// int HLT_HIL1DoubleMuOpen_v1; hlt->SetBranchAddress("HLT_HIL1DoubleMuOpen_v1", &HLT_HIL1DoubleMuOpen_v1);
// int HLT_HIL1DoubleMu10_v1; hlt->SetBranchAddress("HLT_HIL1DoubleMu10_v1", &HLT_HIL1DoubleMu10_v1);
// if(!HLT_HIL1DoubleMuOpen_v1 && !HLT_HIL1DoubleMu10_v1 && !HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1) continue; //
