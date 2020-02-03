#include "dtree.h"

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>

#include "xjjcuti.h"

void skim(std::string inputname, std::string outputname)
{
  //
  TFile* inf = TFile::Open(inputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->mkdir("d0ana_mc_genmatchunswap");
  outf->cd("d0ana_mc_genmatchunswap");
  TTree* nt = new TTree("VertexCompositeNtuple", "");
  ddtree::dtree* idtree = new ddtree::dtree(inf, "d0ana_mc_genmatchunswap/VertexCompositeNtuple");
  ddtree::dtree* odtree = new ddtree::dtree(nt);
  ddtree::gtree* igtree = new ddtree::gtree(inf, "d0ana_mc_genmatchunswap/VertexCompositeNtuple");
  ddtree::gtree* ogtree = new ddtree::gtree(nt);

  //
  std::vector<float> cents = {0, 30, 50, 80};
  std::vector<ddtree::mvas*> mvacents(cents.size()-1);
  for(int l=0; l<cents.size()-1; l++)
    { mvacents[l] = new ddtree::mvas(Form("../mvas/BDTCuts_finePtBins_PrompD0_PbPb_cent%.0fto%.0f_trkPt1GeV.root", cents[l], cents[l+1]), 
                                     Form("cent%.0fto%.0f", cents[l], cents[l+1])); }

  int nentries = idtree->nt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      if(i%10000==0) { xjjc::progressbar(i, nentries); }
      idtree->nt()->GetEntry(i);
      int lcent = -1;
      for(int l=0; l<cents.size(); l++) { if(idtree->centrality < cents[l]) { break; } lcent++; }

      odtree->candSize = 0;
      odtree->fillcands(idtree);
      for(int j=0; j<idtree->candSize; j++)
        {
          if(lcent<cents.size()-1 && !mvacents[lcent]->passmva(idtree->pT[j], idtree->mva[j])) continue;
          odtree->fillcands(idtree, j);
        }
      ogtree->candSize_gen = 0;
      for(int j=0; j<igtree->candSize_gen; j++)
        {
          if(!(abs(igtree->DauID1_gen[j])==321 && abs(igtree->DauID2_gen[j])==211) && 
             !(abs(igtree->DauID1_gen[j])==211 && abs(igtree->DauID2_gen[j])==321)) continue;
          ogtree->fillcands(igtree, j);
        }

      odtree->nt()->Fill();
    }
  xjjc::progressbar_summary(nentries);
  outf->cd("d0ana_mc_genmatchunswap");
  odtree->nt()->Write();
}

int main(int argc, char* argv[])
{
  if(argc==3) { skim(argv[1], argv[2]); return 0; }
  return 1;
}

