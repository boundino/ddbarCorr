#include "xjjrootuti.h"
#include "xjjcuti.h"
#include <TTree.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <TH2D.h>

#ifndef MAX_CAND
#define MAX_CAND 200
#endif

namespace ddtree
{
  class dtree
  {
  public:
    dtree(TTree* nt) : fnt(nt) { buildbranches(); }
    // dtree(TTree* nt, dtree* idtree) : fnt(nt) { buildbranches(); }
    dtree(TFile* inf, std::string treename) { fnt = (TTree*)inf->Get(treename.c_str()); setbranchesaddress(); }
    // branches
    int Ntrkoffline;
    float HFsumETPlus;
    float HFsumETMinus;
    float bestvtxX;
    float bestvtxY;
    float bestvtxZ;
    int centrality;
    int candSize;
    float pT[MAX_CAND];
    float y[MAX_CAND];
    float phi[MAX_CAND];
    float mass[MAX_CAND];
    float flavor[MAX_CAND];
    float dca[MAX_CAND];
    float mva[MAX_CAND];
    TTree* nt() { return fnt; }
    void fillcands(dtree* idgree);
    void fillcands(dtree* idgree, int j);
  private:
    TTree* fnt;
    void buildbranches();
    void setbranchesaddress();
  };
}

void ddtree::dtree::buildbranches()
{
  fnt->Branch("Ntrkoffline", &Ntrkoffline);
  fnt->Branch("HFsumETPlus", &HFsumETPlus);
  fnt->Branch("HFsumETMinus", &HFsumETMinus);
  fnt->Branch("bestvtxX", &bestvtxX);
  fnt->Branch("bestvtxY", &bestvtxY);
  fnt->Branch("bestvtxZ", &bestvtxZ);
  fnt->Branch("centrality", &centrality);
  fnt->Branch("candSize", &candSize);
  fnt->Branch("pT", pT, "pT[candSize]/F");
  fnt->Branch("y", y, "y[candSize]/F");
  fnt->Branch("phi", phi, "phi[candSize]/F");
  fnt->Branch("mass", mass, "mass[candSize]/F");
  fnt->Branch("flavor", flavor, "flavor[candSize]/F");
  fnt->Branch("dca", dca, "dca[candSize]/F");
  fnt->Branch("mva", mva, "mva[candSize]/F");
}

void ddtree::dtree::setbranchesaddress()
{
  fnt->SetBranchAddress("Ntrkoffline", &Ntrkoffline);
  fnt->SetBranchAddress("HFsumETPlus", &HFsumETPlus);
  fnt->SetBranchAddress("HFsumETMinus", &HFsumETMinus);
  fnt->SetBranchAddress("bestvtxX", &bestvtxX);
  fnt->SetBranchAddress("bestvtxY", &bestvtxY);
  fnt->SetBranchAddress("bestvtxZ", &bestvtxZ);
  fnt->SetBranchAddress("centrality", &centrality);
  fnt->SetBranchAddress("candSize", &candSize);
  fnt->SetBranchAddress("pT", pT);
  fnt->SetBranchAddress("y", y);
  fnt->SetBranchAddress("phi", phi);
  fnt->SetBranchAddress("mass", mass);
  fnt->SetBranchAddress("flavor", flavor);
  fnt->SetBranchAddress("dca", dca);
  fnt->SetBranchAddress("mva", mva);
}

void ddtree::dtree::fillcands(dtree* idtree)
{
  Ntrkoffline = idtree->Ntrkoffline;
  HFsumETPlus = idtree->HFsumETPlus;
  HFsumETMinus = idtree->HFsumETMinus;
  bestvtxX = idtree->bestvtxX;
  bestvtxY = idtree->bestvtxY;
  bestvtxZ = idtree->bestvtxZ;
  centrality = idtree->centrality;
}

void ddtree::dtree::fillcands(dtree* idtree, int j)
{
  pT[candSize] = idtree->pT[j];
  y[candSize] = idtree->y[j];
  phi[candSize] = idtree->phi[j];
  mass[candSize] = idtree->mass[j];
  flavor[candSize] = idtree->flavor[j];
  dca[candSize] = idtree->dca[j];
  mva[candSize] = idtree->mva[j];
  candSize++;
}

namespace ddtree
{
  class mvas
  {
  public:
    mvas(std::string inputname, std::string name="");
    bool passmva(float pt, float mva);
  private:
    TH2D* fhmva;
    std::vector<float> fbounds, fvalues;
  };
}

ddtree::mvas::mvas(std::string inputname, std::string name)
{
  TFile* inf = TFile::Open(inputname.c_str());
  fhmva = (TH2D*)inf->Get("hist_bdtcut");
  fhmva->SetName(Form("%s_%s", fhmva->GetName(), name.c_str()));
  int n = fhmva->GetYaxis()->GetNbins();
  for(int i=0;i<n;i++)
    {
      fbounds.push_back(fhmva->GetYaxis()->GetBinLowEdge(i+1));
      if(i==n-1) fbounds.push_back(fhmva->GetYaxis()->GetBinUpEdge(i+1));
      fvalues.push_back(fhmva->GetBinContent(1, i+1));
    }
}

bool ddtree::mvas::passmva(float pt, float mva)
{
  int j=-1;
  for(int i=0; i<fbounds.size()-1; i++)
    {
      if(pt < fbounds[i]) break;
      j++;
    }
  if(j < 0) return false; // !
  return mva > fvalues[j];
}
