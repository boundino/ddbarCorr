#include "TTree.h"
#include "dtree.h"

#include "xjjcuti.h"

#include <vector>
#include <string>
#include <iostream>

class EventMixer
{
  public:
  EventMixer(ddtree::dtree* t, int hiBin_width=10, float vz_width=10., float vz_min=-30., float vz_max=30.);
  ~EventMixer();
  //EventMixer(std::string filename = "mixedTrees.root");
  //void SaveToFile(std::string filename = "mixedTrees.root");
  //void LoadFromFile(std::string filename = "mixedTrees.root");
  void buildTreeBins(ddtree::dtree* t);
  ddtree::dtree* findSimilarEvent(int hiBin, float vz, size_t original_evt);
  private:
  TFile* f;
  Int_t _hiBin_width;
  Float_t _vz_width;
  Float_t _vz_min;
  Float_t _vz_max;
  std::vector<std::vector<ddtree::dtree*>> dtreeBins;
  std::vector<std::vector<TTree*>> trees;
  std::vector<std::vector<std::vector<size_t>>> entrylists;
  std::vector<std::vector<int>> cursors;
};

EventMixer::EventMixer(ddtree::dtree* t, int hiBin_width, float vz_width, float vz_min, float vz_max)
{
  f = new TFile("/data/mjpeters/mixTrees.root","recreate");
  _hiBin_width = hiBin_width;
  _vz_width = vz_width;
  _vz_min = vz_min;
  _vz_max = vz_max;
  buildTreeBins(t);
}

EventMixer::~EventMixer()
{
  f->Close();
  delete f;
}
/*
EventMixer::EventMixer(std::string filename)
{
  LoadFromFile(filename);
}

void EventMixer::LoadFromFile(std::string filename)
{
  TFile* f = TFile::Open(filename.c_str());
  _hiBin_width = (int)f->Get("hiBin_width");
  _vz_width = (float)f->Get("vz_width");
  _vz_min = (float)f->Get("vz_min");
  _vz_max = (float)f->Get("vz_max");
  dtreeBins = (std::vector<std::vector<ddtree::dtree*>>)f->Get("dtreeBins");
  trees = (std::vector<std::vector<TTree*>>)f->Get("trees");
  entrylists = (std::vector<std::vector<std::vector<size_t>>>)f->Get("entrylists");
  cursors = (std::vector<std::vector<int>>)f->Get("cursors");
}

void EventMixer::SaveToFile(std::string filename)
{
  TFile* f = new TFile(filename.c_str());
  _hiBin_width.Write();
  _vz_width.Write();
  _vz_min.Write();
  _vz_max.Write();
  dtreeBins.Write();
  trees.Write();
  entrylists.Write();
  cursors.Write();
}
*/

void EventMixer::buildTreeBins(ddtree::dtree* t)
{
  // setup 2D vectors of TTrees and dtrees
  dtreeBins.clear();
  int nhiBins = 200/_hiBin_width;
  int nvzBins = ceil((_vz_max-_vz_min)/_vz_width);
  trees.resize(nhiBins);
  dtreeBins.resize(nhiBins);
  entrylists.resize(nhiBins);
  cursors.resize(nhiBins);
  for(int i=0;i<nhiBins;i++)
  {
    dtreeBins[i].resize(nvzBins);
    trees[i].resize(nvzBins);
    entrylists[i].resize(nvzBins);
    cursors[i].resize(nvzBins);
  }

  // copy structure for all TTrees, and initialize cursors
  for(int i=0;i<nhiBins;i++)
  {
    for(int j=0;j<nvzBins;j++)
    {
      trees[i][j] = t->nt()->CloneTree(0);
      std::string tname = t->nt()->GetName();
      trees[i][j]->SetName((tname+std::to_string(i)+"_"+std::to_string(j)).c_str());
      cursors[i][j] = 1;
    }
  }

  // loop over original TTree, fill appropriate TTree
  size_t nentries = t->nt()->GetEntries();
  std::cout << "Building event-mixing TTrees...\n";
  for(size_t evt=0;evt<nentries;evt++)
  {
    if(evt % 10000 == 0) xjjc::progressbar(evt,nentries);
    t->nt()->GetEntry(evt);
    int hiBin_index = floor(t->centrality/_hiBin_width);
    int vz_index = floor((t->bestvtxZ-_vz_min)/_vz_width);
    trees[hiBin_index][vz_index]->Fill();
    entrylists[hiBin_index][vz_index].push_back(evt);
  }
  
  // make dtrees out of the TTrees
  for(int i=0;i<nhiBins;i++)
  {
    for(int j=0;j<nvzBins;j++)
    {
      dtreeBins[i][j] = new ddtree::dtree(f,(t->nt()->GetName()+std::to_string(i)+"_"+std::to_string(j)).c_str());
    }
  }
  f->Write();
}

ddtree::dtree* EventMixer::findSimilarEvent(int hiBin, float vz, size_t originalentry)
{
  int nhiBins = ceil(200/_hiBin_width);
  int nvzBins = ceil((_vz_max-_vz_min)/_vz_width);
  int hiBin_index = floor(hiBin/_hiBin_width);
  int vz_index = floor((vz-_vz_min)/_vz_width);
  int cursor = cursors[hiBin_index][vz_index];
  int nentries = dtreeBins[hiBin_index][vz_index]->nt()->GetEntries();
  if(originalentry == entrylists[hiBin_index][vz_index][cursor])
  {
    cursors[hiBin_index][vz_index]++;
    cursor = cursors[hiBin_index][vz_index];
  }
  if(cursor==nentries-1) cursors[hiBin_index][vz_index] = 0;
  else cursors[hiBin_index][vz_index]++;
  dtreeBins[hiBin_index][vz_index]->nt()->GetEntry(cursor);
  return dtreeBins[hiBin_index][vz_index];
}
