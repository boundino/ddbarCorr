#ifndef EVENT_MIXER_H
#define EVENT_MIXER_H

#include "TTree.h"
#include "TParameter.h"
#include "TFile.h"
#include "dtree_pp.h"

#include "xjjcuti.h"

#include <vector>
#include <string>
#include <iostream>
#include <random>

class EventMixer
{
  public:
  EventMixer(ddtree::dtree* t, float vz_width=0.2, float vz_min=-4., float vz_max=4.);
  ~EventMixer();
  //EventMixer(std::string filename = "mixedTrees.root");
  //void SaveToFile(std::string filename = "mixedTrees.root");
  //void LoadFromFile(std::string filename = "mixedTrees.root");
  void buildTreeBins(ddtree::dtree* t);
  ddtree::dtree* findSimilarEvent(float vz, size_t original_evt);
  private:
  TFile* f;
  std::string _fname;
  float _vz_width;
  float _vz_min;
  float _vz_max;
  std::vector<ddtree::dtree*> dtreeBins;
  std::vector<TTree*> trees;
  std::vector<int> cursors;
  std::default_random_engine r;
};

EventMixer::EventMixer(ddtree::dtree* t, float vz_width, float vz_min, float vz_max)
{
  _vz_width = vz_width;
  _vz_min = vz_min;
  _vz_max = vz_max;
  _fname = "/data/mjpeters/mixTrees_pp.root";
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
  // check to see if previous TFile exists
  bool sameFileExists = false;
  int nvzBins = ceil((_vz_max-_vz_min)/_vz_width);
  std::cout << "about to open\n";
  f = TFile::Open(_fname.c_str(),"read");
  std::cout << "opened\n";
  if(!f)
  {
    std::cout << "was zombie\n";
  }
  else
  {
    std::cout << "was not zombie\n";
    TParameter<float>* vzw = (TParameter<float>*)f->Get("vz_width");
    TParameter<float>* vzmn = (TParameter<float>*)f->Get("vz_min");
    TParameter<float>* vzmx = (TParameter<float>*)f->Get("vz_max");
    float vzwidth_old; 
    if(!vzw) vzwidth_old = 1e9;
    else vzwidth_old = vzw->GetVal();
    float vzmin_old; 
    if(!vzmn) vzmin_old = 1e9;
    else vzmin_old = vzmn->GetVal();
    float vzmax_old; 
    if(!vzmx) vzmax_old = 1e9;
    else vzmax_old = vzmx->GetVal();
    std::cout << "vzwidth_old: " << vzwidth_old << "\n";
    std::cout << "vzmin_old: " << vzmin_old << "\n";
    std::cout << "vzmax_old: " << vzmax_old << "\n";
    if(fabs(_vz_width-vzwidth_old)<1e-4 && fabs(_vz_min-vzmin_old)<1e-4 && fabs(_vz_max-vzmax_old)<1e-4)
    {
      std::cout << "file already exists\n";
      sameFileExists = true;
      trees.resize(nvzBins);
      dtreeBins.resize(nvzBins);
      cursors.resize(nvzBins);
      for(int j=0;j<nvzBins;j++)
      {
        std::string tname = t->nt()->GetName();
        std::string jname = tname+"_"+std::to_string(j);
        trees[j] = (TTree*)f->Get(jname.c_str());
        dtreeBins[j] = new ddtree::dtree(f,(std::string(t->nt()->GetName())+"_"+std::to_string(j)).c_str());
      }
    }
  }
  std::cout << "done reading\n";
  // if no file exists or the file that exists has different parameters, rebuild mixer tree 
  if(!sameFileExists)
  {
    f = new TFile(_fname.c_str(),"recreate");
    // setup 2D vectors of TTrees and dtrees
    dtreeBins.clear();
    dtreeBins.resize(nvzBins);
    trees.resize(nvzBins);
    cursors.resize(nvzBins);
    ULong64_t current_entry[nvzBins];
    // copy structure for all TTrees, and initialize cursor
    for(int j=0;j<nvzBins;j++)
    {
      trees[j] = t->nt()->CloneTree(0);
      trees[j]->Branch("original_entry",&(current_entry[j]));
      std::string tname = t->nt()->GetName();
      trees[j]->SetName((tname+"_"+std::to_string(j)).c_str());
      cursors[j] = 1;
    }

    // loop over original TTree, fill appropriate TTree
    size_t nentries = t->nt()->GetEntries();
    std::cout << "Building event-mixing TTrees...\n";
    for(ULong64_t evt=0;evt<nentries;evt++)
    {
      if(evt % 10000 == 0) xjjc::progressbar(evt,nentries);
      t->nt()->GetEntry(evt);
      int vz_index = floor((t->bestvtxZ-_vz_min)/_vz_width);
      if(vz_index<0) vz_index = 0;
      if(vz_index>=nvzBins) vz_index = nvzBins-1;
      current_entry[vz_index] = evt;
      trees[vz_index]->Fill();
    }
  
    // make dtrees out of the TTrees
    for(int j=0;j<nvzBins;j++)
    {
      dtreeBins[j] = new ddtree::dtree(f,(std::string(t->nt()->GetName())+"_"+std::to_string(j)).c_str());
      std::uniform_int_distribution<int> dist(1,trees[j]->GetEntries()-1);
      cursors[j] = dist(r);
    }
    TParameter<float> vzwidth("vz_width",_vz_width);
    TParameter<float> vzmin("vz_min",_vz_min);
    TParameter<float> vzmax("vz_max",_vz_max);
    vzwidth.Write();
    vzmin.Write();
    vzmax.Write();
    f->Write();
  }
  std::cout << "\nBuilt event-mixing TTrees.\n";
}

ddtree::dtree* EventMixer::findSimilarEvent(float vz, size_t originalentry)
{
  int nvzBins = ceil((_vz_max-_vz_min)/_vz_width);
  int vz_index = floor((vz-_vz_min)/_vz_width);
  int cursor = cursors[vz_index];
  int nentries = dtreeBins[vz_index]->nt()->GetEntries();
  ULong64_t current_entry;
  dtreeBins[vz_index]->nt()->SetBranchAddress("original_entry",&current_entry);
  dtreeBins[vz_index]->nt()->GetEntry(cursor);
  if(originalentry == current_entry)
  {
    cursors[vz_index]++;
    cursor = cursors[vz_index];
  }
  if(cursor==nentries-1) cursors[vz_index] = 0;
  else cursors[vz_index]++;
  dtreeBins[vz_index]->nt()->GetEntry(cursor);
  return dtreeBins[vz_index];
}
#endif
