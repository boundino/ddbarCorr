#ifndef EVENTMIXER_H
#define EVENTMIXER_H

#include "TTree.h"
#include "dtree.h"
#include "xjjcuti.h"

#include <iostream>
#include <random>
#include <vector>

class EventMixer
{
public:
  EventMixer(TTree* t, int hiBin_width = 1, float vz_width = 0.5, float vz_min =-30., float vz_max=30.);
  EventMixer(ddtree::dtree* t, int hiBin_width = 1, float vz_width = 0.1, float vz_min=-30.,float vz_max=30.);
  size_t findSimilarEvent(int hiBin, float vz, size_t original_entry);
  ddtree::dtree* mix_tree;
private:
  void em_init(int hiBin_width, float vz_width, float vz_min, float vz_max);
  int _hiBin_width;
  float _vz_width;
  float _vz_min;
  float _vz_max;
  std::default_random_engine r;
  std::vector<std::vector<std::vector<size_t>>> similarEvents;
  std::vector<std::vector<std::uniform_int_distribution<int>>> seDists;
};

EventMixer::EventMixer(TTree* t, int hiBin_width, float vz_width, float vz_min, float vz_max)
{
  mix_tree = new ddtree::dtree(t);
  em_init(hiBin_width,vz_width,vz_min,vz_max);
}

EventMixer::EventMixer(ddtree::dtree* t, int hiBin_width, float vz_width, float vz_min, float vz_max)
{
  mix_tree = t;
  em_init(hiBin_width,vz_width,vz_min,vz_max);
}

void EventMixer::em_init(int hiBin_width, float vz_width, float vz_min, float vz_max)
{
  _hiBin_width = hiBin_width;
  _vz_width = vz_width;
  _vz_min = vz_min;
  _vz_max = vz_max;
  int nHiBins = ceil(200./hiBin_width);
  int nVzBins = ceil((vz_max-vz_min)/vz_width);
  std::cout << "nHiBins: " << nHiBins << "\n";
  std::cout << "nVzBins: " << nVzBins << "\n";
  similarEvents.resize(nHiBins);
  seDists.resize(nHiBins);
  for(int i=0;i<nHiBins;i++)
  {
    similarEvents[i].resize(nVzBins);
    seDists[i].resize(nVzBins);
  }
  size_t nevents = mix_tree->nt()->GetEntries();
  std::cout << "Filling event mixing lookup table...\n";
  for(size_t i=0;i<nevents;i++)
  {
    if(i%10000==0) xjjc::progressbar(i,nevents);
    mix_tree->nt()->GetEntry(i);
    if(mix_tree->bestvtxZ < vz_min || mix_tree->bestvtxZ > vz_max) continue;
    int cbin = mix_tree->centrality / hiBin_width;
    int vzbin = floor((mix_tree->bestvtxZ - vz_min) / vz_width);
    similarEvents[cbin][vzbin].push_back(i);
  }
  std::cout << "done.\n";
  std::cout << "Initializing event mixing random engine...\n";
  for(int c=0;c<nHiBins;c++)
  {
    for(int z=0;z<nVzBins;z++)
    {
      if(similarEvents[c][z].size() != 0)
      {
        seDists[c][z] = std::uniform_int_distribution<int>(0,similarEvents[c][z].size()-1);
      }
    }
  }
  std::cout << "done.\n";
}

size_t EventMixer::findSimilarEvent(int hiBin, float vz, size_t original_entry)
{
  int cbin = hiBin / _hiBin_width;
  int vzbin = floor((vz - _vz_min) / _vz_width);
  if(similarEvents[cbin][vzbin].size() <= 1) return -1;
  size_t evt = original_entry;
  do
  {
    int rand_choice = seDists[cbin][vzbin](r);
    evt = similarEvents[cbin][vzbin][rand_choice];
  } while(evt == original_entry);
  mix_tree->nt()->GetEntry(evt);
  return evt;
}

#endif
