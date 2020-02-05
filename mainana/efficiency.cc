#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TString.h"

#include <string>
#include <iostream>
#include <algorithm>

const std::vector<int> centralities = {0, 30, 50, 80};
const unsigned nCentrality = centralities.size() - 1;
const std::string MVA_FILE = "DntupleRun2018/mvas/BDTCuts_finePtBins_PrompD0_PbPb_cent%dto%d_trkPt1GeV.root";
const std::string MC_TREE = "d0ana_mc_genmatchunswap/VertexCompositeNtuple";

const bool true_mc = false;

const float default_pt2min  = 2;
const float default_pt2max  = 999;
const float default_dy      = 1;
const float default_centmin = 0;
const float default_centmax = 80;

const TString plotdir = "plots";

TString centralityString(unsigned iCent) {
  return TString::Format("cent%dto%d", centralities[iCent], centralities[iCent + 1]);
}

/* Save the plot with some extension and dir
 */
void savePlot(TString name, TCanvas* canvas, unsigned iCent, TString extension = "png") {
  canvas->SaveAs(plotdir + "/" + name + "_" + centralityString(iCent) + "." + extension);
}

unsigned centralityID(int centrality) {
  if (centrality == 0) {
    return 0;
  }
  unsigned id = std::lower_bound(centralities.begin(), centralities.end(),
                                 centrality, [=] (float i, float j) {return i < j / 2;}) - centralities.begin();
  return id - 1;
}

/* Get D0 efficiency as a function of pT from MC reco vs gen,
 * and save it in a separate histogram
 * @param inputmc    Name of the input MC root file
 * @param out_name    Name of the output root file
 */
void ddbar_getEfficiency(std::string inputmc, std::string out_name,
                         float ptmin = default_pt2min, float ptmax = default_pt2max,
                         float dy = default_dy, float centmin = default_centmin, float centmax = default_centmax) {
  TFile outfile(out_name.c_str(), "recreate");
  // Get MVA bin width / pT (should be 0.065) from BDT histogram
  std::vector<TH2D*> hist_bdtcut(nCentrality);
  std::vector<TH1D*> hist_gen(nCentrality);
  std::vector<TH1D*> hist_reco(nCentrality);

  for (unsigned i = 0; i < nCentrality; ++i) {
    TFile* mva = TFile::Open(TString::Format(MVA_FILE.c_str(), centralities[i], centralities[i + 1]));
    // Create the histograms in the output file directory
    outfile.cd();
    hist_bdtcut[i] = (TH2D*) mva->Get("hist_bdtcut");
    hist_gen[i] = (TH1D*) hist_bdtcut[i]->ProjectionY("hgen" + centralityString(i));
    hist_gen[i]->Reset();
    hist_gen[i]->SetTitle("gen D candidates " + centralityString(i));
    hist_reco[i] = (TH1D*) hist_gen[i]->Clone("hreco" + centralityString(i));
    hist_reco[i]->SetTitle("reco D candidates " + centralityString(i));
    mva->Close();
  }


  // Get efficiency
  TFile* infmc = TFile::Open(inputmc.c_str());
  TTreeReader reader(MC_TREE.c_str(), infmc);
  TTreeReaderValue<Int_t> centrality(reader, "centrality");
  TTreeReaderArray<Float_t> pT(reader, "pT");
  TTreeReaderArray<Float_t> pT_gen(reader, "pT_gen");
  TTreeReaderArray<Float_t> y(reader, "y");
  TTreeReaderArray<Float_t> y_gen(reader, "y_gen");
  // TTreeReaderArray<Int_t> dID1(reader, "DauID1_gen");
  // TTreeReaderArray<Int_t> dID2(reader, "DauID2_gen");

  // traverse the events
  while(reader.Next()) {
    if(*centrality > centmax*2 || *centrality < centmin*2) {
      continue;
    }
    // select reco candidates
    for (unsigned iCand = 0; iCand < pT.GetSize(); ++iCand) {
      if (pT[iCand] < ptmin || pT[iCand] > ptmax || fabs(y[iCand]) > dy) {
        continue;
      }
      hist_reco[centralityID(*centrality)]->Fill(pT[iCand]);
    }
    // select gen candidates
    for (unsigned iCand = 0; iCand < pT_gen.GetSize(); ++iCand) {
      if (pT_gen[iCand] < ptmin || pT_gen[iCand] > ptmax || fabs(y_gen[iCand]) > dy) {
        continue;
      }
      // make sure the decay products are K and pi
      // if (true_mc) {
      //   if (!(abs(dID1[iCand]) == 211 && abs(dID2[iCand]) == 321) &&
      //       !(abs(dID1[iCand]) == 321 && abs(dID2[iCand]) == 211)) {
      //     continue;
      //   }
      // }
      hist_gen[centralityID(*centrality)]->Fill(pT_gen[iCand]);
    }
  }
  outfile.cd();
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas("c", "", 1200, 600);
  for (unsigned i = 0; i < nCentrality; ++i) {
    hist_gen[i]->Draw();
    savePlot("gen_pt", can, i);
    hist_reco[i]->Draw();
    savePlot("reco_pt", can, i);
    TEfficiency eff(*hist_reco[i], *hist_gen[i]);
    eff.SetName("efficiency_" + centralityString(i));
    eff.SetTitle("reco D / gen D " + centralityString(i));
    eff.Draw();
    savePlot("eff", can, i);
    eff.Write();
  }
  outfile.Write();
}


int main(int argc, char* argv[])
{
  const int default_argc = 3;
  if (argc == default_argc) {
      ddbar_getEfficiency(argv[1], argv[2]);
    } else if (argc == 8) {
      ddbar_getEfficiency(argv[1], argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
  } else {
    std::cout << "not enough or too much arguments" << "\n";
    return 1;
  }
  return 0;
}
