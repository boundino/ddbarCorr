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
#include "TMath.h"
#include "TRandomGen.h"

#include <string>
#include <iostream>
#include <algorithm>

const std::vector<int> centralities = {0, 30, 50, 80};
const std::vector<double> d_centralities = {0.,30.,50.,80.};
const unsigned nCentrality = centralities.size() - 1;
const std::string MVA_FILE = "DntupleRun2018/mvas/BDTCuts_finePtBins_PrompD0_PbPb_cent%dto%d_trkPt1GeV.root";
const std::string MC_TREE = "d0ana_mc_genmatchunswap/VertexCompositeNtuple";

//const Double_t phibinLowEdges[5] = {-TMath::Pi(),-1.,0.,1.,TMath::Pi()};

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
  std::vector<TH2D*> hist_gen(nCentrality);
  std::vector<TH2D*> hist_reco(nCentrality);
  std::vector<TH2D*> eff(nCentrality);
  TH1F* mc_cent = new TH1F("mc_cent","MC Reco Centrality",nCentrality,d_centralities.data());
  TH1F* mc_hiBin = new TH1F("mc_hiBin","MC hiBin",200,0,200);
  for (unsigned i = 0; i < nCentrality; ++i) {
    TFile* mva = TFile::Open(TString::Format(MVA_FILE.c_str(), centralities[i], centralities[i + 1]));
    // Create the histograms in the output file directory
    outfile.cd();
    hist_bdtcut[i] = (TH2D*) mva->Get("hist_bdtcut");
    hist_bdtcut[i]->Rebin2D(1,3);
    // Get pT bins from MVA histogram
    const int nbins = hist_bdtcut[i]->GetNbinsY();
    Double_t ptbins[nbins];
    hist_bdtcut[i]->GetYaxis()->GetLowEdge(ptbins);
    hist_gen[i] = new TH2D("hgen"+centralityString(i),"gen D candidates "+centralityString(i),nbins-1,ptbins,8,-TMath::Pi(),TMath::Pi());
    hist_reco[i] = (TH2D*) hist_gen[i]->Clone("hreco" + centralityString(i));
    hist_reco[i]->SetTitle("reco D candidates " + centralityString(i));
    eff[i] = new TH2D("efficiency_"+centralityString(i),"reco D / gen D "+centralityString(i),nbins-1,ptbins,8,-TMath::Pi(),TMath::Pi());
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
  TTreeReaderArray<Float_t> phi(reader,"phi");
  // randomly fill gen phi bins (assumption that gen phi is flat)
  TRandomMixMax* r = new TRandomMixMax();
//  TTreeReaderArray<Float_t> phi_gen(reader,"phi_gen");
  // TTreeReaderArray<Int_t> dID1(reader, "DauID1_gen");
  // TTreeReaderArray<Int_t> dID2(reader, "DauID2_gen");

  // traverse the events
  while(reader.Next()) {
    if(*centrality > centmax*2 || *centrality < centmin*2) {
      continue;
    }
    mc_cent->Fill((*centrality)/2);
    mc_hiBin->Fill((*centrality));
    // select reco candidates
    for (unsigned iCand = 0; iCand < pT.GetSize(); ++iCand) {
      if (pT[iCand] < ptmin || pT[iCand] > ptmax || fabs(y[iCand]) > dy) {
        continue;
      }
      hist_reco[centralityID(*centrality)]->Fill(pT[iCand],phi[iCand]);
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
      hist_gen[centralityID(*centrality)]->Fill(pT_gen[iCand],r->Uniform(-TMath::Pi(),TMath::Pi()));
    }
    
  }
  outfile.cd();
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas("c", "", 1200, 600);
  for (unsigned i = 0; i < nCentrality; ++i) {
    hist_gen[i]->Sumw2();
    hist_reco[i]->Sumw2();
    hist_gen[i]->Draw();
    savePlot("gen_pt", can, i);
    hist_reco[i]->Draw();
    savePlot("reco_pt", can, i);
    eff[i]->Divide(hist_reco[i],hist_gen[i]);
    eff[i]->Draw();
    savePlot("eff", can, i);
    eff[i]->Write();
  }
  mc_cent->Sumw2();
  mc_cent->Scale(1./mc_cent->Integral(),"width");
  mc_cent->Write();
  mc_hiBin->Write();
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
