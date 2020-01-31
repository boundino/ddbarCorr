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

const std::string MVA_FILE = "/raid5/data/wangj/DntupleRun2018/mvas/BDTCuts_finePtBins_PrompD0_PbPb_cent0to30_trkPt1GeV.root";
const std::string MC_TREE = "d0ana_mc_genmatchunswap/VertexCompositeNtuple";

const float default_pt2min  = 2;
const float default_pt2max  = 999;
const float default_dy      = 1;
const float default_centmin = 0;
const float default_centmax = 80;

const TString plotdir = "plots";

/* Save the plot with some extension and dir
 */
void savePlot(TString name, TCanvas* canvas, TString extension = "png")
{
  canvas->SaveAs(plotdir + "/" + name + "." + extension);
}

/* Get D0 efficiency as a function of pT from MC reco vs gen,
 * and save it in a separate histogram
 * @param inputmc    Name of the input MC root file
 * @param out_name    Name of the output root file
 */
void ddbar_getEfficiency(std::string inputmc, std::string out_name,
                         float ptmin = default_pt2min, float ptmax = default_pt2max,
                         float dy = default_dy, float centmin = default_centmin, float centmax = default_centmax)
{
  // Get MVA bin width / pT (should be 0.065) from BDT histogram
  TFile* mva = TFile::Open(MVA_FILE.c_str());
  TFile outfile(out_name.c_str(), "recreate");

  TH2D* hist_bdtcut = (TH2D*) mva->Get("hist_bdtcut");
  TH1D* hist_gen = (TH1D*) hist_bdtcut->ProjectionY("hgen");
  hist_gen->Reset();
  hist_gen->SetTitle("gen D candidates");
  TH1D* hist_reco = (TH1D*) hist_gen->Clone("hreco");
  hist_reco->SetTitle("reco D candidates");

  // Get efficiency
  TFile* infmc = TFile::Open(inputmc.c_str());
  TTreeReader reader(MC_TREE.c_str(), infmc);
  TTreeReaderValue<Int_t> centrality(reader, "centrality");
  TTreeReaderArray<Float_t> pT(reader, "pT");
  TTreeReaderArray<Float_t> pT_gen(reader, "pT_gen");
  TTreeReaderArray<Float_t> y(reader, "y");
  TTreeReaderArray<Float_t> y_gen(reader, "y_gen");

  // traverse the events
  while(reader.Next())
    {
      if(*centrality > centmax*2 || *centrality < centmin*2)
        {
          continue;
        }
      // select reco candidates
      for (unsigned iCand = 0; iCand < pT.GetSize(); ++iCand)
        {
          if (pT[iCand] < ptmin || pT[iCand] > ptmax || fabs(y[iCand]) > dy)
            {
              continue;
            }
          hist_reco->Fill(pT[iCand]);
        }
      // select gen candidates
      for (unsigned iCand = 0; iCand < pT_gen.GetSize(); ++iCand)
        {
          if (pT_gen[iCand] < ptmin || pT_gen[iCand] > ptmax || fabs(y_gen[iCand]) > dy)
            {
              continue;
            }
          hist_gen->Fill(pT_gen[iCand]);
        }
    }
  outfile.cd();
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas("c", "", 1200, 600);
  hist_gen->Draw();
  savePlot("gen_pt", can);
  hist_reco->Draw();
  savePlot("reco_pt", can);
  TEfficiency eff(*hist_reco, *hist_gen);
  eff.Draw();
  savePlot("eff", can);
  eff.Write();
  outfile.Write();
}


int main(int argc, char* argv[])
{
  const int default_argc = 3;
  if (argc == default_argc)
    {
      ddbar_getEfficiency(argv[1], argv[2]);
    }
  else
    {
      ddbar_getEfficiency(argv[1], argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
    }
  return 1;
}
