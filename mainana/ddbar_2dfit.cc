#include <TFile.h>
#include <TH1F.h>
#include <string>
#include <vector>
#include <TString.h>
#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "fit2d.hpp"
#include "ddbar.h"

#include <algorithm>

void ddbar_2dfit(std::string inputname, std::string swapname, std::string outputdir)
{
  TFile* inf = TFile::Open(inputname.c_str());
  TFile *infswap = TFile::Open(swapname.c_str());
  // ddbar::kinematics* kinfo = new ddbar::kinematics(inf);
  ddbar::kinematics* kinfo = new ddbar::kinematics(2, 999, 2, 999, 1, 0, 80);
  ddbar::binning* binfo = new ddbar::binning(inf);

  // import D0 mass and weight tree
  TTree* tmass = (TTree*) inf->Get("dmass");
  TTree *swapmass = (TTree*) infswap->Get("swapdp");

  TH1F* hdphi_sub = new TH1F("hdphi_sub", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sub_norm = new TH1F("hdphi_sub_norm", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;1/N dN/d#Delta#phi", binfo->nphibin(), binfo->phibin().data());

  std::vector<float> phi_range;
  for (auto i = 0; i <= binfo->nphibin() / 2; ++i) {
    phi_range.push_back(-1 + i * 0.25);
  }
  std::cout << binfo->nphibin() << "\n";

  std::string outputname = "plots/" + outputdir;
  xjjroot::mkdir(outputname+"/idx/");

  std::vector<TString> label = {
    kinfo->texpt().c_str(),
    kinfo->texyd().c_str(),
    kinfo->texcent().c_str()
  };

  // xjjroot::dfitter* fitter = new xjjroot::dfitter("YL");
  xjjroot::fit2d* fitter = new xjjroot::fit2d("YLR");
  fitter->SetSignalregion(ddbar::signalwidth);
  fitter->SetSidebandL(ddbar::sideband_l);
  fitter->SetSidebandH(ddbar::sideband_h);
  // fitter->SetTexLinespc(0.0);

  // float sidebandscale = fitter->GetSidebandScale();
  fitter->SetOption("Y");
  std::vector<RooFitResult> status;
  for(int k=0; k<binfo->nphibin(); k++)
  // for(int k=0; k<2; k++)
    {
      label.push_back(Form("%s < #Delta#phi/#pi < %s", xjjc::number_remove_zero(binfo->phibin()[k]).c_str(), xjjc::number_remove_zero(binfo->phibin()[k+1]).c_str()));
      auto result = fitter->simfit(tmass, swapmass, k, "PbPb", Form("%s/idx/c%02d", outputname.c_str(), k), label);
      status.push_back(*result);

      hdphi_sub->SetBinContent(k + 1, fitter->GetY());
      // hdphi_sub->SetBinError(k + 1, fitter->GetYE());
      hdphi_sub->SetBinError(k + 1, fitter->GetYE());
      hdphi_sub->Sumw2();

      hdphi_sub_norm->SetBinContent(k+1, hdphi_sub->GetBinContent(k+1) / hdphi_sub_norm->GetBinWidth(k+1));
      hdphi_sub_norm->SetBinError(k+1, hdphi_sub->GetBinError(k+1) / hdphi_sub_norm->GetBinWidth(k+1));

      label.pop_back();
    }

  int idphi = 0;
  for (auto i : status) {
    std::cout << idphi++ << "\n";
    i.Print();
  }
  hdphi_sub_norm->Scale(1./hdphi_sub_norm->Integral("width"));

  hdphi_sub_norm->GetXaxis()->SetNdivisions(-504);
  hdphi_sub_norm->GetYaxis()->SetNdivisions(-505);
  hdphi_sub_norm->SetMinimum(0);
  hdphi_sub_norm->SetMaximum(2.0);
  hdphi_sub->GetXaxis()->SetNdivisions(-504);
  hdphi_sub->SetMinimum(0);
  hdphi_sub->SetMaximum(hdphi_sub->GetMaximum()*1.5);

  xjjroot::sethempty(hdphi_sub, 0, 0.5);
  xjjroot::setthgrstyle(hdphi_sub, kBlack, 21, 1.2, kBlack, 1, 1);
  xjjroot::sethempty(hdphi_sub_norm, 0, 0.5);
  xjjroot::setthgrstyle(hdphi_sub_norm, kBlack, 21, 1.2, kBlack, 1, 1);

  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  hdphi_sub->Draw();
  kinfo->drawlabel(0.60, 0.85);
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2018 PbPb 5.02 TeV");

  c->cd(2);
  hdphi_sub_norm->Draw();
  kinfo->drawlabel(0.60, 0.85);
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2018 PbPb 5.02 TeV");
  c->SaveAs(Form("%s/chdphi.pdf", outputname.c_str()));
  delete c;

  TFile* outf = new TFile(Form("rootfiles/%s/fithist.root", outputdir.c_str()), "recreate");
  kinfo->write();
  hdphi_sub->Write();
  hdphi_sub_norm->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if (argc==4) {
    ddbar_2dfit(argv[1], argv[2], argv[3]);
    return 0;
  } else {
    return 1;
  }
}
