#include <TFile.h>
#include <TH1F.h>
#include <string>
#include <vector>
#include <TString.h>
#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fit2d.hpp"

#ifdef SYSTEM_PBPB
#include "ddbar.h"
#endif

#ifdef SYSTEM_PP
#include "ddbar_pp.h"
#endif

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

void ddbar_fithist(std::string inputname, std::string swapname, std::string outputdir, float pt1min, float pt1max,
  float pt2min, float pt2max, float dy, float centmin, float centmax, std::string lbl="2pi")
{
  TFile* inf = TFile::Open(inputname.c_str());
  TFile *infswap = TFile::Open(swapname.c_str());
  //ddbar::kinematics* kinfo = new ddbar::kinematics(inf);
#ifdef SYSTEM_PBPB
  ddbar::kinematics* kinfo = new ddbar::kinematics(pt1min, pt1max, pt2min, pt2max, dy, centmin, centmax);
#endif
#ifdef SYSTEM_PP
  ddbar::kinematics* kinfo = new ddbar::kinematics(pt1min, pt1max, pt2min, pt2max, dy);
#endif
  //ddbar::binning* binfo = new ddbar::binning(inf);
  ddbar::binning* binfo = new ddbar::binning(lbl);

  // import D0 mass and weight tree
  TTree* tmass = (TTree*)inf->Get("dmass");
  TTree* tmass_sm = (TTree*)inf->Get("dmass_sigmix");
  TTree* tmass_mm = (TTree*)inf->Get("dmass_mixmix");
  TTree* swapmass = (TTree*)infswap->Get("swapdp");

/*
  TChain* tmass = new TChain("dmass");
  TChain* tmass_sm = new TChain("dmass_sigmix");
  TChain* tmass_mm = new TChain("dmass_mixmix");
  TChain* swapmass = new TChain("swapdp");
  std::fstream listfile_dp;
  std::fstream listfile_swap;
  listfile_dp.open("filelist_dptree.txt",std::ios::in);
  listfile_swap.open("filelist_swaptree.txt",std::ios::in);
  std::string f;
  while(std::getline(listfile_dp,f))
  {
    tmass->Add(f.c_str());
    tmass_sm->Add(f.c_str());
    tmass_mm->Add(f.c_str());
  }
  while(std::getline(listfile_swap,f)){
    swapmass->Add(f.c_str());
  }
  tmass->LoadTree(0);
  tmass_sm->LoadTree(0);
  tmass_mm->LoadTree(0);
  swapmass->LoadTree(0);
*/
  TH1F* hdphi_incl = new TH1F("hdphi_incl",";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sigmix = new TH1F("hdphi_sigmix",";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_mixmix = new TH1F("hdphi_mixmix",";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
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
    kinfo->texpt1().c_str(),
    kinfo->texyd().c_str(),
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
      std::cout << k << "\n";
      label.push_back(Form("%s < #Delta#phi/#pi < %s", xjjc::number_remove_zero(binfo->phibin()[k]).c_str(), xjjc::number_remove_zero(binfo->phibin()[k+1]).c_str()));
      auto result = fitter->simfit(tmass, swapmass, k, "PbPb", Form("%s/idx/c%02d", outputname.c_str(), k), label);
      std::cout << "signal fit done\n";
      status.push_back(*result);

      hdphi_incl->SetBinContent(k + 1, fitter->GetY());
      // hdphi_sub->SetBinError(k + 1, fitter->GetYE());
      hdphi_incl->SetBinError(k + 1, fitter->GetYE());
      hdphi_incl->Sumw2();

      auto result_sm = fitter->simfit(tmass_sm, swapmass, k, "PbPb", Form("%s/idx/c%02d",outputname.c_str(),k), label);
      std::cout << "signal-mixed fit done\n";
      status.push_back(*result_sm);

      hdphi_sigmix->SetBinContent(k+1,fitter->GetY());
      hdphi_sigmix->SetBinError(k+1,fitter->GetYE());
      hdphi_sigmix->Sumw2();

      auto result_mm = fitter->simfit(tmass_mm, swapmass, k, "PbPb", Form("%s/idx/c%02d",outputname.c_str(),k),label);
      std::cout << "mixed-mixed fit done\n";
      status.push_back(*result_mm);

      hdphi_mixmix->SetBinContent(k+1,fitter->GetY());
      hdphi_mixmix->SetBinError(k+1,fitter->GetYE());
      hdphi_mixmix->Sumw2();

      hdphi_sub->SetBinContent(k+1,hdphi_incl->GetBinContent(k+1)-hdphi_sigmix->GetBinContent(k+1)+hdphi_mixmix->GetBinContent(k+1));
      hdphi_sub->SetBinError(k+1,sqrt(fabs(hdphi_incl->GetBinError(k+1)*hdphi_incl->GetBinError(k+1)+hdphi_sigmix->GetBinError(k+1)*hdphi_sigmix->GetBinError(k+1)+hdphi_mixmix->GetBinError(k+1)*hdphi_mixmix->GetBinError(k+1))));

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
  if (argc==12) {
    ddbar_fithist(argv[1], argv[2], argv[3], std::stof(argv[4]),
                  std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]),
                  std::stof(argv[8]), std::stof(argv[9]), std::stof(argv[10]),
                  argv[11]);
    return 0;
  } else {
    return 1;
  }
}
