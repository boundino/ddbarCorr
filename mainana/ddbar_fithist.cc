#include <TFile.h>
#include <TH1F.h>
#include <string>
#include <vector>
#include <TString.h>
#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "dfitter.h"
#include "ddbar.h"

void ddbar_fithist(std::string inputname, std::string outputdir)
{
  TFile* inf = TFile::Open(inputname.c_str());
  ddbar::kinematics* kinfo = new ddbar::kinematics(inf);
  ddbar::binning* binfo = new ddbar::binning(inf);
  
  TH1F* hmassmc_sgl = (TH1F*)inf->Get("hmassmc_sgl");
  TH1F* hmassmc_swp = (TH1F*)inf->Get("hmassmc_swp");
  std::vector<TH1F*> hmass_incl(binfo->nphibin(), 0), hmass_sdbd(binfo->nphibin(), 0);
  for(int k=0; k<binfo->nphibin(); k++)
    {
      hmass_incl[k] = (TH1F*)inf->Get(Form("hmass_dphi%d_incl", k));
      hmass_sdbd[k] = (TH1F*)inf->Get(Form("hmass_dphi%d_sdbd", k));
    }
  TH1F* hmass_trig = (TH1F*)inf->Get("hmass_trig");

  TH1F* hdphi_incl = new TH1F("hdphi_incl", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;dN/d#Delta#phi", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sdbd = new TH1F("hdphi_sdbd", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;dN/d#Delta#phi", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sdbd_noscale = new TH1F("hdphi_sdbd_noscale", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;dN/d#Delta#phi", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sub = new TH1F("hdphi_sub", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;dN/d#Delta#phi", binfo->nphibin(), binfo->phibin().data());

  std::string outputname = "plots/" + outputdir;
  xjjroot::mkdir(outputname+"/idx/");

  std::vector<std::string> label = {
    kinfo->texpt1(),
    kinfo->texyd(),
    kinfo->texcent()
  };

  xjjroot::dfitter* fitter = new xjjroot::dfitter("YL");
  fitter->SetSignalregion(ddbar::signalwidth);
  fitter->SetSidebandL(ddbar::sideband_l);
  fitter->SetSidebandH(ddbar::sideband_h);
  // fitter->SetTexLinespc(0.0);

  fitter->fit(hmass_trig, hmassmc_sgl, hmassmc_swp, "PbPb", Form("%s/c%s", outputname.c_str(), hmass_trig->GetName()), label);
  // float sidebandscale = ddbar::signalwidth/(ddbar::sideband_h-ddbar::sideband_l);
  TF1* f_not_mass = fitter->GetFun_not_mass();
  float sidebandscale = f_not_mass->Integral(MASS_DZERO-ddbar::signalwidth, MASS_DZERO+ddbar::signalwidth) / 
    (f_not_mass->Integral(MASS_DZERO-ddbar::sideband_h, MASS_DZERO-ddbar::sideband_l) + f_not_mass->Integral(MASS_DZERO+ddbar::sideband_l, MASS_DZERO+ddbar::sideband_h));
  label.push_back(kinfo->texpt2());
  fitter->SetOption("Y");
  for(int k=0; k<binfo->nphibin(); k++)
    {
      label.push_back(Form("%s < #Delta#phi/#pi < %s", xjjc::number_remove_zero(binfo->phibin()[k]).c_str(), xjjc::number_remove_zero(binfo->phibin()[k+1]).c_str()));
      fitter->fit(hmass_incl[k], hmassmc_sgl, hmassmc_swp, "PbPb", Form("%s/idx/c%s", outputname.c_str(), hmass_incl[k]->GetName()), label);
      hdphi_incl->SetBinContent(k+1, fitter->GetY() / hdphi_incl->GetBinWidth(k+1));
      hdphi_incl->SetBinError(k+1, fitter->GetYE() / hdphi_incl->GetBinWidth(k+1));

      fitter->fit(hmass_sdbd[k], hmassmc_sgl, hmassmc_swp, "PbPb", Form("%s/idx/c%s", outputname.c_str(), hmass_sdbd[k]->GetName()), label);
      hdphi_sdbd_noscale->SetBinContent(k+1, fitter->GetY() / hdphi_sdbd_noscale->GetBinWidth(k+1));
      hdphi_sdbd_noscale->SetBinError(k+1, fitter->GetYE() / hdphi_sdbd_noscale->GetBinWidth(k+1));
      hdphi_sdbd->SetBinContent(k+1, fitter->GetY()*sidebandscale / hdphi_sdbd->GetBinWidth(k+1));
      hdphi_sdbd->SetBinError(k+1, fitter->GetYE()*sidebandscale / hdphi_sdbd->GetBinWidth(k+1));

      hdphi_sub->SetBinContent(k+1, hdphi_incl->GetBinContent(k+1)-hdphi_sdbd->GetBinContent(k+1));
      hdphi_sub->SetBinError(k+1, sqrt(fabs(hdphi_incl->GetBinError(k+1)*hdphi_incl->GetBinError(k+1)+hdphi_sdbd->GetBinError(k+1)*hdphi_sdbd->GetBinError(k+1))));

      label.pop_back();
    }
  // hdphi_sdbd->Scale(ddbar::signalwidth/(ddbar::sideband_h-ddbar::sideband_l));
  hdphi_sub->Scale(1./hdphi_sub->Integral("width"));

  hdphi_incl->GetXaxis()->SetNdivisions(-504);
  hdphi_incl->SetMinimum(0);
  hdphi_incl->SetMaximum(hdphi_incl->GetMaximum()*1.2);
  hdphi_sub->GetXaxis()->SetNdivisions(-504);  
  // if(hdphi_sub->GetMinimum() > 0) hdphi_sub->SetMinimum(0);
  // hdphi_sub->SetMaximum(hdphi_sub->GetMaximum()*1.4);
  hdphi_sub->SetMinimum(0);
  hdphi_sub->SetMaximum(2.0);

  xjjroot::sethempty(hdphi_incl, 0, 0.7);
  xjjroot::setthgrstyle(hdphi_incl, xjjroot::mycolor_satmiddle["azure"], 21, 1.2, xjjroot::mycolor_satmiddle["azure"], 1, 1);
  xjjroot::sethempty(hdphi_sdbd, 0, 0.7);
  xjjroot::setthgrstyle(hdphi_sdbd, xjjroot::mycolor_satmiddle["red"], 21, 1.2, xjjroot::mycolor_satmiddle["red"], 1, 1);
  xjjroot::sethempty(hdphi_sdbd_noscale, 0, 0.7);
  xjjroot::setthgrstyle(hdphi_sdbd_noscale, xjjroot::mycolor_satmiddle["green"], 21, 1.2, xjjroot::mycolor_satmiddle["green"], 1, 1);
  xjjroot::sethempty(hdphi_sub, 0, 0.7);
  xjjroot::setthgrstyle(hdphi_sub, kBlack, 21, 1.2, kBlack, 1, 1);

  TLegend* leg1 = new TLegend(0.60, 0.60-0.04*3, 0.88, 0.60);
  xjjroot::setleg(leg1, 0.04);
  leg1->AddEntry(hdphi_incl, "Inclusive", "pl");
  leg1->AddEntry(hdphi_sdbd_noscale, "Sideband", "pl");
  leg1->AddEntry(hdphi_sdbd, "Sideband scaled", "pl");
  TLegend* leg2 = new TLegend(0.64, 0.60-0.04*1, 0.90, 0.60);
  xjjroot::setleg(leg2, 0.04);
  leg2->AddEntry(hdphi_sub, "Signal", "pl");

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  hdphi_incl->Draw();
  hdphi_sdbd_noscale->Draw("same");
  hdphi_sdbd->Draw("same");
  kinfo->drawlabel(0.60, 0.85);
  leg1->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2018 PbPb 5.02 TeV");
  c->cd(2);
  hdphi_sub->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg2->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2018 PbPb 5.02 TeV");
  c->SaveAs(Form("%s/chdphi.pdf", outputname.c_str()));

  TFile* outf = new TFile(Form("rootfiles/%s/fithist.root", outputdir.c_str()), "recreate");
  kinfo->write();
  hdphi_incl->Write();
  hdphi_sdbd_noscale->Write();
  hdphi_sdbd->Write();
  hdphi_sub->Write();
  hmass_trig->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==3) { ddbar_fithist(argv[1], argv[2]); return 0; }
  return 1;
}
