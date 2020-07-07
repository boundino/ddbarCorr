#include <TFile.h>
#include <TH1F.h>
#include <string>
#include <vector>
#include <TString.h>
#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "dfitter.h"
#include "ddbar_pp.h"

#include <algorithm>

void ddbar_fithist(std::string inputname, std::string outputdir)
{
  TFile* inf = TFile::Open(inputname.c_str());
  ddbar::kinematics* kinfo = new ddbar::kinematics(inf);
  ddbar::binning* binfo = new ddbar::binning(inf);
  
  TH1F* hmassmc_sgl = (TH1F*)inf->Get("hmassmc_sgl");
  TH1F* hmassmc_swp = (TH1F*)inf->Get("hmassmc_swp");
  std::vector<TH1F*> hmass_incl(binfo->nphibin(), 0), hmass_sdbd(binfo->nphibin(), 0);
  std::vector<TH1F*> hmass_incl_det(binfo->nphibin(), 0);
  for(int k=0; k<binfo->nphibin(); k++)
    {
      hmass_incl[k] = (TH1F*)inf->Get(Form("hmass_dphi%d_incl", k));
      hmass_sdbd[k] = (TH1F*)inf->Get(Form("hmass_dphi%d_sdbd", k));
      hmass_incl_det[k] = (TH1F*)inf->Get(Form("hmass_phi%d_incl_det", k));
    }
  TH1F* hmass_trig = (TH1F*)inf->Get("hmass_trig");

  TH1F* hdphi_incl = new TH1F("hdphi_incl", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sdbd = new TH1F("hdphi_sdbd", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sdbd_noscale = new TH1F("hdphi_sdbd_noscale", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sub = new TH1F("hdphi_sub", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin(), binfo->phibin().data());
  TH1F* hdphi_sub_norm = new TH1F("hdphi_sub_norm", ";#Delta#phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;1/N dN/d#Delta#phi", binfo->nphibin(), binfo->phibin().data());

  std::vector<float> phi_range;
  for (auto i = 0; i <= binfo->nphibin() / 2; ++i) {
    phi_range.push_back(-1 + i * 0.25);
  }
  std::cout << binfo->nphibin() << "\n";
  TH1F* hphi_incl = new TH1F("hphi_incl", ";detector #phi(D^{0}, #bar{D^{#lower[0.2]{0}}}) / #pi;N(D^{0}-#bar{D^{#lower[0.2]{0}}})", binfo->nphibin() / 2, phi_range.data());

  std::string outputname = "plots/" + outputdir;
  xjjroot::mkdir(outputname+"/idx/");

  std::vector<TString> label = {
    kinfo->texpt1().c_str(),
    kinfo->texyd().c_str(),
  };

  xjjroot::dfitter* fitter = new xjjroot::dfitter("YL");
  fitter->SetSignalregion(ddbar::signalwidth);
  fitter->SetSidebandL(ddbar::sideband_l);
  fitter->SetSidebandH(ddbar::sideband_h);
  // fitter->SetTexLinespc(0.0);

  fitter->fit(hmass_trig, hmassmc_sgl, hmassmc_swp, "pp", Form("%s/c%s", outputname.c_str(), hmass_trig->GetName()), label);
  // float sidebandscale = ddbar::signalwidth/(ddbar::sideband_h-ddbar::sideband_l);
  TF1* f_not_mass = fitter->GetFun_not_mass();
  float sidebandscale = f_not_mass->Integral(MASS_DZERO-ddbar::signalwidth, MASS_DZERO+ddbar::signalwidth) / 
    (f_not_mass->Integral(MASS_DZERO-ddbar::sideband_h, MASS_DZERO-ddbar::sideband_l) + f_not_mass->Integral(MASS_DZERO+ddbar::sideband_l, MASS_DZERO+ddbar::sideband_h));
  label.push_back(kinfo->texpt2());
  fitter->SetOption("Y");
  for(int k=0; k<binfo->nphibin(); k++)
    {
      label.push_back(Form("%s < #Delta#phi/#pi < %s", xjjc::number_remove_zero(binfo->phibin()[k]).c_str(), xjjc::number_remove_zero(binfo->phibin()[k+1]).c_str()));
      fitter->fit(hmass_incl[k], hmassmc_sgl, hmassmc_swp, "pp", Form("%s/idx/c%s", outputname.c_str(), hmass_incl[k]->GetName()), label);
      hdphi_incl->SetBinContent(k+1, fitter->GetY());
      hdphi_incl->SetBinError(k+1, fitter->GetYE());
      hdphi_incl->Sumw2();

      if (k < binfo->nphibin() / 2) {
        fitter->fit(hmass_incl_det[k], hmassmc_sgl, hmassmc_swp, "pp", Form("%s/idx/c%s", outputname.c_str(), hmass_incl_det[k]->GetName()), label);
        hphi_incl->SetBinContent(k+1, fitter->GetY());
        hphi_incl->SetBinError(k+1, fitter->GetYE());
        hphi_incl->Sumw2();
      }

      fitter->fit(hmass_sdbd[k], hmassmc_sgl, hmassmc_swp, "pp", Form("%s/idx/c%s", outputname.c_str(), hmass_sdbd[k]->GetName()), label);
      hdphi_sdbd_noscale->SetBinContent(k+1, fitter->GetY());
      hdphi_sdbd_noscale->SetBinError(k+1, fitter->GetYE());
      hdphi_sdbd->SetBinContent(k+1, fitter->GetY()*sidebandscale);
      hdphi_sdbd->SetBinError(k+1, fitter->GetYE()*sidebandscale);
      hdphi_sdbd->Sumw2();

      hdphi_sub->SetBinContent(k+1, hdphi_incl->GetBinContent(k+1)-hdphi_sdbd->GetBinContent(k+1));
      hdphi_sub->SetBinError(k+1, sqrt(fabs(hdphi_incl->GetBinError(k+1)*hdphi_incl->GetBinError(k+1)+hdphi_sdbd->GetBinError(k+1)*hdphi_sdbd->GetBinError(k+1))));

      hdphi_sub_norm->SetBinContent(k+1, hdphi_sub->GetBinContent(k+1) / hdphi_sub_norm->GetBinWidth(k+1));
      hdphi_sub_norm->SetBinError(k+1, hdphi_sub->GetBinError(k+1) / hdphi_sub_norm->GetBinWidth(k+1));

      label.pop_back();
    }
  hdphi_sub_norm->Scale(1./hdphi_sub_norm->Integral("width"));

  hdphi_incl->GetXaxis()->SetNdivisions(-504);
  hdphi_incl->SetMinimum(0);
  hdphi_incl->SetMaximum(hdphi_incl->GetMaximum()*1.2);
  hdphi_sub_norm->GetXaxis()->SetNdivisions(-504);  
  hdphi_sub_norm->GetYaxis()->SetNdivisions(-505);  
  hdphi_sub_norm->SetMinimum(0);
  hdphi_sub_norm->SetMaximum(2.0);
  hdphi_sub->GetXaxis()->SetNdivisions(-504);  
  hdphi_sub->SetMinimum(0);
  hdphi_sub->SetMaximum(hdphi_sub->GetMaximum()*1.5);

  xjjroot::sethempty(hdphi_incl, 0, 0.5);
  xjjroot::setthgrstyle(hdphi_incl, xjjroot::mycolor_satmiddle["azure"], 21, 1.2, xjjroot::mycolor_satmiddle["azure"], 1, 1);
  xjjroot::sethempty(hdphi_sdbd, 0, 0.5);
  xjjroot::setthgrstyle(hdphi_sdbd, xjjroot::mycolor_satmiddle["red"], 21, 1.2, xjjroot::mycolor_satmiddle["red"], 1, 1);
  xjjroot::sethempty(hdphi_sdbd_noscale, 0, 0.5);
  xjjroot::setthgrstyle(hdphi_sdbd_noscale, xjjroot::mycolor_satmiddle["green"], 21, 1.2, xjjroot::mycolor_satmiddle["green"], 1, 1);
  xjjroot::sethempty(hdphi_sub, 0, 0.5);
  xjjroot::setthgrstyle(hdphi_sub, kBlack, 21, 1.2, kBlack, 1, 1);
  xjjroot::sethempty(hdphi_sub_norm, 0, 0.5);
  xjjroot::setthgrstyle(hdphi_sub_norm, kBlack, 21, 1.2, kBlack, 1, 1);

  TLegend* leg1 = new TLegend(0.60, 0.60-0.04*3, 0.88, 0.60);
  xjjroot::setleg(leg1, 0.04);
  leg1->AddEntry(hdphi_incl, "Inclusive", "pl");
  leg1->AddEntry(hdphi_sdbd_noscale, "Sideband", "pl");
  leg1->AddEntry(hdphi_sdbd, "Sideband scaled", "pl");
  TLegend* leg11 = new TLegend(0.60, 0.60-0.04*2, 0.88, 0.60);
  xjjroot::setleg(leg11, 0.04);
  leg11->AddEntry(hdphi_incl, "Inclusive", "pl");
  leg11->AddEntry(hdphi_sdbd, "Background", "pl");
  TLegend* leg2 = new TLegend(0.64, 0.60-0.04*1, 0.90, 0.60);
  xjjroot::setleg(leg2, 0.04);
  leg2->AddEntry(hdphi_sub, "Signal", "pl");

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  hdphi_incl->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg1->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  // c->SaveAs(Form("%s/chdphi_step0.pdf", outputname.c_str()));
  c->cd(1);
  hdphi_sdbd_noscale->Draw("same");
  // c->SaveAs(Form("%s/chdphi_step1.pdf", outputname.c_str()));
  c->cd(1);
  hdphi_sdbd->Draw("same");
  c->cd(2);
  hdphi_sub->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg2->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  // c->SaveAs(Form("%s/chdphi_step2.pdf", outputname.c_str()));

  delete c;
  c = new TCanvas("c", "", 1800, 600);
  c->Divide(3, 1);
  c->cd(1);
  hdphi_incl->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg1->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  // c->SaveAs(Form("%s/chdphi_step0.pdf", outputname.c_str()));
  c->cd(1);
  hdphi_sdbd_noscale->Draw("same");
  // c->SaveAs(Form("%s/chdphi_step1.pdf", outputname.c_str()));
  c->cd(1);
  hdphi_sdbd->Draw("same");
  c->cd(2);
  hdphi_sub->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg2->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  c->cd(3);
  hdphi_sub_norm->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg2->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  c->SaveAs(Form("%s/chdphi_3p.pdf", outputname.c_str()));

  delete c;
  c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  hdphi_incl->Draw();
  hdphi_sdbd->Draw("same");
  kinfo->drawlabel(0.60, 0.85);
  leg11->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  c->cd(2);
  hdphi_sub_norm->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg2->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  c->SaveAs(Form("%s/chdphi.pdf", outputname.c_str()));

  delete c;
  c = new TCanvas("c", "", 600, 600);
  hdphi_incl->Draw();
  hdphi_sdbd->Draw("same");
  kinfo->drawlabel(0.60, 0.85);
  leg11->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  // c->SaveAs(Form("%s/chdphi_panel1.pdf", outputname.c_str()));

  delete c;
  c = new TCanvas("c", "", 600, 600);
  hdphi_sub_norm->Draw();
  kinfo->drawlabel(0.60, 0.85);
  leg2->Draw();
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright("2017 pp 5.02 TeV");
  // c->SaveAs(Form("%s/chdphi_panel2.pdf", outputname.c_str()));

  TFile* outf = new TFile(Form("rootfiles/%s/fithist.root", outputdir.c_str()), "recreate");
  kinfo->write();
  hdphi_incl->Write();
  hphi_incl->Write();
  hdphi_sdbd_noscale->Write();
  hdphi_sdbd->Write();
  hdphi_sub->Write();
  hdphi_sub_norm->Write();
  hmass_trig->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==3) { ddbar_fithist(argv[1], argv[2]); return 0; }
  return 1;
}
