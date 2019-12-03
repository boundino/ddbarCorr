#include <vector>
#include <string>
#include <map>

#include <xjjrootuti.h>

#ifndef MASS_DZERO
#define MASS_DZERO 1.8648
#endif

namespace ddbar
{
  class binning
  {
  public:
    binning(std::string option) : foption(option) { fphibin = bins[foption]; fnphibin = fphibin.size() - 1; }
    binning(TFile* inf);
    int nphibin() { return fnphibin; }
    std::vector<float> phibin() { return fphibin; }
    std::string option() { return foption; }
    float getdphi(float phi1, float phi2, int& idphi);
    void write();
  private:
    int fnphibin;
    std::vector<float> fphibin;
    std::string foption;
    std::map<std::string, std::vector<float>> bins = 
      {
        std::pair<std::string, std::vector<float>>("2pi", std::vector<float>({-1./2, -1./4, 0, 1./4, 1./2, 3./4, 1., 5./4, 3./2})),
        std::pair<std::string, std::vector<float>>("4pi", std::vector<float>({-2., -7./4, -6./4, -5./4, -1, -3./4, -2./4, -1./4, 0, 1./4, 1./2, 3./4, 1., 5./4, 3./2, 7./4, 2.})),
        // std::pair<std::string, std::vector<float>>("pi", std::vector<float>({0, 1./6, 1./3, 1./2, 2./3, 5./6, 1.})),
        std::pair<std::string, std::vector<float>>("pi", std::vector<float>({0, 1./8, 2./8, 3./8, 4./8, 5./8, 6./8, 7./8, 1.})),
      };
  };

  const float signalwidth = 0.05, sideband_l = 0.05, sideband_h = 0.10;

  class kinematics
  {
  public:
    kinematics(float pt1min, float pt1max, float pt2min, float pt2max, float yd, float centmin, float centmax) : 
      fpt1min(pt1min), fpt1max(pt1max), fpt2min(pt2min), fpt2max(pt2max), fyd(yd), fcentmin(centmin), fcentmax(centmax) { ; }
    kinematics(TFile* inf);
    void write();
    float pt1min() { return fpt1min; }
    float pt1max() { return fpt1max; }
    float pt2min() { return fpt2min; }
    float pt2max() { return fpt2max; }
    float yd() { return fyd; }
    float centmin() { return fcentmin; }
    float centmax() { return fcentmax; }
    std::string texpt1();
    std::string texpt2();
    std::string texyd() { return Form("|y^{D}| < %s", xjjc::number_remove_zero(fyd).c_str()); }
    std::string texcent() { return Form("Centrality %.0f-%.0f%s", fcentmin, fcentmax, "%"); }
    void drawlabel(float xx, float yy);
  private:
    float fpt1min, fpt1max, fpt2min, fpt2max, fyd, fcentmin, fcentmax;
    bool infpt1() { return fpt1max>900; }
    bool infpt2() { return fpt2max>900; }
  };

}

float ddbar::binning::getdphi(float phi1, float phi2, int& idphi)
{
  float dphi = phi2 - phi1;
  if(foption == "2pi")
    {
      if(dphi < -M_PI) dphi += M_PI*2;
      else if(dphi > M_PI) dphi -= M_PI*2;
      if(dphi < -M_PI/2.) dphi += M_PI*2;
    }
  if(foption == "pi")
    {
      if(dphi < -M_PI) dphi += M_PI*2;
      else if(dphi > M_PI) dphi -= M_PI*2;
      dphi = fabs(dphi);
    }
  dphi /= M_PI;

  if(dphi >= fphibin[0])
    {
      for(idphi=0; idphi<fnphibin; idphi++)
        {
          if(dphi <= fphibin[idphi+1]) break;
        }
    }
  else idphi = -1;

  return dphi;
}

ddbar::binning::binning(TFile* inf)
{
  std::string *getoption = 0;
  TTree* binfo = (TTree*)inf->Get("binfo");
  binfo->SetBranchAddress("option", &getoption);
  binfo->GetEntry(0);
  foption = *getoption;
  fphibin = bins[foption];
  fnphibin = fphibin.size() - 1;
}

void ddbar::binning::write()
{
  TTree* binfo = new TTree("binfo", "binning");
  binfo->Branch("option", &foption);
  binfo->Fill();
  binfo->Write();
}

ddbar::kinematics::kinematics(TFile* inf)
{
  TTree* kinfo = (TTree*)inf->Get("kinfo");
  kinfo->SetBranchAddress("pt1min", &fpt1min);
  kinfo->SetBranchAddress("pt1max", &fpt1max);
  kinfo->SetBranchAddress("pt2min", &fpt2min);
  kinfo->SetBranchAddress("pt2max", &fpt2max);
  kinfo->SetBranchAddress("centmin", &fcentmin);
  kinfo->SetBranchAddress("centmax", &fcentmax);
  kinfo->SetBranchAddress("yd", &fyd);
  kinfo->GetEntry(0);
}

void ddbar::kinematics::write()
{
  TTree* kinfo = new TTree("kinfo", "kinematics");
  kinfo->Branch("pt1min", &fpt1min);
  kinfo->Branch("pt1max", &fpt1max);
  kinfo->Branch("pt2min", &fpt2min);
  kinfo->Branch("pt2max", &fpt2max);
  kinfo->Branch("centmin", &fcentmin);
  kinfo->Branch("centmax", &fcentmax);
  kinfo->Branch("yd", &fyd);
  kinfo->Fill();
  kinfo->Write();
}

std::string ddbar::kinematics::texpt1()
{
  std::string tt(Form("%s < p_{T}^{trig D} < %s GeV/c", xjjc::number_remove_zero(fpt1min).c_str(), xjjc::number_remove_zero(fpt1max).c_str()));
  if(infpt1()) tt = Form("p_{T}^{trig D} > %s GeV/c", xjjc::number_remove_zero(fpt1min).c_str());
  return tt;
}

std::string ddbar::kinematics::texpt2()
{
  std::string tt(Form("%s < p_{T}^{asso D} < %s GeV/c", xjjc::number_remove_zero(fpt2min).c_str(), xjjc::number_remove_zero(fpt2max).c_str()));
  if(infpt2()) tt = Form("p_{T}^{asso D} > %s GeV/c", xjjc::number_remove_zero(fpt2min).c_str());
  return tt;
}

void ddbar::kinematics::drawlabel(float xx, float yy)
{
  std::vector<std::string> label = {
    texpt1(),
    texpt2(),
    texyd(),
    texcent(),
  };
  float deltay = 0.06; yy += deltay;
  for(auto& tt : label)
    xjjroot::drawtex(xx, yy-=deltay, tt.c_str(), 0.038, 13);
}

