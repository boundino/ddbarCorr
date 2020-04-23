#ifndef _XJJROOT_DFITTER_H_
#define _XJJROOT_DFITTER_H_

// Created by Jing Wang (7/26/2017)

/******************************************************************************************
 * Class : xjjroot::dfitter                                                               * 
 * This class provides tools fitting Dzero invariant mass spectra.                        * 
 * The object to be used can be declared by                                               * 
 *                                                                                        * 
 *    xjjroot::dfitter obj(options);                                                      * 
 *                                                                                        * 
 * Options supported are listed below                                                     * 
 *                                                                                        * 
 *   "3" : Using 3-Gaussian function to model signal (default is 2-Gaussian function)     * 
 *   "Y" : Draw yield info                                                                * 
 *   "C" : Draw chi2 info                                                                 * 
 *   "L" : Draw lines at signal region and side bands                                     * 
 *   "S" : Draw significance info and lines at signal region                              * 
 *   "D" : Draw all details                                                               * 
 *   "V" : Switch off Quiet mode of fitting                                               * 
 *   "X" : Do not save plots                                                              * 
 *                                                                                        * 
 * The core function of this class is                                                     * 
 *                                                                                        * 
 *    TF1* fit(TH1*, TH1*, TH1*, TString, TString, std::vector<TString>)                  * 
 *                                                                                        * 
 ******************************************************************************************/

#include <TF1.h>
#include <TH1.h>
#include <TString.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TFitResult.h>

namespace xjjroot
{
  class dfitter
  {
  public:
    dfitter(Option_t* option="") {foption = option; resolveoption(); init();}
    ~dfitter() {};

    TF1* fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped, TString collisionsyst="", TString outputname="cmass", const std::vector<TString>& vtex=std::vector<TString>());
    Bool_t isFitted() const {return fparamfuns;}

    Double_t GetS() const {return S;}
    Double_t GetB() const {return B;}
    Double_t GetSig() const {return Sig;}
    Double_t GetY() const {return yield;}
    Double_t GetYE() const {return yieldErr;}
    Double_t GetChi2() const {return Chi2;}
    Double_t GetNDF() const {return NDF;}
    Double_t GetChi2Prob() const {return Chi2Prob;}

    TF1* GetFun_f(TString name="Fun_f") const {if(!fparamfuns){return 0;} return clonefun(fun_f, name);}
    TF1* GetFun_mass(TString name="Fun_mass") const {if(!fparamfuns){return 0;} return clonefun(fun_mass, name);}
    TF1* GetFun_swap(TString name="Fun_swap") const {if(!fparamfuns){return 0;} return clonefun(fun_swap, name);}
    TF1* GetFun_background(TString name="Fun_background") const {if(!fparamfuns){return 0;} return clonefun(fun_background, name);}
    TF1* GetFun_not_mass(TString name="Fun_not_mass") const {if(!fparamfuns){return 0;} return clonefun(fun_not_mass, name);}

    void SetOption(Option_t* option="") {foption = option; resolveoption();}
    void SetSignalregion(Double_t d_mass_signal_) {d_mass_signal =  d_mass_signal_; calvar();}
    void SetSidebandL(Double_t d_mass_sideband_l_) {d_mass_sideband_l = d_mass_sideband_l_; calvar();}
    void SetSidebandH(Double_t d_mass_sideband_h_) {d_mass_sideband_h = d_mass_sideband_h_; calvar();}
    void SetTexLinespc(Double_t linespc=0) {texlinespc = linespc;}
    
    Double_t GetMassL() const {return min_hist_dzero;}
    Double_t GetMassH() const {return max_hist_dzero;}

  private:
    Double_t S;
    Double_t B;
    Double_t Sig;
    Double_t yield;
    Double_t yieldErr;
    Double_t Chi2;
    Double_t NDF;
    Double_t Chi2Prob;

    TF1* fun_f;
    TF1* fun_mass;
    TF1* fun_swap;
    TF1* fun_background;
    TF1* fun_not_mass;

    TFitResultPtr r;

    TString foption;
    Bool_t fparamfun_f;
    Bool_t fparamfuns;
    Bool_t f3gaus;
    Bool_t fdrawyield;
    Bool_t fdrawchi2;
    Bool_t fdrawsignif;
    Bool_t fdrawsigsid;
    Bool_t ffitverbose;
    Bool_t fdrawdetail;
    Bool_t fsaveplot;

    const Double_t setparam0 = 100.;
    const Double_t setparam1 = 1.865;
    const Double_t setparam2 = 0.03;
    const Double_t setparam10 = 0.005;
    const Double_t setparam13 = 0.002;
    const Double_t setparam8 = 0.1;
    const Double_t setparam9 = 0.1;
    const Double_t setparam12 = 0.5;
    const Double_t fixparam1 = 1.865;
    Double_t fixparam7;
    
    const Double_t mass_dzero = 1.8649;
    Double_t d_mass_signal = 0.045;
    Double_t d_mass_sideband_l = 0.07;
    Double_t d_mass_sideband_h = 0.12;

    Double_t mass_dzero_signal_l = mass_dzero - d_mass_signal;
    Double_t mass_dzero_signal_h = mass_dzero + d_mass_signal;
    Double_t mass_dzero_sideband_l_p = d_mass_signal + d_mass_sideband_l;
    Double_t mass_dzero_sideband_h_p = d_mass_signal + d_mass_sideband_h;
    Double_t mass_dzero_sideband_l_n = d_mass_signal - d_mass_sideband_l;
    Double_t mass_dzero_sideband_h_n = d_mass_signal - d_mass_sideband_h;
    
    const Double_t n_hist_dzero = 56;
    const Double_t min_hist_dzero = 1.72;
    const Double_t max_hist_dzero = 2.0;
    const Double_t binwid_hist_dzero = (max_hist_dzero-min_hist_dzero)/n_hist_dzero;

    Double_t texlinespc = 0;

    void init();
    void reset();
    void createfun();
    void deletefun();
    void clearvar();
    void calvar();
    void resolveoption();
    void setfunparameters();
    void setfunstyle();

    TF1* clonefun(const TF1* fun, TString fun_name) const;
    void sethist(TH1* h) const;
    void drawCMS(TString collision, TString snn="5.02") const;
    void drawtex(Double_t x, Double_t y, const char* text, Float_t tsize=0.04, Short_t align=12) const;
    void drawleg(TH1* h) const;
    void drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor=kBlack, Style_t lstyle=1, Width_t lwidth=2) const;
    void setgstyle() const;
  };
}

void xjjroot::dfitter::resolveoption()
{
  f3gaus = false;
  if(foption.Contains("3")) f3gaus = true;
  fdrawyield = false;
  if(foption.Contains("Y")) fdrawyield = true;
  fdrawchi2 = false;
  if(foption.Contains("C")) fdrawchi2 = true;
  fdrawsigsid = false;
  if(foption.Contains("L")) fdrawsigsid = true;
  fdrawsignif = false;
  if(foption.Contains("S")) fdrawsignif = true;
  fdrawdetail = false;
  if(foption.Contains("D")) {fdrawdetail = true; fdrawyield = true; fdrawchi2 = true; fdrawsigsid = true; fdrawsignif = true;}
  ffitverbose = false;
  if(foption.Contains("V")) ffitverbose = true;
  fsaveplot = true;
  if(foption.Contains("X")) fsaveplot = false;
}

TF1* xjjroot::dfitter::fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped, TString collisionsyst/*=""*/, TString outputname/*="cmass"*/, const std::vector<TString> &vtex/*=std::vector<TString>()*/)
{
  reset();
  init();

  TH1* h = (TH1*)hmass->Clone("h");
  TH1* hMCSignal = (TH1*)hmassMCSignal->Clone("hMCSignal");
  TH1* hMCSwapped = (TH1*)hmassMCSwapped->Clone("hMCSwapped");
  sethist(h);
  sethist(hMCSignal);
  sethist(hMCSwapped);

  TString fitoption = ffitverbose?"L m":"L m q";
  setgstyle();
  TCanvas* c = new TCanvas("c", "" , 600, 600);
  
  fun_f->SetParLimits(4,  -1000, 1000);
  fun_f->SetParLimits(10, 0.001, 0.05);
  fun_f->SetParLimits(2,  0.01,  0.5);
  fun_f->SetParLimits(8,  0.02,  0.2);
  if(f3gaus) fun_f->SetParLimits(13, 0.002, 0.1);
  fun_f->SetParLimits(7,  0,     1);
  fun_f->SetParLimits(9,  0,     1);
  if(f3gaus) fun_f->SetParLimits(12, 0, 1);
    
  // -- fit MC
  fun_f->FixParameter(3, 0);
  fun_f->FixParameter(4, 0);
  fun_f->FixParameter(5, 0);
  fun_f->FixParameter(6, 0);

  //  - fit signal 
  fun_f->FixParameter(11, 0);
  fun_f->FixParameter(1,  fixparam1);
  fun_f->FixParameter(7,  1);
  fun_f->FixParameter(8,  1);
  fun_f->SetParameter(0,  setparam0);
  fun_f->SetParameter(1,  setparam1);
  fun_f->SetParameter(2,  setparam2);
  fun_f->SetParameter(10, setparam10);
  if(f3gaus) fun_f->SetParameter(13, setparam13);
  fun_f->SetParameter(9,  setparam9);
  if(f3gaus) fun_f->SetParameter(12, setparam12);

  hMCSignal->Fit("fun_f", "q", "", min_hist_dzero, max_hist_dzero);
  hMCSignal->Fit("fun_f", "q", "", min_hist_dzero, max_hist_dzero);
  fun_f->ReleaseParameter(1);
  hMCSignal->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  hMCSignal->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  hMCSignal->Fit("fun_f", fitoption, "", min_hist_dzero, max_hist_dzero);

  fun_f->FixParameter(1, fun_f->GetParameter(1));
  fun_f->FixParameter(2, fun_f->GetParameter(2));
  fun_f->FixParameter(10, fun_f->GetParameter(10));
  if(f3gaus) fun_f->FixParameter(13, fun_f->GetParameter(13));
  fun_f->FixParameter(9, fun_f->GetParameter(9));
  if(f3gaus) fun_f->FixParameter(12, fun_f->GetParameter(12));

  fixparam7 = fun_f->GetParameter(0);

  //   - fit swapped
  fun_f->FixParameter(7,0);
  fun_f->ReleaseParameter(8);
  fun_f->SetParameter(8,setparam8);
  
  hMCSwapped->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  hMCSwapped->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  hMCSwapped->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  hMCSwapped->Fit("fun_f",fitoption,"", min_hist_dzero, max_hist_dzero);
  
  fixparam7 = fixparam7/(fun_f->GetParameter(0)+fixparam7);
  fun_f->FixParameter(7, fixparam7);
  fun_f->FixParameter(8,fun_f->GetParameter(8));

  //  -- fit data
  fun_f->ReleaseParameter(3);
  fun_f->ReleaseParameter(4);
  fun_f->ReleaseParameter(5);
  fun_f->ReleaseParameter(6);
  
  h->Fit("fun_f", "q", "", min_hist_dzero, max_hist_dzero);
  h->Fit("fun_f", "q", "", min_hist_dzero, max_hist_dzero);
  fun_f->ReleaseParameter(1);
  fun_f->SetParLimits(1, 1.86, 1.87);
  fun_f->ReleaseParameter(11);
  fun_f->SetParLimits(11, -0.5, 0.5);
  h->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  h->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  h->Fit("fun_f", "L q", "", min_hist_dzero, max_hist_dzero);
  r = h->Fit("fun_f", Form("%s S",fitoption.Data()),"", min_hist_dzero, max_hist_dzero);

  fparamfun_f = true;
  setfunparameters();
  calvar();
  
  h->Draw("e");
  fun_background->Draw("same");   
  fun_mass->Draw("same");
  fun_swap->Draw("same");
  if(fdrawsigsid || fdrawsignif)
    {
      fun_not_mass->SetRange(mass_dzero_signal_l,mass_dzero_signal_h);
      fun_not_mass->Draw("same");
      drawline(mass_dzero_signal_l, 0, mass_dzero_signal_l, fun_f->Eval(mass_dzero_signal_l), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_signal_h, 0, mass_dzero_signal_h, fun_f->Eval(mass_dzero_signal_h), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    }
  if(fdrawsigsid)
    {
      drawline(mass_dzero_sideband_l_p, 0, mass_dzero_sideband_l_p, fun_f->Eval(mass_dzero_sideband_l_p), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_sideband_h_p, 0, mass_dzero_sideband_h_p, fun_f->Eval(mass_dzero_sideband_h_p), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_sideband_l_n, 0, mass_dzero_sideband_l_n, fun_f->Eval(mass_dzero_sideband_l_n), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_sideband_h_n, 0, mass_dzero_sideband_h_n, fun_f->Eval(mass_dzero_sideband_h_n), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    }
  fun_f->Draw("same");

  drawleg(h);
  drawCMS(collisionsyst);

  Float_t texxpos = 0.22, texypos = 0.90, texdypos = 0.053;
  if(!vtex.empty())
    {
      texypos+=texlinespc;
      for(std::vector<TString>::const_iterator it=vtex.begin(); it!=vtex.end(); it++) 
        drawtex(texxpos, texypos=(texypos-texdypos-texlinespc), *it);
    }
  if(fdrawyield) drawtex(texxpos, texypos=(texypos-texdypos), Form("N = %.0f #pm %.0f",yield,yieldErr));
  if(fdrawchi2) 
    {
      drawtex(texxpos, texypos=(texypos-texdypos), Form("#chi^{2} / ndf = %.1f / %.0f",Chi2,NDF));
      drawtex(texxpos, texypos=(texypos-texdypos), Form("Prob = %.2f%s",Chi2Prob*100,"%"));
    }
  if(fdrawsignif)
    {
      drawtex(texxpos, texypos=(texypos-texdypos), Form("S = %.0f, B = %.0f",S,B));
      drawtex(texxpos, (texypos=(texypos-texdypos)) + 0.01, Form("S/#sqrt{S+B} = %.1f",S/TMath::Sqrt(S+B)));
    }
  if(fdrawdetail)
    {
      drawtex(texxpos, texypos=(texypos-texdypos), Form("N#scale[0.6]{#lower[0.7]{sig}}/(N#scale[0.6]{#lower[0.7]{sig}}+N#scale[0.6]{#lower[0.7]{swap}}) = %.2f",fun_f->GetParameter(7)));
    }
  
  if(fsaveplot) c->SaveAs(Form("%s.pdf",outputname.Data()));

  delete c;
  delete h;
  delete hMCSignal;
  delete hMCSwapped;

  return clonefun(fun_mass, "Fun_mass");
}

void xjjroot::dfitter::reset()
{
  clearvar();
  deletefun();
}

void xjjroot::dfitter::init()
{
  clearvar();
  createfun();
  setfunstyle();
}

void xjjroot::dfitter::clearvar()
{
  S = 0;
  B = 0;
  Sig = 0;
  yield = 0;
  yieldErr = 0;
  Chi2 = 0;
  NDF = 0;
  Chi2Prob = 0;
}

void xjjroot::dfitter::calvar()
{
  mass_dzero_signal_l = mass_dzero - d_mass_signal;
  mass_dzero_signal_h = mass_dzero + d_mass_signal;
  mass_dzero_sideband_l_p = mass_dzero + d_mass_sideband_l;
  mass_dzero_sideband_h_p = mass_dzero + d_mass_sideband_h;
  mass_dzero_sideband_l_n = mass_dzero - d_mass_sideband_l;
  mass_dzero_sideband_h_n = mass_dzero - d_mass_sideband_h;

  if(!fparamfuns) return;
  S = fun_mass->Integral(mass_dzero_signal_l,mass_dzero_signal_h)/binwid_hist_dzero;
  B = fun_background->Integral(mass_dzero_signal_l,mass_dzero_signal_h)/binwid_hist_dzero + fun_swap->Integral(mass_dzero_signal_h,mass_dzero_signal_h)/binwid_hist_dzero;
  Sig = S/TMath::Sqrt(S+B);
  yield = fun_mass->Integral(min_hist_dzero,max_hist_dzero)/binwid_hist_dzero;
  yieldErr = fun_mass->Integral(min_hist_dzero,max_hist_dzero)/binwid_hist_dzero*fun_mass->GetParError(0)/fun_mass->GetParameter(0);
  Chi2 = 2.*r->MinFcnValue();
  NDF = fun_f->GetNDF();
  Chi2Prob = TMath::Prob(Chi2,NDF);
}

void xjjroot::dfitter::createfun()
{
  TString str_fun_f = f3gaus?
    "[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*([12]*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11]))+(1-[12])*TMath::Gaus(x,[1],[13]*(1+[11]))/(sqrt(2*3.14159)*[13]*(1+[11]))))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x":
    "[0]*([7]*([9]*TMath::Gaus(x,[1],[2]*(1+[11]))/(sqrt(2*3.14159)*[2]*(1+[11]))+(1-[9])*TMath::Gaus(x,[1],[10]*(1+[11]))/(sqrt(2*3.14159)*[10]*(1+[11])))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[11]))/(sqrt(2*3.14159)*[8]*(1+[11])))+[3]+[4]*x+[5]*x*x+[6]*x*x*x";
  TString str_fun_mass = f3gaus?
    "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*([7]*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))+(1-[7])*TMath::Gaus(x,[1],[8]*(1+[6]))/(sqrt(2*3.14159)*[8]*(1+[6])))))":
    "[0]*([3]*([4]*TMath::Gaus(x,[1],[2]*(1+[6]))/(sqrt(2*3.14159)*[2]*(1+[6]))+(1-[4])*TMath::Gaus(x,[1],[5]*(1+[6]))/(sqrt(2*3.14159)*[5]*(1+[6]))))";
  fun_f = new TF1("fun_f", str_fun_f.Data(), min_hist_dzero, max_hist_dzero);  
  fun_mass = new TF1("fun_mass", str_fun_mass.Data(), min_hist_dzero, max_hist_dzero);
  fun_background = new TF1("fun_background", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", min_hist_dzero, max_hist_dzero);
  fun_swap = new TF1("fun_swap", "[0]*(1-[2])*TMath::Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))", min_hist_dzero, max_hist_dzero);
  fun_not_mass = new TF1("fun_not_mass", "[0]*(1-[2])*TMath::Gaus(x,[1],[3]*(1+[4]))/(sqrt(2*3.14159)*[3]*(1+[4]))+[5]+[6]*x+[7]*x*x+[8]*x*x*x", min_hist_dzero, max_hist_dzero);

  fparamfun_f = false;
  fparamfuns = false;
}

void xjjroot::dfitter::deletefun()
{
  delete fun_f;
  delete fun_background;
  delete fun_mass;
  delete fun_swap;
  delete fun_not_mass;

  fparamfun_f = false;
  fparamfuns = false;
}

void xjjroot::dfitter::setfunstyle()
{
  fun_f->SetNpx(2000);
  fun_f->SetLineColor(2);
  fun_f->SetLineStyle(1);
  fun_f->SetLineWidth(3);

  fun_background->SetNpx(2000);
  fun_background->SetLineColor(4);
  fun_background->SetLineStyle(2);
  fun_background->SetLineWidth(3);

  fun_mass->SetNpx(2000);
  fun_mass->SetLineColor(kOrange-3);
  fun_mass->SetLineStyle(2);
  fun_mass->SetLineWidth(3);
  fun_mass->SetFillColor(kOrange-3);
  fun_mass->SetFillStyle(3002);

  fun_swap->SetNpx(2000);
  fun_swap->SetLineColor(kGreen+4);
  fun_swap->SetLineStyle(1);
  fun_swap->SetLineWidth(3);
  fun_swap->SetFillColor(kGreen+4);
  fun_swap->SetFillStyle(3005);

  fun_not_mass->SetLineColor(12);
  fun_not_mass->SetLineStyle(2);
  fun_not_mass->SetLineWidth(3);
}

void xjjroot::dfitter::setfunparameters()
{
  fun_background->SetParameters(fun_f->GetParameter(3), fun_f->GetParameter(4), fun_f->GetParameter(5), fun_f->GetParameter(6));
  fun_background->SetParError(0,fun_f->GetParError(3));
  fun_background->SetParError(1,fun_f->GetParError(4));
  fun_background->SetParError(2,fun_f->GetParError(5));
  fun_background->SetParError(3,fun_f->GetParError(6));
  
  fun_mass->SetParameters(fun_f->GetParameter(0),fun_f->GetParameter(1),fun_f->GetParameter(2),fun_f->GetParameter(7),fun_f->GetParameter(9),fun_f->GetParameter(10),fun_f->GetParameter(11));
  fun_mass->SetParError(0,fun_f->GetParError(0));
  fun_mass->SetParError(1,fun_f->GetParError(1));
  fun_mass->SetParError(2,fun_f->GetParError(2));
  fun_mass->SetParError(3,fun_f->GetParError(7));
  fun_mass->SetParError(4,fun_f->GetParError(9));
  fun_mass->SetParError(5,fun_f->GetParError(10));
  fun_mass->SetParError(6,fun_f->GetParError(11));
  if(f3gaus)
    {
      fun_mass->SetParameter(7, fun_f->GetParameter(12));
      fun_mass->SetParameter(8, fun_f->GetParameter(13));
      fun_mass->SetParError(7, fun_f->GetParError(12));
      fun_mass->SetParError(8, fun_f->GetParError(13));
    }   
 
  fun_swap->SetParameters(fun_f->GetParameter(0),fun_f->GetParameter(1),fun_f->GetParameter(7),fun_f->GetParameter(8),fun_f->GetParameter(11));
  fun_swap->SetParError(0,fun_f->GetParError(0));
  fun_swap->SetParError(1,fun_f->GetParError(1));
  fun_swap->SetParError(2,fun_f->GetParError(7));
  fun_swap->SetParError(3,fun_f->GetParError(8));
  fun_swap->SetParError(4,fun_f->GetParError(11));

  fun_not_mass->SetParameters(fun_swap->GetParameter(0),fun_swap->GetParameter(1),fun_swap->GetParameter(2),fun_swap->GetParameter(3),fun_swap->GetParameter(4),
                              fun_background->GetParameter(0),fun_background->GetParameter(1),fun_background->GetParameter(2),fun_background->GetParameter(3));
  fun_not_mass->SetParError(0,fun_swap->GetParError(0));
  fun_not_mass->SetParError(1,fun_swap->GetParError(1));
  fun_not_mass->SetParError(2,fun_swap->GetParError(2));
  fun_not_mass->SetParError(3,fun_swap->GetParError(3));
  fun_not_mass->SetParError(4,fun_swap->GetParError(4));
  fun_not_mass->SetParError(5,fun_background->GetParError(0));
  fun_not_mass->SetParError(6,fun_background->GetParError(1));
  fun_not_mass->SetParError(7,fun_background->GetParError(2));
  fun_not_mass->SetParError(8,fun_background->GetParError(3));

  fparamfuns = true;
}

TF1* xjjroot::dfitter::clonefun(const TF1* fun, TString fun_name) const
{
  TF1* newfun = new TF1(*fun);
  newfun->SetName(fun_name);
  return newfun;
}

void xjjroot::dfitter::sethist(TH1* h) const
{
  h->SetMaximum(-1111);
  h->SetXTitle("m_{#piK} (GeV/c^{2})");
  h->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", binwid_hist_dzero*1.e+3));
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.8);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  h->SetStats(0);
}

void xjjroot::dfitter::drawleg(TH1* h) const
{
  TLegend* leg = new TLegend(0.65, 0.58, 0.85, 0.88, NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  leg->AddEntry(h,"Data", "pl");
  leg->AddEntry(fun_f,"Fit", "l");
  leg->AddEntry(fun_mass,"D^{0}+#bar{D^{#lower[0.2]{0}}} Signal", "f");
  leg->AddEntry(fun_swap,"K-#pi swapped", "f");
  leg->AddEntry(fun_background,"Combinatorial", "l");

  leg->Draw("same");
}

void xjjroot::dfitter::drawCMS(TString collision, TString snn/*="5.02"*/) const
{
  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  TLatex* texCol = new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = %s TeV", collision.Data(), snn.Data()));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);

  texCms->Draw();
  texCol->Draw();
}

void xjjroot::dfitter::drawtex(Double_t x, Double_t y, const char* text, Float_t tsize/*=0.04*/, Short_t align/*=12*/) const
{
  TLatex* tex = new TLatex(x, y, text);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
  tex->Draw();
}

void xjjroot::dfitter::drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=2*/) const
{
  TLine* l = new TLine(x1, y1, x2, y2);
  l->SetLineColor(lcolor);
  l->SetLineStyle(lstyle);
  l->SetLineWidth(lwidth);
  l->Draw();
}

void xjjroot::dfitter::setgstyle() const
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
}

#endif
