#ifndef _XJJROOT_DROOFIT_H_
#define _XJJROOT_DROOFIT_H_

// Adapted from dfitter.h

/******************************************************************************************
 * Class : xjjroot::droofit
 * This class provides tools fitting Dzero invariant mass spectra.                        *
 * The object to be used can be declared by                                               *
 *                                                                                        *
 *    xjjroot::droofit obj(options);                                                      *
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
 *   "R" : Use roofit for fitting                                                         *
 *                                                                                        *
 * The core function of this class is                                                     *
 *                                                                                        *
 *    TF1* fit(TTree*, int, int, TString, TString, std::vector<TString>)                  *
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


// RooFit libs
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooDataHist.h"

using namespace RooFit;

namespace xjjroot
{
  class droofit
  {
  public:
    droofit(Option_t* option="") {foption = option; resolveoption(); init();}
    ~droofit() {};

    TF1* fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped, TString collisionsyst="", TString outputname="cmass", const std::vector<TString>& vtex=std::vector<TString>());
    Bool_t isFitted() const {return fparamfuns;}

    void fit(TTree* tree, const int isSig, const int iBin,
             TString collisionsyst = "",
             TString outputname = "cmass",
             const std::vector<TString> &vtex = std::vector<TString>());

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

    Double_t GetSidebandScale() const {return sidebandScale;}

  private:
    Double_t S;
    Double_t B;
    Double_t Sig;
    Double_t yield;
    Double_t yieldErr;
    Double_t Chi2;
    Double_t NDF;
    Double_t Chi2Prob;
    Double_t sidebandScale;

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
    Bool_t froofit;

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

    const Double_t n_hist_dzero = 60;
    const Double_t min_hist_dzero = 1.7;
    const Double_t max_hist_dzero = 2.0;
    const Double_t binwid_hist_dzero = (max_hist_dzero-min_hist_dzero)/n_hist_dzero;

    const Double_t min_weight = 0;
    const Double_t max_weight = 10000;

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

    float Nintegral(RooAddPdf& function, RooRealVar& integrationVar, double lowerLimit, double upperLimit) {
      integrationVar.setRange("integralRange", lowerLimit, upperLimit);
      RooAbsReal* integral = function.createIntegral(integrationVar, NormSet(integrationVar), Range("integralRange")) ;
      double normalizedIntegralValue = integral->getVal();
      return normalizedIntegralValue;
    }
  };
}

void xjjroot::droofit::resolveoption()
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
  froofit = false;
  if(foption.Contains("R")) froofit = true;
}

/*
  Do chi2 fit
*/
TF1* xjjroot::droofit::fit(const TH1* hmass, const TH1* hmassMCSignal, const TH1* hmassMCSwapped, TString collisionsyst/*=""*/, TString outputname/*="cmass"*/, const std::vector<TString> &vtex/*=std::vector<TString>()*/)
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

  // C r e a t e   o b s e r v a b l e   a n d   u n w e i g h t e d   d a t a s e t
  // -------------------------------------------------------------------------------

  // Declare observable
  RooRealVar x("x", "m_{#pi K} / GeV", min_hist_dzero, max_hist_dzero);
  x.setBins(n_hist_dzero);

  // double gaussian as signal p.d.f.
  // RooRealVar mean("mean", "mean of gaussians", 1.865, 1.85, 1.88);
  RooRealVar mean("mean", "mean of gaussians", mass_dzero, mass_dzero, mass_dzero);
  RooRealVar sigma1("sigma1", "width of gaussians", 0.02, 0.01, 0.03);
  RooRealVar sigma2("sigma2", "width of gaussians", 0.005, 0.001, 0.01);

  RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
  RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);

  // Build Chebychev polynomial p.d.f. as combinatorial background
  RooRealVar a1("a1", "a1", -0.4, -1., 1.);
  RooRealVar a2("a2", "a2", 0.0, -1., 1.);
  RooRealVar a3("a3", "a3", 0.0, -1., 1.);
  RooChebychev bkg("bkg", "Background", x, RooArgSet(a1, a2, a3));

  // gaussian as swapped k pi mass
  // RooRealVar means("means", "mean of swapped mass", 1.86, 1.8, 1.92);
  RooRealVar sigmas("sigmas", "width of swapped mass", 0.1, 0.05, 2);
  RooGaussian swapped("swapped", "swapped k pi mass", x, mean, sigmas);

  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.5, 0., 1.);
  RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);

  const long num_max = 400000000000;
  RooRealVar nsig("nsig", "number of signal events", 600000, 0., num_max);
  RooRealVar nbkg("nbkg", "number of background events", 6000000, 0, num_max);
  // RooRealVar nswapped("nswapped", "number of swapped events", 500000, 0, 10000000);

  // // Associated nsig/nbkg as expected number of events with sig/bkg
  // RooExtendPdf esig("esig", "extended signal p.d.f", sig, nsig);
  // RooExtendPdf ebkg("ebkg", "extended background p.d.f", bkg, nbkg);
  // RooExtendPdf eswapped("eswapped", "extended swapped p.d.f", swapped, nswapped);

  RooAddPdf model("model", "composite pdf", RooArgList(bkg, swapped, sig), RooArgList(nbkg, nsig, nsig));

  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  RooDataHist dh("dh", "dh", x, Import(*h));

  // set range
  // x.setRange("D0range",  min_hist_dzero, max_hist_dzero);
  // x.setRange("D0range", 1.74, 1.99);
  x.setRange("D0range", 1.73, max_hist_dzero);

  // U n b i n n e d   M L   f i t   t o   w e i g h t e d   d a t a
  // A first order correction to estimated parameter errors
  // RooFitResult *r_ml_wgt_corr = model.fitTo(dh, Save(), SumW2Error(kTRUE), Range("D0range"), Extended(kTRUE));

  RooFitResult *r_ml_wgt_corr = model.fitTo(dh, Save(), Range("D0range"), Extended(kTRUE));

  // RooFitResult *r_ml_wgt_corr = model.chi2FitTo(dh, Save(), Range("D0range"), Extended(kTRUE));
  r_ml_wgt_corr->Print();

  // Calculate yields and error
  // yield = Nintegral(sig, x, mass_dzero_signal_l, mass_dzero_signal_h);
  // yieldErr = yield * nsig.getError() / nsig.getValV();
  yield = nsig.getValV();
  yieldErr = nsig.getError();

  // calculate background scale
  RooAddPdf not_mass("not_mass", "backgrounds and combinatorial", RooArgList(bkg, swapped), RooArgList(nbkg, nsig));
  sidebandScale = Nintegral(not_mass, x, mass_dzero_signal_l, mass_dzero_signal_h) /
    (Nintegral(not_mass, x, mass_dzero_sideband_h_n, mass_dzero_sideband_l_n) +
     Nintegral(not_mass, x, mass_dzero_sideband_l_p, mass_dzero_sideband_h_p));


  auto canv = new TCanvas("Canvas", "Canvas", 800, 600);

  RooPlot * frame = x.frame(Title("Full range fitted"));
  dh.plotOn(frame, DataError(RooAbsData::SumW2));
  model.plotOn(frame, Range("D0range"));
  model.plotOn(frame, Range("D0range"), VisualizeError(*r_ml_wgt_corr));
  model.plotOn(frame, Range("D0range"));
  dh.plotOn(frame);

  model.plotOn(frame, Components(sig), LineStyle(kDotted), Range("D0range"));
  model.plotOn(frame, Components(bkg), LineStyle(kDashed), Range("D0range"));
  model.plotOn(frame, Components(swapped), LineStyle(kDashed), LineColor(kRed), Range("D0range"));

  frame->Draw();
  // c->Draw();
  // canv->Draw();

  // c->cd();
  // h->Draw("e");
  // fun_background->Draw("same");
  // fun_mass->Draw("same");
  // fun_swap->Draw("same");

  if(fdrawsigsid || fdrawsignif)
    {
      // fun_not_mass->SetRange(mass_dzero_signal_l,mass_dzero_signal_h);
      // fun_not_mass->Draw("same");
      // drawline(mass_dzero_signal_l, 0, mass_dzero_signal_l, fun_f->Eval(mass_dzero_signal_l), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      // drawline(mass_dzero_signal_h, 0, mass_dzero_signal_h, fun_f->Eval(mass_dzero_signal_h), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_signal_l, 0, mass_dzero_signal_l, 1000000, fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_signal_h, 0, mass_dzero_signal_h, 1000000, fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    }
  if(fdrawsigsid)
    {
      std::cout << mass_dzero_sideband_l_n << ", " << mass_dzero_sideband_l_p <<
        mass_dzero_sideband_h_n << ", " << mass_dzero_sideband_h_p << "\n";
      // drawline(mass_dzero_sideband_l_p, 0, mass_dzero_sideband_l_p, fun_f->Eval(mass_dzero_sideband_l_p), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      // drawline(mass_dzero_sideband_h_p, 0, mass_dzero_sideband_h_p, fun_f->Eval(mass_dzero_sideband_h_p), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      // drawline(mass_dzero_sideband_l_n, 0, mass_dzero_sideband_l_n, fun_f->Eval(mass_dzero_sideband_l_n), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      // drawline(mass_dzero_sideband_h_n, 0, mass_dzero_sideband_h_n, fun_f->Eval(mass_dzero_sideband_h_n), fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_sideband_l_p, 0, mass_dzero_sideband_l_p, 1000000, fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_sideband_h_p, 0, mass_dzero_sideband_h_p, 1000000, fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_sideband_l_n, 0, mass_dzero_sideband_l_n, 1000000, fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
      drawline(mass_dzero_sideband_h_n, 0, mass_dzero_sideband_h_n, 1000000, fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    }
  // fun_f->Draw("same");

  // drawleg(h);
  drawCMS(collisionsyst);

  Float_t texxpos = 0.22, texypos = 0.90, texdypos = 0.053;
  if(!vtex.empty())
    {
      texypos+=texlinespc;
      for(std::vector<TString>::const_iterator it=vtex.begin(); it!=vtex.end(); it++)
        drawtex(texxpos, texypos=(texypos-texdypos-texlinespc), *it);
    }
  if(fdrawyield) drawtex(texxpos, texypos=(texypos-texdypos), Form("N = %.0f #pm %.0f", yield, yieldErr));
  // if(fdrawchi2)
  //   {
  //     drawtex(texxpos, texypos=(texypos-texdypos), Form("#chi^{2} / ndf = %.1f / %.0f",Chi2,NDF));
  //     drawtex(texxpos, texypos=(texypos-texdypos), Form("Prob = %.2f%s",Chi2Prob*100,"%"));
  //   }
  // if(fdrawsignif)
  //   {
  //     drawtex(texxpos, texypos=(texypos-texdypos), Form("S = %.0f, B = %.0f",S,B));
  //     drawtex(texxpos, (texypos=(texypos-texdypos)) + 0.01, Form("S/#sqrt{S+B} = %.1f",S/TMath::Sqrt(S+B)));
  //   }
  // if(fdrawdetail)
  //   {
  //     drawtex(texxpos, texypos=(texypos-texdypos), Form("N#scale[0.6]{#lower[0.7]{sig}}/(N#scale[0.6]{#lower[0.7]{sig}}+N#scale[0.6]{#lower[0.7]{swap}}) = %.2f",fun_f->GetParameter(7)));
  //   }

  if(fsaveplot)
    {
      // c->SaveAs(Form("%s.pdf",outputname.Data()));
      canv->SaveAs(Form("%s.pdf", outputname.Data()));
    }

  delete c;
  delete h;
  delete hMCSignal;
  delete hMCSwapped;

  return clonefun(fun_mass, "Fun_mass");
}

/*
  Do an unbinned extended ML fit on D0 mass in a certain dphi bin
  Calculate the signal yield and error
  @param tree    TTree*  tree containing data
  @param isSig   int     1:signal 0:sideband
  @param iBin    int     ID of delta phi bin
*/
void xjjroot::droofit::fit(
    TTree *tree, const int isSig, const int iBin, TString collisionsyst /*=""*/,
    TString outputname /*="cmass"*/,
    const std::vector<TString> &vtex /*=std::vector<TString>()*/) {
  reset();
  init();

  TString fitoption = ffitverbose ? "L m" : "L m q";
  setgstyle();

  // Declare an observable for D0 mass
  RooRealVar x("mass", "m_{#pi K} / GeV", min_hist_dzero, max_hist_dzero);

  // weight of the D0 mass
  RooRealVar w("weight", "weight", min_weight, max_weight);

  // set range for D0 mass
  x.setRange("D0range", 1.73, max_hist_dzero);

  RooRealVar iPhi("iPhi", "dphi bin ID", 0, 20);
  RooRealVar sigcat("signal", "signal or sideband", 0, 1);
  // Construct unbinned dataset importing tree branches mass and weight matching
  // between branches and RooRealVars
  RooDataSet imported("ds", "ds", RooArgSet(x, w, iPhi, sigcat), Import(*tree),
                      WeightVar(w));

  // Select the dphi bin in question
  RooDataSet *ds = (RooDataSet *)imported.reduce(RooArgSet(x, w),
      TString::Format("signal == %i && iPhi == %i", isSig, iBin));

  // Double gaussian as signal p.d.f.
  // Fix the mean at D0 mass
  RooRealVar mean("mean", "mean of gaussians", mass_dzero, mass_dzero, mass_dzero);
  RooRealVar sigma1("sigma1", "width of gaussians", 0.02, 0.01, 0.03);
  RooRealVar sigma2("sigma2", "width of gaussians", 0.005, 0.001, 0.01);
  RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
  RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);

  // Build 3rd-order Chebychev polynomial p.d.f. as combinatorial background
  RooRealVar a1("a1", "a1", -0.4, -1., 1.);
  RooRealVar a2("a2", "a2", 0.0, -1., 1.);
  RooRealVar a3("a3", "a3", 0.0, -1., 1.);
  RooChebychev bkg("bkg", "Background", x, RooArgSet(a1, a2, a3));

  // Gaussian as swapped K pi mass
  RooRealVar sigmas("sigmas", "width of swapped mass", 0.1, 0.05, 2);
  RooGaussian swapped("swapped", "swapped k pi mass", x, mean, sigmas);

  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.5, 0., 1.);
  RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);

  const long num_max = 400000000000;
  // Initially, Nsig is set to 3% of the total entries, and Nswap = Nsig.
  RooRealVar nsig("nsig", "number of signal entries",
                  0.03 * ds->sumEntries("", "D0range"), 0., num_max);
  RooRealVar nbkg("nbkg", "number of background entries",
                  0.95 * ds->sumEntries("", "D0range"), 0, num_max);

  // // Associated nsig/nbkg as expected number of events with sig/bkg
  // RooExtendPdf esig("esig", "extended signal p.d.f", sig, nsig);
  // RooExtendPdf ebkg("ebkg", "extended background p.d.f", bkg, nbkg);
  // RooExtendPdf eswapped("eswapped", "extended swapped p.d.f", swapped,
  // nswapped);

  RooAddPdf model("model", "composite pdf", RooArgList(bkg, swapped, sig),
                  RooArgList(nbkg, nsig, nsig));

  // U n b i n n e d   M L   f i t   t o   w e i g h t e d   d a t a
  // Calculate the covariance matrix as V' = V C-1 V
  // Extended to get Nsig directly
  RooFitResult *r_ml_wgt_corr =
      model.fitTo(*ds, Save(), SumW2Error(kTRUE), Range("D0range"),
                  Extended(kTRUE), NumCPU(6));
  r_ml_wgt_corr->Print();

  // Calculate yields and error
  // yield = Nintegral(sig, x, mass_dzero_signal_l, mass_dzero_signal_h);
  // yieldErr = yield * nsig.getError() / nsig.getValV();
  yield = nsig.getValV();
  yieldErr = nsig.getError();

  // Calculate background scale
  RooAddPdf not_mass("not_mass", "backgrounds and combinatorial",
                     RooArgList(bkg, swapped), RooArgList(nbkg, nsig));
  sidebandScale =
      Nintegral(not_mass, x, mass_dzero_signal_l, mass_dzero_signal_h) /
      (Nintegral(not_mass, x, mass_dzero_sideband_h_n,
                 mass_dzero_sideband_l_n) +
       Nintegral(not_mass, x, mass_dzero_sideband_l_p,
                 mass_dzero_sideband_h_p));

  // Plot the result
  auto canv = new TCanvas("Canvas", "Canvas", 800, 600);

  RooPlot *frame = x.frame(1.73, max_hist_dzero);
  ds->plotOn(frame, Name("pdat"), DataError(RooAbsData::SumW2));
  model.plotOn(frame, VisualizeError(*r_ml_wgt_corr));
  model.plotOn(frame, Name("pfit"));
  model.plotOn(frame, Name("psig"), Components(sig), DrawOption("LF"),
               FillStyle(3002), FillColor(kOrange - 3), LineStyle(2),
               LineColor(kOrange - 3));
  model.plotOn(frame, Name("pswp"), Components(swapped), DrawOption("LF"),
               FillStyle(3005), FillColor(kGreen + 4), LineStyle(1),
               LineColor(kGreen + 4));
  model.plotOn(frame, Name("pbkg"), Components(bkg), LineStyle(2),
               LineColor(4));

  frame->Draw();
  if (fdrawsigsid || fdrawsignif) {
    // fun_not_mass->SetRange(mass_dzero_signal_l,mass_dzero_signal_h);
    // fun_not_mass->Draw("same");
    // drawline(mass_dzero_signal_l, 0, mass_dzero_signal_l,
    // fun_f->Eval(mass_dzero_signal_l), fun_not_mass->GetLineColor(),
    // fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    // drawline(mass_dzero_signal_h, 0, mass_dzero_signal_h,
    // fun_f->Eval(mass_dzero_signal_h), fun_not_mass->GetLineColor(),
    // fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    drawline(mass_dzero_signal_l, 0, mass_dzero_signal_l, 1000000,
             fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(),
             fun_not_mass->GetLineWidth());
    drawline(mass_dzero_signal_h, 0, mass_dzero_signal_h, 1000000,
             fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(),
             fun_not_mass->GetLineWidth());
  }
  if (fdrawsigsid) {
    std::cout << mass_dzero_sideband_l_n << ", " << mass_dzero_sideband_l_p
              << mass_dzero_sideband_h_n << ", " << mass_dzero_sideband_h_p
              << "\n";
    // drawline(mass_dzero_sideband_l_p, 0, mass_dzero_sideband_l_p,
    // fun_f->Eval(mass_dzero_sideband_l_p), fun_not_mass->GetLineColor(),
    // fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    // drawline(mass_dzero_sideband_h_p, 0, mass_dzero_sideband_h_p,
    // fun_f->Eval(mass_dzero_sideband_h_p), fun_not_mass->GetLineColor(),
    // fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    // drawline(mass_dzero_sideband_l_n, 0, mass_dzero_sideband_l_n,
    // fun_f->Eval(mass_dzero_sideband_l_n), fun_not_mass->GetLineColor(),
    // fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    // drawline(mass_dzero_sideband_h_n, 0, mass_dzero_sideband_h_n,
    // fun_f->Eval(mass_dzero_sideband_h_n), fun_not_mass->GetLineColor(),
    // fun_not_mass->GetLineStyle(), fun_not_mass->GetLineWidth());
    drawline(mass_dzero_sideband_l_p, 0, mass_dzero_sideband_l_p, 1000000,
             fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(),
             fun_not_mass->GetLineWidth());
    drawline(mass_dzero_sideband_h_p, 0, mass_dzero_sideband_h_p, 1000000,
             fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(),
             fun_not_mass->GetLineWidth());
    drawline(mass_dzero_sideband_l_n, 0, mass_dzero_sideband_l_n, 1000000,
             fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(),
             fun_not_mass->GetLineWidth());
    drawline(mass_dzero_sideband_h_n, 0, mass_dzero_sideband_h_n, 1000000,
             fun_not_mass->GetLineColor(), fun_not_mass->GetLineStyle(),
             fun_not_mass->GetLineWidth());
  }

  TLegend *leg = new TLegend(0.65, 0.58, 0.85, 0.88, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  leg->AddEntry("pdat", "Data", "pl");
  leg->AddEntry("pfit", "Fit", "l");
  leg->AddEntry("psig", "D^{0}+#bar{D^{#lower[0.2]{0}}} Signal", "f");
  leg->AddEntry("pswp", "K-#pi swapped", "f");
  leg->AddEntry("pbkg", "Combinatorial", "l");

  leg->Draw("same");

  drawCMS(collisionsyst);

  Float_t texxpos = 0.22, texypos = 0.90, texdypos = 0.053;
  if (!vtex.empty()) {
    texypos += texlinespc;
    for (std::vector<TString>::const_iterator it = vtex.begin();
         it != vtex.end(); it++)
      drawtex(texxpos, texypos = (texypos - texdypos - texlinespc), *it);
  }
  if (fdrawyield) {
    drawtex(texxpos, texypos = (texypos - texdypos),
            Form("N = %.0f #pm %.0f", yield, yieldErr));
  }

  if (fsaveplot) {
    canv->SaveAs(Form("%s.pdf", outputname.Data()));
  }
  return;
}

void xjjroot::droofit::reset()
{
  clearvar();
  deletefun();
}

void xjjroot::droofit::init()
{
  clearvar();
  createfun();
  setfunstyle();
}

void xjjroot::droofit::clearvar()
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

void xjjroot::droofit::calvar()
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
  // yield = fun_mass->Integral(min_hist_dzero,max_hist_dzero)/binwid_hist_dzero;
  // yieldErr = fun_mass->Integral(min_hist_dzero,max_hist_dzero)/binwid_hist_dzero*fun_mass->GetParError(0)/fun_mass->GetParameter(0);
  Chi2 = 2.*r->MinFcnValue();
  NDF = fun_f->GetNDF();
  Chi2Prob = TMath::Prob(Chi2,NDF);
}

void xjjroot::droofit::createfun()
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

void xjjroot::droofit::deletefun()
{
  delete fun_f;
  delete fun_background;
  delete fun_mass;
  delete fun_swap;
  delete fun_not_mass;

  fparamfun_f = false;
  fparamfuns = false;
}

void xjjroot::droofit::setfunstyle()
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

void xjjroot::droofit::setfunparameters()
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

TF1* xjjroot::droofit::clonefun(const TF1* fun, TString fun_name) const
{
  TF1* newfun = new TF1(*fun);
  newfun->SetName(fun_name);
  return newfun;
}

void xjjroot::droofit::sethist(TH1* h) const
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

void xjjroot::droofit::drawleg(TH1* h) const
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

void xjjroot::droofit::drawCMS(TString collision, TString snn/*="5.02"*/) const
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

void xjjroot::droofit::drawtex(Double_t x, Double_t y, const char* text, Float_t tsize/*=0.04*/, Short_t align/*=12*/) const
{
  TLatex* tex = new TLatex(x, y, text);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
  tex->Draw();
}

void xjjroot::droofit::drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=2*/) const
{
  TLine* l = new TLine(x1, y1, x2, y2);
  l->SetLineColor(lcolor);
  l->SetLineStyle(lstyle);
  l->SetLineWidth(lwidth);
  l->Draw();
}

void xjjroot::droofit::setgstyle() const
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
