#ifndef _XJJROOT_FIT2D_H_
#define _XJJROOT_FIT2D_H_

// Adapted from dfitter.h

/******************************************************************************************
 * Class : xjjroot::fit2d
 * This class provides tools fitting Dzero invariant mass spectra.                        *
 * The object to be used can be declared by                                               *
 *                                                                                        *
 *    xjjroot::fit2d obj(options);                                                      *
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

#include "TH1.h"
#include "TString.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFitResult.h"

// RooFit libs
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"

using namespace RooFit;

namespace xjjroot
{
  class fit2d
  {
  public:
    fit2d(TTree *sigtree, TTree *swaptree, Option_t *option = "")
        : sigtree(sigtree), swaptree(swaptree) {
      foption = option;
      resolveoption();
      init();
    }
    ~fit2d(){};

    RooFitResult *
    fit(const int iBin,
        TString collisionsyst = "", TString outputname = "cmass",
        const std::vector<TString> &vtex = std::vector<TString>());

    RooFitResult *
    simfit(TTree *sigtree, TTree *swaptree, const int iBin,
        TString collisionsyst = "", TString outputname = "cmass",
        const std::vector<TString> &vtex = std::vector<TString>());

    Double_t GetS() const {return S;}
    Double_t GetB() const {return B;}
    Double_t GetSig() const {return Sig;}
    Double_t GetY() const {return yield;}
    Double_t GetYE() const {return yieldErr;}
    Double_t GetChi2() const {return Chi2;}
    Double_t GetNDF() const {return NDF;}
    Double_t GetChi2Prob() const {return Chi2Prob;}

    void SetOption(Option_t* option="") {foption = option; resolveoption();}
    void SetSignalregion(Double_t d_mass_signal_) {d_mass_signal =  d_mass_signal_; calvar();}
    void SetSidebandL(Double_t d_mass_sideband_l_) {d_mass_sideband_l = d_mass_sideband_l_; calvar();}
    void SetSidebandH(Double_t d_mass_sideband_h_) {d_mass_sideband_h = d_mass_sideband_h_; calvar();}
    void SetTexLinespc(Double_t linespc=0) {texlinespc = linespc;}

    Double_t GetMassL() const {return min_hist_dzero;}
    Double_t GetMassH() const {return max_hist_dzero;}

    Double_t GetSidebandScale() const {return sidebandScale;}

  private:
    TTree *sigtree;
    TTree *swaptree;
    Double_t S;
    Double_t B;
    Double_t Sig;
    Double_t yield;
    Double_t yieldErr;
    Double_t Chi2;
    Double_t NDF;
    Double_t Chi2Prob;
    Double_t sidebandScale;

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
    const Double_t max_weight = (int) 1e7;

    Double_t texlinespc = 0;

    const Double_t min_fit_dzero = 1.73;

    // Declare an observable for D0 mass
    RooRealVar m1 = {"m1", "m_{D} / GeV", min_fit_dzero, max_hist_dzero};
    RooRealVar m2 = {"m2", "m_{#bar{D}} / GeV", min_fit_dzero, max_hist_dzero};

    RooRealVar iPhi = {"iPhi", "dphi bin ID", 0, 20};
    RooRealVar isSwap1 = {"isSwap1", "whether m1 is swapped", 0, 1};
    RooRealVar isSwap2 = {"isSwap2", "whether m2 is swapped", 0, 1};
    RooRealVar pt1 = {"pt1", "pt of D", 0, 60};
    RooRealVar pt2 = {"pt2", "pt of Dbar", 0, 60};
    RooRealVar pT1 = {"pT1", "pt of D", 0, 60};
    RooRealVar pT2 = {"pT2", "pt of Dbar", 0, 60};

    // weight of the D0 mass
    // RooRealVar w("weight", "weight", min_weight, max_weight);
    // RooRealVar werr("weightErr", "error of weight", min_weight, max_weight);

    // Construct an unbinned dataset from tree branches
    RooDataSet total_ds = {"ds", "ds", RooArgSet(m1, m2, iPhi, pT1, pT2),
                           Import(*sigtree)};

    RooDataSet swapds = {"swapds", "swapped dataset",
                         RooArgSet(m1, isSwap1, pt1, m2, isSwap2, pt2),
                         Import(*swaptree)};

    // Reference pt bins = {1.0,  2.0,  2.5,  3.0,  3.5,  4.0, 4.5,
    // 5.5, 6.0, 6.5, 7.0, 8.0, 10.0,
    // 12.5, 15.0, 20.0, 30.0, 40.0, 60.0};
    // Use wider ranges to ensure enough events in each pt bin
    const std::vector<double> ptbins = {1.0, 2.0, 2.5, 3.0, 3.5,
                                        4.0, 5.5, 6.5, 60.0};
    // Only go through the following (idphi, xpt, ypt) bins
    std::vector<std::array<double, 3> > skiptobins;

    void init();
    void reset();
    void createfun();
    void deletefun();
    void clearvar();
    void calvar();
    void resolveoption();
    void setfunparameters();
    void setfunstyle();

    void sethist(TH1* h) const;
    void drawCMS(TString collision, TString snn="5.02") const;
    void drawtex(Double_t x, Double_t y, const char* text, Float_t tsize=0.04, Short_t align=12, int color=kBlack) const;
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

#endif
