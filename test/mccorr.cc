#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRandomizeParamMCSModule.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"

#include <algorithm>

// whether to set the initial values of the fitting parameters to the answer
constexpr bool cheat = true;
// whether to include the swap components
constexpr bool use_swap = false;
// whether to use data instead of MC (Not implemented yet)
constexpr bool use_data = false;
// whether to fix the signal shape during the fit
constexpr bool fix_sig = false;



const double nSigSig = 6;
const double nSigBkg = 167;
const double nBkgSig = 184;
const double nBkgBkg = 5642;

using namespace RooFit;

namespace color {
  int night0 = TColor::GetColor("#2e3440");
  int night1 = TColor::GetColor("#3b4252");
  int night2 = TColor::GetColor("#434c5e");
  int night3 = TColor::GetColor("#4c566a");
  int snow0 = TColor::GetColor("#d8dee9");
  int snow1 = TColor::GetColor("#e5e9f0");
  int snow2 = TColor::GetColor("#eceff4");
  int frost0 = TColor::GetColor("#8fbcbb");
  int frost1 = TColor::GetColor("#88c0d0");
  int frost2 = TColor::GetColor("#81a1c1");
  int frost3 = TColor::GetColor("#5e81ac");
  int aurora0 = TColor::GetColor("#bf616a");
  int aurora1 = TColor::GetColor("#d08770");
  int aurora2 = TColor::GetColor("#ebcb8b");
  int aurora3 = TColor::GetColor("#a3be8c");
  int aurora4 = TColor::GetColor("#b48ead");
};

struct vars
{
  Int_t status;
  double nss;
  double nsserr;
  double nbb;
  double nbberr;
  double nsb;
  double nsberr;
  double nbs;
  double nbserr;
};

void pull(std::vector<vars>);

void mccorr(UInt_t nsamples = 200, bool exit_on_corr = false) {

  const Double_t mass_dzero = 1.8649;

  const Double_t min_fit_dzero = 1.73;
  const Double_t max_hist_dzero = 2.0;

  RooRealVar m1 = {"m1", "m_{D} / GeV", min_fit_dzero, max_hist_dzero};
  RooRealVar m2 = {"m2", "m_{#bar{D}} / GeV", min_fit_dzero, max_hist_dzero};

  RooRealVar mean("mean", "mean of gaussians", mass_dzero, mass_dzero - 0.01,
                  mass_dzero + 0.01);
  RooRealVar sigma1x("sigma1x", "width of gaussians", 0.04, 0.011, 0.15);
  RooRealVar sigma2x("sigma2x", "width of gaussians", 0.01, 0.001, 0.011);
  RooGaussian sig1x("sig1x", "Signal component 1", m1, mean, sigma1x);
  RooGaussian sig2x("sig2x", "Signal component 2", m1, mean, sigma2x);

  RooRealVar sigma1y("sigma1y", "width of gaussians", 0.04, 0.011, 0.15);
  RooRealVar sigma2y("sigma2y", "width of gaussians", 0.01, 0.001, 0.011);
  RooGaussian sig1y("sig1y", "Signal component 1", m2, mean, sigma1y);
  RooGaussian sig2y("sig2y", "Signal component 2", m2, mean, sigma2y);

  // Build 3rd-order Chebychev polynomial p.d.f. as combinatorial background
  RooRealVar a1("a1", "a1", -0.2, -.8, .8);
  RooRealVar a2("a2", "a2", 0.02, -.1, .1);
  RooRealVar a3("a3", "a3", 0.0, -.1, .1);
  RooChebychev bkgx("bkgx", "Background", m1, RooArgSet(a1, a2, a3));

  RooRealVar b1("b1", "b1", -0.2, -.8, .8);
  RooRealVar b2("b2", "b2", 0.02, -.1, .1);
  RooRealVar b3("b3", "b3", 0.0, -.1, .1);
  RooChebychev bkgy("bkgy", "Background", m2, RooArgSet(b1, b2, b3));

  // Gaussian as swapped K pi mass
  RooRealVar sigmasx("sigmasx", "width of swapped mass", 0.1, 0.03, 2);
  RooRealVar sigmasy("sigmasy", "width of swapped mass", 0.1, 0.03, 2);

  RooGaussian swpx("swpx", "swapped k pi mass", m1, mean, sigmasx);
  RooGaussian swpy("swpy", "swapped k pi mass", m2, mean, sigmasy);

  RooProdPdf swapped("swapped", "background PDF", RooArgSet(swpx, swpy));

  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sigfracx("sigfracx", "fraction of component 1 in signal", 0.5,
                      0.001, 0.999);
  RooAddPdf sigx("sigx", "Signal", RooArgList(sig1x, sig2x), sigfracx);
  RooRealVar sigfracy("sigfracy", "fraction of component 1 in signal", 0.5,
                      0.001, 0.999);
  RooAddPdf sigy("sigy", "Signal", RooArgList(sig1y, sig2y), sigfracy);

  RooProdPdf sigsig("sigsig", "PDF", RooArgSet(sigx, sigy));
  RooProdPdf sigswp("sigswp", "PDF", RooArgSet(sigx, swpy));
  RooProdPdf sigbkg("sigbkg", "PDF", RooArgSet(sigx, bkgy));
  RooProdPdf swpsig("swpsig", "PDF", RooArgSet(swpx, sigy));
  RooProdPdf swpswp("swpswp", "PDF", RooArgSet(swpx, swpy));
  RooProdPdf swpbkg("swpbkg", "PDF", RooArgSet(swpx, bkgy));
  RooProdPdf bkgsig("bkgsig", "PDF", RooArgSet(bkgx, sigy));
  RooProdPdf bkgswp("bkgswp", "PDF", RooArgSet(bkgx, swpy));
  RooProdPdf bkgbkg("bkgbkg", "PDF", RooArgSet(bkgx, bkgy));

  int num_max = 5999;
  RooRealVar nss("nss", "number of signal entries", nSigSig, -0.02 * num_max,
                 0.02 * num_max);
  RooRealVar nsb("nsb", "number of swapped entries", nSigBkg,
                 0.001 * num_max, 0.3 * num_max);
  RooRealVar nbs("nbs", "number of swapped entries", nBkgSig,
                 0.001 * num_max, 0.3 * num_max);
  RooRealVar nbb("nbb", "number of background entries", nBkgBkg,
                 0 * num_max, 1.1 * num_max);

  // 2D PDFs
  RooAddPdf modelsig("modelsig", "composite pdf", RooArgList(sigsig), RooArgList(nss));
  RooAddPdf modelsigbkg("modelsigbkg", "composite pdf", RooArgList(sigbkg), RooArgList(nsb));
  RooAddPdf modelbkgsig("modelbkgsig", "composite pdf", RooArgList(bkgsig), RooArgList(nbs));
  RooAddPdf modelbkg("modelbkg", "composite pdf", RooArgList(bkgbkg), RooArgList(nbb));

  RooAddPdf modelswpswp("modelswpswp", "composite pdf", RooArgList(swpswp), RooArgList(nss));
  RooAddPdf modelsigswp("modelsigswp", "composite pdf", RooArgList(sigswp), RooArgList(nss));
  RooAddPdf modelswpsig("modelswpsig", "composite pdf", RooArgList(swpsig), RooArgList(nss));
  RooAddPdf modelswpbkg("modelswpbkg", "composite pdf", RooArgList(swpbkg), RooArgList(nsb));
  RooAddPdf modelbkgswp("modelbkgswp", "composite pdf", RooArgList(bkgswp), RooArgList(nbs));

  RooAddPdf modelsw("modelsw", "composite pdf of signal and swap",
                    RooArgList(sigsig, sigswp, swpsig, swpswp), RooArgList(nss, nss, nss, nss));

  RooAddPdf modelsb(
      "modelsb", "composite pdf",
      RooArgList(sigsig, sigbkg, bkgsig, bkgbkg), RooArgList(nss, nsb, nbs, nbb));

  // model including the swap comonent, using the same event number as signal
  RooAddPdf modelsbw(
                     "modelsbw", "composite pdf",
                     RooArgList(sigsig, sigswp, sigbkg, swpsig, swpswp, swpbkg,
                                bkgsig, bkgswp, bkgbkg),
                     RooArgList(nss, nss, nsb, nss, nss, nsb, nbs, nbs, nbb));

  RooAddPdf* modelptr;
  if (!use_swap) {
    modelptr = new RooAddPdf(modelsb, "model");
      } else {
    modelptr = new RooAddPdf(modelsbw, "model");
  }
  RooAddPdf model = *modelptr;


  auto mset = RooArgSet(m1, m2);
  RooDataSet* dss = modelsig.generate(mset, nSigSig);
  RooDataSet* dsb = modelsigbkg.generate(mset, nSigBkg);
  RooDataSet* dbs = modelbkgsig.generate(mset, nBkgSig);
  RooDataSet* dbb = modelbkg.generate(mset, nBkgBkg);
  RooDataSet* sum = new RooDataSet("sum", "sum", mset);
  RooDataSet* dww = modelswpswp.generate(mset, nSigSig);
  RooDataSet* dsw = modelsigswp.generate(mset, nSigSig);
  RooDataSet* dws = modelswpsig.generate(mset, nSigSig);
  RooDataSet* dwb = modelswpbkg.generate(mset, nSigBkg);
  RooDataSet* dbw = modelbkgswp.generate(mset, nBkgSig);

  vars fitted;
  TTree* res = new TTree("res", "res");
  res->Branch("status", &fitted.status, "status/I");
  res->Branch("nss", &fitted.nss, "nss/D");
  res->Branch("nsserr", &fitted.nsserr, "nsserr/D");
  res->Branch("nbb", &fitted.nbb, "nbb/D");
  res->Branch("nbberr", &fitted.nbberr, "nbberr/D");
  res->Branch("nsb", &fitted.nsb, "nsb/D");
  res->Branch("nsberr", &fitted.nsberr, "nsberr/D");
  res->Branch("nbs", &fitted.nbs, "nbs/D");
  res->Branch("nbserr", &fitted.nbserr, "nbserr/D");
  res->Branch("dss", &dss);
  res->Branch("dsb", &dsb);
  res->Branch("dbs", &dbs);
  res->Branch("dbb", &dbb);
  res->Branch("sum", &sum);

  // Suppress most messages
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().setSilentMode(true);

  // plotting
  auto canv = new TCanvas("ptcanvas", "Canvas", 800, 600);
  canv->SetFillColor(color::snow2);

  std::vector<vars> fitresults;
  ROOT::EnableThreadSafety();
  for (unsigned i = 0; i < nsamples; ++i) {
      // Generate each component and sum them up
      sum = new RooDataSet("sum", "sum", mset);
      dss = modelsig.generate(mset, nSigSig);
      dsb = modelsigbkg.generate(mset, nSigBkg);
      dbs = modelbkgsig.generate(mset, nBkgSig);
      dbb = modelbkg.generate(mset, nBkgBkg);
      sum->append(*dss);
      sum->append(*dsb);
      sum->append(*dbs);
      sum->append(*dbb);

      dww = modelswpswp.generate(mset, nSigSig);
      dsw = modelsigswp.generate(mset, nSigSig);
      dws = modelswpsig.generate(mset, nSigSig);
      dwb = modelswpbkg.generate(mset, nSigBkg);
      dbw = modelbkgswp.generate(mset, nBkgSig);


      if (use_swap) {
        sum->append(*dww);
        sum->append(*dsw);
        sum->append(*dws);
        sum->append(*dbw);
        sum->append(*dwb);
      }

      // Use gen signals and swaps to fit the shapes of signal and swap PDF
      // Fix the signal shape during fitting
      if (fix_sig) {
        RooDataSet* dnobkg = new RooDataSet("dnobkg", "PDF about signals and swaps", mset);

        auto dssmc = modelsig.generate(mset, 30 * nSigSig);
        auto dwwmc = modelswpswp.generate(mset, 30 * nSigSig);
        auto dswmc = modelsigswp.generate(mset, 30 * nSigSig);
        auto dwsmc = modelswpsig.generate(mset, 30 * nSigSig);

        dnobkg->append(*dssmc);
        dnobkg->append(*dswmc);
        dnobkg->append(*dwsmc);
        dnobkg->append(*dwwmc);

        sigfracx.setConstant(false);
        sigma1x.setConstant(false);
        sigma2x.setConstant(false);
        sigfracy.setConstant(false);
        sigma1y.setConstant(false);
        sigma2y.setConstant(false);
        sigmasx.setConstant(false);
        sigmasy.setConstant(false);
        RooFitResult* result;
        result = modelsw.fitTo(*dnobkg, Extended(), Save(), NumCPU(4));

        sigfracx.setConstant();
        sigma1x.setConstant();
        sigma2x.setConstant();
        sigfracy.setConstant();
        sigma1y.setConstant();
        sigma2y.setConstant();
        sigmasx.setConstant();
        sigmasy.setConstant();

        // result->Print();
      }

      // reset fitted variables to the "answer"
      if (cheat) {
        nss.setVal(nSigSig);
        nsb.setVal(nSigBkg);
        nbs.setVal(nBkgSig);
        nbb.setVal(nBkgBkg);
        mean.setVal(mass_dzero);
        a1.setVal(-0.2);
        a2.setVal(0.02);
        a3.setVal(0.0);
        b1.setVal(-0.2);
        b2.setVal(0.02);
        b3.setVal(0.0);
      }

      // fit on the combined PDF
      RooFitResult* result;
      result = model.fitTo(*sum, Extended(), Save(), NumCPU(6));
      fitted.status = result->status();

      std::cout << "Fitting result: " <<
        "nss: " << nss.getValV() <<
        ", nsb: " << nsb.getValV() <<
        ", nbs: " << nbs.getValV() <<
        ", nbb: " << nbb.getValV() << "\n";

      fitted.nss = nss.getValV();
      fitted.nsserr = nss.getError();
      fitted.nbb = nbb.getValV();
      fitted.nbberr = nbb.getError();
      fitted.nbs = nbs.getValV();
      fitted.nbserr = nbs.getError();
      fitted.nsb = nsb.getValV();
      fitted.nsberr = nsb.getError();
      res->Fill();

      fitresults.push_back(fitted);
      if (fitted.status) {
        std::cout << "=== Fitting failed for event " << i << "! ===" << "\n";
        std::cout << "Nss global correlation coeff.: " << result->globalCorr("nss") << "\n";
        std::cout << "Nbb global correlation coeff.: " << result->globalCorr("nbb") << "\n";
        std::cout << "Nsb global correlation coeff.: " << result->globalCorr("nsb") << "\n";
        std::cout << "Nbs global correlation coeff.: " << result->globalCorr("nbs") << "\n";
        std::cout << "mean global correlation coeff.: " << result->globalCorr("mean") << "\n";
        std::cout << "a1 global correlation coeff.: " << result->globalCorr("a1") << "\n";
        std::cout << "a2 global correlation coeff.: " << result->globalCorr("a2") << "\n";
        std::cout << "a3 global correlation coeff.: " << result->globalCorr("a3") << "\n";
        std::cout << "b1 global correlation coeff.: " << result->globalCorr("b1") << "\n";
        std::cout << "b2 global correlation coeff.: " << result->globalCorr("b2") << "\n";
        std::cout << "b3 global correlation coeff.: " << result->globalCorr("b3") << "\n";
        result->Print();
        result->correlationMatrix().Print();
        std::cout << std::endl;


        // Produce string including multiple components with correct order
        auto compstr = [](TString name, TString thisdim,
                          std::initializer_list<TString> thatdim) {
          TString str;
          for (auto i : thatdim) {
            if (name == "m1") {
              str += thisdim + i + ",";
            } else {
              str += i + thisdim + ",";
            }
          }
          str.Remove(str.Length() - 1);
          return str;
        };

        std::vector<RooRealVar> dim = {m1, m2};
        for (auto m : dim) {
          TString name(m.GetName());
          RooPlot *plot = m.frame(Title(name + "_plot"), Bins(60));
          sum->plotOn(plot, Name("pdat" + name));
          modelsb.plotOn(plot, VisualizeError(*result));
          sum->plotOn(plot, Name("pdat" + name),
                      DataError(RooAbsData::SumW2));
          model.plotOn(plot, Name("pfit" + name));
          model.plotOn(plot, Name("psig" + name), Components("sigsig"),
                         DrawOption("LF"), FillStyle(3002), FillColor(color::aurora0),
                         LineStyle(2), LineColor(color::aurora0));
          model.plotOn(plot, Name("pbkg" + name), Components("bkgbkg"),
                         DrawOption("LF"), LineStyle(2), LineColor(color::frost1));
          model.plotOn(plot, Name("psig2" + name),
                         Components(compstr(name, "sig", {"bkg"})),
                         DrawOption("LF"), FillStyle(3002), FillColor(color::aurora2),
                         LineStyle(2), LineColor(color::aurora2));
          model.plotOn(plot, Name("pbkg2" + name),
                         Components(compstr(name, "bkg", {"sig"})),
                         DrawOption("LF"), LineStyle(2), LineColor(color::frost3));

          if (use_swap) {
            model.plotOn(plot, Name("psig3" + name),
                         ProjectionRange("signal_box_" + name),
                         // Components(compstr(name, "sig", {"sig", "swp", "bkg"})),
                         Components(compstr(name, "sig", {"swp"})), DrawOption("LF"),
                         FillStyle(3002), FillColor(color::aurora4), LineStyle(2),
                         LineColor(color::aurora4));
            model.plotOn(plot, Name("pbkg3" + name),
                         ProjectionRange("signal_box_" + name),
                         Components(compstr(name, "bkg", {"swp"})), LineStyle(2),
                         LineColor(color::frost2));
            model.plotOn(plot, Name("pswp3" + name),
                         ProjectionRange("signal_box_" + name),
                         Components(compstr(name, "swp", {"swp"})), DrawOption("LF"),
                         FillStyle(3005), FillColor(color::night0), LineStyle(1),
                         LineColor(color::night1));
            model.plotOn(plot, Name("pswp2" + name),
                         ProjectionRange("signal_box_" + name),
                         Components(compstr(name, "swp", {"bkg"})), DrawOption("LF"),
                         FillStyle(3005), FillColor(color::night1), LineStyle(1),
                         LineColor(color::night2));
            model.plotOn(plot, Name("pswp" + name),
                         ProjectionRange("signal_box_" + name),
                         Components(compstr(name, "swp", {"sig"})), DrawOption("LF"),
                         FillStyle(3005), FillColor(color::night2), LineStyle(1),
                         LineColor(color::night3));
          }

          plot->Draw();

          TLegend *leg = new TLegend(0.11, 0.18, 0.31, 0.78, NULL, "brNDC");
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextFont(42);
          leg->SetTextSize(0.04);

          leg->AddEntry("pdat" + name, "Data", "pl");
          leg->AddEntry("pfit" + name, "Fit", "l");

          std::vector<TString> dzerotex = {"D^{0}", "#bar{D^{#lower[0.2]{0}}}"};
          auto make_legend = [dzerotex](TString name, TString thisdim,
                                        TString thatdim) {
            TString str;
            if (name == "m1") {
              str = dzerotex[0] + thisdim + "+" + dzerotex[1] + thatdim;
            } else {
              str = dzerotex[0] + thatdim + "+" + dzerotex[1] + thisdim;
            }
            return str;
          };
          leg->AddEntry("psig2" + name, make_legend(name, "sig", "bkg"), "f");
          leg->AddEntry("psig" + name, make_legend(name, "sig", "sig"), "f");
          leg->AddEntry("pbkg" + name, make_legend(name, "bkg", "bkg"), "l");
          leg->AddEntry("pbkg2" + name, make_legend(name, "bkg", "sig"), "l");

          if (use_swap) {
            leg->AddEntry("pswp3" + name, make_legend(name, "swap", "swap"), "f");
            leg->AddEntry("pswp2" + name, make_legend(name, "swap", "bkg"), "f");
            leg->AddEntry("pswp" + name, make_legend(name, "swp", "sig"), "f");
            leg->AddEntry("pbkg3" + name, make_legend(name, "bkg", "swap"), "l");
            leg->AddEntry("psig3" + name, make_legend(name, "sig", "swap"), "f");
          }
          leg->Draw("same");

          canv->SaveAs(TString::Format("evt_%i_", i) + name + ".png");
        }

        if (exit_on_corr) break;
      }
    }

  gDirectory->Add(res);
  TFile fout("mccorr.root", "RECREATE");
  res->Write();
  fout.Close();
  pull(fitresults);
}

void pull(std::vector<vars> results) {
  TH1D hSigpull("sigpull", "pull of signal-signal", 60, -5, 5);
  TH1D hBkgpull("bkgpull", "pull of bkg-bkg", 60, -5, 5);
  for (auto evt : results) {
    double sigPull = (evt.nss - nSigSig) / evt.nsserr;
    double bkgPull = (evt.nbb - nBkgBkg) / evt.nbberr;
    hSigpull.Fill(sigPull);
    hBkgpull.Fill(bkgPull);
  }
  auto canv = new TCanvas("ptcanvas", "Canvas", 800, 600);
  canv->SetFillColor(color::snow2);
  hSigpull.Draw();
  canv->SaveAs("sigpull.png");
  hBkgpull.Draw();
  canv->SaveAs("bkgpull.png");
  return;
}
