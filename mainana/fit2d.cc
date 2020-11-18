#include "fit2d.hpp"
#include "RooCategory.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "RooSimultaneous.h"
#include "TFile.h"

#include "color.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <numeric>

#include <fstream>

using namespace RooFit;

void xjjroot::fit2d::resolveoption()
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
  Do an unbinned extended ML fit on D0 mass in a certain dphi bin
  Calculate the signal yield and error
  @param sigtree    TTree*  tree containing data
  @param swaptree   TTree*  tree containing D0 mass and weight from swapped K pi
  @param iBin    int     ID of delta phi bin
*/
void xjjroot::fit2d::fit(
    const int iBin, TString collisionsyst /*=""*/,
    TString outputname /*="cmass"*/,
    const std::vector<TString> &vtex /*=std::vector<TString>()*/) {
  reset();
  init();

  // Suppress RooFit messages
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);

  TString fitoption = ffitverbose ? "L m" : "L m q";
  setgstyle();

  // set range for D0 mass
  m1.setRange("D0range", min_fit_dzero, max_hist_dzero);
  m2.setRange("D0range", min_fit_dzero, max_hist_dzero);

  // select the dphi bin in question
  RooDataSet ds = *(RooDataSet *)total_ds.reduce(RooArgSet(m1, m2, pT1, pT2),
                                                 Form("iPhi == %i", iBin));

  // Divide the data into pt bins
  using dsvec = std::vector<RooDataSet *>;
  std::vector<dsvec> dsets(ptbins.size() - 1, dsvec(ptbins.size() - 1));
  for (unsigned i = 0; i < ptbins.size() - 1; ++i) {
    for (unsigned j = 0; j < ptbins.size() - 1; ++j) {
      std::array<double, 3> ptvec = {(double)iBin, ptbins[i], ptbins[j]};
      if (skiptobins.size() > 0 &&
          std::find(skiptobins.cbegin(), skiptobins.cend(), ptvec) ==
              skiptobins.end()) {
        continue;
      }
      dsets[i][j] = (RooDataSet *)ds.reduce(
          RooArgSet(m1, m2),
          Form("pT1 > %f && pT1 < %f && pT2 > %f && pT2 < %f", ptbins[i],
               ptbins[i + 1], ptbins[j], ptbins[j + 1]));
    }
  }

  dsvec sigdsm1(ptbins.size() - 1);
  dsvec sigdsm2(ptbins.size() - 1);
  dsvec swapdsm1(ptbins.size() - 1);
  dsvec swapdsm2(ptbins.size() - 1);
  for (unsigned i = 0; i < ptbins.size() - 1; ++i) {
    sigdsm1[i] = (RooDataSet *)swapds.reduce(
        RooArgSet(m1),
        Form("!isSwap1 && pt1 > %f && pt1 < %f", ptbins[i], ptbins[i + 1]));
    sigdsm2[i] = (RooDataSet *)swapds.reduce(
        RooArgSet(m2),
        Form("!isSwap2 && pt2 > %f && pt2 < %f", ptbins[i], ptbins[i + 1]));
    swapdsm1[i] = (RooDataSet *)swapds.reduce(
        RooArgSet(m1),
        Form("isSwap1 && pt1 > %f && pt1 < %f", ptbins[i], ptbins[i + 1]));
    swapdsm2[i] = (RooDataSet *)swapds.reduce(
        RooArgSet(m2),
        Form("isSwap2 && pt2 > %f && pt2 < %f", ptbins[i], ptbins[i + 1]));
  }
  // Create a binned dataset (only when there is performance issue)
  // RooDataHist *swapdh = swapds.binnedClone();

  // Represent data in dh as pdf in x
  // RooHistPdf swapped("histpdfsw", "histpdfsw", x, *swapdh, 2);

  // Double gaussian as signal p.d.f.
  // Fix the mean at D0 mass
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

  auto canv = new TCanvas("Canvas", "Canvas", 800, 600);
  canv->SetFillColor(color::snow2);
  const int nPrefitBins = 60;
  m1.setRange("m1_3sigma", mass_dzero - 0.035, mass_dzero + 0.035);
  m2.setRange("m2_3sigma", mass_dzero - 0.035, mass_dzero + 0.035);
  // 1D fit
  for (unsigned ipt = 0; ipt < ptbins.size() - 1; ++ipt) {
    std::cout << ptbins[ipt] << " < pt < " << ptbins[ipt + 1] << std::endl;
    sigx.fitTo(*sigdsm1[ipt], Range("m1_3sigma"), Save());
    RooPlot *framex = m1.frame(Title("m1_plot"), Bins(nPrefitBins));
    sigdsm1[ipt]->plotOn(framex, DataError(RooAbsData::SumW2));
    sigx.plotOn(framex);
    sigx.paramOn(framex, Layout(0.6), Format("NEU", AutoPrecision(1)),
                 ShowConstants(false));
    framex->Draw();
    TLatex *pttex = new TLatex(
        0.4, 0.7, Form("%.1f < p_{T} < %.1f", ptbins[ipt], ptbins[ipt + 1]));
    pttex->SetNDC();
    pttex->SetTextAlign(32);
    pttex->SetTextSize(0.04);
    pttex->SetTextFont(42);
    pttex->Draw();
    canv->SaveAs(outputname +
                 Form("_sigfitm1_%.1f_%.1f.pdf", ptbins[ipt], ptbins[ipt + 1]));

    sigy.fitTo(*sigdsm2[ipt], Range("m2_3sigma"), Save());
    RooPlot *framey = m2.frame(Title("m2_plot"), Bins(nPrefitBins));
    sigdsm2[ipt]->plotOn(framey, DataError(RooAbsData::SumW2));
    sigy.plotOn(framey);
    sigy.paramOn(framey, Layout(0.6), Format("NEU", AutoPrecision(1)),
                 ShowConstants(false));
    framey->Draw();
    pttex->Draw();
    canv->SaveAs(outputname +
                 Form("_sigfitm2_%.1f_%.1f.pdf", ptbins[ipt], ptbins[ipt + 1]));

    swpy.fitTo(*swapdsm2[ipt], Save());
    RooPlot *frameswpy = m2.frame(Title("m2_plot"), Bins(nPrefitBins));
    swapdsm2[ipt]->plotOn(frameswpy, DataError(RooAbsData::SumW2));
    swpy.plotOn(frameswpy);
    swpy.paramOn(frameswpy, Layout(0.6), Format("NEU", AutoPrecision(1)),
                 ShowConstants(false));
    frameswpy->Draw();
    pttex->Draw();
    canv->SaveAs(outputname + Form("_swapfitm2_%.1f_%.1f.pdf", ptbins[ipt],
                                   ptbins[ipt + 1]));
  }

  std::vector<RooRealVar *> fixed_pars = {&sigma1x,  &sigma2x, &sigfracx,
                                          &sigmasx,  &sigma1y, &sigma2y,
                                          &sigfracy, &sigmasy};

  using vec2d = std::vector<std::vector<double>>;
  vec2d yields(ptbins.size() - 1);
  vec2d yieldErrs(ptbins.size() - 1);
  vec2d notsigs(ptbins.size() - 1);
  std::vector<std::array<double, 2>> failed;
  RooFitResult *r_ml_wgt_corr;
  // 2D fit
  for (unsigned xpt = 0; xpt < ptbins.size() - 1; ++xpt) {
    for (auto par : fixed_pars)
      par->setConstant(false);
    sigx.fitTo(*sigdsm1[xpt], Range("m1_3sigma"), Save());
    swpx.fitTo(*swapdsm1[xpt], Save());
    for (unsigned ypt = 0; ypt < ptbins.size() - 1; ++ypt) {
      for (auto par : fixed_pars)
        par->setConstant(false);
      sigy.fitTo(*sigdsm2[ypt], Range("m2_3sigma"), Save());
      swpy.fitTo(*swapdsm2[ypt], Save());
      for (auto par : fixed_pars)
        par->setConstant(true);

      std::array<double, 3> ptvec = {(double)iBin, ptbins[xpt], ptbins[ypt]};
      if (skiptobins.size() > 0 &&
          std::find(skiptobins.cbegin(), skiptobins.cend(), ptvec) ==
              skiptobins.end()) {
        yields[xpt].push_back(0);
        yieldErrs[xpt].push_back(0);
        notsigs[xpt].push_back(0);
        continue;
      }

      auto dataset = *dsets[xpt][ypt];
      const long num_max = dataset.sumEntries("", "D0range");

      if (num_max == 0) {
        yields[xpt].push_back(0);
        yieldErrs[xpt].push_back(0);
        notsigs[xpt].push_back(0);
        continue;
      }
      // Initially, Nsig is set to 3% of the total entries, and Nswap = Nsig.
      RooRealVar nss("nss", "number of signal entries", 0.03 * num_max,
                     0 * num_max, 0.2 * num_max);
      RooRealVar nsb("nsb", "number of swapped entries", 0.06 * num_max,
                     0.001 * num_max, 0.3 * num_max);
      RooRealVar nbs("nbs", "number of swapped entries", 0.06 * num_max,
                     0.001 * num_max, 0.3 * num_max);
      RooRealVar nbb("nbb", "number of background entries", 0.9 * num_max,
                     0.5 * num_max, 1.1 * num_max);

      // // Associated nsig/nbkg as expected number of events with sig/bkg
      // RooExtendPdf esig("esig", "extended signal p.d.f", sig, nsig);
      // RooExtendPdf ebkg("ebkg", "extended background p.d.f", bkg, nbkg);
      // RooExtendPdf eswapped("eswapped", "extended swapped p.d.f", swapped,
      // nsig); RooAddPdf model("model", "composite pdf", RooArgList(esig, ebkg,
      // eswapped));

      RooAddPdf model("model", "composite pdf",
                      RooArgList(sigsig, sigswp, sigbkg, swpsig, swpswp, swpbkg,
                                 bkgsig, bkgswp, bkgbkg),
                      RooArgList(nss, nss, nsb, nss, nss, nsb, nbs, nbs, nbb));

      // Unbinned ML fit
      // Extended to get Nsig directly
      r_ml_wgt_corr = model.fitTo(dataset, Save(), Extended(kTRUE), NumCPU(6));
      r_ml_wgt_corr->Print();

      // Record the pt of failed fit
      if (r_ml_wgt_corr->status() != 0) {
        failed.push_back({ptbins[xpt], ptbins[ypt]});
      }

      // Record yields and error
      yields[xpt].push_back(nss.getValV());
      yieldErrs[xpt].push_back(nss.getError());
      notsigs[xpt].push_back(nsb.getValV() + nbs.getValV() + nbb.getValV());

      projectionPlot(model, r_ml_wgt_corr, dataset, xpt, ypt, nss, outputname ,collisionsyst, vtex);
    }
  }
  // Fit the whole dataset without dividing pt
  sigdsm1.push_back((RooDataSet *)swapds.reduce(RooArgSet(m1), "!isSwap1"));
  sigdsm2.push_back((RooDataSet *)swapds.reduce(RooArgSet(m2), "!isSwap2"));
  swapdsm1.push_back((RooDataSet *)swapds.reduce(RooArgSet(m1), "isSwap1"));
  swapdsm2.push_back((RooDataSet *)swapds.reduce(RooArgSet(m2), "isSwap2"));
  for (auto par : fixed_pars)
    par->setConstant(false);
  sigx.fitTo(*sigdsm1.back(), Range("m1_3sigma"), Save());
  swpx.fitTo(*swapdsm1.back(), Save());
  sigy.fitTo(*sigdsm2.back(), Range("m2_3sigma"), Save());
  swpy.fitTo(*swapdsm2.back(), Save());
  for (auto par : fixed_pars)
    par->setConstant(true);

  const long num_max = ds.sumEntries("", "D0range");
  // Initially, Nsig is set to 3% of the total entries, and Nswap = Nsig.
  RooRealVar nss("nss", "number of signal entries", 0.03 * num_max, 0 * num_max,
                 0.2 * num_max);
  RooRealVar nsb("nsb", "number of swapped entries", 0.06 * num_max,
                 0.001 * num_max, 0.3 * num_max);
  RooRealVar nbs("nbs", "number of swapped entries", 0.06 * num_max,
                 0.001 * num_max, 0.3 * num_max);
  RooRealVar nbb("nbb", "number of background entries", 0.9 * num_max,
                 0.5 * num_max, 1.1 * num_max);

  RooAddPdf model("model", "composite pdf",
                  RooArgList(sigsig, sigswp, sigbkg, swpsig, swpswp, swpbkg,
                             bkgsig, bkgswp, bkgbkg),
                  RooArgList(nss, nss, nsb, nss, nss, nsb, nbs, nbs, nbb));

  // Unbinned ML fit
  // Extended to get Nsig directly
  r_ml_wgt_corr = model.fitTo(ds, Save(), Extended(kTRUE), NumCPU(8));
  r_ml_wgt_corr->Print();
  projectionPlot(model, r_ml_wgt_corr, ds, 0, 0, nss, outputname,
                 collisionsyst, vtex);

  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.18);
  auto canv2d = new TCanvas("Canvas2D", "Canvas", 800, 600);
  int nbins = 27;
  double mass_low = 1.73;
  double mass_upp = 2.0;
  TString ddname = Form("dd%d", iBin);
  TH2D *dd = new TH2D(ddname, "mass; m_{D}; m_{#bar{D}}", nbins, mass_low,
                      mass_upp, nbins, mass_low, mass_upp);
  sigtree->Draw("m2:m1 >> " + ddname, Form("iPhi == %d", iBin));
  dd->Draw("colorz");
  TString dir = outputname(0, outputname.Length() - 3);
  canv2d->SaveAs(Form("%s" + ddname + ".pdf", dir.Data()));

  yield = std::accumulate(
      yields.cbegin(), yields.cend(), 0, [](auto lhs, const auto &rhs) {
        return std::accumulate(rhs.cbegin(), rhs.cend(), lhs);
      });

  yieldErr = 0;
  for (auto verr : yieldErrs) {
    for (auto err : verr) {
      yieldErr += err * err;
    }
    }
    yieldErr = std::sqrt(yieldErr);

    S = yield;
    B = std::accumulate(notsigs.cbegin(), notsigs.cend(), 0,
                        [](auto lhs, const auto &rhs) {
                          return std::accumulate(rhs.cbegin(), rhs.cend(), lhs);
                        });
    std::cout << "++++++++ failed pt bins for dphi " << iBin << "++++++++++"
              << "\n";
    for (auto fit : failed) {
      std::cout << "pt1 " << fit[0] << ", pt2 " << fit[1] << "\n";
    }
    std::cout << std::flush;
    return;
  }

void xjjroot::fit2d::projectionPlot(RooAddPdf model, RooFitResult *result,
                                    RooDataSet dataset, int xpt, int ypt,
                                    RooRealVar nss, TString outputname,
                                    TString collisionsyst,
                                    const std::vector<TString> &vtex) {
  // Cut ranges for N - 1 plots
  m1.setRange("signal_box_m1", min_fit_dzero, max_hist_dzero);
  m2.setRange("signal_box_m1", mass_dzero_signal_l, mass_dzero_signal_h);

  m1.setRange("signal_box_m2", mass_dzero_signal_l, mass_dzero_signal_h);
  m2.setRange("signal_box_m2", min_fit_dzero, max_hist_dzero);

  std::vector<RooRealVar> dim = {m1, m2};
  std::vector<TString> dzerotex = {"D^{0}", "#bar{D^{#lower[0.2]{0}}}"};

  auto canv = new TCanvas("ptcanvas", "Canvas", 800, 600);
  canv->SetFillColor(color::snow2);

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

  // Plot data and fitting for D and Dbar
  for (auto m : dim) {
    TString name(m.GetName());
    RooPlot *plot = m.frame(Title(name + "_plot"), Bins(60));
    dataset.plotOn(plot, Name("pdat" + name), CutRange("signal_box_" + name));
    model.plotOn(plot, ProjectionRange("signal_box_" + name),
                 VisualizeError(*result));
    dataset.plotOn(plot, Name("pdat" + name), CutRange("signal_box_" + name),
                   DataError(RooAbsData::SumW2));
    model.plotOn(plot, Name("pfit" + name),
                 ProjectionRange("signal_box_" + name));
    model.plotOn(plot, Name("psig3" + name),
                 ProjectionRange("signal_box_" + name),
                 // Components(compstr(name, "sig", {"sig", "swp", "bkg"})),
                 Components(compstr(name, "sig", {"bkg"})), DrawOption("LF"),
                 FillStyle(3002), FillColor(color::aurora4), LineStyle(2),
                 LineColor(color::aurora4));
    model.plotOn(plot, Name("psig2" + name),
                 ProjectionRange("signal_box_" + name),
                 Components(compstr(name, "sig", {"swp"})), DrawOption("LF"),
                 FillStyle(3002), FillColor(color::aurora2), LineStyle(2),
                 LineColor(color::aurora2));
    model.plotOn(plot, Name("psig" + name),
                 ProjectionRange("signal_box_" + name), Components("sigsig"),
                 DrawOption("LF"), FillStyle(3002), FillColor(color::aurora0),
                 LineStyle(2), LineColor(color::aurora0));
    model.plotOn(plot, Name("pbkg3" + name),
                 ProjectionRange("signal_box_" + name),
                 Components(compstr(name, "bkg", {"bkg"})), LineStyle(2),
                 LineColor(color::frost3));
    model.plotOn(plot, Name("pbkg2" + name),
                 ProjectionRange("signal_box_" + name),
                 Components(compstr(name, "bkg", {"swp"})), LineStyle(2),
                 LineColor(color::frost2));
    model.plotOn(plot, Name("pbkg" + name),
                 ProjectionRange("signal_box_" + name),
                 Components(compstr(name, "bkg", {"sig"})), LineStyle(2),
                 LineColor(color::frost0));
    model.plotOn(plot, Name("pswp3" + name),
                 ProjectionRange("signal_box_" + name),
                 Components(compstr(name, "swp", {"bkg"})), DrawOption("LF"),
                 FillStyle(3005), FillColor(color::night0), LineStyle(1),
                 LineColor(color::night1));
    model.plotOn(plot, Name("pswp2" + name),
                 ProjectionRange("signal_box_" + name),
                 Components(compstr(name, "swp", {"swp"})), DrawOption("LF"),
                 FillStyle(3005), FillColor(color::night1), LineStyle(1),
                 LineColor(color::night2));
    model.plotOn(plot, Name("pswp" + name),
                 ProjectionRange("signal_box_" + name),
                 Components(compstr(name, "swp", {"sig"})), DrawOption("LF"),
                 FillStyle(3005), FillColor(color::night2), LineStyle(1),
                 LineColor(color::night3));
    plot->Draw();

    TLegend *leg = new TLegend(0.69, 0.18, 0.89, 0.78, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);

    leg->AddEntry("pdat" + name, "Data", "pl");
    leg->AddEntry("pfit" + name, "Fit", "l");

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

    leg->AddEntry("psig3" + name, make_legend(name, "sig", "bkg"), "f");
    leg->AddEntry("psig2" + name, make_legend(name, "sig", "swap"), "f");
    leg->AddEntry("psig" + name, make_legend(name, "sig", "sig"), "f");
    leg->AddEntry("pswp3" + name, make_legend(name, "swap", "bkg"), "f");
    leg->AddEntry("pswp2" + name, make_legend(name, "swap", "swap"), "f");
    leg->AddEntry("pswp" + name, make_legend(name, "swp", "sig"), "f");
    leg->AddEntry("pbkg3" + name, make_legend(name, "bkg", "bkg"), "l");
    leg->AddEntry("pbkg2" + name, make_legend(name, "bkg", "swap"), "l");
    leg->AddEntry("pbkg" + name, make_legend(name, "bkg", "sig"), "l");

    leg->Draw("same");

    drawCMS(collisionsyst);

    TString other = (name == m1.GetName()) ? "D" : "#bar{D}";
    Float_t texxpos = 0.22, texypos = 0.55, texdypos = 0.053;
    if (!vtex.empty()) {
      texypos += texlinespc;
      for (std::vector<TString>::const_iterator it = vtex.begin();
           it != vtex.end(); it++)
        drawtex(texxpos, texypos = (texypos - texdypos - texlinespc), *it);
    }
    drawtex(texxpos, texypos = (texypos - texdypos),
            Form("(%.3f < m_{%s} < %.3f)", mass_dzero_signal_l, other.Data(),
                 mass_dzero_signal_h));
    if (fdrawyield) {
      drawtex(texxpos, texypos = (texypos - texdypos),
              Form("N = %.0f #pm %.0f", nss.getValV(), nss.getError()));
    }
    // Add pt range onto the plots
    if (xpt > 0 && ypt > 0) {
      drawtex(0.75, 0.85, Form("%.1f < p_{T, D} < %.1f", ptbins[xpt], ptbins[xpt + 1]));
      drawtex(0.75, 0.78, Form("%.1f < p_{T, #bar{D}} < %.1f", ptbins[ypt], ptbins[ypt + 1]));
      canv->SaveAs(outputname + Form("_x%.1f_y%.1f_", ptbins[xpt], ptbins[ypt]) + name + ".pdf");
    } else {
      drawtex(0.75, 0.85, "p_{T} inclusive");
      canv->SaveAs(outputname + "_inclusive_" + name + ".pdf");
    }

  }
}

/*
  Do an unbinned extended simultaneous ML fit on D0 mass in a certain dphi bin
  Calculate the signal yield and error
  @param sigtree    TTree*  tree containing data
  @param swaptree   TTree*  tree containing D0 mass and weight from swapped K pi
  @param isSig   int     1:signal 0:sideband
  @param iBin    int     ID of delta phi bin
*/
void xjjroot::fit2d::simfit(
    TTree *sigtree, TTree *swaptree, const int iBin,
    TString collisionsyst /*=""*/, TString outputname /*="cmass"*/,
    const std::vector<TString> &vtex /*=std::vector<TString>()*/) {
  reset();
  init();

  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);

  TString fitoption = ffitverbose ? "L m" : "L m q";
  setgstyle();

  const Double_t min_fit_dzero = 1.73;
  // Declare an observable for D0 mass
  RooRealVar m1("m1", "m_{#pi K}^{1} / GeV", min_fit_dzero, max_hist_dzero);
  RooRealVar m2("m2", "m_{#pi K}^{2} / GeV", min_fit_dzero, max_hist_dzero);

  // weight of the D0 mass
  RooRealVar w("weight", "weight", min_weight, max_weight);

  // set range for D0 mass
  m1.setRange("D0range", min_fit_dzero, max_hist_dzero);
  m2.setRange("D0range", min_fit_dzero, max_hist_dzero);

  RooRealVar iPhi("iPhi", "dphi bin ID", 0, 20);
  RooRealVar werr("weightErr", "error of weight", min_weight, max_weight);
  RooRealVar isSwap1("isSwap1", "whether m1 is swapped", 0, 1);
  RooRealVar isSwap2("isSwap2", "whether m2 is swapped", 0, 1);
  RooRealVar pt1("pt1", "pt of D", 0, 60);
  RooRealVar pt2("pt2", "pt of Dbar", 0, 60);
  RooRealVar pT1("pT1", "pt of D", 0, 60);
  RooRealVar pT2("pT2", "pt of Dbar", 0, 60);

  // Define category to distinguish trigger and associate samples events
  RooCategory sample("sample", "sample");
  sample.defineType("x");
  sample.defineType("y");

  const std::vector<double> ptbins = {1.0,  2.0,  2.5,  3.0,  3.5,  4.0, 4.5,
                                      5.0,  5.5,  6.0,  6.5,  7.0,  8.0, 10.0,
                                      12.5, 15.0, 20.0, 30.0, 40.0, 60.0};

  // Construct unbinned dataset importing tree branches mass matching
  // between branches and RooRealVars
  // Select the dphi bin in question
  RooDataSet ds("ds", "ds", RooArgSet(m1, m2, iPhi, pT1, pT2),
                Cut(TString::Format("iPhi == %i", iBin)), Import(*sigtree));

  RooDataSet dsm1 = *(RooDataSet *)ds.reduce(RooArgSet(m1, pT1));
  RooDataSet dsm2 = *(RooDataSet *)ds.reduce(RooArgSet(m2, pT2));
  // RooDataSet ds_asso("dsasso", "dsasso", RooArgSet(m1), Import(y));
  std::vector<RooDataSet *> dsetsm1(ptbins.size() - 1);
  std::vector<RooDataSet *> dsetsm2(ptbins.size() - 1);
  for (unsigned i = 1; i < ptbins.size() - 1; ++i) {
    dsetsm1[i] = (RooDataSet *)dsm1.reduce(
        RooArgSet(m1), Form("pT1 > %f && pT1 < %f", ptbins[i], ptbins[i + 1]));
    dsetsm2[i] = (RooDataSet *)dsm2.reduce(
        RooArgSet(m2), Form("pT2 > %f && pT2 < %f", ptbins[i], ptbins[i + 1]));
  }

  RooDataSet swapds("swapds", "swapped dataset",
                    RooArgSet(m1, isSwap1, pt1, m2, isSwap2, pt2),
                    Import(*swaptree));
  RooDataSet swapdsm1 = *(RooDataSet *)swapds.reduce(RooArgSet(m1), "isSwap1");
  RooDataSet swapdsm2 = *(RooDataSet *)swapds.reduce(RooArgSet(m2), "isSwap2");

  RooDataSet sigdsm1set =
      *(RooDataSet *)swapds.reduce(RooArgSet(m1, pt1), "!isSwap1");
  RooDataSet sigdsm2set =
      *(RooDataSet *)swapds.reduce(RooArgSet(m2, pt2), "!isSwap2");
  std::vector<RooDataSet *> sigdsm1(ptbins.size() - 1);
  std::vector<RooDataSet *> sigdsm2(ptbins.size() - 1);
  for (unsigned i = 1; i < ptbins.size() - 1; ++i) {
    sigdsm1[i] = (RooDataSet *)sigdsm1set.reduce(
        RooArgSet(m1), Form("pt1 > %f && pt1 < %f", ptbins[i], ptbins[i + 1]));
    sigdsm2[i] = (RooDataSet *)sigdsm2set.reduce(
        RooArgSet(m2), Form("pt2 > %f && pt2 < %f", ptbins[i], ptbins[i + 1]));
  }

  // Construct combined dataset in (x,sample)
  std::vector<RooDataSet *> dsets2d(ptbins.size() - 1);
  for (unsigned i = 1; i < ptbins.size() - 1; ++i) {
    dsets2d[i] = new RooDataSet(
        Form("combData%i", i), "combined data", RooArgSet(m1, m2),
        Index(sample), Import("x", *dsetsm1[i]), Import("y", *dsetsm2[i]));
  }

  // Double gaussian as signal p.d.f.
  // Don't fix the mean at nominal D0 mass
  RooRealVar mean("mean", "mean of gaussians", mass_dzero, mass_dzero - 0.01,
                  mass_dzero + 0.01);
  RooRealVar sigma1x("sigma1x", "width of gaussians", 0.04, 0.02, 0.1);
  RooRealVar sigma2x("sigma2x", "width of gaussians", 0.01, 0.001, 0.02);
  // RooRealVar sigma3x("sigma3x", "width of gaussians", 0.005, 0.001, 0.01);
  RooGaussian sig1x("sig1x", "Signal component 1", m1, mean, sigma1x);
  RooGaussian sig2x("sig2x", "Signal component 2", m1, mean, sigma2x);
  // RooGaussian sig3x("sig3x", "Signal component 3", m1, mean, sigma3x);
  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sigfracx("sigfracx", "fraction of component 1 in signal", 0.1,
                      0.001, 0.2);
  // RooRealVar sigfracx2("sigfracx2", "fraction of component 2 in signal",
  // 0.5, 0.01, 0.8);
  RooAddPdf sigx("sigx", "Signal", RooArgList(sig1x, sig2x), sigfracx);
  // RooAddPdf sigx("sigx", "Signal", RooArgList(sig1x, sig2x, sig3x),
  // RooArgList(sigfracx, sigfracx2));
  // RooGaussian sigx(sig2x, "sigx");

  RooGaussian sig1y("sig1y", "Signal component 1", m2, mean, sigma1x);
  RooGaussian sig2y("sig2y", "Signal component 2", m2, mean, sigma2x);
  // RooGaussian sig3y("sig3y", "Signal component 3", m2, mean, sigma3x);
  // RooAddPdf sigy("sigy", "Signal", RooArgList(sig1y, sig2y, sig3y),
  //                RooArgList(sigfracx, sigfracx2));
  RooAddPdf sigy("sigy", "Signal", RooArgList(sig1y, sig2y), sigfracx);

  // Build 3rd-order Chebychev polynomial p.d.f. as combinatorial background
  RooRealVar a1("a1", "a1", -0.2, -.8, .8);
  // RooRealVar a2("a2", "a2", 0.02, -.1, .1);
  // RooRealVar a3("a3", "a3", 0.0, -.1, .1);
  // RooChebychev bkgx("bkgx", "Background", m1, RooArgSet(a1, a2, a3));
  // RooChebychev bkgx("bkgx", "Background", m1, RooArgSet(a1, a3));
  RooChebychev bkgx("bkgx", "Background", m1, RooArgSet(a1));

  RooRealVar b1("b1", "b1", -0.2, -.8, .8);
  // RooRealVar b2("b2", "b2", 0.02, -.1, .1);
  // RooRealVar b3("b3", "b3", 0.0, -.1, .1);
  // RooChebychev bkgy("bkgy", "Background", m2, RooArgSet(b1, b2, b3));
  RooChebychev bkgy("bkgy", "Background", m2, RooArgSet(b1));

  // Gaussian as swapped K pi mass
  // RooRealVar swmean("swmean", "mean of gaussians", mass_dzero, mass_dzero -
  // 0.05, mass_dzero + 0.05);
  RooRealVar sigmas1x("sigmas1x", "width of swapped mass", 0.1, 0.05, 2);
  // RooRealVar sigmas1y("sigmas1y", "width of swapped mass", 0.1, 0.05, 2);

  RooGaussian swpx("swpx", "swapped k pi mass", m1, mean, sigmas1x);
  RooGaussian swpy("swpy", "swapped k pi mass", m2, mean, sigmas1x);

  const int nPrefitBins = 60;
  swpy.fitTo(swapdsm2, Save());
  RooPlot *frame3 = m2.frame(Title("swappedm2_plot"), Bins(nPrefitBins));
  swapdsm2.plotOn(frame3);
  swpy.plotOn(frame3);
  swpy.paramOn(frame3);

  // Plot the result
  auto canv = new TCanvas("Canvas", "Canvas", 800, 600);
  canv->SetFillColor(color::snow2);
  frame3->Draw();
  canv->SaveAs(outputname + "_swapfitm2" + ".pdf");

  swpx.fitTo(swapdsm1, Save());
  RooPlot *frame2 = m1.frame(Title("swappedm1_plot"), Bins(nPrefitBins));
  swapdsm1.plotOn(frame2);
  swpx.plotOn(frame2);
  swpx.paramOn(frame2);
  sigmas1x.setConstant(true);
  frame2->Draw();
  canv->SaveAs(outputname + "_swapfitm1" + ".pdf");

  sigy.fitTo(*sigdsm2[1], Save());
  RooPlot *framesig2 = m2.frame(Title("sigm2_plot"), Bins(nPrefitBins));
  sigdsm2[1]->plotOn(framesig2);
  sigy.plotOn(framesig2);
  sigy.paramOn(framesig2, Layout(0.6), Format("NEU", AutoPrecision(1)),
               ShowConstants(false));
  framesig2->Draw();
  canv->SaveAs(outputname + "_sigfitm2" + ".pdf");

  RooFitResult *r_ml_wgt_corr;

  RooPlot *framesig1 = m1.frame(Title("sigm1_plot"), Bins(nPrefitBins));
  m1.setRange("3sigma", mass_dzero - 0.04, mass_dzero + 0.04);
  for (unsigned ipt = 0; ipt < ptbins.size() - 1; ++ipt) {
    sigx.fitTo(*sigdsm1[ipt], Range("3sigma"), Save());
    framesig1->remove();
    framesig1->remove();
    framesig1->remove();
    framesig1->Clear();
    canv->Clear();
    canv->SetLogy(0);
    sigdsm1[ipt]->plotOn(framesig1);
    sigx.plotOn(framesig1);
    sigx.paramOn(framesig1, Layout(0.6), Format("NEU", AutoPrecision(1)),
                 ShowConstants(false));

    TLatex *pttex = new TLatex(
        0.4, 0.7, Form("%.1f < p_{T} < %.1f", ptbins[ipt], ptbins[ipt + 1]));
    pttex->SetNDC();
    pttex->SetTextAlign(32);
    pttex->SetTextSize(0.04);
    pttex->SetTextFont(42);
    pttex->Draw();
    canv->SaveAs(outputname +
                 Form("_sigfitm1_%.1f_%.1f", ptbins[ipt], ptbins[ipt + 1]) +
                 ".pdf");
    canv->SetLogy();
    canv->SaveAs(outputname +
                 Form("_sigfitm1_%.1f_%.1f_log", ptbins[ipt], ptbins[ipt + 1]) +
                 ".pdf");
    canv->SetLogy(0);
    // sigma1x.setConstant(true);
    // sigma2x.setConstant(true);
    // sigma3x.setConstant(true);
    // sigfracx.setConstant(true);
    // sigfracx2.setConstant(true);

    // mean.setConstant(true);
    // swmean.setConstant(true);
    // RooRealVar rsigswap("rsigswap", "coeff of signal", 1. / 3.);
    // RooAddPdf sig("sig", "Signal", RooArgList(s

    const long num_max = dsets2d[ipt]->sumEntries("", "D0range");

    // Initially, Nsig is set to 3% of the total
    RooRealVar nss("nss", "number of signal-signal entries", 0.02 * num_max,
                   0.001 * num_max, 0.1 * num_max);
    RooRealVar nbb("nbb", "number of bg-bg entries", 0.92 * num_max,
                   0.700 * num_max, 1.1 * num_max);
    RooRealVar nww("nww", "number of swap-swap entries", 0.02 * num_max,
                   0.001 * num_max, 0.1 * num_max);
    RooRealVar nsb("nsb", "number of signal-bg entries",
                   // 0.02 * num_max, 0.001 * num_max, 0.1 * num_max);
                   0.02 * num_max, 0.001 * num_max, 0.1 * num_max);
    RooRealVar nbs("nbs", "number of bg-signal entries", 0.02 * num_max,
                   0.001 * num_max, 0.1 * num_max);
    RooRealVar nsw("nsw", "number of signal-swap entries", 0.01 * num_max,
                   0.001 * num_max, 0.1 * num_max);
    RooRealVar nws("nws", "number of swap-signal entries", 0.01 * num_max,
                   0.001 * num_max, 0.1 * num_max);
    RooRealVar nbw("nbw", "number of bg-swap entries", 0.01 * num_max,
                   0.001 * num_max, 0.1 * num_max);
    RooRealVar nwb("nwb", "number of swap-bg entries", 0.01 * num_max,
                   0.001 * num_max, 0.1 * num_max);
    // RooArgList nums(nss, nbb, nww, nsb, nbs, nsw, nws, nbw, nwb, "nums");
    // RooArgList nums(nss, nbb, nww, nsb, nsb, nsw, nsw, nbw, nbw, "nums");
    // RooArgList nums(nss, nbb, nss, nsb, nsb, nss, nss, nsb, nsb, "nums");
    // RooArgList nums(nss, nbb, nww, "nums");
    RooArgList nums(nss, nbb, nss, "nums");

    // copy the PDF so they can be plotted separately
    RooAddPdf sbx(sigx, "sbx");
    RooAddPdf swx(sigx, "swx");
    RooAddPdf sby(sigy, "sby");
    RooAddPdf swy(sigy, "swy");

    // Model for D0
    RooAddPdf modelx("modelx", "composite pdf", RooArgList(sigx, bkgx, swpx), nums);
    // RooArgList(sigx, bkgx, swpx, sbx, bkgx, swx, swpx, bkgx, swpx), nums);
    // Model for D0bar
    RooAddPdf modely("modely", "composite pdf", RooArgList(sigy, bkgy, swpy), nums);
    // RooArgList(sigy, bkgy, swpy, bkgy, sby, swpy, swy, swpy, bkgy), nums);

    // Construct a simultaneous pdf using category sample as index
    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);

    // Associate model with the physics state and model_ctl with the control
    // state
    simPdf.addPdf(modelx, "x");
    simPdf.addPdf(modely, "y");

    // 1D fit
    modelx.fitTo(dsm1, Save(), Extended(kTRUE), NumCPU(10));
    RooPlot *framex = m1.frame(Title("m1_plot"), Bins(nPrefitBins));
    dsm1.plotOn(framex, DataError(RooAbsData::SumW2));
    modelx.plotOn(framex);
    framex->Draw();
    canv->SaveAs(outputname + "_fitm1.pdf");

    modely.fitTo(dsm2, Save(), Extended(kTRUE), NumCPU(10));
    RooPlot *framey = m2.frame(Title("m2_plot"), Bins(nPrefitBins));
    dsm2.plotOn(framey, DataError(RooAbsData::SumW2));
    modely.plotOn(framey);
    framey->Draw();
    canv->SaveAs(outputname + "_fitm2.pdf");

    RooDataSet *combData = dsets2d[ipt];
    // U n b i n n e d   M L   f i t
    // Extended to get Nsig directly
    r_ml_wgt_corr =
        simPdf.fitTo(*combData, Save(), Extended(kTRUE), NumCPU(10));

    r_ml_wgt_corr->Print();
    // Calculate yields and error
    yield = nss.getValV();
    yieldErr = nss.getError();
    // yieldErr = nsig.getAsymErrorHi();

    // std::vector<RooRealVar> nlist = {nss, nbb, nww, nsb, nbs,
    //                                  nsw, nws, nbw, nwb};
    std::vector<RooRealVar> nlist = {nss, nbb};

    int total = 0;
    for (auto &i : nlist) {
      i.printValue(std::cout);
      std::cout << ", " << i.getValV() / num_max << "\n";
      total += i.getValV();
    }
    std::cout << "total: " << num_max << ", sum: " << total << "\n";

    std::vector<TString> dim = {"x", "y"};
    std::vector<TString> dzerotex = {"D^{0}", "#bar{D^{#lower[0.2]{0}}}"};
    std::vector<RooRealVar> vars = {m1, m2};

    // // lazy lambda to change x in strings to y
    // auto compstr = [](TString str, TString com) {
    //   return str.ReplaceAll("x", com);
    // };

    for (unsigned i = 0; i < dim.size(); ++i) {
      auto m = dim[i];
      auto var = vars[i];
      TString name = m;
      RooPlot *plot = var.frame(Title(name + "_plot"), Bins(60));
      combData->plotOn(plot, Name("pdat" + name),
                       Cut("sample==sample::" + name));
      // simPdf.plotOn(plot, Slice(sample, m), ProjWData(sample, *combData),
      //               LineColor(color::frost2), VisualizeError(*r_ml_wgt_corr));
      combData->plotOn(plot, Name("pdat" + name),
                       Cut("sample==sample::" + name),
                       DataError(RooAbsData::SumW2));
      // simPdf.plotOn(plot, Name("pfit" + name), LineColor(color::frost2),
      //               Slice(sample, m), ProjWData(sample, *combData));
      // simPdf.plotOn(plot, Name("psw" + name), Slice(sample, m),
      //               ProjWData(sample, *combData),
      //               Components(compstr("sigx,sbx,swx", m)), DrawOption("LF"),
      //               FillStyle(1001), FillColor(color::aurora4),
      //               LineColor(color::aurora4));
      // simPdf.plotOn(plot, Name("psb" + name), Slice(sample, m),
      //               ProjWData(sample, *combData),
      //               Components(compstr("sigx,sbx", m)), DrawOption("LF"),
      //               FillStyle(1001), FillColor(color::aurora2),
      //               LineColor(color::aurora2));
      // simPdf.plotOn(plot, Name("psig" + name), Slice(sample, m),
      //               ProjWData(sample, *combData),
      //               Components(compstr("sigx", m)), DrawOption("LF"),
      //               FillStyle(1001), FillColor(color::aurora0),
      //               LineColor(color::aurora0));
      // simPdf.plotOn(plot, Name("pbkg" + name), Slice(sample, m),
      //               ProjWData(sample, *combData),
      //               Components(compstr("bkgx", m)), LineStyle(2),
      //               LineColor(color::frost2));
      // simPdf.plotOn(plot, Name("pswp" + name), Slice(sample, m),
      //               ProjWData(sample, *combData),
      //               Components(compstr("swpx", m)), LineStyle(1),
      //               DrawOption("LF"), FillStyle(3005), FillColor(color::frost0),
      //               LineColor(color::frost0));
      plot->Draw();

      float x1 = 0.72, x2 = 0.91, y1 = 0.7, y2 = 0.89;
      TLegend *legtop = new TLegend(x1, y1, x2, y2, NULL, "brNDC");
      legtop->SetBorderSize(0);

      legtop->SetTextFont(42);
      legtop->SetTextSize(0.04);
      legtop->SetFillStyle(0);
      legtop->AddEntry("pdat" + name, "Data", "pl");
      legtop->AddEntry("pfit" + name, "Fit", "l");
      legtop->AddEntry("pbkg" + name, "Combinatorial", "l");

      TLegend *leg = new TLegend(x1, 0.20, x2, 0.50, NULL, "brNDC");
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.04);

      TString thisD = (name == dim[0]) ? dzerotex[0] : dzerotex[1];
      TString otherD = (name == dim[0]) ? dzerotex[1] : dzerotex[0];
      leg->AddEntry("psw" + name, thisD + " Sig - " + otherD + " Swap", "f");
      leg->AddEntry("psb" + name, thisD + " Sig - " + otherD + " Bkg", "f");
      leg->AddEntry("psig" + name, thisD + " Sig - " + otherD + " Sig", "f");
      leg->AddEntry("pswp" + name, "K-#pi swapped", "f");

      // legtop->Draw("same");
      // leg->Draw("same");

      drawCMS(collisionsyst);

      Float_t texxpos = 0.22, texypos = 0.55, texdypos = 0.058;
      if (!vtex.empty()) {
        texypos += texlinespc;
        for (std::vector<TString>::const_iterator it = vtex.begin();
             it != vtex.end(); it++)
          drawtex(texxpos, texypos = (texypos - texdypos - texlinespc), *it);
      }
      // drawtex(texxpos, texypos = (texypos - texdypos),
      //         Form("(%.3f < m_{%s} < %.3f)", mass_dzero_signal_l,
      //         other.Data(),
      //              mass_dzero_signal_h));
      if (fdrawyield) {
        drawtex(texxpos, texypos = (texypos - texdypos),
                Form("N = %.0f #pm %.0f", yield, yieldErr), 0.04, 12,
                color::aurora0);
      }

      if (fsaveplot) {
        canv->SaveAs(Form("%s_" + name + "%.1f-%.1f.pdf", outputname.Data(),
                          ptbins[ipt], ptbins[ipt + 1]));
      }
    }
  }

  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.18);
  auto canv2d = new TCanvas("Canvas2D", "Canvas", 800, 600);
  int nbins = 27;
  double mass_low = 1.73;
  double mass_upp = 2.0;
  TString ddname = Form("dd%d", iBin);
  TH2D *dd = new TH2D(ddname, "mass; m_{D}; m_{#bar{D}}", nbins, mass_low,
                      mass_upp, nbins, mass_low, mass_upp);
  sigtree->Draw("m2:m1 >> " + ddname, Form("iPhi == %d", iBin));
  dd->Draw("colorz");
  TString dir = outputname(0, outputname.Length() - 3);
  canv2d->SaveAs(Form("%s" + ddname + ".pdf", dir.Data()));

  return;
}
void xjjroot::fit2d::reset()
{
  clearvar();
}

void xjjroot::fit2d::init()
{
  clearvar();
  // read bins to skip from the config file
  std::ifstream fin("skip_to_bins.txt");
  double idphi;
  double xpt;
  double ypt;
  while (fin >> idphi >> xpt >> ypt) {
    skiptobins.push_back({idphi, xpt, ypt});
  }
}

void xjjroot::fit2d::clearvar()
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

void xjjroot::fit2d::calvar()
{
  mass_dzero_signal_l = mass_dzero - d_mass_signal;
  mass_dzero_signal_h = mass_dzero + d_mass_signal;
  mass_dzero_sideband_l_p = mass_dzero + d_mass_sideband_l;
  mass_dzero_sideband_h_p = mass_dzero + d_mass_sideband_h;
  mass_dzero_sideband_l_n = mass_dzero - d_mass_sideband_l;
  mass_dzero_sideband_h_n = mass_dzero - d_mass_sideband_h;

  if(!fparamfuns) return;
  Sig = S/TMath::Sqrt(S+B);
  Chi2 = 2.*r->MinFcnValue();
  Chi2Prob = TMath::Prob(Chi2,NDF);
}


void xjjroot::fit2d::sethist(TH1* h) const
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

void xjjroot::fit2d::drawCMS(TString collision, TString snn/*="5.02"*/) const
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

void xjjroot::fit2d::drawtex(Double_t x, Double_t y, const char* text, Float_t tsize/*=0.04*/, Short_t align/*=12*/, int color) const
{
  TLatex* tex = new TLatex(x, y, text);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
  tex->SetTextColor(color);
  tex->Draw();
}

void xjjroot::fit2d::drawline(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Color_t lcolor/*=kBlack*/, Style_t lstyle/*=1*/, Width_t lwidth/*=2*/) const
{
  TLine* l = new TLine(x1, y1, x2, y2);
  l->SetLineColor(lcolor);
  l->SetLineStyle(lstyle);
  l->SetLineWidth(lwidth);
  l->Draw();
}

void xjjroot::fit2d::setgstyle() const
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
