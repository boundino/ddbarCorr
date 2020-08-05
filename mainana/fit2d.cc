#include "fit2d.hpp"
#include "RooCategory.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "RooSimultaneous.h"
#include "TFile.h"

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
  @param isSig   int     1:signal 0:sideband
  @param iBin    int     ID of delta phi bin
*/
RooFitResult* xjjroot::fit2d::fit(
    TTree* sigtree, TTree* swaptree, const int iBin,
    TString collisionsyst /*=""*/, TString outputname /*="cmass"*/,
    const std::vector<TString> &vtex /*=std::vector<TString>()*/) {
  reset();
  init();

  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);

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
  // Construct unbinned dataset importing tree branches mass matching
  // between branches and RooRealVars
  // Select the dphi bin in question
  RooDataSet ds("ds", "ds", RooArgSet(m1, m2, iPhi),
                Cut(TString::Format("iPhi == %i", iBin)), Import(*sigtree));

  RooDataSet swapds("swapds", "swapped dataset", RooArgSet(m1, m2), Import(*swaptree));
  // Create a binned dataset
  // RooDataHist *swapdh = swapds.binnedClone();

  // Represent data in dh as pdf in x
  // RooHistPdf swapped("histpdfsw", "histpdfsw", x, *swapdh, 2);

  // Double gaussian as signal p.d.f.
  // Fix the mean at D0 mass
  RooRealVar mean("mean", "mean of gaussians", mass_dzero, mass_dzero - 0.01, mass_dzero + 0.01);
  mean.setConstant(true);
  RooRealVar sigma1x("sigma1x", "width of gaussians", 0.02, 0.005, 0.04);
  RooRealVar sigma2x("sigma2x", "width of gaussians", 0.005, 0.001, 0.02);
  RooGaussian sig1x("sig1x", "Signal component 1", m1, mean, sigma1x);
  RooGaussian sig2x("sig2x", "Signal component 2", m1, mean, sigma2x);

  RooRealVar sigma1y("sigma1y", "width of gaussians", 0.02, 0.005, 0.04);
  RooRealVar sigma2y("sigma2y", "width of gaussians", 0.005, 0.001, 0.02);
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

  // 2-Gaussian as swapped K pi mass
  RooRealVar sigmas1x("sigmas1x", "width of swapped mass", 0.1, 0.05, 2);
  RooRealVar sigmas2x("sigmas2x", "width of swapped mass", 0.02, 0.01, 0.1);
  RooGaussian swapped1x("swapped1x", "swapped k pi mass", m1, mean, sigmas1x);
  RooGaussian swapped2x("swapped2x", "swapped k pi mass", m1, mean, sigmas2x);

  RooRealVar sigmas1y("sigmas1y", "width of swapped mass", 0.1, 0.05, 2);
  RooRealVar sigmas2y("sigmas2y", "width of swapped mass", 0.02, 0.01, 0.1);
  RooGaussian swapped1y("swapped1y", "swapped k pi mass", m2, mean, sigmas1y);
  RooGaussian swapped2y("swapped2y", "swapped k pi mass", m2, mean, sigmas2y);

  RooRealVar swpfracx("swpfracx", "fraction of component 1 in swapped", 0.5, 0.001, 0.999);
  RooRealVar swpfracy("swpfracy", "fraction of component 1 in swapped", 0.5, 0.001, 0.999);
  RooAddPdf swappedx("swappedx", "swp", RooArgList(swapped1x, swapped2x), swpfracx);
  RooAddPdf swappedy("swappedy", "swp", RooArgList(swapped1y, swapped2y), swpfracy);
  RooProdPdf swapped("swapped", "background PDF", RooArgSet(swappedx, swappedy));

  swapped.fitTo(swapds, Save());
  RooPlot *frame2 = m1.frame(Title("swapped_plot"), Bins(30));
  swapds.plotOn(frame2);
  swapped.plotOn(frame2);

  sigmas1x.setConstant(true);
  sigmas2x.setConstant(true);
  swpfracx.setConstant(true);
  sigmas1y.setConstant(true);
  sigmas2y.setConstant(true);
  swpfracy.setConstant(true);

  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sigfracx("sigfracx", "fraction of component 1 in signal", 0.5, 0.001, 0.999);
  RooAddPdf sigx("sigx", "Signal", RooArgList(sig1x, sig2x), sigfracx);
  RooRealVar sigfracy("sigfracy", "fraction of component 1 in signal", 0.5, 0.001, 0.999);
  RooAddPdf sigy("sigy", "Signal", RooArgList(sig1y, sig2y), sigfracy);

  RooProdPdf sigonly("sigonly", "signal PDF", RooArgSet(sigx, sigy));
  RooRealVar rsigswap("rsigswap", "coeff of signal", 1. / 3.);
  RooAddPdf sig("sig", "Signal", RooArgList(sigonly, swapped), rsigswap);
  RooProdPdf bkg("bkg", "background PDF", RooArgSet(bkgx, bkgy));

  const long num_max = ds.sumEntries("", "D0range");

  // Initially, Nsig is set to 3% of the total entries, and Nswap = Nsig.
  RooRealVar nsig("nsig", "number of signal entries",
                  0.03 * num_max, 0.001 * num_max, 0.1 * num_max);
  RooRealVar nbkg("nbkg", "number of background entries",
                  0.95 * num_max, 0.6 * num_max, 1.1 * num_max);
  RooRealVar nswap("nswap", "number of swapped entries", 0.06 * num_max,
                  0.01 * num_max, 0.2 * num_max);

  // // Associated nsig/nbkg as expected number of events with sig/bkg
  // RooExtendPdf esig("esig", "extended signal p.d.f", sig, nsig);
  // RooExtendPdf ebkg("ebkg", "extended background p.d.f", bkg, nbkg);
  // RooExtendPdf eswapped("eswapped", "extended swapped p.d.f", swapped, nsig);
  // RooAddPdf model("model", "composite pdf", RooArgList(esig, ebkg,
  // eswapped));

  // RooAddPdf model("model", "composite pdf", RooArgList(bkg, swapped, sig),
  //                 RooArgList(nbkg, nsig, nsig));
  RooAddPdf model("model", "composite pdf", RooArgList(bkg, swapped, sigonly), RooArgList(nbkg, nswap, nsig));

  // U n b i n n e d   M L   f i t   t o   w e i g h t e d   d a t a
  // Extended to get Nsig directly
  RooFitResult *r_ml_wgt_corr =
      model.fitTo(ds, Save(), Extended(kTRUE), NumCPU(6));

  r_ml_wgt_corr->Print();

  // Calculate yields and error
  yield = nsig.getValV();
  yieldErr = nsig.getError();
  // yieldErr = nsig.getAsymErrorHi();

  // Plot the result
  auto canv = new TCanvas("Canvas", "Canvas", 800, 600);
  frame2->Draw();
  canv->SaveAs(outputname + "_swapfit" +".pdf");

  // Cut ranges for N - 1 plots
  m1.setRange("signal_box_m1", min_fit_dzero, max_hist_dzero);
  m2.setRange("signal_box_m1", mass_dzero_signal_l, mass_dzero_signal_h);

  m1.setRange("signal_box_m2", mass_dzero_signal_l, mass_dzero_signal_h);
  m2.setRange("signal_box_m2", min_fit_dzero, max_hist_dzero);

  std::vector<RooRealVar> dim = {m1, m2};

  for (auto m : dim) {
    TString name(m.GetName());
    RooPlot *plot = m.frame(Title(name + "_plot"), Bins(60));
    ds.plotOn(plot, Name("pdat" + name), CutRange("signal_box_" + name));
    model.plotOn(plot, ProjectionRange("signal_box_" + name),
                 VisualizeError(*r_ml_wgt_corr));
    ds.plotOn(plot, Name("pdat" + name), CutRange("signal_box_" + name),
              DataError(RooAbsData::SumW2));
    model.plotOn(plot, Name("pfit" + name),
                 ProjectionRange("signal_box_" + name));
    model.plotOn(plot, Name("psig" + name),
                 ProjectionRange("signal_box_" + name), Components(sigonly),
                 DrawOption("LF"), FillStyle(3002), FillColor(kOrange - 3),
                 LineStyle(2), LineColor(kOrange - 3));
    model.plotOn(plot, Name("pbkg" + name),
                 ProjectionRange("signal_box_" + name), Components(bkg),
                 LineStyle(2), LineColor(4));
    model.plotOn(plot, Name("pswp" + name),
                 ProjectionRange("signal_box_" + name), Components(swapped),
                 DrawOption("LF"), FillStyle(3005), FillColor(kGreen + 4),
                 LineStyle(1), LineColor(kGreen + 4));
    plot->Draw();

    TLegend *leg = new TLegend(0.69, 0.18, 0.89, 0.48, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);

    leg->AddEntry("pdat" + name, "Data", "pl");
    leg->AddEntry("pfit" + name, "Fit", "l");
    leg->AddEntry("psig" + name, "D^{0}+#bar{D^{#lower[0.2]{0}}} Signal", "f");
    leg->AddEntry("pswp" + name, "K-#pi swapped", "f");
    leg->AddEntry("pbkg" + name, "Combinatorial", "l");

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
              Form("N = %.0f #pm %.0f", yield, yieldErr));
    }

    if (fsaveplot) {
      canv->SaveAs(Form("%s_" + name + ".pdf", outputname.Data()));
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

  return r_ml_wgt_corr;
}

/*
  Do an unbinned extended simultaneous ML fit on D0 mass in a certain dphi bin
  Calculate the signal yield and error
  @param sigtree    TTree*  tree containing data
  @param swaptree   TTree*  tree containing D0 mass and weight from swapped K pi
  @param isSig   int     1:signal 0:sideband
  @param iBin    int     ID of delta phi bin
*/
RooFitResult* xjjroot::fit2d::simfit(
    TTree* sigtree, TTree* swaptree, const int iBin,
    TString collisionsyst /*=""*/, TString outputname /*="cmass"*/,
    const std::vector<TString> &vtex /*=std::vector<TString>()*/) {
  reset();
  init();

  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);

  TString fitoption = ffitverbose ? "L m" : "L m q";
  setgstyle();

  const Double_t min_fit_dzero = 1.73;
  // Declare an observable for D0 mass
  RooRealVar m1("m1", "m_{#pi K}^{1} / GeV", min_fit_dzero, max_hist_dzero);
  // RooRealVar m2("m2", "m_{#pi K}^{2} / GeV", min_fit_dzero, max_hist_dzero);

  // weight of the D0 mass
  RooRealVar w("weight", "weight", min_weight, max_weight);

  // set range for D0 mass
  m1.setRange("D0range", min_fit_dzero, max_hist_dzero);
  // m2.setRange("D0range", min_fit_dzero, max_hist_dzero);

  RooRealVar iPhi("iPhi", "dphi bin ID", 0, 20);
  RooRealVar werr("weightErr", "error of weight", min_weight, max_weight);

  // Define category to distinguish trigger and associate samples events
  RooCategory sample("sample", "sample");
  sample.defineType("trig");
  sample.defineType("asso");

  // tree containing associate D
  TFile ftmp("tmp.root", "recreate");
  TTree asso("ass", "ass");
  TTreeReader reader(sigtree);
  TTreeReaderValue<Float_t> sigmass(reader, "m2");
  TTreeReaderValue<Int_t> sigiphi(reader, "iPhi");
  Float_t smass;
  asso.Branch("m1", &smass, "m1/F");
  while (reader.Next()) {
    if (*sigiphi != iBin) continue;
    smass = *sigmass;
    asso.Fill();
  }

  // Construct unbinned dataset importing tree branches mass matching
  // between branches and RooRealVars
  // Select the dphi bin in question
  RooDataSet ds("ds", "ds", RooArgSet(m1, iPhi),
                Cut(TString::Format("iPhi == %i", iBin)), Import(*sigtree));

  RooDataSet ds_asso("dsasso", "dsasso", RooArgSet(m1), Import(asso));

  RooDataSet swapds("swapds", "swapped dataset", RooArgSet(m1), Import(*swaptree));

  // Construct combined dataset in (x,sample)
  RooDataSet combData("combData", "combined data", m1, Index(sample),
                      Import("trig", ds), Import("asso", ds_asso));

  // Double gaussian as signal p.d.f.
  // Fix the mean at D0 mass
  RooRealVar mean("mean", "mean of gaussians", mass_dzero, mass_dzero - 0.01, mass_dzero + 0.01);
  mean.setConstant(true);
  RooRealVar sigma1x("sigma1x", "width of gaussians", 0.02, 0.005, 0.04);
  RooRealVar sigma2x("sigma2x", "width of gaussians", 0.005, 0.001, 0.02);
  RooGaussian sig1x("sig1x", "Signal component 1", m1, mean, sigma1x);
  RooGaussian sig2x("sig2x", "Signal component 2", m1, mean, sigma2x);
  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sigfracx("sigfracx", "fraction of component 1 in signal", 0.5, 0.001, 0.999);
  RooAddPdf sigx("sigx", "Signal", RooArgList(sig1x, sig2x), sigfracx);

  RooRealVar sigma1y("sigma1y", "width of gaussians", 0.02, 0.005, 0.04);
  RooRealVar sigma2y("sigma2y", "width of gaussians", 0.005, 0.001, 0.02);
  RooGaussian sig1y("sig1y", "Signal component 1", m1, mean, sigma1y);
  RooGaussian sig2y("sig2y", "Signal component 2", m1, mean, sigma2y);
  // Sum the signal components into a composite signal p.d.f.
  RooRealVar sigfracy("sigfracy", "fraction of component 1 in signal", 0.5, 0.001, 0.999);
  RooAddPdf sigy("sigy", "Signal", RooArgList(sig1y, sig2y), sigfracy);

  // Build 3rd-order Chebychev polynomial p.d.f. as combinatorial background
  RooRealVar a1("a1", "a1", -0.2, -.8, .8);
  RooRealVar a2("a2", "a2", 0.02, -.1, .1);
  RooRealVar a3("a3", "a3", 0.0, -.1, .1);
  RooChebychev bkgx("bkgx", "Background", m1, RooArgSet(a1, a2, a3));

  RooRealVar b1("b1", "b1", -0.2, -.8, .8);
  RooRealVar b2("b2", "b2", 0.02, -.1, .1);
  RooRealVar b3("b3", "b3", 0.0, -.1, .1);
  RooChebychev bkgy("bkgy", "Background", m1, RooArgSet(b1, b2, b3));

  // 2-Gaussian as swapped K pi mass
  RooRealVar sigmas1x("sigmas1x", "width of swapped mass", 0.1, 0.05, 2);
  RooRealVar sigmas2x("sigmas2x", "width of swapped mass", 0.02, 0.01, 0.1);
  RooGaussian swapped1x("swapped1x", "swapped k pi mass", m1, mean, sigmas1x);
  RooGaussian swapped2x("swapped2x", "swapped k pi mass", m1, mean, sigmas2x);
  RooRealVar swpfracx("swpfracx", "fraction of component 1 in swapped", 0.5, 0.001, 0.999);
  RooAddPdf swappedx("swappedx", "swp", RooArgList(swapped1x, swapped2x), swpfracx);

  RooRealVar sigmas1y("sigmas1y", "width of swapped mass", 0.1, 0.05, 2);
  RooRealVar sigmas2y("sigmas2y", "width of swapped mass", 0.02, 0.01, 0.1);
  RooGaussian swapped1y("swapped1y", "swapped k pi mass", m1, mean, sigmas1y);
  RooGaussian swapped2y("swapped2y", "swapped k pi mass", m1, mean, sigmas2y);
  RooRealVar swpfracy("swpfracy", "fraction of component 1 in swapped", 0.5, 0.001, 0.999);
  RooAddPdf swappedy("swappedy", "swp", RooArgList(swapped1y, swapped2y), swpfracy);

  swappedx.fitTo(swapds, Save());
  RooPlot *frame2 = m1.frame(Title("swapped_plot"), Bins(30));
  swapds.plotOn(frame2);
  swappedx.plotOn(frame2);

  sigmas1x.setConstant(true);
  sigmas2x.setConstant(true);
  swpfracx.setConstant(true);

  // RooRealVar rsigswap("rsigswap", "coeff of signal", 1. / 3.);
  // RooAddPdf sig("sig", "Signal", RooArgList(sigx, swapped), rsigswap);

  const long num_max = ds.sumEntries("", "D0range");

  // Initially, Nsig is set to 3% of the total entries, and Nswap = Nsig.
  RooRealVar nss("nss", "number of signal entries",
                  0.03 * num_max, 0.001 * num_max, 0.1 * num_max);
  RooRealVar nbb("nbb", "number of background entries",
                 0.95 * num_max, 0.6 * num_max, 1.1 * num_max);
  RooRealVar nsb("nsb", "number of background entries", 0.03 * num_max,
                 0.001 * num_max, 0.1 * num_max);
  RooRealVar nbs("nbs", "number of background entries", 0.03 * num_max,
                 0.001 * num_max, 0.1 * num_max);
  // RooRealVar nsw("nsw", "number of swapped entries", 0.06 * num_max,
  //                 0.01 * num_max, 0.2 * num_max);
  // RooRealVar nws("nws", "number of swapped entries", 0.06 * num_max,
  //                0.01 * num_max, 0.2 * num_max);

  //
  RooAddPdf modelx("modelx", "composite pdf", RooArgList(sigx, bkgx, sigx, bkgx), RooArgList(nss, nbb, nsb, nbs));
  // RooAddPdf modely("modely", "composite pdf", RooArgList(sigy, bkgy, bkgy, sigy), RooArgList(nss, nbb, nsb, nbs));
  RooAddPdf modely("modely", "composite pdf", RooArgList(sigx, bkgx, bkgx, sigx), RooArgList(nss, nbb, nsb, nbs));

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);

  // Associate model with the physics state and model_ctl with the control state
  simPdf.addPdf(modelx, "trig");
  simPdf.addPdf(modely, "asso");

  // U n b i n n e d   M L   f i t
  // Extended to get Nsig directly
  RooFitResult *r_ml_wgt_corr =
      simPdf.fitTo(combData, Save(), Extended(kTRUE), NumCPU(10));

  r_ml_wgt_corr->Print();

  // Calculate yields and error
  yield = nss.getValV();
  yieldErr = nss.getError();
  // yieldErr = nsig.getAsymErrorHi();

  // Plot the result
  auto canv = new TCanvas("Canvas", "Canvas", 800, 600);
  frame2->Draw();
  canv->SaveAs(outputname + "_swapfit" +".pdf");

  std::vector<TString> dim = {"trig", "asso"};

  for (auto m : dim) {
    TString name = m;
    RooPlot *plot = m1.frame(Title(name + "_plot"), Bins(60));
    combData.plotOn(plot, Name("pdat" + name), Cut("sample==sample::" + name));
    simPdf.plotOn(plot, Slice(sample, m), ProjWData(sample, combData),
                  VisualizeError(*r_ml_wgt_corr));
    combData.plotOn(plot, Name("pdat" + name), Cut("sample==sample::" + name),
              DataError(RooAbsData::SumW2));
    simPdf.plotOn(plot, Name("pfit" + name),
                  Slice(sample, m), ProjWData(sample, combData));
    simPdf.plotOn(plot, Name("psig" + name), Slice(sample, m),
                  ProjWData(sample, combData), Components(sigx),
                 DrawOption("LF"), FillStyle(3002), FillColor(kOrange - 3),
                 LineStyle(2), LineColor(kOrange - 3));
    simPdf.plotOn(plot, Name("pbkg" + name), Slice(sample, m),
                  ProjWData(sample, combData), Components(bkgx),
                 LineStyle(2), LineColor(4));
    // simPdf.plotOn(plot, Name("pswp" + name), Components(swapped),
    //               Slice(sample, m), ProjWData(sample, combData),
    //               DrawOption("LF"), FillStyle(3005), FillColor(kGreen + 4),
    //               LineStyle(1), LineColor(kGreen + 4));
    plot->Draw();

    TLegend *leg = new TLegend(0.69, 0.18, 0.89, 0.48, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);

    leg->AddEntry("pdat" + name, "Data", "pl");
    leg->AddEntry("pfit" + name, "Fit", "l");
    leg->AddEntry("psig" + name, "D^{0}+#bar{D^{#lower[0.2]{0}}} Signal", "f");
    // leg->AddEntry("pswp" + name, "K-#pi swapped", "f");
    leg->AddEntry("pbkg" + name, "Combinatorial", "l");

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
    // drawtex(texxpos, texypos = (texypos - texdypos),
    //         Form("(%.3f < m_{%s} < %.3f)", mass_dzero_signal_l, other.Data(),
    //              mass_dzero_signal_h));
    if (fdrawyield) {
      drawtex(texxpos, texypos = (texypos - texdypos),
              Form("N = %.0f #pm %.0f", yield, yieldErr));
    }

    if (fsaveplot) {
      canv->SaveAs(Form("%s_" + name + ".pdf", outputname.Data()));
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

  return r_ml_wgt_corr;
}
void xjjroot::fit2d::reset()
{
  clearvar();
}

void xjjroot::fit2d::init()
{
  clearvar();
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
  S = fun_mass->Integral(mass_dzero_signal_l,mass_dzero_signal_h)/binwid_hist_dzero;
  B = fun_background->Integral(mass_dzero_signal_l,mass_dzero_signal_h)/binwid_hist_dzero + fun_swap->Integral(mass_dzero_signal_h,mass_dzero_signal_h)/binwid_hist_dzero;
  Sig = S/TMath::Sqrt(S+B);
  // yield = fun_mass->Integral(min_hist_dzero,max_hist_dzero)/binwid_hist_dzero;
  // yieldErr = fun_mass->Integral(min_hist_dzero,max_hist_dzero)/binwid_hist_dzero*fun_mass->GetParError(0)/fun_mass->GetParameter(0);
  Chi2 = 2.*r->MinFcnValue();
  NDF = fun_f->GetNDF();
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

void xjjroot::fit2d::drawleg(TH1* h) const
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

void xjjroot::fit2d::drawtex(Double_t x, Double_t y, const char* text, Float_t tsize/*=0.04*/, Short_t align/*=12*/) const
{
  TLatex* tex = new TLatex(x, y, text);
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
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
