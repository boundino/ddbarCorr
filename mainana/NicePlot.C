#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TColor.h"

void NicePlot(){
	gStyle->SetOptStat(0);
	TFile* f = TFile::Open("savehist.root");
	TH1F* phi_incl = (TH1F*)f->Get("phi_incl");
	TH1F* phi_sdbd = (TH1F*)f->Get("phi_sdbd");
	TH1F* phi_sgl = (TH1F*)f->Get("phi_sgl");
	TH1F* y_sgl = (TH1F*)f->Get("y_sgl");
	TCanvas* c = new TCanvas("c","c",800,800);
	c->DrawFrame(-TMath::Pi(),TMath::Pi(),0.,15e6);
	phi_incl->SetLineColor(kBlack);
        phi_sdbd->SetLineColor(kBlue);
        phi_sgl->SetLineColor(kRed);
	y_sgl->SetLineColor(kRed);
	phi_incl->Draw("E SAME");
	phi_sdbd->Draw("E SAME");
	phi_sgl->Draw("E SAME");
	gPad->BuildLegend();
	TCanvas* d = new TCanvas("d","d",800,800);
	d->cd();
	TFile* g = TFile::Open("/data/wangj/DntupleRun2018/skim_mc_analysisTree_promptD0_official.root");
	TTree* t = (TTree*)g->Get("d0ana_mc_genmatchunswap/VertexCompositeNtuple");
	phi_sgl->Scale(1./phi_sgl->Integral(),"width");
	phi_sgl->Draw("E SAME");
	TH1F* recophi = (TH1F*)phi_sgl->Clone("recophi");
	recophi->Reset();
	recophi->SetTitle("MC reco #phi");
	recophi->SetLineColor(kBlue);
	t->Draw("phi>>+recophi","abs(y)<1 && pT>2 && pT<8 && centrality/2 < 80","goff");
	//TH1F* recophi = (TH1F*)gDirectory->Get("recophi");
	recophi->Scale(1./recophi->Integral(),"width");
	recophi->Draw("E SAME");
	gPad->BuildLegend();
	TCanvas* k = new TCanvas("k","k",800,800);
	k->cd();
	y_sgl->Scale(1./y_sgl->Integral(),"width");
	y_sgl->Draw();
	TH1F* recoy = (TH1F*)y_sgl->Clone("recoy");
	recoy->Reset();
	recoy->SetTitle("MC reco y");
	recoy->SetLineColor(kBlue);
	t->Draw("y>>+recoy","abs(y)<1 && pT>2 && pT<8 && centrality/2 < 80","goff");
	recoy->Scale(1./recoy->Integral(),"width");
	recoy->Draw("E SAME");
	gPad->BuildLegend();
}
