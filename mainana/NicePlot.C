#include "TFile.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TH1.h"
#include "TColor.h"
#include "dtree.h"
#include "ddbar.h"
#include "xjjcuti.h"
#include "xjjrootuti.h"

const std::vector<int> centralities = {0, 30, 50, 80};
const unsigned nCentrality = centralities.size() - 1;

unsigned centralityID(int centrality) {
  if (centrality == 0) {
    return 0;
  }
  unsigned id = std::lower_bound(centralities.begin(), centralities.end(),
                                 centrality, [=] (float i, float j) {return i < j / 2;}) - centralities.begin();
  return id - 1;
}

void NicePlot(){
	gStyle->SetOptStat(0);
std::cout << "0\n";
	TFile* f = TFile::Open("savehist.root");
	TH1F* phi_incl = (TH1F*)f->Get("phi_incl");
	TH1F* phi_sdbd = (TH1F*)f->Get("phi_sdbd");
	TH1F* phi_sgl = (TH1F*)f->Get("phi_sgl");
	TH1F* y_sgl = (TH1F*)f->Get("y_sgl");
	TH1F* pt_sgl = (TH1F*)f->Get("pt_sgl");
	TH1F* dca_sgl = (TH1F*)f->Get("dca_sgl");
std::cout << "1\n";
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
std::cout << "2\n";
	TFile* g = TFile::Open("/data/wangj/DntupleRun2018/skim_mc_analysisTree_promptD0_official.root");
	TTree* t = (TTree*)g->Get("d0ana_mc_genmatchunswap/VertexCompositeNtuple");
	ddtree::dtree* dt = new ddtree::dtree(g,"d0ana_mc_genmatchunswap/VertexCompositeNtuple");
std::cout << "3\n";
	TCanvas* d = new TCanvas("d","d",800,800);
	d->cd();
	phi_sgl->Scale(1./phi_sgl->Integral(),"width");
	phi_sgl->Draw("E SAME");
	TH1F* recophi = (TH1F*)phi_sgl->Clone("recophi");
	recophi->Reset();
	recophi->SetTitle("MC reco #phi");
	recophi->SetLineColor(kBlue);
std::cout << "about to get\n";
	TH1F* data_cent = (TH1F*)f->Get("cent");
	TH1F* mc_cent = (TH1F*)f->Get("mc_cent");
	TH1F* centweight = (TH1F*)data_cent->Clone("centweight");
std::cout << "got\n";
	centweight->Reset();
	centweight->Sumw2();
	centweight->Divide(data_cent,mc_cent);
	int nentries = dt->nt()->GetEntries();
std::cout << "4\n";
	for(int i=0;i<nentries;i++)
	{
		if(i&100000==0) xjjc::progressbar(i,nentries);
		dt->nt()->GetEntry(i);
		if(dt->centrality>160) continue;
		float scale = centweight->GetBinContent(centralityID(dt->centrality));
		for(int j=0;j<dt->candSize;j++) recophi->Fill(dt->phi[j],scale);
	}	
std::cout << "5\n";
	//t->Draw("phi>>+recophi","abs(y)<1 && pT>2 && pT<8 && centrality/2 < 80","goff");
	//TH1F* recophi = (TH1F*)gDirectory->Get("recophi");
	recophi->Scale(1./recophi->Integral(),"width");
	recophi->Draw("E SAME");
	gPad->BuildLegend();
	gPad->Print("phi_weight.png");
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
	gPad->Print("y.png");
std::cout << "6\n";
	TCanvas* canpt = new TCanvas("canpt", "canpt", 800, 1200);
	canpt->cd();
	canpt->SetLogy();
	pt_sgl->Scale(1./pt_sgl->Integral(),"width");
	pt_sgl->SetLineColor(kRed);
	// pt_sgl->Draw();
	pt_sgl->Sumw2();
	TH1F* recopt = (TH1F*)pt_sgl->Clone("recopt");
	recopt->Reset();
	recopt->SetTitle("MC reco pt");
	recopt->SetLineColor(kBlue);
	t->Draw("pT>>+recopt", "abs(y)<1 && centrality/2 < 80","goff");
	recopt->Scale(1./recopt->Integral(),"width");
	// recopt->Draw("E SAME");
	TRatioPlot* rpt = new TRatioPlot(pt_sgl, recopt);
	rpt->SetH1DrawOpt("e");
	rpt->Draw();
	gPad->BuildLegend();
std::cout << "7\n";
	TCanvas* candca = new TCanvas("candca", "candca", 800, 1200);
	candca->SetLogy();
	// candca->SetLogx();
	dca_sgl->Scale(1./dca_sgl->Integral(),"width");
	dca_sgl->SetLineColor(kRed);
	dca_sgl->Sumw2();
	TH1F* recodca = (TH1F*)dca_sgl->Clone("recodca");
	recodca->Reset();
	recodca->SetTitle("MC reco DCA");
	recodca->SetLineColor(kBlue);
	t->Draw("dca>>recodca", "abs(y)<1 && pT>2 && pT<8 && centrality/2 < 80","goff");
	recodca->Scale(1./recodca->Integral(),"width");
	recodca->Draw();
	TRatioPlot* rdca = new TRatioPlot(dca_sgl, recodca);
	rdca->SetH1DrawOpt("e");
	rdca->Draw();
	gPad->BuildLegend();
}

int main(){
	NicePlot();
	return 0;
}
