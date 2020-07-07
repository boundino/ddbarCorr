#include "TFile.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TMath.h"

#include <vector>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "dtree_pp.h"
#include "ddbar_pp.h"

const int event_cutoff = (int) 1e10;
// const int event_cutoff = (int) 1e6;
const float dmass_min = 1.72;
const float dmass_max = 2.0;
/*
const std::vector<int> centralities = {0, 30, 50, 80};
const std::vector<double> d_centralities = {0.,30.,50.,80.};
const unsigned nCentrality = centralities.size() - 1;

TString centralityString(unsigned iCent) {
  return TString::Format("cent%dto%d", centralities[iCent], centralities[iCent + 1]);
}

unsigned centralityID(int centrality) {
  if (centrality == 0) {
    return 0;
  }
  unsigned id = std::lower_bound(centralities.begin(), centralities.end(),
                                 centrality, [=] (float i, float j) {return i < j / 2;}) - centralities.begin();
  return id - 1;
}
*/
void ddbar_savehist(std::string inputdata, std::string inputmc, std::string output,
                    float pt1min, float pt1max, float pt2min, float pt2max,
                    float yd,
                    std::string option="2pi", std::string inputeff ="efficiency.root")
{
  ddbar::kinematics* kinfo = new ddbar::kinematics(pt1min, pt1max, pt2min, pt2max, yd);
  ddbar::binning* binfo = new ddbar::binning(option);
  
  TFile* infmc = TFile::Open(inputmc.c_str());
  TH1F* hmassmc_sgl = (TH1F*)infmc->Get("hHistoRMassSignal_pt_0_dr_0");
  hmassmc_sgl->SetName("hmassmc_sgl");
  // swapped K pi mass
  TH1F* hmassmc_swp = (TH1F*)infmc->Get("hHistoRMassSwapped_pt_0_dr_0");
  hmassmc_swp->SetName("hmassmc_swp");
/*
  // Get D0 efficiency from MC
  TFile* ineff = TFile::Open(inputeff.c_str());
  std::vector<TH2D*> eff(3);
  for (unsigned i = 0; i < nCentrality; ++i)
    {
      eff[i] = (TH2D*) ineff->Get("efficiency_" + centralityString(i));
    }
*/
  TFile* infdata = TFile::Open(inputdata.c_str());
  ddtree::dtree* dnt = new ddtree::dtree(infdata, "ntDkpi");

  std::vector<TH1F*> hmass_incl(binfo->nphibin(), 0), hmass_sdbd(binfo->nphibin(), 0);
  std::vector<TH1F*> hmass_incl_scaled(binfo->nphibin(), 0), hmass_sdbd_scaled(binfo->nphibin(), 0), hmass_incl_det(binfo->nphibin(), 0);
//  TH1F* cent = new TH1F("cent","Centrality",nCentrality,d_centralities.data());
  TH1F* phi = new TH1F("phi", "#phi", 36, -TMath::Pi(), TMath::Pi());
  TH1F* phi_sc = new TH1F("phi_sc", "#phi eff. scaled", 36, -TMath::Pi(), TMath::Pi());
  TH1F* phi_incl = new TH1F("phi_incl","Inclusive #phi eff. scaled",36,-TMath::Pi(),TMath::Pi());
  TH1F* phi_sdbd = new TH1F("phi_sdbd","Sideband #phi eff. scaled",36,-TMath::Pi(),TMath::Pi());
  TH2F* phipt = new TH2F("phipt", "#phi / #p_{T}", 36, -TMath::Pi(), TMath::Pi(), 60, 2, 8);
  TH2F* phipt_sc = new TH2F("phipt_sc", "#phi / #p_{T} eff. scaled", 36, -TMath::Pi(), TMath::Pi(), 60, 2, 8);
  TH1F* y_incl = new TH1F("y_incl","Inclusive y eff. scaled",40,-1.,1.);
  TH1F* y_sdbd = new TH1F("y_sdbd","Sideband y eff. scaled",40,-1.,1.);
  TH1F* pt_incl = new TH1F("pt_incl", "Inclusive pT eff. scaled", 40, 2, 16);
  TH1F* pt_sdbd = new TH1F("pt_sdbd", "Sideband pT eff. scaled", 40, 2, 16);
  TH1F* dca_incl = new TH1F("dca_incl", "Inclusive dca eff. scaled", 100, 0, 1);
  TH1F* dca_sdbd = new TH1F("dca_sdbd", "Sideband dca eff. scaled", 100, 0, 1);
//  TH1F* hiBin = new TH1F("hiBin","Data hiBin",200,0,200);
  for(int k=0; k<binfo->nphibin(); k++)
    {
      // inclusive or signal region
      hmass_incl[k] = new TH1F(Form("hmass_dphi%d_incl", k), ";m_{K#pi} (GeV/c^{2});Entries", 60, dmass_min, dmass_max);
      // side-band region
      hmass_sdbd[k] = new TH1F(Form("hmass_dphi%d_sdbd", k), ";m_{K#pi} (GeV/c^{2});Entries", 60, dmass_min, dmass_max);
      // scaled mass
      hmass_incl_scaled[k] = new TH1F(Form("hmass_dphi%d_incl_sc", k), ";m_{K#pi} (GeV/c^{2});Entries", 60, dmass_min, dmass_max);
      hmass_sdbd_scaled[k] = new TH1F(Form("hmass_dphi%d_sdbd_sc", k), ";m_{K#pi} (GeV/c^{2});Entries", 60, dmass_min, dmass_max);
      hmass_incl_det[k] = new TH1F(Form("hmass_phi%d_incl_det", k), ";m_{K#pi} (GeV/c^{2});Entries", 60, dmass_min, dmass_max);
    }
  TH1F* hmass_trig = new TH1F("hmass_trig", ";m_{K#pi} (GeV/c^{2});Entries", 60, dmass_min, dmass_max);
//  TH1F* mc_cent = (TH1F*)ineff->Get("mc_cent");
  int nentries = dnt->nt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      if(i >= event_cutoff) break;
      if(i%100000==0) xjjc::progressbar(i, nentries);
      dnt->nt()->GetEntry(i);
//      if(dnt->centrality > centmax*2 || dnt->centrality < centmin*2) continue;
//      cent->Fill(dnt->centrality/2);
//      hiBin->Fill(dnt->centrality);
//      unsigned centID = centralityID(dnt->centrality);
      for(int j=0; j<dnt->candSize; j++)
        {
          if(dnt->pT[j] < pt1min || dnt->pT[j] > pt1max) continue;
          if(fabs(dnt->y[j]) > yd) continue;
          hmass_trig->Fill(dnt->mass[j]);
          // int fillid = -1;
          int isinclusive = 0, issideband = 0;
          if(fabs(dnt->mass[j] - MASS_DZERO) < ddbar::signalwidth) isinclusive = 1;
          if(fabs(dnt->mass[j] - MASS_DZERO) > ddbar::sideband_l && fabs(dnt->mass[j] - MASS_DZERO) < ddbar::sideband_h) issideband = 1;
          if(!(isinclusive+issideband)) continue;
          float primary_pT = dnt->pT[j];
          float primary_y = dnt->y[j];
          // Get phi distribution
          phi->Fill(dnt->phi[j]);
          phipt->Fill(dnt->phi[j], dnt->pT[j]);
          double scale = 1.;///eff[centID]->GetBinContent(eff[centID]->FindBin(primary_pT));
          phi_sc->Fill(dnt->phi[j], scale);
          if(isinclusive) {
            phi_incl->Fill(dnt->phi[j], scale);
            y_incl->Fill(dnt->y[j], scale);
            pt_incl->Fill(dnt->pT[j], scale);
            dca_incl->Fill(dnt->dca[j], scale);
          }
          if(issideband) {
            y_sdbd->Fill(dnt->y[j], scale);
            phi_sdbd->Fill(dnt->phi[j], scale);
            pt_sdbd->Fill(dnt->pT[j], scale);
            dca_sdbd->Fill(dnt->dca[j], scale);
          }
          phipt_sc->Fill(dnt->phi[j], dnt->pT[j], scale);

          std::vector<float> phis;
          for (auto i = 0; i <= 8; ++i) {
            phis.push_back((-1 + i * 0.25) * TMath::Pi());
          }
          unsigned iphi = std::lower_bound(phis.begin(), phis.end(), dnt->phi[j]) - phis.begin() - 1;
          hmass_incl_det[iphi]->Fill(dnt->mass[j]);
          for(int l=0; l<dnt->candSize; l++)
            {
              float associate_pT = dnt->pT[l];
              float associate_y = dnt->y[l];
              if(dnt->pT[l] < pt2min || dnt->pT[l] > pt2max) continue;
              if(fabs(dnt->y[l]) > yd) continue;
              if(dnt->pT[l] >= dnt->pT[j] || fabs(dnt->pT[l]-dnt->pT[j])<0.0001 || dnt->flavor[l]*dnt->flavor[j] != 2) continue;
              // Calculate the scale factor from MC efficiency
              double scale_primary = 1.;///eff[centID]->GetBinContent(eff[centID]->FindBin(primary_pT));
              double scale_associate = 1.;///eff[centID]->GetBinContent(eff[centID]->FindBin(associate_pT));
              // Fill the other D0 mass in the specific phi region
              int idphi = -1;
              float dphi = binfo->getdphi(dnt->phi[j], dnt->phi[l], idphi);
              if(idphi < 0) { std::cout<<__FUNCTION__<<": error: invalid dphi calculated."<<std::endl; }
              if(isinclusive)
                {
                  hmass_incl[idphi]->Fill(dnt->mass[l]);
                  hmass_incl_scaled[idphi]->Fill(dnt->mass[l], scale_primary * scale_associate);
                }
              if(issideband)
                {
                  hmass_sdbd[idphi]->Fill(dnt->mass[l]);
                  hmass_sdbd_scaled[idphi]->Fill(dnt->mass[l], scale_primary * scale_associate);
                }
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  TH1F* phi_sgl = (TH1F*)phi_incl->Clone("phi_sgl");
  TH1F* y_sgl = (TH1F*)y_incl->Clone("y_sgl");
  TH1F* pt_sgl = (TH1F*)pt_incl->Clone("pt_sgl");
  TH1F* dca_sgl = (TH1F*)dca_incl->Clone("dca_sgl");
  phi_incl->Sumw2();
  phi_sdbd->Sumw2();
  y_incl->Sumw2();
  y_sdbd->Sumw2();
  pt_incl->Sumw2();
  pt_sdbd->Sumw2();
  phi_sgl->SetTitle("Signal #phi");
  y_sgl->SetTitle("Signal y");
  pt_sgl->SetTitle("Signal pT");
  dca_sgl->SetTitle("Signal dca");
  phi_sgl->Sumw2();
  y_sgl->Sumw2();
  pt_sgl->Sumw2();
  dca_sgl->Sumw2();
  phi_sgl->Add(phi_sdbd, -1);
  y_sgl->Add(y_sdbd, -1);
  pt_sgl->Add(pt_sdbd, -1);
  pt_sgl->Add(pt_sdbd, -1);
  std::string outputname = "rootfiles/" + output + "/savehist.root";
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  hmassmc_sgl->Write();
  hmassmc_swp->Write();
  hmass_trig->Write();
  for(auto& hh : hmass_incl)
    {
      hh->Sumw2();
      hh->Write();
    }
  for(auto& hh : hmass_sdbd)
    {
      hh->Sumw2();
      hh->Write();
    }
  for(auto& hh : hmass_incl_scaled)
    {
      hh->Sumw2();
      hh->Write();
    }
  for(auto& hh : hmass_sdbd_scaled)
    {
      hh->Sumw2();
      hh->Write();
    }
  for(auto& hh : hmass_incl_det)
    {
      hh->Sumw2();
      hh->Write();
    }
//  std::cout << cent->GetEntries() << "\n";
//  cent->Sumw2();
//  cent->Scale(1./cent->Integral(),"width");
//  cent->Write();
//  mc_cent->Write();
  phi->Write();
  phi_sc->Sumw2();
  phi_sc->Write();
  phi_incl->Write();
  phi_sgl->Write();
  phi_sdbd->Write();
  y_incl->Write();
  y_sdbd->Write();
  y_sgl->Write();
  pt_incl->Write();
  pt_sdbd->Write();
  pt_sgl->Write();
  dca_incl->Write();
  dca_sdbd->Write();
  dca_sgl->Write();
  phipt->Write();
  phipt_sc->Write();
//  hiBin->Write();
  kinfo->write();
  binfo->write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==11) { ddbar_savehist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), argv[9], argv[10]); return 0; }
  if(argc==10) { ddbar_savehist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), argv[9]); return 0; }
  if(argc==9) { ddbar_savehist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8])); return 0; }
  return 1;
}

