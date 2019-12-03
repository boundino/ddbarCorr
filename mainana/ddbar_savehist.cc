#include <TFile.h>
#include <TH1F.h>

#include <vector>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "dtree.h"
#include "ddbar.h"

void ddbar_savehist(std::string inputdata, std::string inputmc, std::string output, 
                    float pt1min, float pt1max, float pt2min, float pt2max, float yd, float centmin, float centmax, std::string option="2pi")
{
  ddbar::kinematics* kinfo = new ddbar::kinematics(pt1min, pt1max, pt2min, pt2max, yd, centmin, centmax);
  ddbar::binning* binfo = new ddbar::binning(option);

  TFile* infmc = TFile::Open(inputmc.c_str());
  TH1F* hmassmc_sgl = (TH1F*)infmc->Get("hHistoRMassSignal_pt_0_dr_0");
  hmassmc_sgl->SetName("hmassmc_sgl");
  TH1F* hmassmc_swp = (TH1F*)infmc->Get("hHistoRMassSwapped_pt_0_dr_0");
  hmassmc_swp->SetName("hmassmc_swp");

  TFile* infdata = TFile::Open(inputdata.c_str());
  ddtree::dtree* dnt = new ddtree::dtree(infdata, "d0ana/VertexCompositeNtuple");

  std::vector<TH1F*> hmass_incl(binfo->nphibin(), 0), hmass_sdbd(binfo->nphibin(), 0);
  for(int k=0; k<binfo->nphibin(); k++) 
    {
      hmass_incl[k] = new TH1F(Form("hmass_dphi%d_incl", k), ";m_{K#pi} (GeV/c^{2});Entries", 60, 1.7, 2.0);
      hmass_sdbd[k] = new TH1F(Form("hmass_dphi%d_sdbd", k), ";m_{K#pi} (GeV/c^{2});Entries", 60, 1.7, 2.0);
    }
  TH1F* hmass_trig = new TH1F("hmass_trig", ";m_{K#pi} (GeV/c^{2});Entries", 60, 1.7, 2.0);

  int nentries = dnt->nt()->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      if(i%100000==0) xjjc::progressbar(i, nentries);
      dnt->nt()->GetEntry(i);
      if(dnt->centrality > centmax*2 || dnt->centrality < centmin*2) continue;

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

          for(int l=0; l<dnt->candSize; l++)
            {
              if(dnt->pT[l] < pt2min || dnt->pT[l] > pt2max) continue;
              if(fabs(dnt->y[l]) > yd) continue;
              if(dnt->pT[l] >= dnt->pT[j] || fabs(dnt->pT[l]-dnt->pT[j])<0.0001 || dnt->flavor[l]*dnt->flavor[j] > 0) continue;

              int idphi = -1;
              float dphi = binfo->getdphi(dnt->phi[j], dnt->phi[l], idphi);
              if(idphi < 0) { std::cout<<__FUNCTION__<<": error: invalid dphi calculated."<<std::endl; }
              if(isinclusive) hmass_incl[idphi]->Fill(dnt->mass[l]);
              if(issideband) hmass_sdbd[idphi]->Fill(dnt->mass[l]);
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  
  std::string outputname = "rootfiles/" + output + "/savehist.root";
  xjjroot::mkdir(outputname.c_str());
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  hmassmc_sgl->Write();
  hmassmc_swp->Write();
  hmass_trig->Write();
  for(auto& hh : hmass_incl) hh->Write();
  for(auto& hh : hmass_sdbd) hh->Write();
  kinfo->write();
  binfo->write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==12) { ddbar_savehist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]), argv[11]); return 0; }
  if(argc==11) { ddbar_savehist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10])); return 0; }
  return 1;
}

