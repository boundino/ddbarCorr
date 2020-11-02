#include "TEfficiency.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <vector>
#include <string>
#include <algorithm>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#ifdef SYSTEM_PBPB

#include "dtree.h"
#include "ddbar.h"
#include "event_mixer.h"

#endif 
#ifdef SYSTEM_PP

#include "dtree_pp.h"
#include "ddbar_pp.h"
#include "event_mixer_pp.h"

#endif

const int event_cutoff = (int) 1e10;
const float dmass_min = 1.72;
const float dmass_max = 2.0;
const std::vector<int> centralities = {0, 30, 50, 80};
const std::vector<double> d_centralities = {0.,30.,50.,80.};
const unsigned nCentrality = centralities.size() - 1;
const std::string MC_TREE = "d0ana_mc_genmatch/VertexCompositeNtuple";

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

void ddbar_saveDpair(std::string inputdata, std::string inputmc, std::string inputmixdata,
                     std::string output, std::string inputeff, float pt1min,
                     float pt1max, float pt2min, float pt2max, float dy,
                     float centmin, float centmax, std::string option = "2pi",
                     std::string treename = "d0ana/VertexCompositeNtuple") {
std::cout << "am alive\n";
  ddbar::binning *binfo = new ddbar::binning(option);
  // TFile *infmc = TFile::Open(inputmc.c_str());

  // Get D0 efficiency from MC
  TFile *ineff = TFile::Open(inputeff.c_str());
  std::vector<TH2D *> eff(nCentrality);
  for (unsigned i = 0; i < nCentrality; ++i) {
    ineff->GetObject<TH2D>("efficiency_" + centralityString(i), eff[i]);
  }

  TFile* infdata = TFile::Open(inputdata.c_str());
  ddtree::dtree *dnt = new ddtree::dtree(infdata,treename.c_str());

  TFile* infmix = TFile::Open(inputmixdata.c_str());
  ddtree::dtree *mdt = new ddtree::dtree(infmix,treename.c_str());
  EventMixer mixer(mdt);

  std::string outputname = "rootfiles/" + output + "/dptree.root";
  xjjroot::mkdir(outputname.c_str());
  Int_t iPhi;
  Float_t m1; // mass of D0
  Float_t m2; // mass of D0bar
  Float_t y1; // mass of D0
  Float_t y2; // mass of D0bar
  Float_t weight, weightErr;
  Float_t dphi;

  Int_t iPhi_sm;
  Float_t m1_sm;
  Float_t m2_sm;
  Float_t y1_sm;
  Float_t y2_sm;
  Float_t weight_sm, weightErr_sm;
  Float_t dphi_sm;

  Int_t iPhi_mm;
  Float_t m1_mm;
  Float_t m2_mm;
  Float_t y1_mm;
  Float_t y2_mm;
  Float_t weight_mm, weightErr_mm;
  Float_t dphi_mm;

  TFile *outf = new TFile(outputname.c_str(), "RECREATE");

  TTree* dmptree = new TTree("dmass", "D0 mass pairs (signal+signal)");
  dmptree->Branch("m1", &m1, "m1/F");
  dmptree->Branch("m2", &m2, "m2/F");
  dmptree->Branch("weight", &weight, "weight/F");
  dmptree->Branch("weightErr", &weightErr, "weightErr/F");
  dmptree->Branch("iPhi", &iPhi, "iPhi/I");
  dmptree->Branch("dphi", &dphi, "dphi/F");
  dmptree->Branch("y1", &y1, "y1/F");
  dmptree->Branch("y2", &y2, "y2/F");

  TTree* dmptree_sigmix = new TTree("dmass_sigmix", "D0 mass (signal+mixed)");
  dmptree_sigmix->SetDirectory(outf);
  dmptree_sigmix->Branch("m1", &m1_sm, "m1/F");
  dmptree_sigmix->Branch("m2", &m2_sm, "m2/F");
  dmptree_sigmix->Branch("weight", &weight_sm, "weight/F");
  dmptree_sigmix->Branch("weightErr", &weightErr_sm, "weightErr/F");
  dmptree_sigmix->Branch("iPhi", &iPhi_sm, "iPhi/I");
  dmptree_sigmix->Branch("dphi", &dphi_sm, "dphi/F");
  dmptree_sigmix->Branch("y1", &y1_sm, "y1/F");
  dmptree_sigmix->Branch("y2", &y2_sm, "y2/F");
  
  TTree* dmptree_mixmix = new TTree("dmass_mixmix", "D0 mass (mixed+mixed)");
  dmptree_mixmix->SetDirectory(outf);
  dmptree_mixmix->Branch("m1", &m1_mm, "m1/F");
  dmptree_mixmix->Branch("m2", &m2_mm, "m2/F");
  dmptree_mixmix->Branch("weight", &weight_mm, "weight/F");
  dmptree_mixmix->Branch("weightErr", &weightErr_mm, "weightErr/F");
  dmptree_mixmix->Branch("iPhi", &iPhi_mm, "iPhi/I");
  dmptree_mixmix->Branch("dphi", &dphi_mm, "dphi/F");
  dmptree_mixmix->Branch("y1", &y1_mm, "y1/F");
  dmptree_mixmix->Branch("y2", &y2_mm, "y2/F");
  
  std::vector<float> cents = {0, 30, 50, 80};
  std::vector<ddtree::mvas*> mvacents(cents.size()-1);
  for(int l=0; l<cents.size()-1; l++)
    { mvacents[l] = new ddtree::mvas(Form("DntupleRun2018/mvas/BDTCuts_finePtBins_PrompD0_PbPb_cent%.0fto%.0f_trkPt1GeV.root", cents[l], cents[l+1]),Form("cent%.0fto%.0f", cents[l], cents[l+1])); }

  // Event loop
  //int nthreads = 1;
  //ROOT::EnableImplicitMT(nthreads);
  int nentries = dnt->nt()->GetEntries();
  if (nentries >= event_cutoff) {
    std::cout << "[Warning] Limited to " << event_cutoff << "events"
              << "\n";
  }
  
  for (int iEvt = 0; iEvt < nentries; iEvt++) {
    if (iEvt >= event_cutoff)
      break;
    if (iEvt % 10000 == 0)
      xjjc::progressbar(iEvt, nentries);
    dnt->nt()->GetEntry(iEvt);
#ifdef SYSTEM_PBPB
    if (dnt->centrality > centmax * 2 || dnt->centrality < centmin * 2)
      continue;

    ddtree::dtree *mixtree = mixer.findSimilarEvent(dnt->centrality,dnt->bestvtxZ,iEvt);

    unsigned centID = centralityID(dnt->centrality);
    unsigned mixcentID = centralityID(mixtree->centrality);
#endif

#ifdef SYSTEM_PP
    ddtree::dtree *mixtree = mixer.findSimilarEvent(dnt->bestvtxZ,iEvt);
#endif
    for (int j = 0; j < dnt->candSize; j++) {
      float pT1 = dnt->pT[j];
      if (pT1 < pt1min || pT1 > pt1max || fabs(dnt->y[j]) > dy) continue;
      for (int l = 0; l < dnt->candSize; l++) {
        float pT2 = dnt->pT[l];
        if (pT2 < pt2min || pT2 > pt2max || fabs(dnt->y[l]) > dy) continue;
        if (pT2 > pt1min && l > j) continue;
        if (dnt->flavor[l] * dnt->flavor[j] > 0) continue;
        if (fabs(pT1 - pT2) < 1e-4 && fabs(dnt->phi[j] - dnt->phi[l]) < 1e-4) continue;

        int im1;
        int im2;
        iPhi = -1;
        if (dnt->flavor[j] > 0) {
          im1 = j;
          im2 = l;
        } else {
          im1 = l;
          im2 = j;
        }
        dphi = binfo->getdphi(dnt->phi[im1], dnt->phi[im2], iPhi);
        m1 = dnt->mass[im1];
        m2 = dnt->mass[im2];
        y1 = dnt->y[im1];
        y2 = dnt->y[im2];
        if (iPhi < 0) {
          std::cout << __FUNCTION__ << ": error: invalid dphi calculated."
                    << std::endl;
        }
#ifdef SYSTEM_PBPB
         //double scale1 =
         //    1. / eff[centID]->GetBinContent(eff[centID]->FindBin(pT1, phi1));
         //double scale2 =
         //    1. / eff[centID]->GetBinContent(eff[centID]->FindBin(pT2, phi2));
         //double error1 =
         //    eff[centID]->GetBinError(eff[centID]->FindBin(pT1, phi1)) * scale1 *
         //    scale1;
         //double error2 =
         //    eff[centID]->GetBinError(eff[centID]->FindBin(pT2, phi2)) * scale2 *
         //    scale2;
         //weight = scale1 * scale2;
         //weightErr = error1 * error2;
#endif
        dmptree->Fill();
      }
      for(int m = 0; m < mixtree->candSize; m++) {
        float pT2 = dnt->pT[m];
        if (pT2 < pt2min || pT2 > pt2max || fabs(mixtree->y[m]) > dy) continue;
        if (pT2 > pt1min && m > j) continue;
        if (mixtree->flavor[m] * dnt->flavor[j] > 0) continue;
        if (fabs(pT1 - pT2) < 1e-4 && fabs(dnt->phi[j] - mixtree->phi[m]) < 1e-4) continue;
#ifdef SYSTEM_PBPB
        if(mixcentID<cents.size()-1 && !mvacents[mixcentID]->passmva(mixtree->pT[m], mixtree->mva[m])) continue;
#endif

        iPhi_sm = -1;
        if (dnt->flavor[j] > 0) {
          dphi_sm = binfo->getdphi(dnt->phi[j],mixtree->phi[m],iPhi_sm);
          m1_sm = dnt->mass[j];
          y1_sm = dnt->y[j];
          m2_sm = mixtree->mass[m];
          y2_sm = mixtree->y[m];
        } else {
          dphi_sm = binfo->getdphi(mixtree->phi[m],dnt->phi[j],iPhi_sm);
          m1_sm = mixtree->mass[m];
          y1_sm = mixtree->y[m];
          m2_sm = dnt->mass[j];
          y2_sm = dnt->y[j];
        }
        if (iPhi_sm < 0) {
          std::cout << __FUNCTION__ << ": error: invalid dphi calculated."
                    << std::endl;
        }
        dmptree_sigmix->Fill();
      }
    }
    for(int n = 0; n < mixtree->candSize; n++) {
      float pT1 = mixtree->pT[n];
      if (pT1 < pt1min || pT1 > pt1max || fabs(mixtree->y[n]) > dy) continue;
#ifdef SYSTEM_PBPB
      if(mixcentID<cents.size()-1 && !mvacents[mixcentID]->passmva(mixtree->pT[n], mixtree->mva[n])) continue;
#endif
      for (int m = 0; m < mixtree->candSize; m++) {
        float pT2 = mixtree->pT[m];
        if (pT2 < pt2min || pT2 > pt2max || fabs(mixtree->y[m]) > dy) continue;
        if (pT2 > pt1min && m > n) continue;
        if (mixtree->flavor[n] * mixtree->flavor[m] > 0) continue;
        if (fabs(pT1 - pT2) < 1e-4 && fabs(mixtree->phi[n] - mixtree->phi[m]) < 1e-4) continue;
#ifdef SYSTEM_PBPB
        if(mixcentID<cents.size()-1 && !mvacents[mixcentID]->passmva(mixtree->pT[m], mixtree->mva[m])) continue;
#endif
        int im1;
        int im2;
        iPhi_mm = -1;
        if (mixtree->flavor[n] > 0) {
          im1 = n;
          im2 = m;
        } else {
          im1 = m;
          im2 = n;
        }
        dphi_mm = binfo->getdphi(dnt->phi[im1], dnt->phi[im2], iPhi_mm);
        m1_mm = dnt->mass[im1];
        m2_mm = dnt->mass[im2];
        y1_mm = dnt->y[im1];
        y2_mm = dnt->y[im2];
        if (iPhi_mm < 0) {
          std::cout << __FUNCTION__ << ": error: invalid dphi calculated."
                    << std::endl;
        }
        dmptree_mixmix->Fill();
      }
    }
  }
  outf->cd();
  dmptree->Write();
  dmptree_sigmix->Write();
  dmptree_mixmix->Write();
  outf->Close();
}

bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
   if (value.GetSetupStatus() < 0) {
      std::cerr << "Error " << value.GetSetupStatus()
                << "setting up reader for " << value.GetBranchName() << '\n';
      return false;
   }
   return true;
}

void ddbar_swapmc(std::string swapmc, std::string output, std::string inputeff,
                  float ptmin, float ptmax,
                  float dy, float centmin, float centmax, std::string option = "2pi")
{
  std::cout << "entering swapmc\n";
  std::cout << "using file " << swapmc << "\n";
  std::cout << "centmin: " << centmin << "\n";
  std::cout << "centmax: " << centmax << "\n";
  ddbar::binning *binfo = new ddbar::binning(option);

  TFile *infmc = TFile::Open(swapmc.c_str());
  TTreeReader reader(MC_TREE.c_str(), infmc);
  TTreeReaderValue<Int_t> centrality(reader, "centrality");
  TTreeReaderArray<Float_t> pT(reader, "pT");
  TTreeReaderArray<Float_t> y(reader, "y");
  TTreeReaderArray<Float_t> phi(reader, "phi");
  TTreeReaderArray<Float_t> mass_d0(reader, "mass");
  TTreeReaderArray<Float_t> flavor(reader, "flavor");
  TTreeReaderArray<Bool_t> isSwap(reader, "isSwap");
  TTreeReaderArray<Int_t> da1(reader, "DauID1_gen");
  TTreeReaderArray<Int_t> da2(reader, "DauID2_gen");
  TTreeReaderArray<Float_t> ptd1(reader, "pTD1");
  TTreeReaderArray<Float_t> ptd2(reader, "pTD2");
  TTreeReaderArray<Float_t> etad1(reader, "EtaD1");
  TTreeReaderArray<Float_t> etad2(reader, "EtaD2");
  std::cout << "got readers\n";
  // Get D0 efficiency from MC
  TFile *ineff = TFile::Open(inputeff.c_str());
  std::vector<TH2D *> eff(nCentrality);
  for (unsigned i = 0; i < nCentrality; ++i) {
    eff[i] = (TH2D *)ineff->Get("efficiency_" + centralityString(i));
  }
  std::cout << "got eff\n";
  std::string outputname = "rootfiles/" + output + "/swaptree.root";
  xjjroot::mkdir(outputname.c_str());
  TFile *outf = new TFile(outputname.c_str(), "recreate");
  Float_t m1, m2, weight;
  Float_t y1, y2;
  Float_t dphi;
  Int_t iPhi;
  Bool_t signal;
  Bool_t isSwap1, isSwap2;
  Bool_t sameK;
  Float_t pt1da1, pt1da2;
  Float_t pt2da1, pt2da2;
  Float_t eta1da1, eta1da2;
  Float_t eta2da1, eta2da2;
  TTree* swaptree = new TTree("swapdp", "swapped D0 pair");
  swaptree->Branch("m1", &m1, "m1/F");
  swaptree->Branch("m2", &m2, "m2/F");
  swaptree->Branch("dphi", &dphi, "dphi/F");
  swaptree->Branch("iPhi", &iPhi, "iPhi/I");
  swaptree->Branch("weight", &weight, "weight/F");
  swaptree->Branch("signal", &signal, "signal/B");
  swaptree->Branch("y1", &y1, "y1/F");
  swaptree->Branch("y2", &y2, "y2/F");
  swaptree->Branch("isSwap1", &isSwap1, "isSwap1/B");
  swaptree->Branch("isSwap2", &isSwap2, "isSwap2/B");
  swaptree->Branch("sameK", &sameK, "sameK/B");
  swaptree->Branch("pt1da1", &pt1da1, "pt1da1/F");
  swaptree->Branch("pt2da1", &pt2da1, "pt2da1/F");
  swaptree->Branch("eta1da1", &eta1da1, "eta1da1/F");
  swaptree->Branch("eta2da1", &eta2da1, "eta2da1/F");
  swaptree->Branch("pt1da2", &pt1da2, "pt1da2/F");
  swaptree->Branch("pt2da2", &pt2da2, "pt2da2/F");
  swaptree->Branch("eta1da2", &eta1da2, "eta1da2/F");
  swaptree->Branch("eta2da2", &eta2da2, "eta2da2/F");
  TH1I* cat = new TH1I("cat", "category; D#bar{D} type; Events", 3, 0, 3);
  TH1I *nCand = new TH1I(
      "nCand", "nCand;# of gen D(#bar{D}) candidates; Events", 10, 0, 10);
  
  std::cout << "started loop\n";
  std::cout << "entries: " << reader.GetEntries() << "\n";
// traverse the events
//  while (reader.Next()) {
  for(size_t i=0;i<reader.GetEntries();i++){
    reader.Next();
    reader.GetTree()->GetEntry(i);
    if(i%10000==0) std::cout << "Entry " << i << "/" << reader.GetEntries() << "\n";
    if (*centrality > centmax * 2 || *centrality < centmin * 2) {
      continue;
    }
    int size = 0;
    for (unsigned iGen = 0; iGen < da1.GetSize(); ++iGen) {
      if (da1[iGen] == 321 && da2[iGen] == 211) {
        size++;
      }
    }
    nCand->Fill(size);

    for (unsigned iCand = 0; iCand < pT.GetSize(); ++iCand) {
      if (pT[iCand] < ptmin || pT[iCand] > ptmax || fabs(y[iCand]) > dy) {
        continue;
      }
       float pT1 = pT[iCand];
       float phi1 = phi[iCand];
      for (unsigned jCand = 0; jCand < iCand; ++jCand) {
        float pT2 = pT[jCand];
        float phi2 = phi[jCand];
        sameK = false;
        if (pT[jCand] < ptmin || pT[jCand] > ptmax || fabs(y[jCand]) > dy) {
          continue;
          // They should have different flavor
        } else if (flavor[iCand] * flavor[jCand] > 0) {
          continue;
          // and should not be 2 candidates from the same 2 tracks
        // } else if (fabs(pT1 - pT2) < 1e-4 && fabs(phi1 - phi2) < 1e-4) {
        } else if (ptd1[iCand] == ptd1[jCand] && ptd2[iCand] == ptd2[jCand] &&
                       etad1[iCand] == etad1[jCand] && etad2[iCand] == etad2[jCand]) {
          // continue;
          sameK = true;
        } else if (ptd1[iCand] == ptd2[jCand] && ptd2[iCand] == ptd1[jCand] &&
                   etad1[iCand] == etad2[jCand] &&
                   etad2[iCand] == etad1[jCand]) {
        }
        iPhi = -1;
        unsigned i1, i2;
        if (flavor[iCand] > 0) {
          i1 = iCand;
          i2 = jCand;
        } else {
          i1 = jCand;
          i2 = iCand;
        }
        dphi = binfo->getdphi(phi[i1], phi[i2], iPhi);
        m1 = mass_d0[i1];
        m2 = mass_d0[i2];
        y1 = y[i1];
        y2 = y[i2];
        isSwap1 = isSwap[i1];
        isSwap2 = isSwap[i2];
        pt1da1 = ptd1[i1];
        pt1da2 = ptd2[i1];
        pt2da1 = ptd1[i2];
        pt2da2 = ptd2[i2];
        eta1da1 = etad1[i1];
        eta1da2 = etad2[i1];
        eta2da1 = etad1[i2];
        eta2da2 = etad2[i2];
#ifdef SYSTEM_PBPB
        // weight = scale;
#endif
        if (iPhi < 0) {
          std::cout << __FUNCTION__ << ": error: invalid dphi calculated."
                    << std::endl;
        }
        // only keep the swap event
        if (isSwap[iCand] ^ isSwap[jCand]) {
          signal = false;
          cat->Fill(1);
        } else if (!isSwap[iCand] & !isSwap[jCand]) {
          signal = true;
          cat->Fill(0);
        } else {
          cat->Fill(2);
        }
        swaptree->Fill();
      }
#ifdef SYSTEM_PBPB
      unsigned centID = centralityID(*centrality);
      // Calculate the scale factor from MC efficiency
      double scale = 1. /
        eff[centID]->GetBinContent(eff[centID]->FindBin(pT[iCand],
          phi[iCand]));
#endif
    }
  }
  outf->cd();
  outf->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if (argc < 14) {
    std::cerr << "Not enough number of arguments! (" << argc << ")" << "\n";
    return 1;
  }
  ddbar_saveDpair(argv[1], argv[2], argv[3], argv[4], argv[13], std::stof(argv[5]),
                  std::stof(argv[6]), std::stof(argv[7]), std::stof(argv[8]),
                  std::stof(argv[9]), std::stof(argv[10]), std::stof(argv[11]),
                  argv[12]);
  ddbar_swapmc(argv[14], argv[4], argv[13], std::stof(argv[7]),
               std::stof(argv[8]), std::stof(argv[9]), std::stof(argv[10]),
               std::stof(argv[11]), argv[12]);
  std::cout << "exited swapmc\n";
  return 0;
}

