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

#include "dtree.h"
#include "ddbar.h"

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

void ddbar_saveDpair(std::string inputdata, std::string inputmc,
                     std::string output, std::string inputeff, float pt1min,
                     float pt1max, float pt2min, float pt2max, float dy,
                     float centmin, float centmax, std::string option = "2pi",
                     std::string treename = "d0ana/VertexCompositeNtuple") {

  ddbar::binning *binfo = new ddbar::binning(option);
  // TFile *infmc = TFile::Open(inputmc.c_str());

  // Get D0 efficiency from MC
  TFile *ineff = TFile::Open(inputeff.c_str());
  std::vector<TH2D *> eff(nCentrality);
  for (unsigned i = 0; i < nCentrality; ++i) {
    ineff->GetObject<TH2D>("efficiency_" + centralityString(i), eff[i]);
  }

  TFile *infdata = TFile::Open(inputdata.c_str());
  ddtree::dtree *dnt = new ddtree::dtree(infdata, treename);

  std::string outputname = "rootfiles/" + output + "/dptree.root";
  xjjroot::mkdir(outputname.c_str());
  TFile *outf = new TFile(outputname.c_str(), "recreate");
  Int_t iPhi;
  Float_t m1; // mass of D0
  Float_t m2; // mass of D0bar
  Float_t y1; // mass of D0
  Float_t y2; // mass of D0bar
  Float_t phi1, phi2; // phi angles
  Float_t pT1, pT2;
  Float_t weight, weightErr;
  Float_t dphi;

  TTree* dmptree = new TTree("dmass", "D0 mass pairs");
  dmptree->Branch("m1", &m1, "m1/F");
  dmptree->Branch("m2", &m2, "m2/F");
  dmptree->Branch("weight", &weight, "weight/F");
  dmptree->Branch("weightErr", &weightErr, "weightErr/F");
  dmptree->Branch("iPhi", &iPhi, "iPhi/I");
  dmptree->Branch("dphi", &dphi, "dphi/F");
  dmptree->Branch("y1", &y1, "y1/F");
  dmptree->Branch("y2", &y2, "y2/F");
  dmptree->Branch("phi1", &phi1, "phi1/F");
  dmptree->Branch("phi2", &phi2, "phi2/F");
  dmptree->Branch("pT1", &pT1, "pT1/F");
  dmptree->Branch("pT2", &pT2, "pT2/F");

  // Event loop
  // int nthreads = 1;
  // ROOT::EnableImplicitMT(nthreads);
  int nentries = dnt->nt()->GetEntries();
  if (nentries >= event_cutoff) {
    std::cout << "[Warning] Limited to " << event_cutoff << "events"
              << "\n";
  }
  for (int iEvt = 0; iEvt < nentries; iEvt++) {
    if (iEvt >= event_cutoff)
      break;
    if (iEvt % 100000 == 0)
      xjjc::progressbar(iEvt, nentries);
    dnt->nt()->GetEntry(iEvt);
    if (dnt->centrality > centmax * 2 || dnt->centrality < centmin * 2)
      continue;

    // unsigned centID = centralityID(dnt->centrality);
    for (int j = 0; j < dnt->candSize; j++) {
      pT1 = dnt->pT[j];
      if (pT1 < pt1min || pT1 > pt1max || fabs(dnt->y[j]) > dy) continue;
      for (int l = 0; l < j; l++) {
         pT2 = dnt->pT[l];
        if (pT2 < pt1min || pT2 > pt1max || fabs(dnt->y[l]) > dy) continue;
        if (dnt->flavor[l] * dnt->flavor[j] > 0) {
          continue;
        } else if (fabs(dnt->pT[j] - dnt->pT[l]) < 1e-4 && fabs(dnt->phi[j] - dnt->phi[l]) < 1e-4) {
          continue;
        }

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
        pT1 = dnt->pT[im1];
        pT2 = dnt->pT[im2];
        phi1 = dnt->phi[im1];
        phi2 = dnt->phi[im2];
        if (iPhi < 0) {
          std::cout << __FUNCTION__ << ": error: invalid dphi calculated."
                    << std::endl;
        }
        // double scale1 =
        //     1. / eff[centID]->GetBinContent(eff[centID]->FindBin(pT1, phi1));
        // double scale2 =
        //     1. / eff[centID]->GetBinContent(eff[centID]->FindBin(pT2, phi2));
        // double error1 =
        //     eff[centID]->GetBinError(eff[centID]->FindBin(pT1, phi1)) * scale1 *
        //     scale1;
        // double error2 =
        //     eff[centID]->GetBinError(eff[centID]->FindBin(pT2, phi2)) * scale2 *
        //     scale2;
        // weight = scale1 * scale2;
        // weightErr = error1 * error2;
        dmptree->Fill();
      }
    }
  }
  outf->cd();
  dmptree->Write();
  binfo->write();
  outf->Close();
}


void ddbar_swapmc(std::string swapmc, std::string output, std::string inputeff,
                  float ptmin, float ptmax,
                  float dy, float centmin, float centmax, std::string option = "2pi")
{
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
  TTreeReaderValue<Int_t> candSize_gen(reader, "candSize_gen");
  TTreeReaderArray<Float_t> pT_gen(reader, "pT_gen");
  TTreeReaderArray<Float_t> eta_gen(reader, "eta_gen");
  TTreeReaderArray<Float_t> y_gen(reader, "y_gen");
  TTreeReaderArray<Int_t> da1(reader, "DauID1_gen");
  TTreeReaderArray<Int_t> da2(reader, "DauID2_gen");
  TTreeReaderArray<Float_t> ptd1(reader, "pTD1");
  TTreeReaderArray<Float_t> ptd2(reader, "pTD2");
  TTreeReaderArray<Float_t> pterrd1(reader, "pTerrD1");
  TTreeReaderArray<Float_t> pterrd2(reader, "pTerrD2");
  TTreeReaderArray<Float_t> etad1(reader, "EtaD1");
  TTreeReaderArray<Float_t> etad2(reader, "EtaD2");

  // Get D0 efficiency from MC
  TFile *ineff = TFile::Open(inputeff.c_str());
  std::vector<TH2D *> eff(nCentrality);
  for (unsigned i = 0; i < nCentrality; ++i) {
    eff[i] = (TH2D *)ineff->Get("efficiency_" + centralityString(i));
  }

  std::string outputname = "rootfiles/" + output + "/swaptree.root";
  xjjroot::mkdir(outputname.c_str());
  TFile *outf = new TFile(outputname.c_str(), "recreate");
  Float_t m1, m2, weight;
  Float_t y1, y2;
  Float_t pt1, pt2;
  Float_t dphi;
  Int_t iPhi;
  Bool_t signal;
  Bool_t isSwap1, isSwap2;
  Bool_t sameK;
  Float_t pt1da1, pt1da2;
  Float_t pt2da1, pt2da2;
  Float_t pterr1da1, pterr1da2;
  Float_t pterr2da1, pterr2da2;
  Float_t eta1da1, eta1da2;
  Float_t eta2da1, eta2da2;
  Float_t dgen1, dgen2;
  Int_t lund1da1, lund1da2;
  Int_t lund2da1, lund2da2;
  TTree* swaptree = new TTree("swapdp", "swapped D0 pair");
  swaptree->Branch("m1", &m1, "m1/F");
  swaptree->Branch("m2", &m2, "m2/F");
  swaptree->Branch("dphi", &dphi, "dphi/F");
  swaptree->Branch("iPhi", &iPhi, "iPhi/I");
  swaptree->Branch("weight", &weight, "weight/F");
  swaptree->Branch("signal", &signal, "signal/B");
  swaptree->Branch("y1", &y1, "y1/F");
  swaptree->Branch("y2", &y2, "y2/F");
  swaptree->Branch("pt1", &pt1, "pt1/F");
  swaptree->Branch("pt2", &pt2, "pt2/F");
  swaptree->Branch("pterr1da1", &pterr1da1, "pterr1da1/F");
  swaptree->Branch("pterr1da2", &pterr1da2, "pterr1da2/F");
  swaptree->Branch("pterr2da1", &pterr2da1, "pterr2da1/F");
  swaptree->Branch("pterr2da2", &pterr2da2, "pterr2da2/F");
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

  // gen level
  swaptree->Branch("dgen1", &dgen1, "dgen1/F");
  swaptree->Branch("dgen2", &dgen2, "dgen2/F");
  swaptree->Branch("lund1da1", &lund1da1, "lund1da1/I");
  swaptree->Branch("lund2da1", &lund2da1, "lund2da1/I");
  swaptree->Branch("lund1da2", &lund1da2, "lund1da2/I");
  swaptree->Branch("lund2da2", &lund2da2, "lund2da2/I");
  TH1I* cat = new TH1I("cat", "category; D#bar{D} type; Events", 3, 0, 3);
  TH1I *nCand = new TH1I(
      "nCand", "nCand;# of gen D(#bar{D}) candidates; Events", 10, 0, 10);
  

// traverse the events
  while (reader.Next()) {
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
      // float pT1 = pT[iCand];
      // float phi1 = phi[iCand];
      for (unsigned jCand = 0; jCand < iCand; ++jCand) {
        // float pT2 = pT[jCand];
        // float phi2 = phi[jCand];
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
        pt1 = pT[i1];
        pt2 = pT[i2];
        pt1da1 = ptd1[i1];
        pt1da2 = ptd2[i1];
        pt2da1 = ptd1[i2];
        pt2da2 = ptd2[i2];
        pterr1da1 = pterrd1[i1];
        pterr1da2 = pterrd2[i1];
        pterr2da1 = pterrd1[i2];
        pterr2da2 = pterrd2[i2];
        eta1da1 = etad1[i1];
        eta1da2 = etad2[i1];
        eta2da1 = etad1[i2];
        eta2da2 = etad2[i2];
        // weight = scale;
        // match gen particles
        lund1da1 = -1;
        lund1da2 = -1;
        lund2da1 = -1;
        lund2da2 = -1;

        auto distance = [](float pt, float ptgen, float y, float ygen) {
          return std::pow(fabs(pt - ptgen) / 0.03, 2) +
            std::pow(fabs(y - ygen) / 0.002, 2);
        };
        std::vector<float> dist1;
        std::vector<float> dist2;
        for (int igen = 0; igen < *candSize_gen; ++igen) {
          dist1.push_back(distance(pT[i1], pT_gen[igen], y[i1], y_gen[igen]));
          dist2.push_back(distance(pT[i2], pT_gen[igen], y[i2], y_gen[igen]));
        }
        std::vector<float>::iterator min1 = std::min_element(dist1.begin(), dist1.end());
        int igen1 = std::distance(dist1.begin(), min1);
        std::vector<float>::iterator min2 = std::min_element(dist2.begin(), dist2.end());
        int igen2 = std::distance(dist2.begin(), min2);
        lund1da1 = da1[igen1];
        lund1da2 = da2[igen1];
        lund2da1 = da1[igen2];
        lund2da2 = da2[igen2];
        dgen1 = dist1[igen1];
        dgen2 = dist2[igen2];

        if (lund1da1 == -1 || lund2da1 == -1) {
          std::cout << "Warning: no matching gen particle!" << "\n";
        }


        if (iPhi < 0) {
          std::cout << __FUNCTION__ << ": error: invalid dphi calculated."
                    << std::endl;
        }
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
      // unsigned centID = centralityID(*centrality);
      // Calculate the scale factor from MC efficiency
      // double scale = 1. /
      // eff[centID]->GetBinContent(eff[centID]->FindBin(pT[iCand],
      // phi[iCand]));
    }
  }
  outf->cd();
  outf->Write();
  binfo->write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if (argc < 14) {
    std::cerr << "Not enough number of arguments! (" << argc << ")" << "\n";
    return 1;
  }
  ddbar_saveDpair(argv[1], argv[2], argv[3], argv[12], std::stof(argv[4]),
                  std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]),
                  std::stof(argv[8]), std::stof(argv[9]), std::stof(argv[10]),
                  argv[11]);
  ddbar_swapmc(argv[13], argv[3], argv[12], std::stof(argv[6]),
               std::stof(argv[7]), std::stof(argv[8]), std::stof(argv[9]),
               std::stof(argv[10]), argv[11]);
  return 0;
}

