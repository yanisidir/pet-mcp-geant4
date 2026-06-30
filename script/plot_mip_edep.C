#include "TAxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  const int kMaxMcpBranchesToScan = 10;
  const int kNumberOfBins = 120;

  struct BranchStats
  {
    Long64_t entries;
    Long64_t nonZeroEntries;
    double sumKeV;
    double sum2KeV;
    double minKeV;
    double maxKeV;

    BranchStats()
      : entries(0),
        nonZeroEntries(0),
        sumKeV(0.0),
        sum2KeV(0.0),
        minKeV(std::numeric_limits<double>::max()),
        maxKeV(-std::numeric_limits<double>::max())
    {
    }
  };

  bool HasBranch(TTree* tree, const std::string& branchName)
  {
    return tree && tree->GetBranch(branchName.c_str());
  }

  int DetectNumberOfMCPs(TFile* file, TTree* tree)
  {
    TTree* plateTree = nullptr;
    if (file) {
      file->GetObject("McpPlateStatsTree", plateTree);
    }

    if (plateTree && plateTree->GetBranch("mcpIndex")) {
      int mcpIndex = -1;
      int maxMcpIndex = -1;
      plateTree->SetBranchStatus("*", 0);
      plateTree->SetBranchStatus("mcpIndex", 1);
      plateTree->SetBranchAddress("mcpIndex", &mcpIndex);

      const Long64_t entries = plateTree->GetEntries();
      for (Long64_t entry = 0; entry < entries; ++entry) {
        plateTree->GetEntry(entry);
        maxMcpIndex = std::max(maxMcpIndex, mcpIndex);
      }

      plateTree->ResetBranchAddresses();
      plateTree->SetBranchStatus("*", 1);

      if (maxMcpIndex >= 0) {
        return maxMcpIndex + 1;
      }
    }

    int numberOfMCPs = 0;
    for (int i = 0; i < kMaxMcpBranchesToScan; ++i) {
      std::ostringstream branchName;
      branchName << "edepPlusMcp" << i << "_keV";
      if (HasBranch(tree, branchName.str())) {
        numberOfMCPs = i + 1;
      }
    }
    return numberOfMCPs;
  }

  BranchStats ComputeBranchStats(TTree* tree,
                                 const std::string& branchName)
  {
    BranchStats stats;
    if (!HasBranch(tree, branchName)) {
      std::cerr << "Missing branch: " << branchName << std::endl;
      return stats;
    }

    double valueKeV = 0.0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(branchName.c_str(), 1);
    tree->SetBranchAddress(branchName.c_str(), &valueKeV);

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
      tree->GetEntry(entry);
      ++stats.entries;
      stats.sumKeV += valueKeV;
      stats.sum2KeV += valueKeV*valueKeV;
      if (valueKeV > 0.0) {
        ++stats.nonZeroEntries;
      }
      stats.minKeV = std::min(stats.minKeV, valueKeV);
      stats.maxKeV = std::max(stats.maxKeV, valueKeV);
    }

    tree->ResetBranchAddresses();
    tree->SetBranchStatus("*", 1);

    if (stats.entries == 0) {
      stats.minKeV = 0.0;
      stats.maxKeV = 0.0;
    }
    return stats;
  }

  double MeanMeV(const BranchStats& stats)
  {
    return stats.entries > 0 ? stats.sumKeV/stats.entries/1000.0 : 0.0;
  }

  double RmsMeV(const BranchStats& stats)
  {
    if (stats.entries <= 0) {
      return 0.0;
    }

    const double meanKeV = stats.sumKeV/stats.entries;
    const double mean2KeV = stats.sum2KeV/stats.entries;
    const double variance = std::max(0.0, mean2KeV - meanKeV*meanKeV);
    return std::sqrt(variance)/1000.0;
  }

  void FillHistogramFromBranch(TTree* tree,
                               const std::string& branchName,
                               TH1D* histogram,
                               double scale)
  {
    if (!HasBranch(tree, branchName) || !histogram) {
      return;
    }

    double value = 0.0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(branchName.c_str(), 1);
    tree->SetBranchAddress(branchName.c_str(), &value);

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
      tree->GetEntry(entry);
      if (value > 0.0) {
        histogram->Fill(value*scale);
      }
    }

    tree->ResetBranchAddresses();
    tree->SetBranchStatus("*", 1);
  }
}

void plot_mip_edep(const char* fileName = "build/mcp_output.root")
{
  gSystem->mkdir("Fig/current", kTRUE);

  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open file: " << fileName << std::endl;
    return;
  }

  TTree* tree = nullptr;
  file->GetObject("EventSummaryTree", tree);
  if (!tree) {
    std::cerr << "EventSummaryTree not found." << std::endl;
    file->Close();
    return;
  }

  gStyle->SetOptStat(1110);

  const int numberOfMCPs = DetectNumberOfMCPs(file, tree);
  if (numberOfMCPs <= 0) {
    std::cerr << "No edepPlusMcp*_keV branches found." << std::endl;
    file->Close();
    delete file;
    return;
  }

  const BranchStats totalStats =
    ComputeBranchStats(tree, "edepTotal_keV");
  const BranchStats plusStats =
    ComputeBranchStats(tree, "edepPlus_keV");
  const BranchStats minusStats =
    ComputeBranchStats(tree, "edepMinus_keV");

  const double maxTotalMeV =
    totalStats.maxKeV > 0.0 ? 1.05*totalStats.maxKeV/1000.0 : 1.0;

  TCanvas* c1 = new TCanvas("c1", "MIP deposited energy", 900, 700);
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SetGridy();

  TH1D* hTotal = new TH1D("hTotal",
                          "MIP deposited energy per event;E_{dep} [MeV];Events",
                          kNumberOfBins,
                          0.0,
                          maxTotalMeV);
  hTotal->SetDirectory(nullptr);
  FillHistogramFromBranch(tree, "edepTotal_keV", hTotal, 1.0/1000.0);
  hTotal->SetLineColor(kBlack);
  hTotal->SetLineWidth(3);
  hTotal->Draw("HIST");

  c1->SaveAs("Fig/current/mip_edep_total.png");

  gStyle->SetOptStat(0);
  TCanvas* c2 = new TCanvas("c2", "Mean deposited energy per side", 900, 700);
  c2->SetLeftMargin(0.13);
  c2->SetBottomMargin(0.12);
  c2->SetGridy();

  TH1D* hSide = new TH1D("hSide",
                         "Mean MIP deposited energy per side;Detector side;Mean E_{dep} [MeV]",
                         2,
                         0,
                         2);

  hSide->GetXaxis()->SetBinLabel(1, "+z stack");
  hSide->GetXaxis()->SetBinLabel(2, "-z stack");

  hSide->SetBinContent(1, MeanMeV(plusStats));
  hSide->SetBinContent(2, MeanMeV(minusStats));

  hSide->SetBarWidth(0.6);
  hSide->SetBarOffset(0.2);
  hSide->SetFillColor(kBlue - 9);
  hSide->Draw("BAR");

  c2->SaveAs("Fig/current/mip_edep_per_side.png");

  TCanvas* c3 = new TCanvas("c3", "Mean deposited energy per MCP", 1000, 700);
  c3->SetLeftMargin(0.12);
  c3->SetBottomMargin(0.12);
  c3->SetGridy();

  TH1D* hPlus = new TH1D("hPlus",
                         "Mean MIP deposited energy per MCP;MCP index;Mean E_{dep} [MeV]",
                         numberOfMCPs,
                         -0.5,
                         numberOfMCPs - 0.5);

  TH1D* hMinus = new TH1D("hMinus",
                          "Mean MIP deposited energy per MCP;MCP index;Mean E_{dep} [MeV]",
                          numberOfMCPs,
                          -0.5,
                          numberOfMCPs - 0.5);

  std::vector<BranchStats> plusByMcp(numberOfMCPs);
  std::vector<BranchStats> minusByMcp(numberOfMCPs);
  for (int i = 0; i < numberOfMCPs; ++i) {
    std::ostringstream plusName;
    plusName << "edepPlusMcp" << i << "_keV";

    std::ostringstream minusName;
    minusName << "edepMinusMcp" << i << "_keV";

    plusByMcp[i] = ComputeBranchStats(tree, plusName.str());
    minusByMcp[i] = ComputeBranchStats(tree, minusName.str());

    hPlus->SetBinContent(i + 1, MeanMeV(plusByMcp[i]));
    hMinus->SetBinContent(i + 1, MeanMeV(minusByMcp[i]));

    std::ostringstream label;
    label << "MCP " << i;
    hPlus->GetXaxis()->SetBinLabel(i + 1, label.str().c_str());
    hMinus->GetXaxis()->SetBinLabel(i + 1, label.str().c_str());
  }

  hPlus->SetLineColor(kRed + 1);
  hPlus->SetMarkerColor(kRed + 1);
  hPlus->SetMarkerStyle(20);
  hPlus->SetLineWidth(2);

  hMinus->SetLineColor(kBlue + 1);
  hMinus->SetMarkerColor(kBlue + 1);
  hMinus->SetMarkerStyle(21);
  hMinus->SetLineWidth(2);

  hPlus->SetMinimum(0.0);
  hPlus->SetMaximum(1.25*std::max(hPlus->GetMaximum(), hMinus->GetMaximum()));
  hPlus->Draw("HIST P");
  hMinus->Draw("HIST P SAME");

  TLegend* legend = new TLegend(0.72, 0.72, 0.92, 0.88);
  
  legend->AddEntry(hPlus, "+z stack", "lp");
  legend->AddEntry(hMinus, "-z stack", "lp");
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();

  c3->SaveAs("Fig/current/mip_edep_per_mcp.png");

  TCanvas* c4 = new TCanvas("c4", "Edep distribution per MCP", 1000, 700);
  c4->SetLeftMargin(0.12);
  c4->SetBottomMargin(0.12);
  c4->SetGridy();

  const Color_t colors[10] = {
    kRed + 1, kOrange + 7, kBlue + 1, kGreen + 2, kMagenta + 1,
    kCyan + 2, kViolet + 1, kPink + 6, kAzure + 7, kGray + 2
  };
  std::vector<TH1D*> plateHists;
  double largest = 0.0;
  for (int i = 0; i < numberOfMCPs; ++i) {
    std::ostringstream branchName;
    branchName << "edepPlusMcp" << i << "_keV";

    std::ostringstream histName;
    histName << "hPlusMcp" << i;
    TH1D* hist = new TH1D(histName.str().c_str(),
                          "",
                          kNumberOfBins,
                          0.0,
                          maxTotalMeV);
    hist->SetDirectory(nullptr);
    hist->SetLineColor(colors[i % 10]);
    hist->SetLineWidth(2);
    FillHistogramFromBranch(tree, branchName.str(), hist, 1.0/1000.0);
    if (hist->Integral() > 0.0) {
      hist->Scale(1.0/hist->Integral());
      largest = std::max(largest, hist->GetMaximum());
    }
    plateHists.push_back(hist);
  }

  TH1D* frame = new TH1D("mipEdepPlusFrame",
                         "Normalized MIP E_{dep} per +z MCP;E_{dep} [MeV];Normalized entries",
                         kNumberOfBins,
                         0.0,
                         maxTotalMeV);
  frame->SetDirectory(nullptr);
  frame->SetMinimum(0.0);
  frame->SetMaximum(largest > 0.0 ? 1.25*largest : 1.0);
  frame->Draw("HIST");

  TLegend* plateLegend = new TLegend(0.70, 0.62, 0.92, 0.88);
  plateLegend->SetBorderSize(0);
  plateLegend->SetFillStyle(0);
  for (int i = 0; i < numberOfMCPs; ++i) {
    plateHists[i]->Draw("HIST SAME");
    std::ostringstream label;
    label << "+z MCP " << i;
    plateLegend->AddEntry(plateHists[i], label.str().c_str(), "l");
  }
  plateLegend->Draw();
  c4->SaveAs("Fig/current/mip_edep_plus_mcp_shapes.png");

  std::cout << "\n=== MIP deposited energy summary ===" << std::endl;
  std::cout << "Events: " << tree->GetEntries() << std::endl;
  std::cout << "Events with Edep > 0: "
            << totalStats.nonZeroEntries << std::endl;
  std::cout << "Mean total Edep = "
            << MeanMeV(totalStats)
            << " MeV/event, RMS = " << RmsMeV(totalStats)
            << " MeV" << std::endl;
  std::cout << "Mean +z Edep    = "
            << MeanMeV(plusStats)
            << " MeV/event" << std::endl;
  std::cout << "Mean -z Edep    = "
            << MeanMeV(minusStats)
            << " MeV/event" << std::endl;

  std::cout << "\nPer-plate mean deposited energy:" << std::endl;
  std::cout << std::fixed << std::setprecision(5);
  std::cout
    << std::setw(6) << "side"
    << " | " << std::setw(8) << "mcpIndex"
    << " | " << std::setw(13) << "mean [MeV]"
    << " | " << std::setw(12) << "RMS [MeV]"
    << " | " << std::setw(12) << "nonzero"
    << std::endl;
  for (int i = 0; i < numberOfMCPs; ++i) {
    std::cout
      << std::setw(6) << "+z"
      << " | " << std::setw(8) << i
      << " | " << std::setw(13) << MeanMeV(plusByMcp[i])
      << " | " << std::setw(12) << RmsMeV(plusByMcp[i])
      << " | " << std::setw(12) << plusByMcp[i].nonZeroEntries
      << std::endl;
    std::cout
      << std::setw(6) << "-z"
      << " | " << std::setw(8) << i
      << " | " << std::setw(13) << MeanMeV(minusByMcp[i])
      << " | " << std::setw(12) << RmsMeV(minusByMcp[i])
      << " | " << std::setw(12) << minusByMcp[i].nonZeroEntries
      << std::endl;
  }

  std::cout << "\nSaved:" << std::endl;
  std::cout << "  mip_edep_total.png" << std::endl;
  std::cout << "  mip_edep_per_side.png" << std::endl;
  std::cout << "  mip_edep_per_mcp.png" << std::endl;
  std::cout << "  mip_edep_plus_mcp_shapes.png" << std::endl;

  file->Close();
  delete file;
}
