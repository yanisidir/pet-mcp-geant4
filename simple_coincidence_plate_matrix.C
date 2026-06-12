#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  struct SimpleCoincidencePlates
  {
    std::set<Int_t> plusIndices;
    std::set<Int_t> minusIndices;
  };

  std::string PlateLabel(const char* side, Int_t index)
  {
    std::ostringstream label;
    label << side << " MCP " << index;
    return label.str();
  }
}

void simple_coincidence_plate_matrix(
  const char* fileName = "build/mcp_output_t0.root")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* eventTree = nullptr;
  TTree* plateTree = nullptr;
  file->GetObject("EventSummaryTree", eventTree);
  file->GetObject("McpPlateStatsTree", plateTree);

  if (!eventTree || !plateTree) {
    std::cerr << "Missing EventSummaryTree or McpPlateStatsTree in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t eventID = -1;
  Int_t plusCount = 0;
  Int_t minusCount = 0;

  // Modifier uniquement ces deux valeurs pour étudier un autre cas.
  const Int_t kRequiredPlusCount = 1;
  const Int_t kRequiredMinusCount = 1;

  std::ostringstream selectionText;
  selectionText << "N(+z)=" << kRequiredPlusCount
                << ", N(-z)=" << kRequiredMinusCount;

  std::ostringstream outputSuffix;
  outputSuffix << kRequiredPlusCount << "_" << kRequiredMinusCount;

  const std::string countFileName =
    "simple_coincidence_plate_matrix_counts_" +
    outputSuffix.str() + ".png";
  const std::string percentFileName =
    "simple_coincidence_plate_matrix_percent_" +
    outputSuffix.str() + ".png";
  const std::string countTitle =
    "Simple coincidence plate pairs (" +
    selectionText.str() + ");+z MCP;-z MCP";
  const std::string percentTitle =
    "Simple coincidence plate pairs (" +
    selectionText.str() + ") (%);+z MCP;-z MCP";

  eventTree->SetBranchAddress("eventID", &eventID);
  eventTree->SetBranchAddress("electronChannelPlusCount", &plusCount);
  eventTree->SetBranchAddress("electronChannelMinusCount", &minusCount);

  std::set<Int_t> selectedEventIDs;
  const Long64_t eventEntries = eventTree->GetEntries();
  for (Long64_t entry = 0; entry < eventEntries; ++entry) {
    eventTree->GetEntry(entry);
    if (plusCount == kRequiredPlusCount &&
        minusCount == kRequiredMinusCount) {
      selectedEventIDs.insert(eventID);
    }
  }

  Int_t plateEventID = -1;
  Int_t side = 0;
  Int_t mcpIndex = -1;
  Int_t electronChannelCount = 0;

  plateTree->SetBranchAddress("eventID", &plateEventID);
  plateTree->SetBranchAddress("side", &side);
  plateTree->SetBranchAddress("mcpIndex", &mcpIndex);
  plateTree->SetBranchAddress("electronChannelCount",
                              &electronChannelCount);

  std::map<Int_t, SimpleCoincidencePlates> platesByEvent;
  std::set<Int_t> allPlusIndices;
  std::set<Int_t> allMinusIndices;

  const Long64_t plateEntries = plateTree->GetEntries();
  for (Long64_t entry = 0; entry < plateEntries; ++entry) {
    plateTree->GetEntry(entry);

    if (side > 0) {
      allPlusIndices.insert(mcpIndex);
    } else if (side < 0) {
      allMinusIndices.insert(mcpIndex);
    }

    if (selectedEventIDs.count(plateEventID) == 0 ||
        electronChannelCount <= 0) {
      continue;
    }

    if (side > 0) {
      platesByEvent[plateEventID].plusIndices.insert(mcpIndex);
    } else if (side < 0) {
      platesByEvent[plateEventID].minusIndices.insert(mcpIndex);
    }
  }

  if (allPlusIndices.empty() || allMinusIndices.empty()) {
    std::cerr << "No MCP plate indices found on one or both sides."
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  const std::vector<Int_t> plusIndices(allPlusIndices.begin(),
                                       allPlusIndices.end());
  const std::vector<Int_t> minusIndices(allMinusIndices.begin(),
                                        allMinusIndices.end());

  std::map<std::pair<Int_t, Int_t>, Long64_t> pairCounts;
  Long64_t validEventCount = 0;
  Long64_t invalidEventCount = 0;
  Long64_t totalPairContributions = 0;

  for (std::set<Int_t>::const_iterator eventIt =
         selectedEventIDs.begin();
       eventIt != selectedEventIDs.end(); ++eventIt) {
    const std::map<Int_t, SimpleCoincidencePlates>::const_iterator found =
      platesByEvent.find(*eventIt);

    if (found == platesByEvent.end() ||
        found->second.plusIndices.empty() ||
        found->second.minusIndices.empty()) {
      ++invalidEventCount;
      continue;
    }

    // Un même événement peut contribuer à plusieurs cases si plusieurs
    // plaques sont actives sur un côté ou sur les deux côtés.
    for (std::set<Int_t>::const_iterator minusIt =
           found->second.minusIndices.begin();
         minusIt != found->second.minusIndices.end(); ++minusIt) {
      for (std::set<Int_t>::const_iterator plusIt =
             found->second.plusIndices.begin();
           plusIt != found->second.plusIndices.end(); ++plusIt) {
        ++pairCounts[std::make_pair(*minusIt, *plusIt)];
        ++totalPairContributions;
      }
    }

    ++validEventCount;
  }

  TH2D* countHistogram =
    new TH2D("hSimpleCoincidencePlateCounts",
             countTitle.c_str(),
             static_cast<Int_t>(plusIndices.size()),
             0.0,
             static_cast<Double_t>(plusIndices.size()),
             static_cast<Int_t>(minusIndices.size()),
             0.0,
             static_cast<Double_t>(minusIndices.size()));

  for (std::size_t plusBin = 0; plusBin < plusIndices.size(); ++plusBin) {
    const std::string label = PlateLabel("+z", plusIndices[plusBin]);
    countHistogram->GetXaxis()->SetBinLabel(
      static_cast<Int_t>(plusBin + 1), label.c_str());
  }

  for (std::size_t minusBin = 0;
       minusBin < minusIndices.size();
       ++minusBin) {
    const std::string label = PlateLabel("-z", minusIndices[minusBin]);
    countHistogram->GetYaxis()->SetBinLabel(
      static_cast<Int_t>(minusBin + 1), label.c_str());

    for (std::size_t plusBin = 0;
         plusBin < plusIndices.size();
         ++plusBin) {
      countHistogram->SetBinContent(
        static_cast<Int_t>(plusBin + 1),
        static_cast<Int_t>(minusBin + 1),
        pairCounts[
          std::make_pair(minusIndices[minusBin],
                         plusIndices[plusBin])]);
    }
  }

  TH2D* percentHistogram =
    static_cast<TH2D*>(
      countHistogram->Clone("hSimpleCoincidencePlatePercent"));
  percentHistogram->SetTitle(
    percentTitle.c_str());
  if (validEventCount > 0) {
    percentHistogram->Scale(
      100.0/static_cast<Double_t>(validEventCount));
  }

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("g");

  TCanvas* countCanvas =
    new TCanvas("cSimpleCoincidencePlateCounts",
                "Simple coincidence plate counts",
                900,
                750);
  countCanvas->SetLeftMargin(0.15);
  countCanvas->SetBottomMargin(0.15);
  countCanvas->SetRightMargin(0.15);
  countHistogram->SetMarkerSize(1.5);
  countHistogram->Draw("COLZ TEXT");
  countCanvas->SaveAs(countFileName.c_str());

  gStyle->SetPaintTextFormat("4.1f");
  TCanvas* percentCanvas =
    new TCanvas("cSimpleCoincidencePlatePercent",
                "Simple coincidence plate percent",
                900,
                750);
  percentCanvas->SetLeftMargin(0.15);
  percentCanvas->SetBottomMargin(0.15);
  percentCanvas->SetRightMargin(0.15);
  percentHistogram->SetMarkerSize(1.5);
  percentHistogram->Draw("COLZ TEXT");
  percentCanvas->SaveAs(percentFileName.c_str());

  std::cout << "File: " << fileName << std::endl;
  std::cout << "\nTotal events with "
            << selectionText.str() << ": "
            << selectedEventIDs.size() << std::endl;
  std::cout << "Valid events with at least one active plate per side: "
            << validEventCount << std::endl;
  std::cout << "Invalid events missing an active plate on one side: "
            << invalidEventCount << std::endl;
  std::cout << "Total plate-to-plate contributions: "
            << totalPairContributions << std::endl;
  std::cout << "Contributions/event: "
            << (validEventCount > 0
                  ? static_cast<Double_t>(totalPairContributions)/
                    validEventCount
                  : 0.0)
            << std::endl;

  std::cout << "\nContribution of each plate pair:" << std::endl;
  std::cout << std::fixed << std::setprecision(3);
  for (std::vector<Int_t>::const_iterator minusIt =
         minusIndices.begin();
       minusIt != minusIndices.end(); ++minusIt) {
    for (std::vector<Int_t>::const_iterator plusIt =
           plusIndices.begin();
         plusIt != plusIndices.end(); ++plusIt) {
      const Long64_t count =
        pairCounts[std::make_pair(*minusIt, *plusIt)];
      const Double_t percent =
        validEventCount > 0
          ? 100.0*static_cast<Double_t>(count)/validEventCount
          : 0.0;

      std::cout << "  (-z MCP " << *minusIt
                << ", +z MCP " << *plusIt << "): "
                << count << " events, "
                << percent << " %" << std::endl;
    }
  }

  std::cout << "\nSaved: " << countFileName << std::endl;
  std::cout << "Saved: " << percentFileName << std::endl;

  file->Close();
  delete file;
}
