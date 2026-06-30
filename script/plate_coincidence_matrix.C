#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2I.h"
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
  struct EventPlateHits
  {
    std::set<Int_t> minusMcpIndices;
    std::set<Int_t> plusMcpIndices;
  };

  std::string McpLabel(const char* side, Int_t mcpIndex)
  {
    std::ostringstream label;
    label << side << " MCP " << mcpIndex;
    return label.str();
  }
}

void plate_coincidence_matrix(
  const char* fileName = "build/mcp_output.root")
{
  gSystem->mkdir("Fig/current", kTRUE);

  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* plateTree = nullptr;
  TTree* eventTree = nullptr;
  file->GetObject("McpPlateStatsTree", plateTree);
  file->GetObject("EventSummaryTree", eventTree);

  if (!plateTree || !eventTree) {
    std::cerr << "Missing McpPlateStatsTree or EventSummaryTree in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t eventID = -1;
  Int_t side = 0;
  Int_t mcpIndex = -1;
  Int_t electronChannelCount = 0;

  plateTree->SetBranchAddress("eventID", &eventID);
  plateTree->SetBranchAddress("side", &side);
  plateTree->SetBranchAddress("mcpIndex", &mcpIndex);
  plateTree->SetBranchAddress("electronChannelCount",
                              &electronChannelCount);

  std::map<Int_t, EventPlateHits> hitsByEvent;
  std::set<Int_t> minusMcpIndices;
  std::set<Int_t> plusMcpIndices;

  const Long64_t plateEntries = plateTree->GetEntries();
  for (Long64_t entry = 0; entry < plateEntries; ++entry) {
    plateTree->GetEntry(entry);

    if (side < 0) {
      minusMcpIndices.insert(mcpIndex);
      if (electronChannelCount > 0) {
        hitsByEvent[eventID].minusMcpIndices.insert(mcpIndex);
      } else {
        hitsByEvent[eventID];
      }
    } else if (side > 0) {
      plusMcpIndices.insert(mcpIndex);
      if (electronChannelCount > 0) {
        hitsByEvent[eventID].plusMcpIndices.insert(mcpIndex);
      } else {
        hitsByEvent[eventID];
      }
    } else {
      std::cerr << "Ignoring invalid side=0 for event "
                << eventID << ", MCP " << mcpIndex << std::endl;
    }
  }

  if (minusMcpIndices.empty() || plusMcpIndices.empty()) {
    std::cerr << "No MCP indices found on one or both detector sides."
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  const std::vector<Int_t> minusIndices(minusMcpIndices.begin(),
                                         minusMcpIndices.end());
  const std::vector<Int_t> plusIndices(plusMcpIndices.begin(),
                                       plusMcpIndices.end());

  std::map<std::pair<Int_t, Int_t>, Long64_t> coincidenceMatrix;
  std::set<Int_t> plateCoincidenceEvents;

  for (std::map<Int_t, EventPlateHits>::const_iterator eventIt =
         hitsByEvent.begin();
       eventIt != hitsByEvent.end(); ++eventIt) {
    const EventPlateHits& hits = eventIt->second;
    if (hits.minusMcpIndices.empty() ||
        hits.plusMcpIndices.empty()) {
      continue;
    }

    plateCoincidenceEvents.insert(eventIt->first);
    for (std::set<Int_t>::const_iterator minusIt =
           hits.minusMcpIndices.begin();
         minusIt != hits.minusMcpIndices.end(); ++minusIt) {
      for (std::set<Int_t>::const_iterator plusIt =
             hits.plusMcpIndices.begin();
           plusIt != hits.plusMcpIndices.end(); ++plusIt) {
        ++coincidenceMatrix[std::make_pair(*minusIt, *plusIt)];
      }
    }
  }

  Int_t summaryEventID = -1;
  Bool_t isCoincidence = false;
  eventTree->SetBranchAddress("eventID", &summaryEventID);
  eventTree->SetBranchAddress("isCoincidence", &isCoincidence);

  std::set<Int_t> summaryCoincidenceEvents;
  const Long64_t eventEntries = eventTree->GetEntries();
  for (Long64_t entry = 0; entry < eventEntries; ++entry) {
    eventTree->GetEntry(entry);
    if (isCoincidence) {
      summaryCoincidenceEvents.insert(summaryEventID);
    }
  }

  std::cout << "File: " << fileName << std::endl;
  std::cout << "\nPlate-to-plate coincidence matrix" << std::endl;
  std::cout << "Rows: -z MCP, columns: +z MCP" << std::endl;

  const Int_t cellWidth = 12;
  std::cout << std::setw(cellWidth) << "";
  for (std::vector<Int_t>::const_iterator plusIt =
         plusIndices.begin();
       plusIt != plusIndices.end(); ++plusIt) {
    std::ostringstream header;
    header << "+z MCP " << *plusIt;
    std::cout << std::setw(cellWidth) << header.str();
  }
  std::cout << std::endl;

  for (std::vector<Int_t>::const_iterator minusIt =
         minusIndices.begin();
       minusIt != minusIndices.end(); ++minusIt) {
    std::ostringstream rowLabel;
    rowLabel << "-z MCP " << *minusIt;
    std::cout << std::setw(cellWidth) << rowLabel.str();

    for (std::vector<Int_t>::const_iterator plusIt =
           plusIndices.begin();
         plusIt != plusIndices.end(); ++plusIt) {
      std::cout << std::setw(cellWidth)
                << coincidenceMatrix[
                     std::make_pair(*minusIt, *plusIt)];
    }
    std::cout << std::endl;
  }

  std::cout << "\nCoincidences by plate pair:" << std::endl;
  for (std::vector<Int_t>::const_iterator minusIt =
         minusIndices.begin();
       minusIt != minusIndices.end(); ++minusIt) {
    for (std::vector<Int_t>::const_iterator plusIt =
           plusIndices.begin();
         plusIt != plusIndices.end(); ++plusIt) {
      std::cout << "  (-z MCP " << *minusIt
                << ", +z MCP " << *plusIt << "): "
                << coincidenceMatrix[
                     std::make_pair(*minusIt, *plusIt)]
                << std::endl;
    }
  }

  std::cout << "\nEvents with at least one plate coincidence: "
            << plateCoincidenceEvents.size() << std::endl;
  std::cout << "Events with isCoincidence=true: "
            << summaryCoincidenceEvents.size() << std::endl;

  const bool eventSetsMatch =
    plateCoincidenceEvents == summaryCoincidenceEvents;
  std::cout << "Event-by-event consistency: "
            << (eventSetsMatch ? "[OK]" : "[MISMATCH]")
            << std::endl;

  if (!eventSetsMatch) {
    std::cout << "  Only in McpPlateStatsTree:";
    for (std::set<Int_t>::const_iterator it =
           plateCoincidenceEvents.begin();
         it != plateCoincidenceEvents.end(); ++it) {
      if (summaryCoincidenceEvents.count(*it) == 0) {
        std::cout << " " << *it;
      }
    }
    std::cout << std::endl;

    std::cout << "  Only in EventSummaryTree:";
    for (std::set<Int_t>::const_iterator it =
           summaryCoincidenceEvents.begin();
         it != summaryCoincidenceEvents.end(); ++it) {
      if (plateCoincidenceEvents.count(*it) == 0) {
        std::cout << " " << *it;
      }
    }
    std::cout << std::endl;
  }

  TH2I* histogram =
    new TH2I("hMcpPlateCoincidence",
             "MCP plate coincidence matrix;+z MCP;-z MCP",
             static_cast<Int_t>(plusIndices.size()),
             0.0,
             static_cast<Double_t>(plusIndices.size()),
             static_cast<Int_t>(minusIndices.size()),
             0.0,
             static_cast<Double_t>(minusIndices.size()));

  for (std::size_t plusBin = 0;
       plusBin < plusIndices.size();
       ++plusBin) {
    const std::string label =
      McpLabel("+z", plusIndices[plusBin]);
    histogram->GetXaxis()->SetBinLabel(
      static_cast<Int_t>(plusBin + 1), label.c_str());
  }

  for (std::size_t minusBin = 0;
       minusBin < minusIndices.size();
       ++minusBin) {
    const std::string label =
      McpLabel("-z", minusIndices[minusBin]);
    histogram->GetYaxis()->SetBinLabel(
      static_cast<Int_t>(minusBin + 1), label.c_str());

    for (std::size_t plusBin = 0;
         plusBin < plusIndices.size();
         ++plusBin) {
      histogram->SetBinContent(
        static_cast<Int_t>(plusBin + 1),
        static_cast<Int_t>(minusBin + 1),
        coincidenceMatrix[
          std::make_pair(minusIndices[minusBin],
                         plusIndices[plusBin])]);
    }
  }

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("g");
  TCanvas* canvas =
    new TCanvas("cMcpPlateCoincidence",
                "MCP plate coincidences",
                900,
                750);
  canvas->SetLeftMargin(0.15);
  canvas->SetBottomMargin(0.15);
  histogram->SetMarkerSize(1.5);
  histogram->Draw("COLZ TEXT");
  canvas->SaveAs("Fig/current/mcp_plate_coincidence_matrix.png");

  std::cout << "Saved: mcp_plate_coincidence_matrix.png" << std::endl;

  file->Close();
  delete file;
}
