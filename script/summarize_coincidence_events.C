#include "TCanvas.h"
#include "TSystem.h"
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
#include <utility>
#include <vector>

namespace
{
  struct CoincidencePlateHits
  {
    std::set<Int_t> plusIndices;
    std::set<Int_t> minusIndices;
  };

  std::string CoincidencePlateLabel(const char* side, Int_t mcpIndex)
  {
    std::ostringstream label;
    label << side << " MCP " << mcpIndex;
    return label.str();
  }
}

void summarize_coincidence_events(
  const char* fileName = "build/mcp_output.root")
{
  gSystem->mkdir("Fig/current", kTRUE);

  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* eventTree = nullptr;
  TTree* plateTree = nullptr;
  TTree* gammaTree = nullptr;
  file->GetObject("EventSummaryTree", eventTree);
  file->GetObject("McpPlateStatsTree", plateTree);
  file->GetObject("GammaInteractionTree", gammaTree);

  if (!eventTree || !plateTree || !gammaTree) {
    std::cerr << "Missing EventSummaryTree, McpPlateStatsTree or "
              << "GammaInteractionTree in " << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t eventID = -1;
  Bool_t isCoincidence = false;
  eventTree->SetBranchAddress("eventID", &eventID);
  eventTree->SetBranchAddress("isCoincidence", &isCoincidence);

  std::set<Int_t> coincidenceEventIDs;
  const Long64_t totalEventCount = eventTree->GetEntries();
  for (Long64_t entry = 0; entry < totalEventCount; ++entry) {
    eventTree->GetEntry(entry);
    if (isCoincidence) {
      coincidenceEventIDs.insert(eventID);
    }
  }

  const Long64_t coincidenceEventCount =
    static_cast<Long64_t>(coincidenceEventIDs.size());
  const Double_t coincidenceEfficiency =
    totalEventCount > 0
      ? static_cast<Double_t>(coincidenceEventCount)/totalEventCount
      : 0.0;

  Int_t plateEventID = -1;
  Int_t side = 0;
  Int_t mcpIndex = -1;
  Int_t electronChannelCount = 0;
  plateTree->SetBranchAddress("eventID", &plateEventID);
  plateTree->SetBranchAddress("side", &side);
  plateTree->SetBranchAddress("mcpIndex", &mcpIndex);
  plateTree->SetBranchAddress("electronChannelCount",
                              &electronChannelCount);

  std::set<Int_t> allPlusIndices;
  std::set<Int_t> allMinusIndices;
  std::map<Int_t, CoincidencePlateHits> platesByEvent;
  std::map<Int_t, Long64_t> plusEventsByPlate;
  std::map<Int_t, Long64_t> minusEventsByPlate;

  const Long64_t plateEntries = plateTree->GetEntries();
  for (Long64_t entry = 0; entry < plateEntries; ++entry) {
    plateTree->GetEntry(entry);

    if (side > 0) {
      allPlusIndices.insert(mcpIndex);
    } else if (side < 0) {
      allMinusIndices.insert(mcpIndex);
    }

    if (coincidenceEventIDs.count(plateEventID) == 0 ||
        electronChannelCount <= 0) {
      continue;
    }

    if (side > 0) {
      platesByEvent[plateEventID].plusIndices.insert(mcpIndex);
    } else if (side < 0) {
      platesByEvent[plateEventID].minusIndices.insert(mcpIndex);
    }
  }

  std::map<std::pair<Int_t, Int_t>, Long64_t> platePairCounts;
  Long64_t validPlateCoincidenceEvents = 0;
  Long64_t totalPlatePairContributions = 0;

  for (std::set<Int_t>::const_iterator eventIt =
         coincidenceEventIDs.begin();
       eventIt != coincidenceEventIDs.end(); ++eventIt) {
    const std::map<Int_t, CoincidencePlateHits>::const_iterator found =
      platesByEvent.find(*eventIt);
    if (found == platesByEvent.end() ||
        found->second.plusIndices.empty() ||
        found->second.minusIndices.empty()) {
      continue;
    }

    ++validPlateCoincidenceEvents;

    for (std::set<Int_t>::const_iterator plusIt =
           found->second.plusIndices.begin();
         plusIt != found->second.plusIndices.end(); ++plusIt) {
      ++plusEventsByPlate[*plusIt];
    }
    for (std::set<Int_t>::const_iterator minusIt =
           found->second.minusIndices.begin();
         minusIt != found->second.minusIndices.end(); ++minusIt) {
      ++minusEventsByPlate[*minusIt];
    }

    // Un événement peut contribuer à plusieurs cases si plusieurs plaques
    // sont actives sur un côté ou sur les deux côtés.
    for (std::set<Int_t>::const_iterator minusIt =
           found->second.minusIndices.begin();
         minusIt != found->second.minusIndices.end(); ++minusIt) {
      for (std::set<Int_t>::const_iterator plusIt =
             found->second.plusIndices.begin();
           plusIt != found->second.plusIndices.end(); ++plusIt) {
        ++platePairCounts[std::make_pair(*minusIt, *plusIt)];
        ++totalPlatePairContributions;
      }
    }
  }

  Int_t gammaEventID = -1;
  char processName[64] = "";
  gammaTree->SetBranchAddress("eventID", &gammaEventID);
  gammaTree->SetBranchAddress("processName", processName);

  std::map<std::string, Long64_t> gammaInteractionCounts;
  std::set<Int_t> coincidenceEventsWithPhot;
  std::set<Int_t> coincidenceEventsWithCompt;

  const Long64_t gammaEntries = gammaTree->GetEntries();
  for (Long64_t entry = 0; entry < gammaEntries; ++entry) {
    gammaTree->GetEntry(entry);
    if (coincidenceEventIDs.count(gammaEventID) == 0) {
      continue;
    }

    const std::string process(processName);
    ++gammaInteractionCounts[process];
    if (process == "phot") {
      coincidenceEventsWithPhot.insert(gammaEventID);
    } else if (process == "compt") {
      coincidenceEventsWithCompt.insert(gammaEventID);
    }
  }

  std::cout << "File: " << fileName << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\n=== Coincidence summary ===" << std::endl;
  std::cout << "Total events: " << totalEventCount << std::endl;
  std::cout << "Coincidence events: "
            << coincidenceEventCount << std::endl;
  std::cout << "Coincidence efficiency: "
            << coincidenceEfficiency << " ("
            << 100.0*coincidenceEfficiency << " %)"
            << std::endl;

  std::cout << "\nCoincidence events involving each +z plate:"
            << std::endl;
  for (std::set<Int_t>::const_iterator it =
         allPlusIndices.begin();
       it != allPlusIndices.end(); ++it) {
    std::cout << "  +z MCP " << *it << ": "
              << plusEventsByPlate[*it] << std::endl;
  }

  std::cout << "\nCoincidence events involving each -z plate:"
            << std::endl;
  for (std::set<Int_t>::const_iterator it =
         allMinusIndices.begin();
       it != allMinusIndices.end(); ++it) {
    std::cout << "  -z MCP " << *it << ": "
              << minusEventsByPlate[*it] << std::endl;
  }

  std::cout << "\nPlate-to-plate matrix "
            << "(rows=-z MCP, columns=+z MCP):" << std::endl;
  const Int_t columnWidth = 12;
  std::cout << std::setw(columnWidth) << "";
  for (std::set<Int_t>::const_iterator plusIt =
         allPlusIndices.begin();
       plusIt != allPlusIndices.end(); ++plusIt) {
    std::ostringstream label;
    label << "+z MCP " << *plusIt;
    std::cout << std::setw(columnWidth) << label.str();
  }
  std::cout << std::endl;

  for (std::set<Int_t>::const_iterator minusIt =
         allMinusIndices.begin();
       minusIt != allMinusIndices.end(); ++minusIt) {
    std::ostringstream label;
    label << "-z MCP " << *minusIt;
    std::cout << std::setw(columnWidth) << label.str();
    for (std::set<Int_t>::const_iterator plusIt =
           allPlusIndices.begin();
         plusIt != allPlusIndices.end(); ++plusIt) {
      std::cout << std::setw(columnWidth)
                << platePairCounts[
                     std::make_pair(*minusIt, *plusIt)];
    }
    std::cout << std::endl;
  }

  std::cout << "\nGamma interactions in coincidence events:"
            << std::endl;
  std::cout << "  phot: " << gammaInteractionCounts["phot"]
            << std::endl;
  std::cout << "  compt: " << gammaInteractionCounts["compt"]
            << std::endl;
  std::cout << "  Rayl: " << gammaInteractionCounts["Rayl"]
            << std::endl;
  std::cout << "  conv: " << gammaInteractionCounts["conv"]
            << std::endl;
  std::cout << "Coincidence events with at least one phot: "
            << coincidenceEventsWithPhot.size() << std::endl;
  std::cout << "Coincidence events with at least one compt: "
            << coincidenceEventsWithCompt.size() << std::endl;

  std::cout << "\nPlate consistency:" << std::endl;
  std::cout << "  Coincidence events with active plates on both sides: "
            << validPlateCoincidenceEvents << " / "
            << coincidenceEventCount
            << (validPlateCoincidenceEvents == coincidenceEventCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "  Total plate-pair contributions: "
            << totalPlatePairContributions << std::endl;

  if (allPlusIndices.empty() || allMinusIndices.empty()) {
    std::cerr << "Cannot create matrix: MCP indices are missing on "
              << "one or both sides." << std::endl;
    file->Close();
    delete file;
    return;
  }

  const std::vector<Int_t> plusIndices(allPlusIndices.begin(),
                                       allPlusIndices.end());
  const std::vector<Int_t> minusIndices(allMinusIndices.begin(),
                                        allMinusIndices.end());

  TH2D* matrix =
    new TH2D("hCoincidencePlateMatrixSummary",
             "Coincidence events by MCP plate pair;"
             "+z MCP plate;-z MCP plate",
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
      CoincidencePlateLabel("+z", plusIndices[plusBin]);
    matrix->GetXaxis()->SetBinLabel(
      static_cast<Int_t>(plusBin + 1), label.c_str());
  }

  for (std::size_t minusBin = 0;
       minusBin < minusIndices.size();
       ++minusBin) {
    const std::string label =
      CoincidencePlateLabel("-z", minusIndices[minusBin]);
    matrix->GetYaxis()->SetBinLabel(
      static_cast<Int_t>(minusBin + 1), label.c_str());

    for (std::size_t plusBin = 0;
         plusBin < plusIndices.size();
         ++plusBin) {
      matrix->SetBinContent(
        static_cast<Int_t>(plusBin + 1),
        static_cast<Int_t>(minusBin + 1),
        platePairCounts[
          std::make_pair(minusIndices[minusBin],
                         plusIndices[plusBin])]);
    }
  }

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("g");

  TCanvas* canvas =
    new TCanvas("cCoincidencePlateMatrixSummary",
                "Coincidence plate matrix summary",
                900,
                750);
  canvas->SetLeftMargin(0.16);
  canvas->SetBottomMargin(0.16);
  canvas->SetRightMargin(0.15);
  matrix->SetMarkerSize(1.5);
  matrix->Draw("COLZ TEXT");
  canvas->SaveAs("Fig/current/coincidence_plate_matrix_summary.png");

  std::cout << "\nSaved: coincidence_plate_matrix_summary.png"
            << std::endl;

  file->Close();
  delete file;
}
