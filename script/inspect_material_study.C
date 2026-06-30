#include "TFile.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>

void inspect_material_study(
  const char* fileName = "build/mcp_output.root")
{
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
              << "GammaInteractionTree in "
              << fileName << std::endl;
    std::cerr << "Generate a new ROOT file with the updated Geant4 "
              << "executable." << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t electronProducedCount = 0;
  Int_t electronChannelCount = 0;
  Int_t electronChannelPlusCount = 0;
  Int_t photonExitPlusCount = 0;
  Int_t gammaInteractionCount = 0;
  Int_t gammaPhotCount = 0;
  Int_t gammaComptCount = 0;
  Int_t gammaRaylCount = 0;
  Int_t gammaConvCount = 0;

  eventTree->SetBranchAddress("electronProducedCount",
                              &electronProducedCount);
  eventTree->SetBranchAddress("electronChannelCount",
                              &electronChannelCount);
  eventTree->SetBranchAddress("electronChannelPlusCount",
                              &electronChannelPlusCount);
  eventTree->SetBranchAddress("photonExitPlusCount",
                              &photonExitPlusCount);
  eventTree->SetBranchAddress("gammaInteractionCount",
                              &gammaInteractionCount);
  eventTree->SetBranchAddress("gammaPhotCount",
                              &gammaPhotCount);
  eventTree->SetBranchAddress("gammaComptCount",
                              &gammaComptCount);
  eventTree->SetBranchAddress("gammaRaylCount",
                              &gammaRaylCount);
  eventTree->SetBranchAddress("gammaConvCount",
                              &gammaConvCount);

  const Long64_t eventCount = eventTree->GetEntries();
  Long64_t totalElectronsProduced = 0;
  Long64_t totalElectronsInChannels = 0;
  Long64_t totalElectronsInPlusChannels = 0;
  Long64_t totalPhotonsExitingPlus = 0;
  Long64_t totalGammaInteractions = 0;
  Long64_t totalGammaPhot = 0;
  Long64_t totalGammaCompt = 0;
  Long64_t totalGammaRayl = 0;
  Long64_t totalGammaConv = 0;
  Long64_t eventsWithProducedElectron = 0;
  Long64_t eventsWithChannelElectron = 0;
  Long64_t eventsWithPhotonExit = 0;
  Long64_t eventsWithGammaInteraction = 0;
  std::map<Int_t, Long64_t> gammaMultiplicityCounts;

  for (Long64_t entry = 0; entry < eventCount; ++entry) {
    eventTree->GetEntry(entry);

    totalElectronsProduced += electronProducedCount;
    totalElectronsInChannels += electronChannelCount;
    totalElectronsInPlusChannels += electronChannelPlusCount;
    totalPhotonsExitingPlus += photonExitPlusCount;
    totalGammaInteractions += gammaInteractionCount;
    totalGammaPhot += gammaPhotCount;
    totalGammaCompt += gammaComptCount;
    totalGammaRayl += gammaRaylCount;
    totalGammaConv += gammaConvCount;

    if (electronProducedCount > 0) {
      ++eventsWithProducedElectron;
    }
    if (electronChannelCount > 0) {
      ++eventsWithChannelElectron;
    }
    if (photonExitPlusCount > 0) {
      ++eventsWithPhotonExit;
    }
    if (gammaInteractionCount > 0) {
      ++eventsWithGammaInteraction;
    }
    ++gammaMultiplicityCounts[gammaInteractionCount];
  }

  std::cout << "File: " << fileName << std::endl;
  std::cout << "Events: " << eventCount << std::endl;

  if (eventCount > 0) {
    const Double_t events =
      static_cast<Double_t>(eventCount);

    std::cout << std::fixed << std::setprecision(6);
    // std::cout << "\nMean electronProducedCount: "
    //           << totalElectronsProduced/events << std::endl;
    // std::cout << "Mean electronChannelCount: "
    //           << totalElectronsInChannels/events << std::endl;
    // std::cout << "Mean electronChannelPlusCount: "
    //           << totalElectronsInPlusChannels/events << std::endl;
    // std::cout << "Mean photonExitPlusCount: "
    //           << totalPhotonsExitingPlus/events << std::endl;
    // std::cout << "Mean gammaInteractionCount: "
    //           << totalGammaInteractions/events << std::endl;

    std::cout << "\nEvents with produced electrons: "
              << eventsWithProducedElectron << " / " << eventCount
              << " (" << eventsWithProducedElectron/events << ")"
              << std::endl;
    std::cout << "Events with at least one channel electron: "
              << eventsWithChannelElectron << " / " << eventCount
              << " (" << eventsWithChannelElectron/events << ", "
              << 100.0*eventsWithChannelElectron/events << " %)"
              << std::endl;
    std::cout << "Events with a photon exiting the +z MCP: "
              << eventsWithPhotonExit << " / " << eventCount
              << " (" << eventsWithPhotonExit/events << ")"
              << std::endl;
    std::cout << "Events with at least one gamma interaction: "
              << eventsWithGammaInteraction << " / " << eventCount
              << " (" << eventsWithGammaInteraction/events << ", "
              << 100.0*eventsWithGammaInteraction/events << " %)"
              << std::endl;
  }

  std::cout << "\nGamma interactions from EventSummaryTree:"
            << std::endl;
  std::cout << "  total: " << totalGammaInteractions << std::endl;
  std::cout << "  phot: " << totalGammaPhot << std::endl;
  std::cout << "  compt: " << totalGammaCompt << std::endl;
  std::cout << "  Rayl: " << totalGammaRayl << std::endl;
  std::cout << "  conv: " << totalGammaConv << std::endl;

  std::cout << "  multiplicity per event:" << std::endl;
  for (std::map<Int_t, Long64_t>::const_iterator it =
         gammaMultiplicityCounts.begin();
       it != gammaMultiplicityCounts.end(); ++it) {
    std::cout << "    " << it->first
              << " interactions: " << it->second
              << " events" << std::endl;
  }

  Int_t plateSide = 0;
  Int_t plateMcpIndex = -1;
  Int_t plateElectronChannelCount = 0;

  plateTree->SetBranchAddress("side", &plateSide);
  plateTree->SetBranchAddress("mcpIndex", &plateMcpIndex);
  plateTree->SetBranchAddress("electronChannelCount",
                              &plateElectronChannelCount);

  std::map<std::pair<Int_t, Int_t>, Long64_t> plateElectronTotals;
  std::map<std::pair<Int_t, Int_t>, Long64_t> plateActiveEventCounts;

  const Long64_t plateEntries = plateTree->GetEntries();
  for (Long64_t entry = 0; entry < plateEntries; ++entry) {
    plateTree->GetEntry(entry);

    const std::pair<Int_t, Int_t> plateKey =
      std::make_pair(plateSide, plateMcpIndex);
    plateElectronTotals[plateKey] += plateElectronChannelCount;
    if (plateElectronChannelCount > 0) {
      ++plateActiveEventCounts[plateKey];
    }
  }

  std::cout << "\nMCP plate statistics:" << std::endl;
  for (std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
         plateElectronTotals.begin();
       it != plateElectronTotals.end(); ++it) {
    std::cout << "  side=" << it->first.first
              << " mcpIndex=" << it->first.second
              << ": total electrons=" << it->second
              << ", events with electron="
              << plateActiveEventCounts[it->first]
              << std::endl;
  }

  Int_t gammaSide = 0;
  Int_t gammaMcpIndex = -1;
  char gammaProcessName[64] = "";

  gammaTree->SetBranchAddress("side", &gammaSide);
  gammaTree->SetBranchAddress("mcpIndex", &gammaMcpIndex);
  gammaTree->SetBranchAddress("processName", gammaProcessName);

  std::map<std::pair<Int_t, Int_t>, Long64_t>
    gammaInteractionTotalsByPlate;
  std::map<std::pair<Int_t, Int_t>,
           std::map<std::string, Long64_t> >
    gammaProcessTotalsByPlate;

  const Long64_t gammaEntries = gammaTree->GetEntries();
  for (Long64_t entry = 0; entry < gammaEntries; ++entry) {
    gammaTree->GetEntry(entry);

    const std::pair<Int_t, Int_t> plateKey =
      std::make_pair(gammaSide, gammaMcpIndex);
    ++gammaInteractionTotalsByPlate[plateKey];
    ++gammaProcessTotalsByPlate[plateKey][gammaProcessName];
  }

  std::cout << "\nGamma interactions by MCP plate:" << std::endl;
  for (std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
         gammaInteractionTotalsByPlate.begin();
       it != gammaInteractionTotalsByPlate.end(); ++it) {
    const std::map<std::string, Long64_t>& processCounts =
      gammaProcessTotalsByPlate[it->first];

    std::cout << "  side=" << it->first.first
              << " mcpIndex=" << it->first.second
              << ": total=" << it->second;

    const char* processNames[4] = {"phot", "compt", "Rayl", "conv"};
    for (Int_t processIndex = 0; processIndex < 4; ++processIndex) {
      const std::string processName = processNames[processIndex];
      const std::map<std::string, Long64_t>::const_iterator found =
        processCounts.find(processName);
      const Long64_t count =
        found != processCounts.end() ? found->second : 0;
      std::cout << ", " << processName << "=" << count;
    }
    std::cout << std::endl;
  }

  const Long64_t summaryProcessTotal =
    totalGammaPhot + totalGammaCompt +
    totalGammaRayl + totalGammaConv;

  std::cout << "\nGamma consistency checks:" << std::endl;
  std::cout << "  GammaInteractionTree rows=" << gammaEntries
            << ", EventSummaryTree total=" << totalGammaInteractions
            << (gammaEntries == totalGammaInteractions
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "  Process counter sum=" << summaryProcessTotal
            << ", EventSummaryTree total=" << totalGammaInteractions
            << (summaryProcessTotal == totalGammaInteractions
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;

  file->Close();
  delete file;
}
