#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>

void inspect_electrons(const char* fileName = "build/mcp_output_t0.root")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* electronTree = nullptr;
  TTree* photonTree = nullptr;
  TTree* eventTree = nullptr;

  file->GetObject("ElectronChannelHitTree", electronTree);
  file->GetObject("PhotonExitTree", photonTree);
  file->GetObject("EventSummaryTree", eventTree);

  if (!electronTree || !photonTree || !eventTree) {
    std::cerr << "Missing required tree in " << fileName << std::endl;
    std::cerr << "Expected: ElectronChannelHitTree, PhotonExitTree, "
              << "EventSummaryTree" << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << "File: " << fileName << std::endl;
  std::cout << std::fixed << std::setprecision(6);

  // --------------------------------------------------------------------------
  // Electrons reaching MCP_channel
  // --------------------------------------------------------------------------

  Int_t electronEventID = 0;
  Int_t electronTrackID = 0;
  Int_t electronParentID = 0;
  Int_t electronSide = 0;
  Int_t electronMcpIndex = -1;
  Double_t electronEnergy_keV = 0.0;
  Double_t electronTime_ns = 0.0;
  Double_t electronX_cm = 0.0;
  Double_t electronY_cm = 0.0;
  Double_t electronZ_cm = 0.0;
  Double_t electronDirX = 0.0;
  Double_t electronDirY = 0.0;
  Double_t electronDirZ = 0.0;
  char electronVolumeName[64] = "";
  char electronCreatorProcessName[64] = "";

  electronTree->SetBranchAddress("eventID", &electronEventID);
  electronTree->SetBranchAddress("trackID", &electronTrackID);
  electronTree->SetBranchAddress("parentID", &electronParentID);
  electronTree->SetBranchAddress("side", &electronSide);
  electronTree->SetBranchAddress("mcpIndex", &electronMcpIndex);
  electronTree->SetBranchAddress("kineticEnergy_keV", &electronEnergy_keV);
  electronTree->SetBranchAddress("globalTime_ns", &electronTime_ns);
  electronTree->SetBranchAddress("x_cm", &electronX_cm);
  electronTree->SetBranchAddress("y_cm", &electronY_cm);
  electronTree->SetBranchAddress("z_cm", &electronZ_cm);
  electronTree->SetBranchAddress("dirX", &electronDirX);
  electronTree->SetBranchAddress("dirY", &electronDirY);
  electronTree->SetBranchAddress("dirZ", &electronDirZ);
  electronTree->SetBranchAddress("volumeName", electronVolumeName);
  electronTree->SetBranchAddress("creatorProcessName",
                                 electronCreatorProcessName);

  const Long64_t electronEntries = electronTree->GetEntries();
  const Long64_t electronPreview =
    std::min<Long64_t>(electronEntries, 10);

  Double_t electronEnergySum = 0.0;
  Double_t electronEnergyMin = std::numeric_limits<Double_t>::max();
  Double_t electronEnergyMax = -std::numeric_limits<Double_t>::max();
  std::map<std::string, Long64_t> electronCreatorCounts;
  std::map<Int_t, Long64_t> electronSideCounts;
  std::map<std::pair<Int_t, Int_t>, Long64_t> electronMcpCounts;

  std::cout << "\nElectronChannelHitTree" << std::endl;
  std::cout << "  entries: " << electronEntries << std::endl;

  for (Long64_t i = 0; i < electronEntries; ++i) {
    electronTree->GetEntry(i);

    electronEnergySum += electronEnergy_keV;
    electronEnergyMin = std::min(electronEnergyMin, electronEnergy_keV);
    electronEnergyMax = std::max(electronEnergyMax, electronEnergy_keV);
    electronCreatorCounts[electronCreatorProcessName]++;
    electronSideCounts[electronSide]++;
    electronMcpCounts[std::make_pair(electronSide, electronMcpIndex)]++;

    if (i < electronPreview) {
      std::cout << "  [" << i << "]"
                << " event=" << electronEventID
                << " track=" << electronTrackID
                << " parent=" << electronParentID
                << " side=" << electronSide
                << " mcpIndex=" << electronMcpIndex
                << " E=" << electronEnergy_keV << " keV"
                << " t=" << electronTime_ns << " ns"
                << " pos=(" << electronX_cm << ", "
                << electronY_cm << ", " << electronZ_cm << ") cm"
                << " dir=(" << electronDirX << ", "
                << electronDirY << ", " << electronDirZ << ")"
                << " volume=" << electronVolumeName
                << " creator=" << electronCreatorProcessName
                << std::endl;
    }
  }

  if (electronEntries > 0) {
    std::cout << "  energy min/mean/max: "
              << electronEnergyMin << " / "
              << electronEnergySum/electronEntries << " / "
              << electronEnergyMax << " keV" << std::endl;
  }

  std::cout << "  creator processes:" << std::endl;
  for (std::map<std::string, Long64_t>::const_iterator it =
         electronCreatorCounts.begin();
       it != electronCreatorCounts.end(); ++it) {
    std::cout << "    " << it->first << ": " << it->second << std::endl;
  }
  std::cout << "  sides: +z=" << electronSideCounts[1]
            << ", -z=" << electronSideCounts[-1] << std::endl;
  std::cout << "  hits by MCP:" << std::endl;
  for (Int_t mcpIndex = 0; mcpIndex < 3; ++mcpIndex) {
    std::cout << "    MCP " << mcpIndex
              << ": +z="
              << electronMcpCounts[std::make_pair(1, mcpIndex)]
              << ", -z="
              << electronMcpCounts[std::make_pair(-1, mcpIndex)]
              << std::endl;
  }

  // --------------------------------------------------------------------------
  // Gamma photons leaving the MCP
  // --------------------------------------------------------------------------

  Int_t photonEventID = 0;
  Int_t photonTrackID = 0;
  Int_t photonParentID = 0;
  Int_t photonSide = 0;
  Double_t photonEnergy_keV = 0.0;
  Double_t photonTime_ns = 0.0;
  Double_t photonX_cm = 0.0;
  Double_t photonY_cm = 0.0;
  Double_t photonZ_cm = 0.0;
  Double_t photonDirX = 0.0;
  Double_t photonDirY = 0.0;
  Double_t photonDirZ = 0.0;
  char photonVolumeName[64] = "";
  char photonStepProcessName[64] = "";

  photonTree->SetBranchAddress("eventID", &photonEventID);
  photonTree->SetBranchAddress("trackID", &photonTrackID);
  photonTree->SetBranchAddress("parentID", &photonParentID);
  photonTree->SetBranchAddress("side", &photonSide);
  photonTree->SetBranchAddress("kineticEnergy_keV", &photonEnergy_keV);
  photonTree->SetBranchAddress("globalTime_ns", &photonTime_ns);
  photonTree->SetBranchAddress("x_cm", &photonX_cm);
  photonTree->SetBranchAddress("y_cm", &photonY_cm);
  photonTree->SetBranchAddress("z_cm", &photonZ_cm);
  photonTree->SetBranchAddress("dirX", &photonDirX);
  photonTree->SetBranchAddress("dirY", &photonDirY);
  photonTree->SetBranchAddress("dirZ", &photonDirZ);
  photonTree->SetBranchAddress("volumeName", photonVolumeName);
  photonTree->SetBranchAddress("stepProcessName", photonStepProcessName);

  const Long64_t photonEntries = photonTree->GetEntries();
  const Long64_t photonPreview = std::min<Long64_t>(photonEntries, 10);

  Double_t photonEnergySum = 0.0;
  Double_t photonEnergyMin = std::numeric_limits<Double_t>::max();
  Double_t photonEnergyMax = -std::numeric_limits<Double_t>::max();
  std::map<std::string, Long64_t> photonStepProcessCounts;
  std::map<Int_t, Long64_t> photonSideCounts;

  std::cout << "\nPhotonExitTree" << std::endl;
  std::cout << "  entries: " << photonEntries << std::endl;

  for (Long64_t i = 0; i < photonEntries; ++i) {
    photonTree->GetEntry(i);

    photonEnergySum += photonEnergy_keV;
    photonEnergyMin = std::min(photonEnergyMin, photonEnergy_keV);
    photonEnergyMax = std::max(photonEnergyMax, photonEnergy_keV);
    photonStepProcessCounts[photonStepProcessName]++;
    photonSideCounts[photonSide]++;

    if (i < photonPreview) {
      std::cout << "  [" << i << "]"
                << " event=" << photonEventID
                << " track=" << photonTrackID
                << " parent=" << photonParentID
                << " side=" << photonSide
                << " E=" << photonEnergy_keV << " keV"
                << " t=" << photonTime_ns << " ns"
                << " pos=(" << photonX_cm << ", "
                << photonY_cm << ", " << photonZ_cm << ") cm"
                << " dir=(" << photonDirX << ", "
                << photonDirY << ", " << photonDirZ << ")"
                << " exitedVolume=" << photonVolumeName
                << " stepProcess=" << photonStepProcessName
                << std::endl;
    }
  }

  if (photonEntries > 0) {
    std::cout << "  energy min/mean/max: "
              << photonEnergyMin << " / "
              << photonEnergySum/photonEntries << " / "
              << photonEnergyMax << " keV" << std::endl;
  }

  std::cout << "  boundary step processes:" << std::endl;
  for (std::map<std::string, Long64_t>::const_iterator it =
         photonStepProcessCounts.begin();
       it != photonStepProcessCounts.end(); ++it) {
    std::cout << "    " << it->first << ": " << it->second << std::endl;
  }
  std::cout << "  sides: +z=" << photonSideCounts[1]
            << ", -z=" << photonSideCounts[-1] << std::endl;

  // --------------------------------------------------------------------------
  // Event summary and consistency checks
  // --------------------------------------------------------------------------

  Int_t summaryEventID = 0;
  Int_t electronProducedCount = 0;
  Int_t electronChannelCount = 0;
  Int_t electronChannelPlusCount = 0;
  Int_t electronChannelMinusCount = 0;
  Int_t electronChannelPlusMcp0Count = 0;
  Int_t electronChannelPlusMcp1Count = 0;
  Int_t electronChannelPlusMcp2Count = 0;
  Int_t electronChannelMinusMcp0Count = 0;
  Int_t electronChannelMinusMcp1Count = 0;
  Int_t electronChannelMinusMcp2Count = 0;
  Bool_t hasElectronChannelPlus = false;
  Bool_t hasElectronChannelMinus = false;
  Bool_t isCoincidence = false;
  Int_t photonExitCount = 0;
  Int_t photonExitPlusCount = 0;
  Int_t photonExitMinusCount = 0;

  eventTree->SetBranchAddress("eventID", &summaryEventID);
  eventTree->SetBranchAddress("electronProducedCount",
                              &electronProducedCount);
  eventTree->SetBranchAddress("electronChannelCount",
                              &electronChannelCount);
  eventTree->SetBranchAddress("electronChannelPlusCount",
                              &electronChannelPlusCount);
  eventTree->SetBranchAddress("electronChannelMinusCount",
                              &electronChannelMinusCount);
  eventTree->SetBranchAddress("electronChannelPlusMcp0Count",
                              &electronChannelPlusMcp0Count);
  eventTree->SetBranchAddress("electronChannelPlusMcp1Count",
                              &electronChannelPlusMcp1Count);
  eventTree->SetBranchAddress("electronChannelPlusMcp2Count",
                              &electronChannelPlusMcp2Count);
  eventTree->SetBranchAddress("electronChannelMinusMcp0Count",
                              &electronChannelMinusMcp0Count);
  eventTree->SetBranchAddress("electronChannelMinusMcp1Count",
                              &electronChannelMinusMcp1Count);
  eventTree->SetBranchAddress("electronChannelMinusMcp2Count",
                              &electronChannelMinusMcp2Count);
  eventTree->SetBranchAddress("hasElectronChannelPlus",
                              &hasElectronChannelPlus);
  eventTree->SetBranchAddress("hasElectronChannelMinus",
                              &hasElectronChannelMinus);
  eventTree->SetBranchAddress("isCoincidence",
                              &isCoincidence);
  eventTree->SetBranchAddress("photonExitCount", &photonExitCount);
  eventTree->SetBranchAddress("photonExitPlusCount",
                              &photonExitPlusCount);
  eventTree->SetBranchAddress("photonExitMinusCount",
                              &photonExitMinusCount);

  const Long64_t eventEntries = eventTree->GetEntries();
  const Long64_t eventPreview = std::min<Long64_t>(eventEntries, 10);

  Long64_t totalElectronProducedCount = 0;
  Long64_t totalElectronChannelCount = 0;
  Long64_t totalElectronChannelPlusCount = 0;
  Long64_t totalElectronChannelMinusCount = 0;
  Long64_t totalElectronChannelPlusMcpCount[3] = {0, 0, 0};
  Long64_t totalElectronChannelMinusMcpCount[3] = {0, 0, 0};
  Long64_t totalPhotonExitCount = 0;
  Long64_t totalPhotonExitPlusCount = 0;
  Long64_t totalPhotonExitMinusCount = 0;
  Long64_t eventsWithProducedElectron = 0;
  Long64_t eventsWithElectron = 0;
  Long64_t eventsDetectedPlus = 0;
  Long64_t eventsDetectedMinus = 0;
  Long64_t coincidenceEvents = 0;
  Long64_t eventsWithPhotonExit = 0;
  std::map<Int_t, Long64_t> producedMultiplicityCounts;
  std::map<Int_t, Long64_t> multiplicityCounts;

  std::cout << "\nEventSummaryTree" << std::endl;
  std::cout << "  entries: " << eventEntries << std::endl;

  for (Long64_t i = 0; i < eventEntries; ++i) {
    eventTree->GetEntry(i);

    totalElectronProducedCount += electronProducedCount;
    totalElectronChannelCount += electronChannelCount;
    totalElectronChannelPlusCount += electronChannelPlusCount;
    totalElectronChannelMinusCount += electronChannelMinusCount;
    totalElectronChannelPlusMcpCount[0] +=
      electronChannelPlusMcp0Count;
    totalElectronChannelPlusMcpCount[1] +=
      electronChannelPlusMcp1Count;
    totalElectronChannelPlusMcpCount[2] +=
      electronChannelPlusMcp2Count;
    totalElectronChannelMinusMcpCount[0] +=
      electronChannelMinusMcp0Count;
    totalElectronChannelMinusMcpCount[1] +=
      electronChannelMinusMcp1Count;
    totalElectronChannelMinusMcpCount[2] +=
      electronChannelMinusMcp2Count;
    totalPhotonExitCount += photonExitCount;
    totalPhotonExitPlusCount += photonExitPlusCount;
    totalPhotonExitMinusCount += photonExitMinusCount;
    if (electronProducedCount > 0) {
      ++eventsWithProducedElectron;
    }
    if (electronChannelCount > 0) {
      ++eventsWithElectron;
    }
    if (hasElectronChannelPlus) {
      ++eventsDetectedPlus;
    }
    if (hasElectronChannelMinus) {
      ++eventsDetectedMinus;
    }
    if (isCoincidence) {
      ++coincidenceEvents;
    }
    if (photonExitCount > 0) {
      ++eventsWithPhotonExit;
    }
    producedMultiplicityCounts[electronProducedCount]++;
    multiplicityCounts[electronChannelCount]++;

    if (i < eventPreview) {
      std::cout << "  [" << i << "]"
                << " event=" << summaryEventID
                << " electronProducedCount=" << electronProducedCount
                << " electronChannelCount=" << electronChannelCount
                << " (+z=" << electronChannelPlusCount
                << ", -z=" << electronChannelMinusCount << ")"
                << " MCPs +z=("
                << electronChannelPlusMcp0Count << ","
                << electronChannelPlusMcp1Count << ","
                << electronChannelPlusMcp2Count << ")"
                << " -z=("
                << electronChannelMinusMcp0Count << ","
                << electronChannelMinusMcp1Count << ","
                << electronChannelMinusMcp2Count << ")"
                << " coincidence="
                << (isCoincidence ? "yes" : "no")
                << " photonExitCount=" << photonExitCount
                << " (+z=" << photonExitPlusCount
                << ", -z=" << photonExitMinusCount << ")"
                << std::endl;
    }
  }

  if (eventEntries > 0) {
    std::cout << "\nEvent-level results:" << std::endl;
    std::cout << "  events with produced electrons: "
              << eventsWithProducedElectron << " / " << eventEntries
              << " ("
              << static_cast<Double_t>(eventsWithProducedElectron)/
                   eventEntries
              << ")" << std::endl;
    std::cout << "  events with channel electron: "
              << eventsWithElectron << " / " << eventEntries
              << " ("
              << static_cast<Double_t>(eventsWithElectron)/eventEntries
              << ")" << std::endl;
    std::cout << "  events detected on +z: "
              << eventsDetectedPlus << " / " << eventEntries
              << " ("
              << static_cast<Double_t>(eventsDetectedPlus)/eventEntries
              << ")" << std::endl;
    std::cout << "  events detected on -z: "
              << eventsDetectedMinus << " / " << eventEntries
              << " ("
              << static_cast<Double_t>(eventsDetectedMinus)/eventEntries
              << ")" << std::endl;
    std::cout << "  PET coincidence efficiency: "
              << coincidenceEvents << " / " << eventEntries
              << " ("
              << static_cast<Double_t>(coincidenceEvents)/eventEntries
              << ", "
              << 100.0*static_cast<Double_t>(coincidenceEvents)/
                   eventEntries
              << " %)" << std::endl;
    std::cout << "  events with exiting photon: "
              << eventsWithPhotonExit << " / " << eventEntries
              << " ("
              << static_cast<Double_t>(eventsWithPhotonExit)/eventEntries
              << ")" << std::endl;
    std::cout << "  mean produced electrons/event: "
              << static_cast<Double_t>(totalElectronProducedCount)/
                   eventEntries
              << std::endl;
    std::cout << "  mean channel-electron multiplicity/event: "
              << static_cast<Double_t>(totalElectronChannelCount)/
                   eventEntries
              << std::endl;
    std::cout << "  mean photon exits/event: "
              << static_cast<Double_t>(totalPhotonExitCount)/eventEntries
              << std::endl;
    std::cout << "  channel electrons by side: +z="
              << totalElectronChannelPlusCount
              << ", -z=" << totalElectronChannelMinusCount << std::endl;
    std::cout << "  channel electrons by MCP:" << std::endl;
    for (Int_t mcpIndex = 0; mcpIndex < 3; ++mcpIndex) {
      std::cout << "    MCP " << mcpIndex
                << ": +z=" << totalElectronChannelPlusMcpCount[mcpIndex]
                << ", -z="
                << totalElectronChannelMinusMcpCount[mcpIndex]
                << std::endl;
    }
    std::cout << "  photon exits by side: +z="
              << totalPhotonExitPlusCount
              << ", -z=" << totalPhotonExitMinusCount << std::endl;
    std::cout << "  fraction of produced electrons reaching a channel: ";
    if (totalElectronProducedCount > 0) {
      std::cout
        << static_cast<Double_t>(totalElectronChannelCount)/
             totalElectronProducedCount
        << std::endl;
    } else {
      std::cout << "undefined (no produced electrons)" << std::endl;
    }
  }

  std::cout << "  produced-electron multiplicity distribution:"
            << std::endl;
  for (std::map<Int_t, Long64_t>::const_iterator it =
         producedMultiplicityCounts.begin();
       it != producedMultiplicityCounts.end(); ++it) {
    std::cout << "    multiplicity " << it->first
              << ": " << it->second << " events" << std::endl;
  }

  std::cout << "  channel-electron multiplicity distribution:" << std::endl;
  for (std::map<Int_t, Long64_t>::const_iterator it =
         multiplicityCounts.begin();
       it != multiplicityCounts.end(); ++it) {
    std::cout << "    multiplicity " << it->first
              << ": " << it->second << " events" << std::endl;
  }

  std::cout << "\nConsistency checks:" << std::endl;
  std::cout << "  electron rows=" << electronEntries
            << ", sum(event electronChannelCount)="
            << totalElectronChannelCount
            << (electronEntries == totalElectronChannelCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "  photon rows=" << photonEntries
            << ", sum(event photonExitCount)="
            << totalPhotonExitCount
            << (photonEntries == totalPhotonExitCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "  electron side totals="
            << totalElectronChannelPlusCount +
                 totalElectronChannelMinusCount
            << (totalElectronChannelCount ==
                totalElectronChannelPlusCount +
                totalElectronChannelMinusCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  const Long64_t totalElectronPlusByMcp =
    totalElectronChannelPlusMcpCount[0] +
    totalElectronChannelPlusMcpCount[1] +
    totalElectronChannelPlusMcpCount[2];
  const Long64_t totalElectronMinusByMcp =
    totalElectronChannelMinusMcpCount[0] +
    totalElectronChannelMinusMcpCount[1] +
    totalElectronChannelMinusMcpCount[2];
  std::cout << "  +z MCP totals=" << totalElectronPlusByMcp
            << (totalElectronPlusByMcp ==
                totalElectronChannelPlusCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "  -z MCP totals=" << totalElectronMinusByMcp
            << (totalElectronMinusByMcp ==
                totalElectronChannelMinusCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "  photon side totals="
            << totalPhotonExitPlusCount + totalPhotonExitMinusCount
            << (totalPhotonExitCount ==
                totalPhotonExitPlusCount + totalPhotonExitMinusCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;

  file->Close();
  delete file;
}
