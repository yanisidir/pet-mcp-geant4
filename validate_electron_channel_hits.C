#include "TFile.h"
#include "TTree.h"

#include <algorithm>
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
  typedef std::pair<Int_t, Int_t> TrackKey;
  typedef std::pair<Int_t, Int_t> PlateKey;

  struct ElectronHitRecord
  {
    Int_t eventID;
    Int_t trackID;
    Int_t parentID;
    Int_t side;
    Int_t mcpIndex;
    Double_t kineticEnergy_keV;
    std::string creatorProcessName;
    std::string volumeName;
  };

  std::string PlateLabel(const PlateKey& plate)
  {
    std::ostringstream label;
    label << (plate.first > 0 ? "+z" : "-z")
          << " MCP " << plate.second;
    return label.str();
  }

  Long64_t GetCount(const std::map<PlateKey, Long64_t>& counts,
                    const PlateKey& plate)
  {
    const std::map<PlateKey, Long64_t>::const_iterator found =
      counts.find(plate);
    return found != counts.end() ? found->second : 0;
  }
}

void validate_electron_channel_hits(
  const char* fileName = "build/mcp_output_t0.root")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* electronTree = nullptr;
  TTree* plateTree = nullptr;
  TTree* eventTree = nullptr;
  file->GetObject("ElectronChannelHitTree", electronTree);
  file->GetObject("McpPlateStatsTree", plateTree);
  file->GetObject("EventSummaryTree", eventTree);

  if (!electronTree || !plateTree || !eventTree) {
    std::cerr << "Missing ElectronChannelHitTree, McpPlateStatsTree "
              << "or EventSummaryTree in " << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t eventID = -1;
  Int_t trackID = -1;
  Int_t parentID = -1;
  Int_t side = 0;
  Int_t mcpIndex = -1;
  Double_t kineticEnergy_keV = 0.0;
  char creatorProcessName[64] = "";
  char volumeName[64] = "";

  electronTree->SetBranchAddress("eventID", &eventID);
  electronTree->SetBranchAddress("trackID", &trackID);
  electronTree->SetBranchAddress("parentID", &parentID);
  electronTree->SetBranchAddress("side", &side);
  electronTree->SetBranchAddress("mcpIndex", &mcpIndex);
  electronTree->SetBranchAddress("kineticEnergy_keV",
                                 &kineticEnergy_keV);
  electronTree->SetBranchAddress("creatorProcessName",
                                 creatorProcessName);
  electronTree->SetBranchAddress("volumeName", volumeName);

  std::set<TrackKey> uniqueTracks;
  std::map<TrackKey, Long64_t> rowsByTrack;
  std::map<Int_t, Long64_t> electronsBySide;
  std::map<PlateKey, Long64_t> electronsByPlate;
  std::map<PlateKey, std::set<TrackKey> > uniqueTracksByPlate;
  std::map<std::string, Long64_t> creatorProcessCounts;
  std::map<PlateKey, std::map<std::string, Long64_t> >
    creatorProcessCountsByPlate;
  std::map<std::string, Long64_t> volumeCounts;
  std::map<Int_t, std::vector<ElectronHitRecord> > hitsByEvent;
  std::set<PlateKey> allPlates;

  const Long64_t electronEntries = electronTree->GetEntries();
  for (Long64_t entry = 0; entry < electronEntries; ++entry) {
    electronTree->GetEntry(entry);

    const TrackKey track(eventID, trackID);
    const PlateKey plate(side, mcpIndex);
    const std::string creator(creatorProcessName);
    const std::string volume(volumeName);

    uniqueTracks.insert(track);
    ++rowsByTrack[track];
    ++electronsBySide[side];
    ++electronsByPlate[plate];
    uniqueTracksByPlate[plate].insert(track);
    ++creatorProcessCounts[creator];
    ++creatorProcessCountsByPlate[plate][creator];
    ++volumeCounts[volume];
    allPlates.insert(plate);

    ElectronHitRecord hit;
    hit.eventID = eventID;
    hit.trackID = trackID;
    hit.parentID = parentID;
    hit.side = side;
    hit.mcpIndex = mcpIndex;
    hit.kineticEnergy_keV = kineticEnergy_keV;
    hit.creatorProcessName = creator;
    hit.volumeName = volume;
    hitsByEvent[eventID].push_back(hit);
  }

  Long64_t duplicatedTrackCount = 0;
  Long64_t duplicateRowCount = 0;
  for (std::map<TrackKey, Long64_t>::const_iterator it =
         rowsByTrack.begin();
       it != rowsByTrack.end(); ++it) {
    if (it->second > 1) {
      ++duplicatedTrackCount;
      duplicateRowCount += it->second - 1;
    }
  }

  Int_t plateSide = 0;
  Int_t plateMcpIndex = -1;
  Int_t plateElectronChannelCount = 0;
  plateTree->SetBranchAddress("side", &plateSide);
  plateTree->SetBranchAddress("mcpIndex", &plateMcpIndex);
  plateTree->SetBranchAddress("electronChannelCount",
                              &plateElectronChannelCount);

  std::map<PlateKey, Long64_t> plateStatsTotals;
  const Long64_t plateEntries = plateTree->GetEntries();
  for (Long64_t entry = 0; entry < plateEntries; ++entry) {
    plateTree->GetEntry(entry);
    const PlateKey plate(plateSide, plateMcpIndex);
    plateStatsTotals[plate] += plateElectronChannelCount;
    allPlates.insert(plate);
  }

  Int_t summaryElectronChannelCount = 0;
  Int_t summaryElectronChannelPlusCount = 0;
  Int_t summaryElectronChannelMinusCount = 0;
  eventTree->SetBranchAddress("electronChannelCount",
                              &summaryElectronChannelCount);
  eventTree->SetBranchAddress("electronChannelPlusCount",
                              &summaryElectronChannelPlusCount);
  eventTree->SetBranchAddress("electronChannelMinusCount",
                              &summaryElectronChannelMinusCount);

  Long64_t summaryTotal = 0;
  Long64_t summaryPlusTotal = 0;
  Long64_t summaryMinusTotal = 0;
  const Long64_t eventEntries = eventTree->GetEntries();
  for (Long64_t entry = 0; entry < eventEntries; ++entry) {
    eventTree->GetEntry(entry);
    summaryTotal += summaryElectronChannelCount;
    summaryPlusTotal += summaryElectronChannelPlusCount;
    summaryMinusTotal += summaryElectronChannelMinusCount;
  }

  std::cout << "File: " << fileName << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\n=== Electron channel-hit validation ==="
            << std::endl;
  std::cout << "Total ElectronChannelHitTree rows: "
            << electronEntries << std::endl;
  std::cout << "Unique electron tracks (eventID, trackID): "
            << uniqueTracks.size() << std::endl;
  std::cout << "Tracks appearing more than once: "
            << duplicatedTrackCount << std::endl;
  std::cout << "Additional duplicate rows: "
            << duplicateRowCount << std::endl;

  std::cout << "\nElectrons by side:" << std::endl;
  std::cout << "  +z: " << electronsBySide[1] << std::endl;
  std::cout << "  -z: " << electronsBySide[-1] << std::endl;

  std::cout << "\nElectrons by creator process:" << std::endl;
  for (std::map<std::string, Long64_t>::const_iterator it =
         creatorProcessCounts.begin();
       it != creatorProcessCounts.end(); ++it) {
    std::cout << "  " << it->first << ": "
              << it->second << std::endl;
  }

  std::cout << "\nVolume names present in ElectronChannelHitTree:"
            << std::endl;
  bool onlyMcpChannelVolumes = true;
  for (std::map<std::string, Long64_t>::const_iterator it =
         volumeCounts.begin();
       it != volumeCounts.end(); ++it) {
    const bool isMcpChannel =
      it->first.find("MCP_") == 0 &&
      it->first.find("_channel_") != std::string::npos;
    if (!isMcpChannel) {
      onlyMcpChannelVolumes = false;
    }

    std::cout << "  " << it->first << ": " << it->second
              << (isMcpChannel ? " [MCP_channel]" : " [UNEXPECTED]")
              << std::endl;
  }
  std::cout << "MCP_channel-only volume check: "
            << (onlyMcpChannelVolumes ? "PASS" : "FAIL")
            << std::endl;

  std::cout << "\n=== Plate summary ===" << std::endl;
  std::cout << std::left
            << std::setw(12) << "Plate"
            << std::right
            << std::setw(14) << "Tree rows"
            << std::setw(16) << "Unique tracks"
            << std::setw(18) << "PlateStats total"
            << std::endl;

  for (std::set<PlateKey>::const_iterator plateIt =
         allPlates.begin();
       plateIt != allPlates.end(); ++plateIt) {
    const std::map<PlateKey, std::set<TrackKey> >::const_iterator found =
      uniqueTracksByPlate.find(*plateIt);
    const Long64_t uniquePlateTrackCount =
      found != uniqueTracksByPlate.end()
        ? static_cast<Long64_t>(found->second.size())
        : 0;

    std::cout << std::left
              << std::setw(12) << PlateLabel(*plateIt)
              << std::right
              << std::setw(14) << GetCount(electronsByPlate, *plateIt)
              << std::setw(16) << uniquePlateTrackCount
              << std::setw(18) << GetCount(plateStatsTotals, *plateIt)
              << std::endl;
  }

  std::cout << "\nCreator process by MCP plate:" << std::endl;
  for (std::set<PlateKey>::const_iterator plateIt =
         allPlates.begin();
       plateIt != allPlates.end(); ++plateIt) {
    const std::map<std::string, Long64_t>& processCounts =
      creatorProcessCountsByPlate[*plateIt];

    std::cout << "  " << PlateLabel(*plateIt) << ":";
    if (processCounts.empty()) {
      std::cout << " no electron";
    } else {
      for (std::map<std::string, Long64_t>::const_iterator processIt =
             processCounts.begin();
           processIt != processCounts.end(); ++processIt) {
        std::cout << " " << processIt->first
                  << "=" << processIt->second;
      }
    }
    std::cout << std::endl;
  }

  const Long64_t electronPlusTotal = electronsBySide[1];
  const Long64_t electronMinusTotal = electronsBySide[-1];
  Long64_t plateMismatchCount = 0;
  for (std::set<PlateKey>::const_iterator plateIt =
         allPlates.begin();
       plateIt != allPlates.end(); ++plateIt) {
    if (GetCount(electronsByPlate, *plateIt) !=
        GetCount(plateStatsTotals, *plateIt)) {
      ++plateMismatchCount;
    }
  }

  std::cout << "\n=== Consistency checks ===" << std::endl;
  std::cout << "Rows vs unique (eventID, trackID): "
            << electronEntries << " vs " << uniqueTracks.size()
            << (electronEntries ==
                static_cast<Long64_t>(uniqueTracks.size())
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "Tree rows vs sum(EventSummary electronChannelCount): "
            << electronEntries << " vs " << summaryTotal
            << (electronEntries == summaryTotal
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "+z tree vs EventSummary: "
            << electronPlusTotal << " vs " << summaryPlusTotal
            << (electronPlusTotal == summaryPlusTotal
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "-z tree vs EventSummary: "
            << electronMinusTotal << " vs " << summaryMinusTotal
            << (electronMinusTotal == summaryMinusTotal
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "Per-plate comparison with McpPlateStatsTree: "
            << (plateMismatchCount == 0 ? "PASS" : "FAIL")
            << " (" << plateMismatchCount << " mismatched plates)"
            << std::endl;

  const std::size_t maximumExampleEvents = 5;
  std::size_t exampleEventCount = 0;
  std::cout << "\n=== Example events with multiple channel electrons ==="
            << std::endl;
  for (std::map<Int_t, std::vector<ElectronHitRecord> >::const_iterator it =
         hitsByEvent.begin();
       it != hitsByEvent.end() &&
       exampleEventCount < maximumExampleEvents;
       ++it) {
    if (it->second.size() <= 1) {
      continue;
    }

    ++exampleEventCount;
    std::cout << "Event " << it->first << " has "
              << it->second.size() << " electron channel hits:"
              << std::endl;
    for (std::vector<ElectronHitRecord>::const_iterator hitIt =
           it->second.begin();
         hitIt != it->second.end(); ++hitIt) {
      std::cout << "  trackID=" << hitIt->trackID
                << " parentID=" << hitIt->parentID
                << " side=" << hitIt->side
                << " mcpIndex=" << hitIt->mcpIndex
                << " creator=" << hitIt->creatorProcessName
                << " E=" << hitIt->kineticEnergy_keV << " keV"
                << " volume=" << hitIt->volumeName
                << std::endl;
    }
  }
  if (exampleEventCount == 0) {
    std::cout << "No event with multiple channel electrons."
              << std::endl;
  }

  file->Close();
  delete file;
}
