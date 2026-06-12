#include "TFile.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>

namespace
{
  typedef std::pair<Int_t, Int_t> TrackKey;
  typedef std::pair<Int_t, Int_t> PlateKey;

  struct GammaTrackSummary
  {
    GammaTrackSummary()
      : interactionCount(0),
        firstEntry(-1),
        firstTime_ns(std::numeric_limits<Double_t>::max()),
        firstSide(0),
        firstMcpIndex(-1)
    {
    }

    Long64_t interactionCount;
    Long64_t firstEntry;
    Double_t firstTime_ns;
    Int_t firstSide;
    Int_t firstMcpIndex;
    std::string firstProcessName;
    std::set<PlateKey> visitedPlates;
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

void validate_gamma_interactions(
  const char* fileName = "build/mcp_output_t0.root")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* gammaTree = nullptr;
  file->GetObject("GammaInteractionTree", gammaTree);
  if (!gammaTree) {
    std::cerr << "Missing GammaInteractionTree in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t eventID = -1;
  Int_t trackID = -1;
  Int_t side = 0;
  Int_t mcpIndex = -1;
  Double_t globalTime_ns = 0.0;
  char processName[64] = "";
  char volumeName[64] = "";

  gammaTree->SetBranchAddress("eventID", &eventID);
  gammaTree->SetBranchAddress("trackID", &trackID);
  gammaTree->SetBranchAddress("side", &side);
  gammaTree->SetBranchAddress("mcpIndex", &mcpIndex);
  gammaTree->SetBranchAddress("globalTime_ns", &globalTime_ns);
  gammaTree->SetBranchAddress("processName", processName);
  gammaTree->SetBranchAddress("volumeName", volumeName);

  std::map<std::string, Long64_t> processCounts;
  std::map<std::string, Long64_t> volumeCounts;
  std::map<PlateKey, Long64_t> totalByPlate;
  std::map<PlateKey, Long64_t> excludingRaylByPlate;
  std::map<PlateKey, Long64_t> photComptByPlate;
  std::map<TrackKey, GammaTrackSummary> tracks;
  std::set<PlateKey> allPlates;

  const Long64_t totalEntries = gammaTree->GetEntries();
  for (Long64_t entry = 0; entry < totalEntries; ++entry) {
    gammaTree->GetEntry(entry);

    const std::string process(processName);
    const std::string volume(volumeName);
    const PlateKey plate(side, mcpIndex);
    const TrackKey track(eventID, trackID);

    ++processCounts[process];
    ++volumeCounts[volume];
    ++totalByPlate[plate];
    allPlates.insert(plate);

    if (process != "Rayl") {
      ++excludingRaylByPlate[plate];
    }
    if (process == "phot" || process == "compt") {
      ++photComptByPlate[plate];
    }

    GammaTrackSummary& summary = tracks[track];
    ++summary.interactionCount;
    summary.visitedPlates.insert(plate);

    const bool isEarlier =
      globalTime_ns < summary.firstTime_ns ||
      (globalTime_ns == summary.firstTime_ns &&
       (summary.firstEntry < 0 || entry < summary.firstEntry));
    if (isEarlier) {
      summary.firstEntry = entry;
      summary.firstTime_ns = globalTime_ns;
      summary.firstSide = side;
      summary.firstMcpIndex = mcpIndex;
      summary.firstProcessName = process;
    }
  }

  std::map<PlateKey, Long64_t> firstInteractionByPlate;
  std::map<PlateKey, Long64_t> uniqueTracksByPlate;
  std::map<std::string, Long64_t> firstProcessCounts;
  std::map<PlateKey, std::map<std::string, Long64_t> >
    firstProcessCountsByPlate;
  Long64_t tracksWithOneInteraction = 0;
  Long64_t tracksWithTwoInteractions = 0;
  Long64_t tracksWithThreeOrMoreInteractions = 0;

  for (std::map<TrackKey, GammaTrackSummary>::const_iterator it =
         tracks.begin();
       it != tracks.end(); ++it) {
    const GammaTrackSummary& summary = it->second;

    const PlateKey firstPlate(summary.firstSide,
                              summary.firstMcpIndex);
    ++firstInteractionByPlate[firstPlate];
    ++firstProcessCounts[summary.firstProcessName];
    ++firstProcessCountsByPlate[firstPlate][
      summary.firstProcessName];
    allPlates.insert(firstPlate);

    for (std::set<PlateKey>::const_iterator plateIt =
           summary.visitedPlates.begin();
         plateIt != summary.visitedPlates.end(); ++plateIt) {
      ++uniqueTracksByPlate[*plateIt];
      allPlates.insert(*plateIt);
    }

    if (summary.interactionCount == 1) {
      ++tracksWithOneInteraction;
    } else if (summary.interactionCount == 2) {
      ++tracksWithTwoInteractions;
    } else {
      ++tracksWithThreeOrMoreInteractions;
    }
  }

  const Long64_t uniqueTrackCount =
    static_cast<Long64_t>(tracks.size());
  const Long64_t repeatedInteractionCount =
    totalEntries - uniqueTrackCount;

  std::cout << "File: " << fileName << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\n=== Gamma interaction count validation ==="
            << std::endl;
  std::cout << "Total GammaInteractionTree rows: "
            << totalEntries << std::endl;
  std::cout << "Unique interacting gamma tracks (eventID, trackID): "
            << uniqueTrackCount << std::endl;
  std::cout << "Additional rows caused by multiple interactions: "
            << repeatedInteractionCount << std::endl;
  if (uniqueTrackCount > 0) {
    std::cout << "Mean interactions per interacting gamma track: "
              << static_cast<Double_t>(totalEntries)/uniqueTrackCount
              << std::endl;
  }

  std::cout << "\nInteractions by process:" << std::endl;
  for (std::map<std::string, Long64_t>::const_iterator it =
         processCounts.begin();
       it != processCounts.end(); ++it) {
    std::cout << "  " << it->first << ": "
              << it->second << std::endl;
  }

  std::cout << "\nFirst interaction process by unique gamma track:"
            << std::endl;
  for (std::map<std::string, Long64_t>::const_iterator it =
         firstProcessCounts.begin();
       it != firstProcessCounts.end(); ++it) {
    std::cout << "  " << it->first << ": "
              << it->second << std::endl;
  }

  std::cout << "\nFirst interaction process by MCP plate:"
            << std::endl;
  for (std::set<PlateKey>::const_iterator plateIt =
         allPlates.begin();
       plateIt != allPlates.end(); ++plateIt) {
    const std::map<std::string, Long64_t>& processCounts =
      firstProcessCountsByPlate[*plateIt];
    const char* processNames[4] = {"phot", "compt", "Rayl", "conv"};

    std::cout << "  " << PlateLabel(*plateIt) << ":";
    for (Int_t processIndex = 0; processIndex < 4; ++processIndex) {
      const std::string processName = processNames[processIndex];
      const std::map<std::string, Long64_t>::const_iterator found =
        processCounts.find(processName);
      const Long64_t count =
        found != processCounts.end() ? found->second : 0;
      std::cout << " " << processName << "=" << count;
    }
    std::cout << std::endl;
  }

  std::cout << "\nInteracting gamma tracks by multiplicity:"
            << std::endl;
  std::cout << "  1 interaction: "
            << tracksWithOneInteraction << std::endl;
  std::cout << "  2 interactions: "
            << tracksWithTwoInteractions << std::endl;
  std::cout << "  3 interactions or more: "
            << tracksWithThreeOrMoreInteractions << std::endl;

  std::cout << "\nVolume names present in GammaInteractionTree:"
            << std::endl;
  bool onlyMcpBodyVolumes = true;
  for (std::map<std::string, Long64_t>::const_iterator it =
         volumeCounts.begin();
       it != volumeCounts.end(); ++it) {
    const bool isMcpBody =
      it->first.find("MCP_") == 0 &&
      it->first.find("_body_") != std::string::npos;
    if (!isMcpBody) {
      onlyMcpBodyVolumes = false;
    }

    std::cout << "  " << it->first << ": " << it->second
              << (isMcpBody ? " [MCP_body]" : " [UNEXPECTED]")
              << std::endl;
  }
  std::cout << "MCP_body-only volume check: "
            << (onlyMcpBodyVolumes ? "PASS" : "FAIL")
            << std::endl;

  std::cout << "\n=== Plate summary ===" << std::endl;
  std::cout << std::left
            << std::setw(12) << "Plate"
            << std::right
            << std::setw(12) << "Total"
            << std::setw(14) << "No Rayl"
            << std::setw(14) << "phot+compt"
            << std::setw(16) << "First only"
            << std::setw(18) << "Unique tracks"
            << std::endl;

  for (std::set<PlateKey>::const_iterator plateIt =
         allPlates.begin();
       plateIt != allPlates.end(); ++plateIt) {
    std::cout << std::left
              << std::setw(12) << PlateLabel(*plateIt)
              << std::right
              << std::setw(12) << GetCount(totalByPlate, *plateIt)
              << std::setw(14) << GetCount(excludingRaylByPlate, *plateIt)
              << std::setw(14) << GetCount(photComptByPlate, *plateIt)
              << std::setw(16)
              << GetCount(firstInteractionByPlate, *plateIt)
              << std::setw(18)
              << GetCount(uniqueTracksByPlate, *plateIt)
              << std::endl;
  }

  std::cout << "\nDefinitions:" << std::endl;
  std::cout << "  Total: every physical interaction row."
            << std::endl;
  std::cout << "  First only: earliest globalTime_ns for each "
            << "(eventID, trackID)." << std::endl;
  std::cout << "  Unique tracks: each gamma track counted once per "
            << "plate where it interacted." << std::endl;
  std::cout << "  phot+compt: only photoelectric and Compton rows."
            << std::endl;

  file->Close();
  delete file;
}
