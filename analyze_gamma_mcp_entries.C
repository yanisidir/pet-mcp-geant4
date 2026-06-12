#include "TFile.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

namespace
{
  typedef std::pair<Int_t, Int_t> TrackKey;
  typedef std::pair<Int_t, Int_t> PlateKey;
  typedef std::tuple<Int_t, Int_t, Int_t, Int_t> TrackPlateKey;

  struct FirstInteraction
  {
    FirstInteraction()
      : entryNumber(-1),
        globalTime_ns(std::numeric_limits<Double_t>::max()),
        side(0),
        mcpIndex(-1)
    {
    }

    Long64_t entryNumber;
    Double_t globalTime_ns;
    Int_t side;
    Int_t mcpIndex;
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

void analyze_gamma_mcp_entries(
  const char* fileName = "build/mcp_output_t0.root")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* entryTree = nullptr;
  TTree* interactionTree = nullptr;
  file->GetObject("GammaMcpEntryTree", entryTree);
  file->GetObject("GammaInteractionTree", interactionTree);

  if (!entryTree || !interactionTree) {
    std::cerr << "Missing GammaMcpEntryTree or GammaInteractionTree in "
              << fileName << std::endl;
    std::cerr << "Generate a new ROOT file with the updated Geant4 "
              << "executable." << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t entryEventID = -1;
  Int_t entryTrackID = -1;
  Int_t entrySide = 0;
  Int_t entryMcpIndex = -1;
  entryTree->SetBranchAddress("eventID", &entryEventID);
  entryTree->SetBranchAddress("trackID", &entryTrackID);
  entryTree->SetBranchAddress("side", &entrySide);
  entryTree->SetBranchAddress("mcpIndex", &entryMcpIndex);

  std::map<PlateKey, Long64_t> gammaEntriesByPlate;
  std::set<TrackPlateKey> enteredTrackPlates;
  std::set<PlateKey> allPlates;

  const Long64_t entryCount = entryTree->GetEntries();
  for (Long64_t entry = 0; entry < entryCount; ++entry) {
    entryTree->GetEntry(entry);

    const PlateKey plate(entrySide, entryMcpIndex);
    const TrackPlateKey trackPlate(entryEventID,
                                   entryTrackID,
                                   entrySide,
                                   entryMcpIndex);
    ++gammaEntriesByPlate[plate];
    enteredTrackPlates.insert(trackPlate);
    allPlates.insert(plate);
  }
  const Long64_t uniqueEntryKeyCount =
    static_cast<Long64_t>(enteredTrackPlates.size());
  const Long64_t duplicateEntryCount =
    entryCount - uniqueEntryKeyCount;

  Int_t interactionEventID = -1;
  Int_t interactionTrackID = -1;
  Int_t interactionSide = 0;
  Int_t interactionMcpIndex = -1;
  Double_t interactionTime_ns = 0.0;
  interactionTree->SetBranchAddress("eventID", &interactionEventID);
  interactionTree->SetBranchAddress("trackID", &interactionTrackID);
  interactionTree->SetBranchAddress("side", &interactionSide);
  interactionTree->SetBranchAddress("mcpIndex", &interactionMcpIndex);
  interactionTree->SetBranchAddress("globalTime_ns",
                                    &interactionTime_ns);

  std::map<PlateKey, Long64_t> interactionsByPlate;
  std::map<PlateKey, Long64_t> matchedInteractionsByPlate;
  std::set<TrackPlateKey> enteredTracksInteractingByPlate;
  std::map<TrackKey, FirstInteraction> firstInteractions;

  const Long64_t interactionCount = interactionTree->GetEntries();
  for (Long64_t entry = 0; entry < interactionCount; ++entry) {
    interactionTree->GetEntry(entry);

    const PlateKey plate(interactionSide, interactionMcpIndex);
    const TrackKey track(interactionEventID, interactionTrackID);
    const TrackPlateKey trackPlate(interactionEventID,
                                   interactionTrackID,
                                   interactionSide,
                                   interactionMcpIndex);

    ++interactionsByPlate[plate];
    allPlates.insert(plate);

    if (enteredTrackPlates.count(trackPlate) > 0) {
      ++matchedInteractionsByPlate[plate];
      enteredTracksInteractingByPlate.insert(trackPlate);
    }

    FirstInteraction& first = firstInteractions[track];
    const bool isEarlier =
      interactionTime_ns < first.globalTime_ns ||
      (interactionTime_ns == first.globalTime_ns &&
       (first.entryNumber < 0 || entry < first.entryNumber));
    if (isEarlier) {
      first.entryNumber = entry;
      first.globalTime_ns = interactionTime_ns;
      first.side = interactionSide;
      first.mcpIndex = interactionMcpIndex;
    }
  }

  std::map<PlateKey, Long64_t> firstInteractionsByPlate;
  std::map<PlateKey, Long64_t> matchedFirstInteractionsByPlate;
  for (std::map<TrackKey, FirstInteraction>::const_iterator it =
         firstInteractions.begin();
       it != firstInteractions.end(); ++it) {
    const FirstInteraction& first = it->second;
    const PlateKey plate(first.side, first.mcpIndex);
    const TrackPlateKey trackPlate(it->first.first,
                                   it->first.second,
                                   first.side,
                                   first.mcpIndex);

    ++firstInteractionsByPlate[plate];
    if (enteredTrackPlates.count(trackPlate) > 0) {
      ++matchedFirstInteractionsByPlate[plate];
    }
    allPlates.insert(plate);
  }

  std::map<PlateKey, Long64_t> interactingEnteredTracksByPlate;
  for (std::set<TrackPlateKey>::const_iterator it =
         enteredTracksInteractingByPlate.begin();
       it != enteredTracksInteractingByPlate.end(); ++it) {
    const PlateKey plate(std::get<2>(*it), std::get<3>(*it));
    ++interactingEnteredTracksByPlate[plate];
  }

  std::cout << "File: " << fileName << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\n=== Gamma entry and interaction analysis ==="
            << std::endl;
  std::cout << "GammaMcpEntryTree rows: " << entryCount << std::endl;
  std::cout << "Unique (eventID, trackID, side, mcpIndex) keys: "
            << uniqueEntryKeyCount << std::endl;
  std::cout << "Duplicate entry rows: "
            << duplicateEntryCount
            << (duplicateEntryCount == 0 ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "GammaInteractionTree rows: "
            << interactionCount << std::endl;

  std::cout << "\n"
            << std::left
            << std::setw(12) << "Plate"
            << std::right
            << std::setw(12) << "Entries"
            << std::setw(16) << "Interactions"
            << std::setw(14) << "First"
            << std::setw(18) << "Int/entry"
            << std::setw(18) << "First/entry"
            << std::setw(20) << "P(any|entered)"
            << std::endl;

  for (std::set<PlateKey>::const_iterator plateIt =
         allPlates.begin();
       plateIt != allPlates.end(); ++plateIt) {
    const Long64_t entries = GetCount(gammaEntriesByPlate, *plateIt);
    const Long64_t interactions =
      GetCount(interactionsByPlate, *plateIt);
    const Long64_t firstInteractions =
      GetCount(firstInteractionsByPlate, *plateIt);
    const Long64_t interactingEnteredTracks =
      GetCount(interactingEnteredTracksByPlate, *plateIt);

    const Double_t interactionsPerEntry =
      entries > 0
        ? static_cast<Double_t>(interactions)/entries
        : 0.0;
    const Double_t firstPerEntry =
      entries > 0
        ? static_cast<Double_t>(firstInteractions)/entries
        : 0.0;
    const Double_t conditionalProbability =
      entries > 0
        ? static_cast<Double_t>(interactingEnteredTracks)/entries
        : 0.0;

    std::cout << std::left
              << std::setw(12) << PlateLabel(*plateIt)
              << std::right
              << std::setw(12) << entries
              << std::setw(16) << interactions
              << std::setw(14) << firstInteractions
              << std::setw(18) << interactionsPerEntry
              << std::setw(18) << firstPerEntry
              << std::setw(20) << conditionalProbability
              << std::endl;
  }

  Long64_t unmatchedInteractionCount = 0;
  Long64_t unmatchedFirstInteractionCount = 0;
  for (std::set<PlateKey>::const_iterator plateIt =
         allPlates.begin();
       plateIt != allPlates.end(); ++plateIt) {
    unmatchedInteractionCount +=
      GetCount(interactionsByPlate, *plateIt) -
      GetCount(matchedInteractionsByPlate, *plateIt);
    unmatchedFirstInteractionCount +=
      GetCount(firstInteractionsByPlate, *plateIt) -
      GetCount(matchedFirstInteractionsByPlate, *plateIt);
  }

  std::cout << "\nInterpretation:" << std::endl;
  std::cout << "  Int/entry is the requested interactions/entries ratio; "
            << "it is a mean multiplicity and may exceed 1." << std::endl;
  std::cout << "  First/entry uses the first interaction of each gamma "
            << "track, assigned to its first MCP plate." << std::endl;
  std::cout << "  P(any|entered) counts entered gamma tracks with at least "
            << "one interaction in the same plate; it is bounded by 1."
            << std::endl;
  std::cout << "  Interaction rows without a matching entry in the same "
            << "track and plate: " << unmatchedInteractionCount
            << std::endl;
  std::cout << "  First interactions without a matching entry in the same "
            << "track and plate: " << unmatchedFirstInteractionCount
            << std::endl;

  file->Close();
  delete file;
}
