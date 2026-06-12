#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
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

  struct CoincidenceMultiplicity
  {
    Int_t plusCount;
    Int_t minusCount;
  };

  struct ActivePlates
  {
    std::set<Int_t> plusIndices;
    std::set<Int_t> minusIndices;
  };

  struct FirstGammaInteraction
  {
    FirstGammaInteraction()
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

  std::string PlateLabel(const char* side, Int_t mcpIndex)
  {
    std::ostringstream label;
    label << side << " MCP " << mcpIndex;
    return label.str();
  }

  Long64_t PlateCount(const std::map<PlateKey, Long64_t>& counts,
                      const PlateKey& plate)
  {
    const std::map<PlateKey, Long64_t>::const_iterator found =
      counts.find(plate);
    return found != counts.end() ? found->second : 0;
  }
}

void summarize_coincidence_full(
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
  TTree* gammaTree = nullptr;
  TTree* gammaEntryTree = nullptr;
  TTree* electronTree = nullptr;
  file->GetObject("EventSummaryTree", eventTree);
  file->GetObject("McpPlateStatsTree", plateTree);
  file->GetObject("GammaInteractionTree", gammaTree);
  file->GetObject("GammaMcpEntryTree", gammaEntryTree);
  file->GetObject("ElectronChannelHitTree", electronTree);

  if (!eventTree || !plateTree || !gammaTree ||
      !gammaEntryTree || !electronTree) {
    std::cerr << "Missing one or more required trees in "
              << fileName << std::endl;
    std::cerr << "Required: EventSummaryTree, McpPlateStatsTree, "
              << "GammaInteractionTree, GammaMcpEntryTree, "
              << "ElectronChannelHitTree"
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  // --------------------------------------------------------------------------
  // Select coincidence events and collect their electron multiplicities.
  // --------------------------------------------------------------------------

  Int_t eventID = -1;
  Int_t plusCount = 0;
  Int_t minusCount = 0;
  Bool_t isCoincidence = false;

  eventTree->SetBranchAddress("eventID", &eventID);
  eventTree->SetBranchAddress("electronChannelPlusCount", &plusCount);
  eventTree->SetBranchAddress("electronChannelMinusCount", &minusCount);
  eventTree->SetBranchAddress("isCoincidence", &isCoincidence);

  std::set<Int_t> coincidenceEventIDs;
  std::vector<CoincidenceMultiplicity> multiplicities;
  std::map<Int_t, Long64_t> plusMultiplicityDistribution;
  std::map<Int_t, Long64_t> minusMultiplicityDistribution;
  Long64_t plusMultiplicitySum = 0;
  Long64_t minusMultiplicitySum = 0;
  Int_t maximumPlusMultiplicity = 0;
  Int_t maximumMinusMultiplicity = 0;

  const Long64_t totalEventCount = eventTree->GetEntries();
  for (Long64_t entry = 0; entry < totalEventCount; ++entry) {
    eventTree->GetEntry(entry);
    if (!isCoincidence) {
      continue;
    }

    coincidenceEventIDs.insert(eventID);

    CoincidenceMultiplicity multiplicity;
    multiplicity.plusCount = plusCount;
    multiplicity.minusCount = minusCount;
    multiplicities.push_back(multiplicity);

    ++plusMultiplicityDistribution[plusCount];
    ++minusMultiplicityDistribution[minusCount];
    plusMultiplicitySum += plusCount;
    minusMultiplicitySum += minusCount;
    maximumPlusMultiplicity =
      std::max(maximumPlusMultiplicity, plusCount);
    maximumMinusMultiplicity =
      std::max(maximumMinusMultiplicity, minusCount);
  }

  const Long64_t coincidenceEventCount =
    static_cast<Long64_t>(coincidenceEventIDs.size());
  const Double_t coincidenceEfficiency =
    totalEventCount > 0
      ? static_cast<Double_t>(coincidenceEventCount)/totalEventCount
      : 0.0;

  if (coincidenceEventCount == 0) {
    std::cout << "File: " << fileName << std::endl;
    std::cout << "Total events: " << totalEventCount << std::endl;
    std::cout << "Coincidence events: 0" << std::endl;
    std::cout << "Coincidence efficiency: 0" << std::endl;
    file->Close();
    delete file;
    return;
  }

  // --------------------------------------------------------------------------
  // Find active plates in coincidence events and build the plate matrix.
  // --------------------------------------------------------------------------

  Int_t plateEventID = -1;
  Int_t plateSide = 0;
  Int_t plateMcpIndex = -1;
  Int_t plateElectronCount = 0;

  plateTree->SetBranchAddress("eventID", &plateEventID);
  plateTree->SetBranchAddress("side", &plateSide);
  plateTree->SetBranchAddress("mcpIndex", &plateMcpIndex);
  plateTree->SetBranchAddress("electronChannelCount",
                              &plateElectronCount);

  std::set<Int_t> plusMcpIndices;
  std::set<Int_t> minusMcpIndices;
  std::map<Int_t, ActivePlates> activePlatesByEvent;
  std::map<PlateKey, Long64_t> plateStatsElectronCounts;
  std::map<PlateKey, Long64_t> coincidenceEventsByPlate;

  const Long64_t plateEntries = plateTree->GetEntries();
  for (Long64_t entry = 0; entry < plateEntries; ++entry) {
    plateTree->GetEntry(entry);

    if (plateSide > 0) {
      plusMcpIndices.insert(plateMcpIndex);
    } else if (plateSide < 0) {
      minusMcpIndices.insert(plateMcpIndex);
    }

    if (coincidenceEventIDs.count(plateEventID) == 0) {
      continue;
    }

    const PlateKey plate(plateSide, plateMcpIndex);
    plateStatsElectronCounts[plate] += plateElectronCount;

    if (plateElectronCount <= 0) {
      continue;
    }

    ++coincidenceEventsByPlate[plate];
    if (plateSide > 0) {
      activePlatesByEvent[plateEventID].plusIndices.insert(
        plateMcpIndex);
    } else if (plateSide < 0) {
      activePlatesByEvent[plateEventID].minusIndices.insert(
        plateMcpIndex);
    }
  }

  // --------------------------------------------------------------------------
  // Count gamma entries in the same coincidence-event selection.
  // --------------------------------------------------------------------------

  Int_t gammaEntryEventID = -1;
  Int_t gammaEntrySide = 0;
  Int_t gammaEntryMcpIndex = -1;
  gammaEntryTree->SetBranchAddress("eventID", &gammaEntryEventID);
  gammaEntryTree->SetBranchAddress("side", &gammaEntrySide);
  gammaEntryTree->SetBranchAddress("mcpIndex", &gammaEntryMcpIndex);

  std::map<PlateKey, Long64_t> gammaEntriesCoincidence;
  const Long64_t gammaEntryCount = gammaEntryTree->GetEntries();
  for (Long64_t entry = 0; entry < gammaEntryCount; ++entry) {
    gammaEntryTree->GetEntry(entry);
    if (coincidenceEventIDs.count(gammaEntryEventID) == 0) {
      continue;
    }

    ++gammaEntriesCoincidence[
      PlateKey(gammaEntrySide, gammaEntryMcpIndex)];
  }

  std::map<std::pair<Int_t, Int_t>, Long64_t> platePairCounts;
  Long64_t validPlateCoincidenceEvents = 0;
  Long64_t totalPlatePairContributions = 0;

  for (std::set<Int_t>::const_iterator eventIt =
         coincidenceEventIDs.begin();
       eventIt != coincidenceEventIDs.end(); ++eventIt) {
    const std::map<Int_t, ActivePlates>::const_iterator found =
      activePlatesByEvent.find(*eventIt);
    if (found == activePlatesByEvent.end() ||
        found->second.plusIndices.empty() ||
        found->second.minusIndices.empty()) {
      continue;
    }

    ++validPlateCoincidenceEvents;
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

  // --------------------------------------------------------------------------
  // Count detailed electron rows by plate for coincidence events.
  // --------------------------------------------------------------------------

  Int_t electronEventID = -1;
  Int_t electronSide = 0;
  Int_t electronMcpIndex = -1;
  electronTree->SetBranchAddress("eventID", &electronEventID);
  electronTree->SetBranchAddress("side", &electronSide);
  electronTree->SetBranchAddress("mcpIndex", &electronMcpIndex);

  std::map<PlateKey, Long64_t> detectedElectronsByPlate;
  const Long64_t electronEntries = electronTree->GetEntries();
  for (Long64_t entry = 0; entry < electronEntries; ++entry) {
    electronTree->GetEntry(entry);
    if (coincidenceEventIDs.count(electronEventID) == 0) {
      continue;
    }

    ++detectedElectronsByPlate[
      PlateKey(electronSide, electronMcpIndex)];
  }

  // --------------------------------------------------------------------------
  // Count gamma processes and first interactions in coincidence events.
  // --------------------------------------------------------------------------

  Int_t gammaEventID = -1;
  Int_t gammaTrackID = -1;
  Int_t gammaSide = 0;
  Int_t gammaMcpIndex = -1;
  Double_t gammaGlobalTime_ns = 0.0;
  char gammaProcessName[64] = "";

  gammaTree->SetBranchAddress("eventID", &gammaEventID);
  gammaTree->SetBranchAddress("trackID", &gammaTrackID);
  gammaTree->SetBranchAddress("side", &gammaSide);
  gammaTree->SetBranchAddress("mcpIndex", &gammaMcpIndex);
  gammaTree->SetBranchAddress("globalTime_ns",
                              &gammaGlobalTime_ns);
  gammaTree->SetBranchAddress("processName", gammaProcessName);

  std::map<std::string, Long64_t> gammaProcessCounts;
  std::map<std::string, std::set<Int_t> > eventsByGammaProcess;
  std::map<TrackKey, FirstGammaInteraction> firstInteractionByTrack;

  const Long64_t gammaEntries = gammaTree->GetEntries();
  for (Long64_t entry = 0; entry < gammaEntries; ++entry) {
    gammaTree->GetEntry(entry);
    if (coincidenceEventIDs.count(gammaEventID) == 0) {
      continue;
    }

    const std::string process(gammaProcessName);
    ++gammaProcessCounts[process];
    eventsByGammaProcess[process].insert(gammaEventID);

    const TrackKey track(gammaEventID, gammaTrackID);
    FirstGammaInteraction& first = firstInteractionByTrack[track];
    const bool isEarlier =
      gammaGlobalTime_ns < first.globalTime_ns ||
      (gammaGlobalTime_ns == first.globalTime_ns &&
       (first.entryNumber < 0 || entry < first.entryNumber));
    if (isEarlier) {
      first.entryNumber = entry;
      first.globalTime_ns = gammaGlobalTime_ns;
      first.side = gammaSide;
      first.mcpIndex = gammaMcpIndex;
    }
  }

  std::map<PlateKey, Long64_t> firstGammaInteractionsByPlate;
  for (std::map<TrackKey, FirstGammaInteraction>::const_iterator it =
         firstInteractionByTrack.begin();
       it != firstInteractionByTrack.end(); ++it) {
    ++firstGammaInteractionsByPlate[
      PlateKey(it->second.side, it->second.mcpIndex)];
  }

  // --------------------------------------------------------------------------
  // Terminal summary.
  // --------------------------------------------------------------------------

  const Double_t meanPlus =
    static_cast<Double_t>(plusMultiplicitySum)/coincidenceEventCount;
  const Double_t meanMinus =
    static_cast<Double_t>(minusMultiplicitySum)/coincidenceEventCount;

  std::cout << "File: " << fileName << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\n=== Full coincidence summary ===" << std::endl;
  std::cout << "Total events: " << totalEventCount << std::endl;
  std::cout << "Coincidence events: "
            << coincidenceEventCount << std::endl;
  std::cout << "Coincidence efficiency: "
            << coincidenceEfficiency << " ("
            << 100.0*coincidenceEfficiency << " %)"
            << std::endl;
  std::cout << "Mean N(+z): " << meanPlus << std::endl;
  std::cout << "Mean N(-z): " << meanMinus << std::endl;

  std::cout << "\nN(+z) multiplicity distribution:" << std::endl;
  for (std::map<Int_t, Long64_t>::const_iterator it =
         plusMultiplicityDistribution.begin();
       it != plusMultiplicityDistribution.end(); ++it) {
    std::cout << "  N(+z)=" << it->first
              << ": " << it->second << " events" << std::endl;
  }

  std::cout << "\nN(-z) multiplicity distribution:" << std::endl;
  for (std::map<Int_t, Long64_t>::const_iterator it =
         minusMultiplicityDistribution.begin();
       it != minusMultiplicityDistribution.end(); ++it) {
    std::cout << "  N(-z)=" << it->first
              << ": " << it->second << " events" << std::endl;
  }

  const char* processNames[4] = {"phot", "compt", "Rayl", "conv"};
  std::cout << "\nGamma processes in coincidence events:" << std::endl;
  for (Int_t processIndex = 0; processIndex < 4; ++processIndex) {
    const std::string process = processNames[processIndex];
    std::cout << "  " << process
              << ": interactions=" << gammaProcessCounts[process]
              << ", events with process="
              << eventsByGammaProcess[process].size()
              << std::endl;
  }

  std::cout << "\nDetected electrons by MCP plate "
            << "(coincidence events):" << std::endl;
  for (std::set<Int_t>::const_iterator it =
         plusMcpIndices.begin();
       it != plusMcpIndices.end(); ++it) {
    const PlateKey plate(1, *it);
    std::cout << "  " << PlateLabel("+z", *it)
              << ": electrons=" << PlateCount(detectedElectronsByPlate,
                                               plate)
              << ", events=" << PlateCount(coincidenceEventsByPlate,
                                           plate)
              << std::endl;
  }
  for (std::set<Int_t>::const_iterator it =
         minusMcpIndices.begin();
       it != minusMcpIndices.end(); ++it) {
    const PlateKey plate(-1, *it);
    std::cout << "  " << PlateLabel("-z", *it)
              << ": electrons=" << PlateCount(detectedElectronsByPlate,
                                               plate)
              << ", events=" << PlateCount(coincidenceEventsByPlate,
                                           plate)
              << std::endl;
  }

  std::cout << "\nFirst gamma interactions by MCP plate "
            << "(coincidence events):" << std::endl;
  for (std::set<Int_t>::const_iterator it =
         plusMcpIndices.begin();
       it != plusMcpIndices.end(); ++it) {
    std::cout << "  " << PlateLabel("+z", *it) << ": "
              << PlateCount(firstGammaInteractionsByPlate,
                            PlateKey(1, *it))
              << std::endl;
  }
  for (std::set<Int_t>::const_iterator it =
         minusMcpIndices.begin();
       it != minusMcpIndices.end(); ++it) {
    std::cout << "  " << PlateLabel("-z", *it) << ": "
              << PlateCount(firstGammaInteractionsByPlate,
                            PlateKey(-1, *it))
              << std::endl;
  }

  std::cout << "\nPlate detection yield in coincidence events only"
            << std::endl;
  std::cout << std::left
            << std::setw(12) << "Plate"
            << std::right
            << std::setw(18) << "Gamma entries"
            << std::setw(24) << "Events with electron"
            << std::setw(20) << "Channel electrons"
            << std::setw(18) << "Event yield (%)"
            << std::setw(21) << "Electron yield (%)"
            << std::endl;

  std::vector<PlateKey> yieldPlateOrder;
  for (std::set<Int_t>::const_iterator it =
         plusMcpIndices.begin();
       it != plusMcpIndices.end(); ++it) {
    yieldPlateOrder.push_back(PlateKey(1, *it));
  }
  for (std::set<Int_t>::const_iterator it =
         minusMcpIndices.begin();
       it != minusMcpIndices.end(); ++it) {
    yieldPlateOrder.push_back(PlateKey(-1, *it));
  }

  std::vector<Double_t> eventYieldPercentByPlate;
  std::vector<Double_t> electronYieldPercentByPlate;
  for (std::vector<PlateKey>::const_iterator plateIt =
         yieldPlateOrder.begin();
       plateIt != yieldPlateOrder.end(); ++plateIt) {
    const Long64_t gammaEntries =
      PlateCount(gammaEntriesCoincidence, *plateIt);
    const Long64_t eventsWithElectron =
      PlateCount(coincidenceEventsByPlate, *plateIt);
    const Long64_t channelElectrons =
      PlateCount(plateStatsElectronCounts, *plateIt);

    const Double_t eventYieldPercent =
      gammaEntries > 0
        ? 100.0*static_cast<Double_t>(eventsWithElectron)/
            gammaEntries
        : 0.0;
    const Double_t electronYieldPercent =
      gammaEntries > 0
        ? 100.0*static_cast<Double_t>(channelElectrons)/
            gammaEntries
        : 0.0;

    eventYieldPercentByPlate.push_back(eventYieldPercent);
    electronYieldPercentByPlate.push_back(electronYieldPercent);

    std::cout << std::left
              << std::setw(12)
              << PlateLabel(plateIt->first > 0 ? "+z" : "-z",
                            plateIt->second)
              << std::right
              << std::setw(18) << gammaEntries
              << std::setw(24) << eventsWithElectron
              << std::setw(20) << channelElectrons
              << std::setw(18) << eventYieldPercent
              << std::setw(21) << electronYieldPercent
              << std::endl;
  }

  std::cout << "\nConsistency checks:" << std::endl;
  std::cout << "  Coincidence events with active plates on both sides: "
            << validPlateCoincidenceEvents << " / "
            << coincidenceEventCount
            << (validPlateCoincidenceEvents == coincidenceEventCount
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;

  Long64_t detailedElectronTotal = 0;
  Long64_t plateStatsElectronTotal = 0;
  for (std::map<PlateKey, Long64_t>::const_iterator it =
         detectedElectronsByPlate.begin();
       it != detectedElectronsByPlate.end(); ++it) {
    detailedElectronTotal += it->second;
  }
  for (std::map<PlateKey, Long64_t>::const_iterator it =
         plateStatsElectronCounts.begin();
       it != plateStatsElectronCounts.end(); ++it) {
    plateStatsElectronTotal += it->second;
  }
  std::cout << "  ElectronChannelHitTree vs McpPlateStatsTree: "
            << detailedElectronTotal << " vs "
            << plateStatsElectronTotal
            << (detailedElectronTotal == plateStatsElectronTotal
                  ? " [OK]" : " [MISMATCH]")
            << std::endl;
  std::cout << "  Total plate-pair contributions: "
            << totalPlatePairContributions << std::endl;

  // --------------------------------------------------------------------------
  // Figure 1: multiplicity distributions.
  // --------------------------------------------------------------------------

  const Int_t maximumMultiplicity =
    std::max(maximumPlusMultiplicity, maximumMinusMultiplicity);
  TH1D* plusMultiplicityHistogram =
    new TH1D("hFullCoincidencePlusMultiplicity",
             "Electron multiplicity in coincidence events;"
             "Electrons reaching MCP channels;Events",
             maximumMultiplicity + 1,
             -0.5,
             maximumMultiplicity + 0.5);
  TH1D* minusMultiplicityHistogram =
    new TH1D("hFullCoincidenceMinusMultiplicity",
             "Electron multiplicity in coincidence events;"
             "Electrons reaching MCP channels;Events",
             maximumMultiplicity + 1,
             -0.5,
             maximumMultiplicity + 0.5);

  TH2D* plusMinusPercent =
    new TH2D("hFullCoincidencePlusMinusPercent",
             "Electron multiplicity correlation in coincidence events (%);"
             "N(+z);N(-z)",
             maximumPlusMultiplicity + 1,
             -0.5,
             maximumPlusMultiplicity + 0.5,
             maximumMinusMultiplicity + 1,
             -0.5,
             maximumMinusMultiplicity + 0.5);

  for (std::vector<CoincidenceMultiplicity>::const_iterator it =
         multiplicities.begin();
       it != multiplicities.end(); ++it) {
    plusMultiplicityHistogram->Fill(it->plusCount);
    minusMultiplicityHistogram->Fill(it->minusCount);
    plusMinusPercent->Fill(it->plusCount, it->minusCount);
  }
  plusMinusPercent->Scale(
    100.0/static_cast<Double_t>(coincidenceEventCount));

  plusMultiplicityHistogram->SetLineColor(kBlue + 1);
  plusMultiplicityHistogram->SetLineWidth(2);
  minusMultiplicityHistogram->SetLineColor(kRed + 1);
  minusMultiplicityHistogram->SetLineWidth(2);

  gStyle->SetOptStat(0);
  TCanvas* multiplicityCanvas =
    new TCanvas("cFullCoincidenceMultiplicity",
                "Coincidence multiplicity",
                900,
                700);
  plusMultiplicityHistogram->SetMaximum(
    1.15*std::max(plusMultiplicityHistogram->GetMaximum(),
                  minusMultiplicityHistogram->GetMaximum()));
  plusMultiplicityHistogram->Draw("HIST");
  minusMultiplicityHistogram->Draw("HIST SAME");

  TLegend* multiplicityLegend =
    new TLegend(0.68, 0.74, 0.88, 0.88);
  multiplicityLegend->AddEntry(plusMultiplicityHistogram,
                               "+z", "l");
  multiplicityLegend->AddEntry(minusMultiplicityHistogram,
                               "-z", "l");
  multiplicityLegend->Draw();
  multiplicityCanvas->SaveAs("coincidence_multiplicity.png");

  // --------------------------------------------------------------------------
  // Figure 2: N(+z) versus N(-z), normalized to coincidence events.
  // --------------------------------------------------------------------------

  gStyle->SetPaintTextFormat("4.1f");
  TCanvas* correlationCanvas =
    new TCanvas("cFullCoincidencePlusMinusPercent",
                "Coincidence multiplicity correlation",
                900,
                750);
  correlationCanvas->SetRightMargin(0.15);
  plusMinusPercent->SetMarkerSize(1.3);
  plusMinusPercent->Draw("COLZ TEXT");
  correlationCanvas->SaveAs(
    "coincidence_plus_vs_minus_percent.png");

  // --------------------------------------------------------------------------
  // Figure 3: dynamic MCP plate-pair matrix.
  // --------------------------------------------------------------------------

  const std::vector<Int_t> plusIndices(plusMcpIndices.begin(),
                                       plusMcpIndices.end());
  const std::vector<Int_t> minusIndices(minusMcpIndices.begin(),
                                        minusMcpIndices.end());

  if (!plusIndices.empty() && !minusIndices.empty()) {
    TH2D* plateMatrix =
      new TH2D("hFullCoincidencePlateMatrix",
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
        PlateLabel("+z", plusIndices[plusBin]);
      plateMatrix->GetXaxis()->SetBinLabel(
        static_cast<Int_t>(plusBin + 1), label.c_str());
    }

    for (std::size_t minusBin = 0;
         minusBin < minusIndices.size();
         ++minusBin) {
      const std::string label =
        PlateLabel("-z", minusIndices[minusBin]);
      plateMatrix->GetYaxis()->SetBinLabel(
        static_cast<Int_t>(minusBin + 1), label.c_str());

      for (std::size_t plusBin = 0;
           plusBin < plusIndices.size();
           ++plusBin) {
        plateMatrix->SetBinContent(
          static_cast<Int_t>(plusBin + 1),
          static_cast<Int_t>(minusBin + 1),
          platePairCounts[
            std::make_pair(minusIndices[minusBin],
                           plusIndices[plusBin])]);
      }
    }

    gStyle->SetPaintTextFormat("g");
    TCanvas* plateCanvas =
      new TCanvas("cFullCoincidencePlateMatrix",
                  "Coincidence plate matrix",
                  900,
                  750);
    plateCanvas->SetLeftMargin(0.16);
    plateCanvas->SetBottomMargin(0.16);
    plateCanvas->SetRightMargin(0.15);
    plateMatrix->SetMarkerSize(1.5);
    plateMatrix->Draw("COLZ TEXT");
    plateCanvas->SaveAs("coincidence_plate_matrix.png");
  }

  // --------------------------------------------------------------------------
  // Figure 4: first gamma interaction by side and MCP plate.
  // --------------------------------------------------------------------------

  std::vector<PlateKey> gammaPlateOrder;
  for (std::set<Int_t>::const_iterator it =
         plusMcpIndices.begin();
       it != plusMcpIndices.end(); ++it) {
    gammaPlateOrder.push_back(PlateKey(1, *it));
  }
  for (std::set<Int_t>::const_iterator it =
         minusMcpIndices.begin();
       it != minusMcpIndices.end(); ++it) {
    gammaPlateOrder.push_back(PlateKey(-1, *it));
  }

  if (!gammaPlateOrder.empty()) {
    TH1D* firstGammaHistogram =
      new TH1D("hFullFirstGammaInteractionByPlate",
               "First gamma interaction in coincidence events;"
               "MCP plate;First gamma interactions",
               static_cast<Int_t>(gammaPlateOrder.size()),
               0.0,
               static_cast<Double_t>(gammaPlateOrder.size()));

    for (std::size_t bin = 0;
         bin < gammaPlateOrder.size();
         ++bin) {
      const PlateKey plate = gammaPlateOrder[bin];
      const std::string label =
        PlateLabel(plate.first > 0 ? "+z" : "-z",
                   plate.second);
      firstGammaHistogram->GetXaxis()->SetBinLabel(
        static_cast<Int_t>(bin + 1), label.c_str());
      firstGammaHistogram->SetBinContent(
        static_cast<Int_t>(bin + 1),
        PlateCount(firstGammaInteractionsByPlate, plate));
    }

    firstGammaHistogram->SetFillColor(kGreen + 2);
    firstGammaHistogram->SetLineColor(kGreen + 3);
    firstGammaHistogram->SetLineWidth(2);

    TCanvas* firstGammaCanvas =
      new TCanvas("cFullFirstGammaInteractionByPlate",
                  "First gamma interaction by MCP plate",
                  900,
                  700);
    firstGammaCanvas->SetBottomMargin(0.16);
    firstGammaHistogram->Draw("HIST TEXT0");
    firstGammaCanvas->SaveAs(
      "coincidence_gamma_first_interaction_by_plate.png");
  }

  // --------------------------------------------------------------------------
  // Figures 5 and 6: plate yields in coincidence events only.
  // --------------------------------------------------------------------------

  if (!yieldPlateOrder.empty()) {
    TH1D* eventYieldHistogram =
      new TH1D("hCoincidencePlateDetectionYield",
               "Plate detection yield in coincidence events only;"
               "MCP plate;Event yield (%)",
               static_cast<Int_t>(yieldPlateOrder.size()),
               0.0,
               static_cast<Double_t>(yieldPlateOrder.size()));
    TH1D* electronYieldHistogram =
      new TH1D("hCoincidencePlateElectronYield",
               "Plate electron yield in coincidence events only;"
               "MCP plate;Electron yield (%)",
               static_cast<Int_t>(yieldPlateOrder.size()),
               0.0,
               static_cast<Double_t>(yieldPlateOrder.size()));

    for (std::size_t bin = 0;
         bin < yieldPlateOrder.size();
         ++bin) {
      const PlateKey plate = yieldPlateOrder[bin];
      const std::string label =
        PlateLabel(plate.first > 0 ? "+z" : "-z",
                   plate.second);

      eventYieldHistogram->GetXaxis()->SetBinLabel(
        static_cast<Int_t>(bin + 1), label.c_str());
      eventYieldHistogram->SetBinContent(
        static_cast<Int_t>(bin + 1),
        eventYieldPercentByPlate[bin]);

      electronYieldHistogram->GetXaxis()->SetBinLabel(
        static_cast<Int_t>(bin + 1), label.c_str());
      electronYieldHistogram->SetBinContent(
        static_cast<Int_t>(bin + 1),
        electronYieldPercentByPlate[bin]);
    }

    gStyle->SetPaintTextFormat("4.2f");

    eventYieldHistogram->SetFillColor(kBlue + 1);
    eventYieldHistogram->SetLineColor(kBlue + 2);
    eventYieldHistogram->SetLineWidth(2);
    eventYieldHistogram->SetMinimum(0.0);
    eventYieldHistogram->SetMaximum(
      1.20*eventYieldHistogram->GetMaximum());

    TCanvas* eventYieldCanvas =
      new TCanvas("cCoincidencePlateDetectionYield",
                  "Coincidence plate detection yield",
                  900,
                  700);
    eventYieldCanvas->SetBottomMargin(0.16);
    eventYieldHistogram->Draw("HIST TEXT0");
    eventYieldCanvas->SaveAs(
      "coincidence_plate_detection_yield.png");

    electronYieldHistogram->SetFillColor(kGreen + 2);
    electronYieldHistogram->SetLineColor(kGreen + 3);
    electronYieldHistogram->SetLineWidth(2);
    electronYieldHistogram->SetMinimum(0.0);
    electronYieldHistogram->SetMaximum(
      1.20*electronYieldHistogram->GetMaximum());

    TCanvas* electronYieldCanvas =
      new TCanvas("cCoincidencePlateElectronYield",
                  "Coincidence plate electron yield",
                  900,
                  700);
    electronYieldCanvas->SetBottomMargin(0.16);
    electronYieldHistogram->Draw("HIST TEXT0");
    electronYieldCanvas->SaveAs(
      "coincidence_plate_electron_yield.png");
  }

  std::cout << "\nSaved figures:" << std::endl;
  std::cout << "  coincidence_multiplicity.png" << std::endl;
  std::cout << "  coincidence_plus_vs_minus_percent.png" << std::endl;
  std::cout << "  coincidence_plate_matrix.png" << std::endl;
  std::cout << "  coincidence_gamma_first_interaction_by_plate.png"
            << std::endl;
  std::cout << "  coincidence_plate_detection_yield.png"
            << std::endl;
  std::cout << "  coincidence_plate_electron_yield.png"
            << std::endl;

  file->Close();
  delete file;
}
