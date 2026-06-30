#include "TAxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLeaf.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{
typedef std::pair<Int_t, Int_t> GammaKey;

struct DetectedElectronCreation
{
  Double_t xMm;
  Double_t yMm;
  Double_t zMm;
  Double_t timeNs;
  std::string creatorProcessName;
  Int_t parentGammaTrackID;
};

struct PairTransitionStats
{
  Long64_t count;
  Double_t sumAbsDzMm;
  Double_t sumRxyMm;
  Double_t sumDistance3dMm;
  std::vector<Double_t> distance3dMm;

  PairTransitionStats()
    : count(0),
      sumAbsDzMm(0.0),
      sumRxyMm(0.0),
      sumDistance3dMm(0.0)
  {
  }
};

TLeaf* FindLeaf(TTree* tree,
                const std::vector<std::string>& candidates,
                std::string& selectedName)
{
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    if (tree->GetBranch(candidates[i].c_str())) {
      selectedName = candidates[i];
      return tree->GetLeaf(selectedName.c_str());
    }
  }
  selectedName = "";
  return nullptr;
}

Double_t ComputeMean(const std::vector<Double_t>& values)
{
  if (values.empty()) {
    return 0.0;
  }

  Double_t sum = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    sum += values[i];
  }
  return sum/values.size();
}

Double_t ComputeRms(const std::vector<Double_t>& values)
{
  if (values.empty()) {
    return 0.0;
  }

  const Double_t mean = ComputeMean(values);
  Double_t squaredDifferenceSum = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    const Double_t difference = values[i] - mean;
    squaredDifferenceSum += difference*difference;
  }
  return std::sqrt(squaredDifferenceSum/values.size());
}

Double_t ComputeQuantile(std::vector<Double_t> values, Double_t fraction)
{
  if (values.empty()) {
    return 0.0;
  }

  std::sort(values.begin(), values.end());
  if (fraction <= 0.0) {
    return values.front();
  }
  if (fraction >= 1.0) {
    return values.back();
  }

  const Double_t realIndex = fraction*(values.size() - 1);
  const std::size_t lowerIndex =
    static_cast<std::size_t>(std::floor(realIndex));
  const std::size_t upperIndex =
    static_cast<std::size_t>(std::ceil(realIndex));
  const Double_t weight = realIndex - lowerIndex;
  return (1.0 - weight)*values[lowerIndex] +
         weight*values[upperIndex];
}

Double_t ComputeMedian(const std::vector<Double_t>& values)
{
  return ComputeQuantile(values, 0.5);
}

bool EarlierCreation(const DetectedElectronCreation& first,
                     const DetectedElectronCreation& second)
{
  if (first.timeNs != second.timeNs) {
    return first.timeNs < second.timeNs;
  }
  return first.zMm < second.zMm;
}

void PrintStatisticLine(const std::string& quantity,
                        const std::vector<Double_t>& valuesMm)
{
  std::cout
    << std::setw(12) << quantity
    << " | " << std::setw(10) << ComputeMean(valuesMm)
    << " | " << std::setw(10) << ComputeRms(valuesMm)
    << " | " << std::setw(12) << ComputeMedian(valuesMm)
    << " | " << std::setw(9)
    << (valuesMm.empty()
          ? 0.0
          : *std::min_element(valuesMm.begin(), valuesMm.end()))
    << " | " << std::setw(9)
    << (valuesMm.empty()
          ? 0.0
          : *std::max_element(valuesMm.begin(), valuesMm.end()))
    << " | " << std::setw(20) << ComputeQuantile(valuesMm, 0.68)
    << " | " << std::setw(20) << ComputeQuantile(valuesMm, 0.95)
    << std::endl;
}

void DrawHistogram(const std::vector<Double_t>& valuesMm,
                   const std::string& histogramName,
                   const std::string& title,
                   const std::string& xAxisTitle,
                   const std::string& outputName,
                   Color_t color)
{
  const Double_t maximum =
    valuesMm.empty()
      ? 1.0
      : *std::max_element(valuesMm.begin(), valuesMm.end());
  const Double_t histogramMaximum = maximum > 0.0 ? 1.05*maximum : 1.0;

  TH1D* histogram = new TH1D(
    histogramName.c_str(), title.c_str(), 100, 0.0, histogramMaximum);
  histogram->SetDirectory(nullptr);
  histogram->SetLineColor(color);
  histogram->SetFillColorAlpha(color, 0.35);
  histogram->SetLineWidth(2);

  for (std::size_t index = 0; index < valuesMm.size(); ++index) {
    histogram->Fill(valuesMm[index]);
  }

  TCanvas* canvas = new TCanvas(
    (histogramName + "Canvas").c_str(), title.c_str(), 1100, 800);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);
  canvas->SetGridy();

  histogram->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histogram->GetYaxis()->SetTitle("Number of electron pairs");
  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.35);
  histogram->Draw("HIST");
  canvas->SaveAs(outputName.c_str());
}

void AddPairTransition(std::map<std::string, PairTransitionStats>& stats,
                       const std::string& transitionName,
                       Double_t absDzMm,
                       Double_t rxyMm,
                       Double_t distance3dMm)
{
  PairTransitionStats& transitionStats = stats[transitionName];
  ++transitionStats.count;
  transitionStats.sumAbsDzMm += absDzMm;
  transitionStats.sumRxyMm += rxyMm;
  transitionStats.sumDistance3dMm += distance3dMm;
  transitionStats.distance3dMm.push_back(distance3dMm);
}

void PrintPairTransitionStats(
  const std::map<std::string, PairTransitionStats>& stats)
{
  std::cout << "\nSuccessive-pair process transition diagnostic"
            << std::endl;
  std::cout
    << "  Transition |      Pairs | Mean |dz| [mm] | Mean rxy [mm]"
    << " | Mean 3D [mm] | Median 3D [mm]" << std::endl;

  for (std::map<std::string, PairTransitionStats>::const_iterator it =
         stats.begin();
       it != stats.end(); ++it) {
    const PairTransitionStats& transitionStats = it->second;
    const Double_t count =
      static_cast<Double_t>(transitionStats.count);
    std::cout
      << "  " << std::setw(10) << it->first
      << " | " << std::setw(10) << transitionStats.count
      << " | " << std::setw(14)
      << (count > 0.0 ? transitionStats.sumAbsDzMm/count : 0.0)
      << " | " << std::setw(13)
      << (count > 0.0 ? transitionStats.sumRxyMm/count : 0.0)
      << " | " << std::setw(12)
      << (count > 0.0 ? transitionStats.sumDistance3dMm/count : 0.0)
      << " | " << std::setw(14)
      << ComputeMedian(transitionStats.distance3dMm)
      << std::endl;
  }
}
}  // namespace

void AnalyzeDetectedElectronCreationSeparationByPrimaryGamma(
  const char* fileName = "build/mcp_output.root")
{
  gSystem->mkdir("Fig/current", kTRUE);

  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* electronTree = nullptr;
  file->GetObject("ElectronChannelHitTree", electronTree);
  if (!electronTree) {
    std::cerr << "Missing ElectronChannelHitTree in " << fileName
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::string eventName;
  std::string primaryGammaTrackName;
  std::string validCreationName;
  std::string creationXName;
  std::string creationYName;
  std::string creationZName;
  std::string creationTimeName;
  std::string creationSideName;
  std::string creationMcpIndexName;
  std::string creatorProcessName;
  std::string parentGammaTrackName;

  TLeaf* eventLeaf = FindLeaf(
    electronTree, {"eventID", "event"}, eventName);
  TLeaf* primaryGammaTrackLeaf = FindLeaf(
    electronTree,
    {"primaryGammaTrackID"},
    primaryGammaTrackName);
  TLeaf* validCreationLeaf = FindLeaf(
    electronTree,
    {"hasValidCreationInfo"},
    validCreationName);
  TLeaf* creationXLeaf = FindLeaf(
    electronTree,
    {"creationX_mm", "creationX_cm", "creationX"},
    creationXName);
  TLeaf* creationYLeaf = FindLeaf(
    electronTree,
    {"creationY_mm", "creationY_cm", "creationY"},
    creationYName);
  TLeaf* creationZLeaf = FindLeaf(
    electronTree,
    {"creationZ_mm", "creationZ_cm", "creationZ"},
    creationZName);
  TLeaf* creationTimeLeaf = FindLeaf(
    electronTree,
    {"creationTime_ns", "creationTime", "globalTime_ns"},
    creationTimeName);
  TLeaf* creationSideLeaf = FindLeaf(
    electronTree,
    {"creationSide"},
    creationSideName);
  TLeaf* creationMcpIndexLeaf = FindLeaf(
    electronTree,
    {"creationMcpIndex"},
    creationMcpIndexName);
  TLeaf* creatorProcessLeaf = FindLeaf(
    electronTree,
    {"creatorProcessName"},
    creatorProcessName);
  TLeaf* parentGammaTrackLeaf = FindLeaf(
    electronTree,
    {"parentGammaTrackID", "parentID", "parent"},
    parentGammaTrackName);

  if (!eventLeaf || !primaryGammaTrackLeaf || !validCreationLeaf ||
      !creationXLeaf || !creationYLeaf || !creationZLeaf ||
      !creationTimeLeaf || !creationSideLeaf ||
      !creationMcpIndexLeaf || !creatorProcessLeaf ||
      !parentGammaTrackLeaf) {
    std::cerr
      << "ElectronChannelHitTree does not contain the branches "
      << "required for primary-gamma grouping." << std::endl;
    std::cerr << "Detected branches: event=" << eventName
              << ", primaryGammaTrackID=" << primaryGammaTrackName
              << ", hasValidCreationInfo=" << validCreationName
              << ", creationX=" << creationXName
              << ", creationY=" << creationYName
              << ", creationZ=" << creationZName
              << ", creationTime=" << creationTimeName
              << ", creationSide=" << creationSideName
              << ", creationMcpIndex=" << creationMcpIndexName
              << ", creatorProcessName=" << creatorProcessName
              << ", parentGammaTrackID=" << parentGammaTrackName
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  char creatorProcessBuffer[64];
  creatorProcessBuffer[0] = '\0';
  electronTree->SetBranchAddress("creatorProcessName",
                                 creatorProcessBuffer);

  const Double_t coordinateScaleToMm =
    creationXName.find("_cm") != std::string::npos ? 10.0 : 1.0;

  std::cout << "Electron creation branches:" << std::endl;
  std::cout << "  event: " << eventName << std::endl;
  std::cout << "  primary gamma track: " << primaryGammaTrackName
            << std::endl;
  std::cout << "  parent gamma track: " << parentGammaTrackName
            << std::endl;
  std::cout << "  validity flag: " << validCreationName << std::endl;
  std::cout << "  creation coordinates: " << creationXName << ", "
            << creationYName << ", " << creationZName << std::endl;
  std::cout << "  creation time: " << creationTimeName << std::endl;
  std::cout << "  creation plate: " << creationSideName << ", "
            << creationMcpIndexName << std::endl;
  std::cout << "  creator process: " << creatorProcessName << std::endl;
  std::cout << "  coordinate scale to mm: "
            << coordinateScaleToMm << std::endl;
  std::cout << "  This analysis groups detected electrons by "
            << "primaryGammaTrackID, i.e. by the initial PET gamma "
            << "ancestor when known." << std::endl;
  std::cout << "  primaryGammaTrackID == -1 means no known primary "
            << "gamma ancestor, typically Bremsstrahlung gammas from "
            << "secondary electrons." << std::endl;

  std::map<GammaKey, std::vector<DetectedElectronCreation> >
    creationsByPrimaryGamma;

  const Long64_t electronEntries = electronTree->GetEntries();
  Long64_t skippedInvalidCreationRows = 0;
  Long64_t skippedUnknownPrimaryGammaRows = 0;

  // Loop over all detected electrons and group them by (eventID, primaryGammaTrackID)
  for (Long64_t entry = 0; entry < electronEntries; ++entry) {
    electronTree->GetEntry(entry);

    // Skip rows without valid creation information
    const bool hasValidCreationInfo =
      static_cast<Int_t>(validCreationLeaf->GetValue()) != 0;
    if (!hasValidCreationInfo) {
      ++skippedInvalidCreationRows;
      continue;
    }

    // Only consider electrons with a known primary gamma ancestor
    const Int_t primaryGammaTrackID =
      static_cast<Int_t>(primaryGammaTrackLeaf->GetValue());
    if (primaryGammaTrackID <= 0) {
      ++skippedUnknownPrimaryGammaRows;
      continue;
    }

    DetectedElectronCreation creation;
    creation.xMm = coordinateScaleToMm*creationXLeaf->GetValue();
    creation.yMm = coordinateScaleToMm*creationYLeaf->GetValue();
    creation.zMm = coordinateScaleToMm*creationZLeaf->GetValue();
    creation.timeNs = creationTimeLeaf->GetValue();
    creation.creatorProcessName = creatorProcessBuffer;
    creation.parentGammaTrackID =
      static_cast<Int_t>(parentGammaTrackLeaf->GetValue());

    creationsByPrimaryGamma[
      GammaKey(static_cast<Int_t>(eventLeaf->GetValue()),
               primaryGammaTrackID)].push_back(creation);
  }

  std::vector<Double_t> absDzMm;
  std::vector<Double_t> rxyMm;
  std::vector<Double_t> distance3dMm;
  Long64_t groupedGammas = 0;
  Long64_t electronsInMultiElectronGroups = 0;
  Long64_t onlyComptGroups = 0;
  Long64_t onlyPhotGroups = 0;
  Long64_t mixedPhotComptConvGroups = 0;
  Long64_t groupsWithOtherCreatorProcess = 0;
  std::map<Int_t, Long64_t> groupSizeCounts;
  std::map<std::string, Long64_t> creatorProcessCountsInGroups;
  std::map<std::string, PairTransitionStats> pairTransitionStats;
  std::map<std::string, PairTransitionStats> parentTrackPairStats;
  std::map<std::string, PairTransitionStats> photPhotParentTrackPairStats;
  Long64_t pairsDistanceBelow001um = 0;  // distance < 0.001 mm
  Long64_t pairsDistanceBelow01um = 0;   // distance < 0.01 mm
  Long64_t pairsDistanceBelow1um = 0;    // distance < 0.1 mm
  Long64_t groupsWithSingleParentGamma = 0;
  Long64_t groupsWithMultipleParentGammas = 0;
  Long64_t onlyPhotParentGammaDistinctSum = 0;

  for (std::map<GammaKey, std::vector<DetectedElectronCreation> >::iterator
         it = creationsByPrimaryGamma.begin();
       it != creationsByPrimaryGamma.end(); ++it) {
    std::vector<DetectedElectronCreation>& creations = it->second;
    std::sort(creations.begin(), creations.end(), EarlierCreation);

    if (creations.size() < 2) {
      continue;
    }
    groupedGammas++;
    electronsInMultiElectronGroups +=
      static_cast<Long64_t>(creations.size());
    ++groupSizeCounts[static_cast<Int_t>(creations.size())];

    bool hasPhot = false;
    bool hasCompt = false;
    bool hasConv = false;
    bool hasOther = false;
    std::set<Int_t> distinctParentGammaTrackIDs;
    for (std::size_t creationIndex = 0;
         creationIndex < creations.size();
         ++creationIndex) {
      const std::string& processName =
        creations[creationIndex].creatorProcessName;
      ++creatorProcessCountsInGroups[processName];
      if (processName == "phot") {
        hasPhot = true;
      } else if (processName == "compt") {
        hasCompt = true;
      } else if (processName == "conv") {
        hasConv = true;
      } else {
        hasOther = true;
      }
      distinctParentGammaTrackIDs.insert(
        creations[creationIndex].parentGammaTrackID);
    }

    const Int_t knownProcessCount =
      (hasPhot ? 1 : 0) + (hasCompt ? 1 : 0) + (hasConv ? 1 : 0);
    if (hasCompt && !hasPhot && !hasConv && !hasOther) {
      ++onlyComptGroups;
    } else if (hasPhot && !hasCompt && !hasConv && !hasOther) {
      ++onlyPhotGroups;
      onlyPhotParentGammaDistinctSum +=
        static_cast<Long64_t>(distinctParentGammaTrackIDs.size());
    } else if (knownProcessCount > 1) {
      ++mixedPhotComptConvGroups;
    }
    if (hasOther) {
      ++groupsWithOtherCreatorProcess;
    }

    if (distinctParentGammaTrackIDs.size() == 1) {
      ++groupsWithSingleParentGamma;
    } else {
      ++groupsWithMultipleParentGammas;
    }

    for (std::size_t index = 1; index < creations.size(); ++index) {
      const DetectedElectronCreation& previous = creations[index - 1];
      const DetectedElectronCreation& current = creations[index];
      const Double_t dx = current.xMm - previous.xMm;
      const Double_t dy = current.yMm - previous.yMm;
      const Double_t dz = current.zMm - previous.zMm;
      const Double_t rxy = std::sqrt(dx*dx + dy*dy);
      const Double_t distance3d =
        std::sqrt(dx*dx + dy*dy + dz*dz);
      const Double_t absDz = std::abs(dz);
      const std::string transitionName =
        previous.creatorProcessName + "->" +
        current.creatorProcessName;
      const std::string parentTrackCategory =
        previous.parentGammaTrackID == current.parentGammaTrackID
          ? "same parentGammaTrackID"
          : "different parentGammaTrackID";

      absDzMm.push_back(absDz);
      rxyMm.push_back(rxy);
      distance3dMm.push_back(distance3d);
      AddPairTransition(pairTransitionStats,
                        transitionName,
                        absDz,
                        rxy,
                        distance3d);
      AddPairTransition(parentTrackPairStats,
                        parentTrackCategory,
                        absDz,
                        rxy,
                        distance3d);
      if (transitionName == "phot->phot") {
        AddPairTransition(photPhotParentTrackPairStats,
                          parentTrackCategory,
                          absDz,
                          rxy,
                          distance3d);
      }

      if (distance3d < 0.001) {
        ++pairsDistanceBelow001um;
      }
      if (distance3d < 0.01) {
        ++pairsDistanceBelow01um;
      }
      if (distance3d < 0.1) {
        ++pairsDistanceBelow1um;
      }
    }
  }

  std::cout << "\nPrimary-gamma creation-position association summary"
            << std::endl;
  std::cout << "  Total electron entries: " << electronEntries
            << std::endl;
  std::cout << "  Skipped invalid creation rows: "
            << skippedInvalidCreationRows << std::endl;
  std::cout << "  Skipped primaryGammaTrackID == -1 rows: "
            << skippedUnknownPrimaryGammaRows << std::endl;
  std::cout << "  (eventID, primaryGammaTrackID) groups: "
            << creationsByPrimaryGamma.size() << std::endl;
  std::cout << "  Groups with >= 2 detected electrons: "
            << groupedGammas << std::endl;
  std::cout << "  Successive useful pairs: "
            << absDzMm.size() << std::endl;

  std::cout << "\nMulti-electron group creator-process diagnostic"
            << std::endl;
  std::cout << "  Groups analysed: " << groupedGammas << std::endl;
  std::cout << "  Mean electrons per analysed group: "
            << (groupedGammas > 0
                  ? static_cast<Double_t>(electronsInMultiElectronGroups)/
                    groupedGammas
                  : 0.0)
            << std::endl;
  std::cout << "  Groups containing only compt: "
            << onlyComptGroups << std::endl;
  std::cout << "  Groups containing only phot: "
            << onlyPhotGroups << std::endl;
  std::cout << "  Mixed phot/compt/conv groups: "
            << mixedPhotComptConvGroups << std::endl;
  std::cout << "  Groups with another creator process: "
            << groupsWithOtherCreatorProcess << std::endl;

  std::cout << "  Group size distribution:" << std::endl;
  for (std::map<Int_t, Long64_t>::const_iterator it =
         groupSizeCounts.begin();
       it != groupSizeCounts.end(); ++it) {
    std::cout << "    N=" << it->first << ": "
              << it->second << " groups" << std::endl;
  }

  std::cout << "  Electron creatorProcessName distribution in these groups:"
            << std::endl;
  for (std::map<std::string, Long64_t>::const_iterator it =
         creatorProcessCountsInGroups.begin();
       it != creatorProcessCountsInGroups.end(); ++it) {
    std::cout << "    " << it->first << ": "
              << it->second << " electrons" << std::endl;
  }

  std::cout << "  Groups with one distinct parentGammaTrackID: "
            << groupsWithSingleParentGamma << std::endl;
  std::cout << "  Groups with multiple distinct parentGammaTrackID: "
            << groupsWithMultipleParentGammas << std::endl;
  std::cout << "  Mean distinct parentGammaTrackID in only-phot groups: "
            << (onlyPhotGroups > 0
                  ? static_cast<Double_t>(onlyPhotParentGammaDistinctSum)/
                    onlyPhotGroups
                  : 0.0)
            << std::endl;

  PrintPairTransitionStats(pairTransitionStats);

  std::cout << "\nSuccessive-pair parentGammaTrackID diagnostic"
            << std::endl;
  PrintPairTransitionStats(parentTrackPairStats);

  std::cout << "\nphot->phot pairs split by parentGammaTrackID"
            << std::endl;
  PrintPairTransitionStats(photPhotParentTrackPairStats);

  std::cout << "\nQuasi-superposed successive pairs" << std::endl;
  std::cout << "  distance3D < 0.001 mm: "
            << pairsDistanceBelow001um << std::endl;
  std::cout << "  distance3D < 0.01 mm : "
            << pairsDistanceBelow01um << std::endl;
  std::cout << "  distance3D < 0.1 mm  : "
            << pairsDistanceBelow1um << std::endl;

  if (absDzMm.empty()) {
    std::cerr << "\nNo primary gamma group has at least two detected "
              << "electrons with recorded creation positions."
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nDetected-electron creation separation by primary gamma"
            << std::endl;
  std::cout
    << std::setw(12) << "Quantity"
    << " | " << std::setw(10) << "Mean [mm]"
    << " | " << std::setw(10) << "RMS [mm]"
    << " | " << std::setw(12) << "Median [mm]"
    << " | " << std::setw(9) << "Min [mm]"
    << " | " << std::setw(9) << "Max [mm]"
    << " | " << std::setw(20) << "68% containment [mm]"
    << " | " << std::setw(20) << "95% containment [mm]"
    << std::endl;
  PrintStatisticLine("|dz|", absDzMm);
  PrintStatisticLine("rxy", rxyMm);
  PrintStatisticLine("distance3D", distance3dMm);

  DrawHistogram(
    absDzMm,
    "hDetectedElectronCreationAbsDzByPrimaryGamma",
    "Detected-electron creation |#Delta z| by primary gamma",
    "|#Delta z| [mm]",
    "Fig/current/detected_electron_creation_by_primary_gamma_dz.png",
    kBlue + 1);

  DrawHistogram(
    rxyMm,
    "hDetectedElectronCreationRxyByPrimaryGamma",
    "Detected-electron creation transverse separation by primary gamma",
    "r_{xy} [mm]",
    "Fig/current/detected_electron_creation_by_primary_gamma_rxy.png",
    kRed + 1);

  DrawHistogram(
    distance3dMm,
    "hDetectedElectronCreationDistance3dByPrimaryGamma",
    "Detected-electron creation 3D separation by primary gamma",
    "Distance [mm]",
    "Fig/current/detected_electron_creation_by_primary_gamma_distance3d.png",
    kGreen + 2);

  std::cout << "\nSaved:" << std::endl;
  std::cout << "  detected_electron_creation_by_primary_gamma_dz.png"
            << std::endl;
  std::cout << "  detected_electron_creation_by_primary_gamma_rxy.png"
            << std::endl;
  std::cout << "  detected_electron_creation_by_primary_gamma_distance3d.png"
            << std::endl;

  file->Close();
  delete file;
}
