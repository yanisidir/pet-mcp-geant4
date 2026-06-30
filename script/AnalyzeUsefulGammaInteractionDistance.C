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
}  // namespace

void AnalyzeUsefulGammaInteractionDistance(
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
  std::string parentGammaTrackName;
  std::string creationXName;
  std::string creationYName;
  std::string creationZName;
  std::string creationTimeName;
  std::string validCreationName;

  TLeaf* eventLeaf = FindLeaf(
    electronTree, {"eventID", "event"}, eventName);
  TLeaf* parentGammaTrackLeaf = FindLeaf(
    electronTree,
    {"parentGammaTrackID", "parentID", "parent"},
    parentGammaTrackName);
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
  TLeaf* validCreationLeaf = FindLeaf(
    electronTree,
    {"hasValidCreationInfo"},
    validCreationName);

  if (!eventLeaf || !parentGammaTrackLeaf ||
      !creationXLeaf || !creationYLeaf || !creationZLeaf ||
      !creationTimeLeaf || !validCreationLeaf) {
    std::cerr
      << "ElectronChannelHitTree does not contain the exact creation "
      << "branches required by this analysis." << std::endl;
    std::cerr << "Detected branches: event=" << eventName
              << ", parentGammaTrack=" << parentGammaTrackName
              << ", creationX=" << creationXName
              << ", creationY=" << creationYName
              << ", creationZ=" << creationZName
              << ", creationTime=" << creationTimeName
              << ", hasValidCreationInfo=" << validCreationName
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  const Double_t coordinateScaleToMm =
    creationXName.find("_cm") != std::string::npos ? 10.0 : 1.0;

  std::cout << "Exact association based on detected electron creation "
            << "positions" << std::endl;
  std::cout << "  event branch: " << eventName << std::endl;
  std::cout << "  parent gamma branch: " << parentGammaTrackName
            << std::endl;
  std::cout << "  creation coordinates: " << creationXName << ", "
            << creationYName << ", " << creationZName << std::endl;
  std::cout << "  creation time: " << creationTimeName << std::endl;
  std::cout << "  validity flag: " << validCreationName << std::endl;
  std::cout << "  coordinate scale to mm: "
            << coordinateScaleToMm << std::endl;
  std::cout << "  Note: grouping is done by (eventID, parentGammaTrackID)."
            << std::endl;
  std::cout << "  Note: only rows with hasValidCreationInfo == true "
            << "are used." << std::endl;
  std::cout << "  primaryGammaTrackID == -1 means no known primary "
            << "gamma ancestor, typically Bremsstrahlung gammas from "
            << "secondary electrons." << std::endl;

  std::map<GammaKey, std::vector<DetectedElectronCreation> >
    creationsByGamma;

  const Long64_t electronEntries = electronTree->GetEntries();
  Long64_t skippedInvalidCreationRows = 0;
  for (Long64_t entry = 0; entry < electronEntries; ++entry) {
    electronTree->GetEntry(entry);

    const bool hasValidCreationInfo =
      static_cast<Int_t>(validCreationLeaf->GetValue()) != 0;
    if (!hasValidCreationInfo) {
      ++skippedInvalidCreationRows;
      continue;
    }

    const Int_t parentGammaTrackID =
      static_cast<Int_t>(parentGammaTrackLeaf->GetValue());
    if (parentGammaTrackID <= 0) {
      continue;
    }

    DetectedElectronCreation creation;
    creation.xMm = coordinateScaleToMm*creationXLeaf->GetValue();
    creation.yMm = coordinateScaleToMm*creationYLeaf->GetValue();
    creation.zMm = coordinateScaleToMm*creationZLeaf->GetValue();
    creation.timeNs = creationTimeLeaf->GetValue();

    creationsByGamma[
      GammaKey(static_cast<Int_t>(eventLeaf->GetValue()),
               parentGammaTrackID)].push_back(creation);
  }

  std::vector<Double_t> distance3dMm;
  Long64_t groupedGammas = 0;

  for (std::map<GammaKey, std::vector<DetectedElectronCreation> >::iterator
         it = creationsByGamma.begin();
       it != creationsByGamma.end(); ++it) {
    std::vector<DetectedElectronCreation>& creations = it->second;
    std::sort(creations.begin(), creations.end(), EarlierCreation);

    if (creations.size() < 2) {
      continue;
    }
    groupedGammas++;

    for (std::size_t index = 1; index < creations.size(); ++index) {
      const DetectedElectronCreation& previous = creations[index - 1];
      const DetectedElectronCreation& current = creations[index];
      const Double_t dx = current.xMm - previous.xMm;
      const Double_t dy = current.yMm - previous.yMm;
      const Double_t dz = current.zMm - previous.zMm;
      distance3dMm.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
    }
  }

  std::cout << "\nExact creation-separation summary" << std::endl;
  std::cout << "  Electron entries: " << electronEntries << std::endl;
  std::cout << "  Skipped invalid creation rows: "
            << skippedInvalidCreationRows << std::endl;
  std::cout << "  (eventID, parentGammaTrackID) groups: "
            << creationsByGamma.size() << std::endl;
  std::cout << "  Groups with >= 2 detected electrons: "
            << groupedGammas << std::endl;
  std::cout << "  Successive useful pairs: "
            << distance3dMm.size() << std::endl;

  if (distance3dMm.empty()) {
    std::cerr << "\nNo parent gamma produced at least two detected "
              << "electrons. Nothing to plot." << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nDistance statistics [mm]" << std::endl;
  std::cout << "  Mean   : " << ComputeMean(distance3dMm) << std::endl;
  std::cout << "  RMS    : " << ComputeRms(distance3dMm) << std::endl;
  std::cout << "  Median : " << ComputeMedian(distance3dMm) << std::endl;
  std::cout << "  Min    : "
            << *std::min_element(distance3dMm.begin(),
                                 distance3dMm.end()) << std::endl;
  std::cout << "  Max    : "
            << *std::max_element(distance3dMm.begin(),
                                 distance3dMm.end()) << std::endl;
  std::cout << "  68%    : " << ComputeQuantile(distance3dMm, 0.68)
            << std::endl;
  std::cout << "  95%    : " << ComputeQuantile(distance3dMm, 0.95)
            << std::endl;

  const Double_t maximum =
    *std::max_element(distance3dMm.begin(), distance3dMm.end());
  const Double_t histogramMaximum = maximum > 0.0 ? 1.05*maximum : 1.0;

  TH1D* histogram = new TH1D(
    "hUsefulGammaInteractionDistance",
    "Distance between successive useful interaction sites;"
    "Distance [mm];Number of pairs",
    100,
    0.0,
    histogramMaximum);
  histogram->SetDirectory(nullptr);
  histogram->SetLineColor(kBlue + 1);
  histogram->SetFillColorAlpha(kBlue + 1, 0.35);
  histogram->SetLineWidth(2);

  for (std::size_t index = 0; index < distance3dMm.size(); ++index) {
    histogram->Fill(distance3dMm[index]);
  }

  TCanvas* canvas = new TCanvas(
    "cUsefulGammaInteractionDistance",
    "Distance between successive useful interaction sites",
    1100,
    800);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);
  canvas->SetGridy();

  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.35);
  histogram->Draw("HIST");
  canvas->SaveAs("Fig/current/gamma_useful_interaction_distance.png");

  file->Close();
  delete file;
}
