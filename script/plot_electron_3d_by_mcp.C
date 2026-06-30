#include "TAxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
const Int_t kPlusSide = 1;
const Int_t kNumberOfMCPs = 3;
const Int_t kNumberOfBins = 100;
const Int_t kNumberOfVariables = 6;
const Bool_t kUseLogScaleForZoom = false;

enum VariableIndex
{
  kX = 0,
  kY,
  kZ,
  kDirectionX,
  kDirectionY,
  kDirectionZ
};

struct PlateData
{
  Int_t mcpIndex;
  Color_t color;
  std::string label;
  std::vector<Double_t> values[kNumberOfVariables];
};

std::string FindBranch(TTree* tree,
                       const std::vector<std::string>& candidates)
{
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    if (tree->GetBranch(candidates[i].c_str())) {
      return candidates[i];
    }
  }
  return "";
}

TLeaf* FindLeaf(TTree* tree,
                const std::vector<std::string>& candidates,
                std::string& selectedName)
{
  selectedName = FindBranch(tree, candidates);
  if (selectedName.empty()) {
    return nullptr;
  }
  return tree->GetLeaf(selectedName.c_str());
}

PlateData* FindPlate(std::vector<PlateData>& plates, Int_t mcpIndex)
{
  for (std::size_t i = 0; i < plates.size(); ++i) {
    if (plates[i].mcpIndex == mcpIndex) {
      return &plates[i];
    }
  }
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

void ComputeRange(const std::vector<PlateData>& plates,
                  Int_t variableIndex,
                  Double_t& minimum,
                  Double_t& maximum)
{
  minimum = std::numeric_limits<Double_t>::max();
  maximum = -std::numeric_limits<Double_t>::max();

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::vector<Double_t>& values =
      plates[plate].values[variableIndex];
    for (std::size_t value = 0; value < values.size(); ++value) {
      minimum = std::min(minimum, values[value]);
      maximum = std::max(maximum, values[value]);
    }
  }

  if (minimum > maximum) {
    minimum = 0.0;
    maximum = 1.0;
    return;
  }

  Double_t padding = 0.05*(maximum - minimum);
  if (padding == 0.0) {
    padding = std::max(0.05*std::abs(minimum), 0.5);
  }
  minimum -= padding;
  maximum += padding;
}

void ComputePercentileRange(const std::vector<PlateData>& plates,
                            Int_t variableIndex,
                            Double_t& minimum,
                            Double_t& maximum)
{
  std::vector<Double_t> allValues;
  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::vector<Double_t>& values =
      plates[plate].values[variableIndex];
    allValues.insert(allValues.end(), values.begin(), values.end());
  }

  if (allValues.empty()) {
    minimum = 0.0;
    maximum = 1.0;
    return;
  }

  std::sort(allValues.begin(), allValues.end());
  const std::size_t lastIndex = allValues.size() - 1;
  const std::size_t lowerIndex =
    static_cast<std::size_t>(0.01*lastIndex);
  const std::size_t upperIndex =
    static_cast<std::size_t>(0.99*lastIndex);

  minimum = allValues[lowerIndex];
  maximum = allValues[upperIndex];

  if (minimum >= maximum) {
    ComputeRange(plates, variableIndex, minimum, maximum);
  }
}

void DrawNormalizedHistograms(const std::vector<PlateData>& plates,
                              Int_t variableIndex,
                              const std::string& canvasName,
                              const std::string& title,
                              const std::string& xAxisTitle,
                              const std::string& outputName,
                              Bool_t usePercentileRange = false,
                              Bool_t useLogY = false)
{
  Double_t minimum = 0.0;
  Double_t maximum = 1.0;
  if (usePercentileRange) {
    ComputePercentileRange(plates, variableIndex, minimum, maximum);
    std::cout << "  " << outputName << " range (1%-99%): ["
              << minimum << ", " << maximum << "]" << std::endl;
  } else {
    ComputeRange(plates, variableIndex, minimum, maximum);
  }

  std::vector<TH1D*> histograms;
  Double_t largestBinContent = 0.0;
  Double_t smallestPositiveBinContent =
    std::numeric_limits<Double_t>::max();

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::string histogramName =
      canvasName + "_mcp_" + std::to_string(plates[plate].mcpIndex);
    TH1D* histogram = new TH1D(
      histogramName.c_str(), "", kNumberOfBins, minimum, maximum);
    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(plates[plate].color);
    histogram->SetLineWidth(3);

    const std::vector<Double_t>& values =
      plates[plate].values[variableIndex];
    for (std::size_t value = 0; value < values.size(); ++value) {
      histogram->Fill(values[value]);
    }

    if (histogram->Integral() > 0.0) {
      histogram->Scale(1.0/histogram->Integral());
      largestBinContent =
        std::max(largestBinContent, histogram->GetMaximum());
      for (Int_t bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
        const Double_t binContent = histogram->GetBinContent(bin);
        if (binContent > 0.0) {
          smallestPositiveBinContent =
            std::min(smallestPositiveBinContent, binContent);
        }
      }
    }
    histograms.push_back(histogram);
  }

  TCanvas* canvas =
    new TCanvas(canvasName.c_str(), title.c_str(), 1100, 800);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);
  canvas->SetGridy();
  if (useLogY) {
    gPad->SetLogy();
  }

  const std::string frameName = canvasName + "_frame";
  TH1D* frame = new TH1D(
    frameName.c_str(), title.c_str(), kNumberOfBins, minimum, maximum);
  frame->SetDirectory(nullptr);
  frame->SetMinimum(
    useLogY &&
    smallestPositiveBinContent < std::numeric_limits<Double_t>::max()
      ? 0.5*smallestPositiveBinContent
      : 0.0);
  frame->SetMaximum(
    largestBinContent > 0.0 ? 1.20*largestBinContent : 1.0);
  frame->GetXaxis()->SetTitle(xAxisTitle.c_str());
  frame->GetYaxis()->SetTitle("Normalized entries");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->Draw("HIST");

  const Bool_t isDirection =
    variableIndex == kDirectionX ||
    variableIndex == kDirectionY ||
    variableIndex == kDirectionZ;
  TLegend* legend = isDirection
    ? new TLegend(0.15, 0.65, 0.40, 0.88)
    : new TLegend(0.70, 0.68, 0.94, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(isDirection ? 0.032 : 0.035);

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    if (histograms[plate]->GetEntries() > 0.0) {
      histograms[plate]->Draw("HIST SAME");
    }

    const std::string legendLabel = isDirection
      ? plates[plate].label
      : plates[plate].label + " (N=" +
        std::to_string(plates[plate].values[variableIndex].size()) + ")";
    legend->AddEntry(histograms[plate], legendLabel.c_str(), "l");
  }

  legend->Draw();
  canvas->RedrawAxis();
  canvas->SaveAs(outputName.c_str());
}

void DrawElectronCounts(const std::vector<PlateData>& plates,
                        const std::string& outputName)
{
  Long64_t maximumCount = 0;
  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    maximumCount = std::max(
      maximumCount,
      static_cast<Long64_t>(plates[plate].values[kX].size()));
  }

  TCanvas* canvas = new TCanvas(
    "plusStackElectronCountCanvas",
    "Detected electrons in the +z MCP stack",
    1000,
    750);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.14);
  canvas->SetGridy();

  TH1D* frame = new TH1D(
    "plusStackElectronCountFrame",
    "Detected electrons by MCP plate in the +z stack",
    plates.size(),
    0.0,
    static_cast<Double_t>(plates.size()));
  frame->SetDirectory(nullptr);
  frame->SetMinimum(0.0);
  frame->SetMaximum(maximumCount > 0 ? 1.20*maximumCount : 1.0);
  frame->GetYaxis()->SetTitle("Number of detected electrons");
  frame->GetYaxis()->SetTitleOffset(1.45);

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    frame->GetXaxis()->SetBinLabel(
      plate + 1, plates[plate].label.c_str());
  }
  frame->Draw("HIST");

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::string histogramName =
      "plusStackElectronCountMcp_" +
      std::to_string(plates[plate].mcpIndex);
    TH1D* histogram = new TH1D(
      histogramName.c_str(),
      "",
      plates.size(),
      0.0,
      static_cast<Double_t>(plates.size()));
    histogram->SetDirectory(nullptr);
    histogram->SetBinContent(
      plate + 1, plates[plate].values[kX].size());
    histogram->SetFillColor(plates[plate].color);
    histogram->SetLineColor(plates[plate].color);
    histogram->SetBarWidth(0.75);
    histogram->SetBarOffset(0.125);
    histogram->Draw("BAR SAME");
  }

  TH1D* labels = static_cast<TH1D*>(
    frame->Clone("plusStackElectronCountLabels"));
  labels->Reset();
  labels->SetDirectory(nullptr);
  labels->SetMarkerSize(1.2);
  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    labels->SetBinContent(
      plate + 1, plates[plate].values[kX].size());
  }
  gStyle->SetPaintTextFormat(".0f");
  labels->Draw("TEXT0 SAME");

  canvas->RedrawAxis();
  canvas->SaveAs(outputName.c_str());
}

void PrintStatistics(const std::vector<PlateData>& plates,
                     const std::string& directionZName)
{
  std::cout << "\nStatistics for side = +1 only" << std::endl;
  std::cout << std::fixed << std::setprecision(5);
  std::cout
    << std::setw(5) << "MCP"
    << std::setw(12) << "N electrons"
    << std::setw(12) << "mean x"
    << std::setw(12) << "RMS x"
    << std::setw(12) << "mean y"
    << std::setw(12) << "RMS y"
    << std::setw(12) << "mean z"
    << std::setw(12) << "RMS z"
    << std::setw(14) << ("mean " + directionZName)
    << std::setw(14) << ("RMS " + directionZName)
    << std::endl;

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    std::cout
      << std::setw(5) << plates[plate].mcpIndex
      << std::setw(12) << plates[plate].values[kX].size()
      << std::setw(12) << ComputeMean(plates[plate].values[kX])
      << std::setw(12) << ComputeRms(plates[plate].values[kX])
      << std::setw(12) << ComputeMean(plates[plate].values[kY])
      << std::setw(12) << ComputeRms(plates[plate].values[kY])
      << std::setw(12) << ComputeMean(plates[plate].values[kZ])
      << std::setw(12) << ComputeRms(plates[plate].values[kZ])
      << std::setw(14)
      << ComputeMean(plates[plate].values[kDirectionZ])
      << std::setw(14)
      << ComputeRms(plates[plate].values[kDirectionZ])
      << std::endl;
  }
}
}  // namespace

void plot_electron_3d_by_mcp(
  const char* fileName = "build/mcp_output.root")
{
  gSystem->mkdir("Fig/current", kTRUE);

  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* tree = nullptr;
  file->GetObject("ElectronChannelHitTree", tree);
  if (!tree) {
    std::cerr << "ElectronChannelHitTree was not found in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::string branchNames[kNumberOfVariables];
  std::string sideName;
  std::string mcpIndexName;

  TLeaf* leaves[kNumberOfVariables] = {
    FindLeaf(tree, {"x_mm", "x", "posX", "positionX"},
             branchNames[kX]),
    FindLeaf(tree, {"y_mm", "y", "posY", "positionY"},
             branchNames[kY]),
    FindLeaf(tree, {"z_mm", "z", "posZ", "positionZ"},
             branchNames[kZ]),
    FindLeaf(tree, {"dirX", "px", "momentumX", "momX"},
             branchNames[kDirectionX]),
    FindLeaf(tree, {"dirY", "py", "momentumY", "momY"},
             branchNames[kDirectionY]),
    FindLeaf(tree, {"dirZ", "pz", "momentumZ", "momZ"},
             branchNames[kDirectionZ])
  };

  TLeaf* sideLeaf = FindLeaf(tree, {"side"}, sideName);
  TLeaf* mcpIndexLeaf = FindLeaf(
    tree, {"mcpIndex", "McpIndex", "plateIndex"}, mcpIndexName);

  Bool_t branchesAreValid = sideLeaf && mcpIndexLeaf;
  for (Int_t variable = 0; variable < kNumberOfVariables; ++variable) {
    branchesAreValid = branchesAreValid && leaves[variable];
  }

  if (!branchesAreValid) {
    std::cerr << "Required branches could not be detected." << std::endl;
    std::cerr << "  position: " << branchNames[kX] << ", "
              << branchNames[kY] << ", " << branchNames[kZ] << std::endl;
    std::cerr << "  direction/momentum: "
              << branchNames[kDirectionX] << ", "
              << branchNames[kDirectionY] << ", "
              << branchNames[kDirectionZ] << std::endl;
    std::cerr << "  plate: " << sideName << ", " << mcpIndexName
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << "ElectronChannelHitTree branch detection:" << std::endl;
  std::cout << "  position: " << branchNames[kX] << ", "
            << branchNames[kY] << ", " << branchNames[kZ] << std::endl;
  std::cout << "  direction/momentum: "
            << branchNames[kDirectionX] << ", "
            << branchNames[kDirectionY] << ", "
            << branchNames[kDirectionZ] << std::endl;
  std::cout << "  plate: " << sideName << ", " << mcpIndexName
            << std::endl;

  const Color_t colors[kNumberOfMCPs] = {
    kRed + 1, kOrange + 7, kBlue + 1
  };

  std::vector<PlateData> plates;
  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    PlateData plate;
    plate.mcpIndex = mcpIndex;
    plate.color = colors[mcpIndex];
    plate.label = "+z MCP " + std::to_string(mcpIndex);
    plates.push_back(plate);
  }

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);

    if (static_cast<Int_t>(sideLeaf->GetValue()) != kPlusSide) {
      continue;
    }

    PlateData* plate = FindPlate(
      plates, static_cast<Int_t>(mcpIndexLeaf->GetValue()));
    if (!plate) {
      continue;
    }

    for (Int_t variable = 0; variable < kNumberOfVariables; ++variable) {
      plate->values[variable].push_back(leaves[variable]->GetValue());
    }
  }

  PrintStatistics(plates, branchNames[kDirectionZ]);

  gStyle->SetOptStat(0);

  DrawNormalizedHistograms(
    plates,
    kZ,
    "plusStackZCanvas",
    "Electron Z position in the +z MCP stack",
    branchNames[kZ],
    "Fig/current/plus_stack_z_position_by_mcp.png");

  DrawNormalizedHistograms(
    plates,
    kDirectionX,
    "plusStackDirectionXCanvas",
    "Electron X direction or momentum in the +z MCP stack",
    branchNames[kDirectionX],
    "Fig/current/plus_stack_dirx_by_mcp.png");

  DrawNormalizedHistograms(
    plates,
    kDirectionY,
    "plusStackDirectionYCanvas",
    "Electron Y direction or momentum in the +z MCP stack",
    branchNames[kDirectionY],
    "Fig/current/plus_stack_diry_by_mcp.png");

  DrawNormalizedHistograms(
    plates,
    kDirectionZ,
    "plusStackDirectionZCanvas",
    "Electron Z direction or momentum in the +z MCP stack",
    branchNames[kDirectionZ],
    "Fig/current/plus_stack_dirz_by_mcp.png");

  DrawNormalizedHistograms(
    plates,
    kX,
    "plusStackXCanvas",
    "Electron X position in the +z MCP stack",
    branchNames[kX],
    "Fig/current/plus_stack_x_position_by_mcp.png");

  DrawNormalizedHistograms(
    plates,
    kY,
    "plusStackYCanvas",
    "Electron Y position in the +z MCP stack",
    branchNames[kY],
    "Fig/current/plus_stack_y_position_by_mcp.png");

  DrawNormalizedHistograms(
    plates,
    kX,
    "plusStackXZoomCanvas",
    "Electron X position in the +z MCP stack (1%-99%)",
    branchNames[kX],
    "Fig/current/plus_stack_x_position_zoom.png",
    true,
    kUseLogScaleForZoom);

  DrawNormalizedHistograms(
    plates,
    kY,
    "plusStackYZoomCanvas",
    "Electron Y position in the +z MCP stack (1%-99%)",
    branchNames[kY],
    "Fig/current/plus_stack_y_position_zoom.png",
    true,
    kUseLogScaleForZoom);

  DrawElectronCounts(
    plates, "Fig/current/plus_stack_electron_count_by_mcp.png");

  std::cout << "\nSaved:" << std::endl;
  std::cout << "  plus_stack_z_position_by_mcp.png" << std::endl;
  std::cout << "  plus_stack_dirx_by_mcp.png" << std::endl;
  std::cout << "  plus_stack_diry_by_mcp.png" << std::endl;
  std::cout << "  plus_stack_dirz_by_mcp.png" << std::endl;
  std::cout << "  plus_stack_x_position_by_mcp.png" << std::endl;
  std::cout << "  plus_stack_y_position_by_mcp.png" << std::endl;
  std::cout << "  plus_stack_x_position_zoom.png" << std::endl;
  std::cout << "  plus_stack_y_position_zoom.png" << std::endl;
  std::cout << "  plus_stack_electron_count_by_mcp.png" << std::endl;

  file->Close();
  delete file;
}
