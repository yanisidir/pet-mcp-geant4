#include "TAxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
const Int_t kPlusSide = 1;
const Int_t kNumberOfMCPs = 3;
const Int_t kNumberOfBins = 100;
const Int_t kNumberOfCoordinates = 3;

enum CoordinateIndex
{
  kX = 0,
  kY,
  kZ
};

struct PlatePositions
{
  Int_t mcpIndex;
  Color_t color;
  std::string label;
  std::vector<Double_t> values[kNumberOfCoordinates];
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

void ComputeRange(const std::vector<PlatePositions>& plates,
                  Int_t coordinate,
                  Double_t& minimum,
                  Double_t& maximum)
{
  minimum = std::numeric_limits<Double_t>::max();
  maximum = -std::numeric_limits<Double_t>::max();

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::vector<Double_t>& values =
      plates[plate].values[coordinate];
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

void DrawPosition(const std::vector<PlatePositions>& plates,
                  Int_t coordinate,
                  const std::string& canvasName,
                  const std::string& title,
                  const std::string& axisTitle,
                  const std::string& outputName)
{
  Double_t minimum = 0.0;
  Double_t maximum = 1.0;
  ComputeRange(plates, coordinate, minimum, maximum);

  std::vector<TH1D*> histograms;
  Double_t largestBinContent = 0.0;

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::string name =
      canvasName + "_mcp_" + std::to_string(plates[plate].mcpIndex);
    TH1D* histogram = new TH1D(
      name.c_str(), "", kNumberOfBins, minimum, maximum);
    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(plates[plate].color);
    histogram->SetLineWidth(3);

    const std::vector<Double_t>& values =
      plates[plate].values[coordinate];
    for (std::size_t value = 0; value < values.size(); ++value) {
      histogram->Fill(values[value]);
    }

    if (histogram->Integral() > 0.0) {
      histogram->Scale(1.0/histogram->Integral());
      largestBinContent =
        std::max(largestBinContent, histogram->GetMaximum());
    }
    histograms.push_back(histogram);
  }

  TCanvas* canvas =
    new TCanvas(canvasName.c_str(), title.c_str(), 1100, 800);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);
  canvas->SetGridy();

  const std::string frameName = canvasName + "_frame";
  TH1D* frame = new TH1D(
    frameName.c_str(), title.c_str(), kNumberOfBins, minimum, maximum);
  frame->SetDirectory(nullptr);
  frame->SetMinimum(0.0);
  frame->SetMaximum(
    largestBinContent > 0.0 ? 1.20*largestBinContent : 1.0);
  frame->GetXaxis()->SetTitle(axisTitle.c_str());
  frame->GetYaxis()->SetTitle("Normalized entries");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->Draw("HIST");

  TLegend* legend = new TLegend(0.70, 0.68, 0.94, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.035);

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    if (histograms[plate]->GetEntries() > 0.0) {
      histograms[plate]->Draw("HIST SAME");
    }

    const std::string legendLabel =
      plates[plate].label + " (N=" +
      std::to_string(plates[plate].values[coordinate].size()) + ")";
    legend->AddEntry(histograms[plate], legendLabel.c_str(), "l");
  }

  legend->Draw();
  canvas->RedrawAxis();
  canvas->SaveAs(outputName.c_str());
}
}  // namespace

void PlotElectronPositions(
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

  std::string coordinateNames[kNumberOfCoordinates];
  std::string sideName;
  std::string mcpIndexName;

  TLeaf* coordinateLeaves[kNumberOfCoordinates] = {
    FindLeaf(tree, {"x_mm", "x_cm", "x", "posX", "positionX"},
             coordinateNames[kX]),
    FindLeaf(tree, {"y_mm", "y_cm", "y", "posY", "positionY"},
             coordinateNames[kY]),
    FindLeaf(tree, {"z_mm", "z_cm", "z", "posZ", "positionZ"},
             coordinateNames[kZ])
  };
  TLeaf* sideLeaf = FindLeaf(tree, {"side"}, sideName);
  TLeaf* mcpIndexLeaf = FindLeaf(
    tree, {"mcpIndex", "McpIndex", "plateIndex"}, mcpIndexName);

  Bool_t branchesAreValid = sideLeaf && mcpIndexLeaf;
  for (Int_t coordinate = 0;
       coordinate < kNumberOfCoordinates;
       ++coordinate) {
    branchesAreValid =
      branchesAreValid && coordinateLeaves[coordinate];
  }

  if (!branchesAreValid) {
    std::cerr << "Required position branches could not be detected."
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << "Position branches: "
            << coordinateNames[kX] << ", "
            << coordinateNames[kY] << ", "
            << coordinateNames[kZ] << std::endl;
  std::cout << "Filter: " << sideName << " == 1" << std::endl;

  const Color_t colors[kNumberOfMCPs] = {
    kRed + 1, kOrange + 7, kBlue + 1
  };

  std::vector<PlatePositions> plates(kNumberOfMCPs);
  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    plates[mcpIndex].mcpIndex = mcpIndex;
    plates[mcpIndex].color = colors[mcpIndex];
    plates[mcpIndex].label = "+z MCP " + std::to_string(mcpIndex);
  }

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    if (static_cast<Int_t>(sideLeaf->GetValue()) != kPlusSide) {
      continue;
    }

    const Int_t mcpIndex =
      static_cast<Int_t>(mcpIndexLeaf->GetValue());
    if (mcpIndex < 0 || mcpIndex >= kNumberOfMCPs) {
      continue;
    }

    for (Int_t coordinate = 0;
         coordinate < kNumberOfCoordinates;
         ++coordinate) {
      plates[mcpIndex].values[coordinate].push_back(
        coordinateLeaves[coordinate]->GetValue());
    }
  }

  gStyle->SetOptStat(0);

  DrawPosition(
    plates, kX, "plusStackXCanvas",
    "Electron X position in the +z MCP stack",
    coordinateNames[kX], "Fig/current/plus_stack_x_position_by_mcp.png");
  DrawPosition(
    plates, kY, "plusStackYCanvas",
    "Electron Y position in the +z MCP stack",
    coordinateNames[kY], "Fig/current/plus_stack_y_position_by_mcp.png");
  DrawPosition(
    plates, kZ, "plusStackZCanvas",
    "Electron Z position in the +z MCP stack",
    coordinateNames[kZ], "Fig/current/plus_stack_z_position_by_mcp.png");

  std::cout << "Saved position plots for the +z MCP stack."
            << std::endl;

  file->Close();
  delete file;
}
