#include "TAxis.h"
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
const Int_t kNumberOfComponents = 3;

enum ComponentIndex
{
  kDirectionX = 0,
  kDirectionY,
  kDirectionZ
};

struct PlateMomenta
{
  Int_t mcpIndex;
  Color_t color;
  std::string label;
  std::vector<Double_t> values[kNumberOfComponents];
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

void ComputeRange(const std::vector<PlateMomenta>& plates,
                  Int_t component,
                  Double_t& minimum,
                  Double_t& maximum)
{
  minimum = std::numeric_limits<Double_t>::max();
  maximum = -std::numeric_limits<Double_t>::max();

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::vector<Double_t>& values =
      plates[plate].values[component];
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

void DrawMomentum(const std::vector<PlateMomenta>& plates,
                  Int_t component,
                  const std::string& canvasName,
                  const std::string& title,
                  const std::string& axisTitle,
                  const std::string& outputName)
{
  Double_t minimum = 0.0;
  Double_t maximum = 1.0;
  ComputeRange(plates, component, minimum, maximum);

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
      plates[plate].values[component];
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

  TLegend* legend = new TLegend(0.15, 0.65, 0.40, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.032);

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    if (histograms[plate]->GetEntries() > 0.0) {
      histograms[plate]->Draw("HIST SAME");
    }
    legend->AddEntry(
      histograms[plate], plates[plate].label.c_str(), "l");
  }

  legend->Draw();
  canvas->RedrawAxis();
  canvas->SaveAs(outputName.c_str());
}
}  // namespace

void PlotElectronMomenta(
  const char* fileName = "build/mcp_output_t0.root")
{
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

  std::string componentNames[kNumberOfComponents];
  std::string sideName;
  std::string mcpIndexName;

  TLeaf* componentLeaves[kNumberOfComponents] = {
    FindLeaf(tree, {"dirX", "px", "momentumX", "momX"},
             componentNames[kDirectionX]),
    FindLeaf(tree, {"dirY", "py", "momentumY", "momY"},
             componentNames[kDirectionY]),
    FindLeaf(tree, {"dirZ", "pz", "momentumZ", "momZ"},
             componentNames[kDirectionZ])
  };
  TLeaf* sideLeaf = FindLeaf(tree, {"side"}, sideName);
  TLeaf* mcpIndexLeaf = FindLeaf(
    tree, {"mcpIndex", "McpIndex", "plateIndex"}, mcpIndexName);

  Bool_t branchesAreValid = sideLeaf && mcpIndexLeaf;
  for (Int_t component = 0;
       component < kNumberOfComponents;
       ++component) {
    branchesAreValid = branchesAreValid && componentLeaves[component];
  }

  if (!branchesAreValid) {
    std::cerr << "Required direction/momentum branches could not be "
              << "detected." << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << "Direction/momentum branches: "
            << componentNames[kDirectionX] << ", "
            << componentNames[kDirectionY] << ", "
            << componentNames[kDirectionZ] << std::endl;
  std::cout << "Filter: " << sideName << " == 1" << std::endl;

  const Color_t colors[kNumberOfMCPs] = {
    kRed + 1, kOrange + 7, kBlue + 1
  };

  std::vector<PlateMomenta> plates(kNumberOfMCPs);
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

    for (Int_t component = 0;
         component < kNumberOfComponents;
         ++component) {
      plates[mcpIndex].values[component].push_back(
        componentLeaves[component]->GetValue());
    }
  }

  gStyle->SetOptStat(0);

  DrawMomentum(
    plates, kDirectionX, "plusStackDirectionXCanvas",
    "Electron X direction or momentum in the +z MCP stack",
    componentNames[kDirectionX], "plus_stack_dirx_by_mcp.png");
  DrawMomentum(
    plates, kDirectionY, "plusStackDirectionYCanvas",
    "Electron Y direction or momentum in the +z MCP stack",
    componentNames[kDirectionY], "plus_stack_diry_by_mcp.png");
  DrawMomentum(
    plates, kDirectionZ, "plusStackDirectionZCanvas",
    "Electron Z direction or momentum in the +z MCP stack",
    componentNames[kDirectionZ], "plus_stack_dirz_by_mcp.png");

  std::cout << "Saved direction/momentum plots for the +z MCP stack."
            << std::endl;

  file->Close();
  delete file;
}
