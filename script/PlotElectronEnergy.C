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
const Double_t kMinimumEnergyKeV = 0.0;
const Double_t kMaximumEnergyKeV = 5008.0;

struct PlateEnergy
{
  Int_t mcpIndex;
  Color_t color;
  std::string label;
  std::vector<Double_t> energiesKeV;
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

void PrintStatistics(const std::vector<PlateEnergy>& plates)
{
  std::cout << "\nElectron kinetic energy statistics, side = +1"
            << std::endl;
  std::cout << std::fixed << std::setprecision(4);
  std::cout
    << std::setw(10) << "MCP index"
    << " | " << std::setw(11) << "N electrons"
    << " | " << std::setw(14) << "mean E [keV]"
    << " | " << std::setw(13) << "RMS E [keV]"
    << " | " << std::setw(13) << "min E [keV]"
    << " | " << std::setw(13) << "max E [keV]"
    << std::endl;

  for (std::size_t plate = 0; plate < plates.size(); ++plate) {
    const std::vector<Double_t>& energies = plates[plate].energiesKeV;
    const Double_t minimum = energies.empty()
      ? 0.0
      : *std::min_element(energies.begin(), energies.end());
    const Double_t maximum = energies.empty()
      ? 0.0
      : *std::max_element(energies.begin(), energies.end());

    std::cout
      << std::setw(10) << plates[plate].mcpIndex
      << " | " << std::setw(11) << energies.size()
      << " | " << std::setw(14) << ComputeMean(energies)
      << " | " << std::setw(13) << ComputeRms(energies)
      << " | " << std::setw(13) << minimum
      << " | " << std::setw(13) << maximum
      << std::endl;
  }
}
}  // namespace

void PlotElectronEnergy(
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

  std::string energyName;
  std::string sideName;
  std::string mcpIndexName;

  TLeaf* energyLeaf = FindLeaf(
    tree,
    {"kineticEnergy_keV", "energy_keV", "kineticEnergy", "energy"},
    energyName);
  TLeaf* sideLeaf = FindLeaf(tree, {"side"}, sideName);
  TLeaf* mcpIndexLeaf = FindLeaf(
    tree, {"mcpIndex", "McpIndex", "plateIndex"}, mcpIndexName);

  if (!energyLeaf || !sideLeaf || !mcpIndexLeaf) {
    std::cerr << "Required branches could not be detected." << std::endl;
    std::cerr << "  energy: " << energyName << std::endl;
    std::cerr << "  plate: " << sideName << ", " << mcpIndexName
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << "Energy branch: " << energyName << std::endl;
  std::cout << "Filter: " << sideName << " == 1" << std::endl;

  const Color_t colors[kNumberOfMCPs] = {
    kRed + 1, kOrange + 7, kBlue + 1
  };

  // const Color_t colors[kNumberOfMCPs] = {
  //   kRed + 1, kOrange + 7, kBlue + 1, kGreen + 1, kMagenta + 1
  // };

  std::vector<PlateEnergy> plates(kNumberOfMCPs);
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

    plates[mcpIndex].energiesKeV.push_back(energyLeaf->GetValue());
  }

  PrintStatistics(plates);

  std::vector<TH1D*> histograms;
  Double_t largestBinContent = 0.0;

  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    const std::string histogramName =
      "electronEnergyMcp_" + std::to_string(mcpIndex);
    TH1D* histogram = new TH1D(
      histogramName.c_str(),
      "",
      kNumberOfBins,
      kMinimumEnergyKeV,
      kMaximumEnergyKeV);
    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(plates[mcpIndex].color);
    histogram->SetLineWidth(3);

    for (std::size_t value = 0;
         value < plates[mcpIndex].energiesKeV.size();
         ++value) {
      histogram->Fill(plates[mcpIndex].energiesKeV[value]);
    }

    if (histogram->Integral() > 0.0) {
      histogram->Scale(1.0/histogram->Integral());
      largestBinContent =
        std::max(largestBinContent, histogram->GetMaximum());
    }
    histograms.push_back(histogram);
  }

  gStyle->SetOptStat(0);

  TCanvas* canvas = new TCanvas(
    "plusStackElectronEnergyCanvas",
    "Kinetic energy of channel electrons in the +z MCP stack",
    1100,
    800);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);
  canvas->SetGridy();

  TH1D* frame = new TH1D(
    "plusStackElectronEnergyFrame",
    "Kinetic energy of channel electrons in the +z MCP stack",
    kNumberOfBins,
    kMinimumEnergyKeV,
    kMaximumEnergyKeV);
  frame->SetDirectory(nullptr);
  frame->SetMinimum(0.0);
  frame->SetMaximum(
    largestBinContent > 0.0 ? 1.20*largestBinContent : 1.0);
  frame->GetXaxis()->SetTitle("Electron kinetic energy [keV]");
  frame->GetYaxis()->SetTitle("Normalized entries");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->Draw("HIST");

  TLegend* legend = new TLegend(0.75, 0.65, 0.92, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.032);

  for (Int_t mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    if (histograms[mcpIndex]->GetEntries() > 0.0) {
      histograms[mcpIndex]->Draw("HIST SAME");
    }

    const std::string legendLabel =
      plates[mcpIndex].label + " (N=" +
      std::to_string(plates[mcpIndex].energiesKeV.size()) + ")";
    legend->AddEntry(
      histograms[mcpIndex], legendLabel.c_str(), "l");
  }

  legend->Draw();
  canvas->RedrawAxis();
  canvas->SaveAs("Fig/current/plus_stack_electron_energy_by_mcp.png");

  std::cout << "\nSaved: plus_stack_electron_energy_by_mcp.png"
            << std::endl;

  file->Close();
  delete file;
}
