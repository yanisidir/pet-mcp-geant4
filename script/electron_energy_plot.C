#include "TCanvas.h"
#include "TSystem.h"
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
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  const Int_t kNumberOfBins = 100;
  const Double_t kMinimumEnergyKeV = 0.0;
  const Double_t kDefaultMaximumEnergyKeV = 511.0;

  struct PlateKey
  {
    Int_t side;
    Int_t mcpIndex;

    bool operator<(const PlateKey& other) const
    {
      if (side != other.side) {
        return side < other.side;
      }
      return mcpIndex < other.mcpIndex;
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
    Double_t sum = 0.0;

    for (std::size_t i = 0; i < values.size(); ++i) {
      const Double_t diff = values[i] - mean;
      sum += diff*diff;
    }

    return std::sqrt(sum/values.size());
  }

  Double_t MinimumValue(const std::vector<Double_t>& values)
  {
    return values.empty()
      ? 0.0
      : *std::min_element(values.begin(), values.end());
  }

  Double_t MaximumValue(const std::vector<Double_t>& values)
  {
    return values.empty()
      ? 0.0
      : *std::max_element(values.begin(), values.end());
  }

  std::string PlateLabel(const PlateKey& key)
  {
    std::ostringstream label;
    label << (key.side > 0 ? "+z" : "-z")
          << " MCP " << key.mcpIndex;
    return label.str();
  }

  Color_t PlateColor(Int_t side, Int_t mcpIndex)
  {
    const Color_t plusColors[5] = {
      kRed + 1, kOrange + 7, kBlue + 1, kGreen + 2, kMagenta + 1
    };
    const Color_t minusColors[5] = {
      kAzure + 2, kCyan + 2, kViolet + 1, kPink + 6, kGray + 2
    };

    if (side > 0) {
      return plusColors[mcpIndex % 5];
    }
    return minusColors[mcpIndex % 5];
  }

  void PrintEnergyStats(const std::string& label,
                        const std::vector<Double_t>& energies)
  {
    std::cout
      << std::setw(12) << label
      << " | " << std::setw(11) << energies.size()
      << " | " << std::setw(14) << ComputeMean(energies)
      << " | " << std::setw(13) << ComputeRms(energies)
      << " | " << std::setw(13) << MinimumValue(energies)
      << " | " << std::setw(13) << MaximumValue(energies)
      << std::endl;
  }
}

void electron_energy_plot(
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
    tree,
    {"mcpIndex", "McpIndex", "plateIndex"},
    mcpIndexName);

  if (!energyLeaf) {
    std::cerr << "Energy branch could not be detected." << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::cout << "Energy branch: " << energyName << std::endl;
  if (sideLeaf && mcpIndexLeaf) {
    std::cout << "Plate branches: " << sideName
              << ", " << mcpIndexName << std::endl;
  } else {
    std::cout << "Plate branches not found; only global plot will be made."
              << std::endl;
  }

  std::vector<Double_t> energiesKeV;
  std::map<PlateKey, std::vector<Double_t> > energiesByPlate;
  
  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);

    const Double_t energy = energyLeaf->GetValue();
    energiesKeV.push_back(energy);

    if (sideLeaf && mcpIndexLeaf) {
      PlateKey key;
      key.side = static_cast<Int_t>(sideLeaf->GetValue());
      key.mcpIndex = static_cast<Int_t>(mcpIndexLeaf->GetValue());
      energiesByPlate[key].push_back(energy);
    }
  }

  const Double_t minimum = energiesKeV.empty()
    ? 0.0
    : *std::min_element(energiesKeV.begin(), energiesKeV.end());

  const Double_t maximum = energiesKeV.empty()
    ? 0.0
    : *std::max_element(energiesKeV.begin(), energiesKeV.end());
  const Double_t plotMaximumEnergyKeV =
    maximum > kDefaultMaximumEnergyKeV
      ? 1.05*maximum
      : kDefaultMaximumEnergyKeV;

  std::cout << "\nElectron kinetic energy statistics" << std::endl;
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "N electrons   = " << energiesKeV.size() << std::endl;
  std::cout << "Mean E [keV]  = " << ComputeMean(energiesKeV) << std::endl;
  std::cout << "RMS E [keV]   = " << ComputeRms(energiesKeV) << std::endl;
  std::cout << "Min E [keV]   = " << minimum << std::endl;
  std::cout << "Max E [keV]   = " << maximum << std::endl;

  if (!energiesByPlate.empty()) {
    std::cout << "\nPer-plate electron kinetic energy statistics"
              << std::endl;
    std::cout
      << std::setw(12) << "plate"
      << " | " << std::setw(11) << "N electrons"
      << " | " << std::setw(14) << "mean E [keV]"
      << " | " << std::setw(13) << "RMS E [keV]"
      << " | " << std::setw(13) << "min E [keV]"
      << " | " << std::setw(13) << "max E [keV]"
      << std::endl;

    for (std::map<PlateKey, std::vector<Double_t> >::const_iterator it =
           energiesByPlate.begin();
         it != energiesByPlate.end();
         ++it) {
      PrintEnergyStats(PlateLabel(it->first), it->second);
    }
  }

  gStyle->SetOptStat(1110);

  TCanvas* canvas = new TCanvas(
    "electronEnergyCanvas",
    "Electron kinetic energy",
    1000,
    800);

  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);
  canvas->SetGridy();

  TH1D* histogram = new TH1D(
    "electronEnergyAll",
    "Electron kinetic energy;Electron kinetic energy [keV];Entries",
    kNumberOfBins,
    kMinimumEnergyKeV,
    plotMaximumEnergyKeV);

  for (std::size_t i = 0; i < energiesKeV.size(); ++i) {
    histogram->Fill(energiesKeV[i]);
  }

  if (histogram->Integral() > 0.0) {
    histogram->Scale(1.0/histogram->Integral());
  }

  histogram->SetLineColor(kBlue + 1);
  histogram->SetLineWidth(3);
  histogram->Draw("HIST");

  canvas->SaveAs("Fig/current/electron_energy_all_mcp.png");

  if (!energiesByPlate.empty()) {
    gStyle->SetOptStat(0);

    TCanvas* plateCanvas = new TCanvas(
      "electronEnergyByPlateCanvas",
      "Electron kinetic energy by MCP plate",
      1100,
      800);
    plateCanvas->SetLeftMargin(0.12);
    plateCanvas->SetRightMargin(0.05);
    plateCanvas->SetBottomMargin(0.12);
    plateCanvas->SetGridy();

    std::vector<TH1D*> plateHistograms;
    Double_t largestBinContent = 0.0;
    for (std::map<PlateKey, std::vector<Double_t> >::const_iterator it =
           energiesByPlate.begin();
         it != energiesByPlate.end();
         ++it) {
      const PlateKey key = it->first;
      std::ostringstream histName;
      histName << "hElectronEnergy_side" << key.side
               << "_mcp" << key.mcpIndex;

      TH1D* plateHist = new TH1D(
        histName.str().c_str(),
        "",
        kNumberOfBins,
        kMinimumEnergyKeV,
        plotMaximumEnergyKeV);
      plateHist->SetDirectory(nullptr);
      plateHist->SetLineColor(PlateColor(key.side, key.mcpIndex));
      plateHist->SetLineWidth(2);

      for (std::size_t i = 0; i < it->second.size(); ++i) {
        plateHist->Fill(it->second[i]);
      }
      if (plateHist->Integral() > 0.0) {
        plateHist->Scale(1.0/plateHist->Integral());
        largestBinContent =
          std::max(largestBinContent, plateHist->GetMaximum());
      }
      plateHistograms.push_back(plateHist);
    }

    TH1D* frame = new TH1D(
      "electronEnergyByPlateFrame",
      "Electron kinetic energy by MCP plate;Electron kinetic energy [keV];Normalized entries",
      kNumberOfBins,
      kMinimumEnergyKeV,
      plotMaximumEnergyKeV);
    frame->SetDirectory(nullptr);
    frame->SetMinimum(0.0);
    frame->SetMaximum(largestBinContent > 0.0
                        ? 1.25*largestBinContent
                        : 1.0);
    frame->Draw("HIST");

    TLegend* legend = new TLegend(0.68, 0.58, 0.92, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.030);

    std::size_t histIndex = 0;
    for (std::map<PlateKey, std::vector<Double_t> >::const_iterator it =
           energiesByPlate.begin();
         it != energiesByPlate.end();
         ++it, ++histIndex) {
      plateHistograms[histIndex]->Draw("HIST SAME");
      std::ostringstream label;
      label << PlateLabel(it->first) << " (N=" << it->second.size() << ")";
      legend->AddEntry(plateHistograms[histIndex], label.str().c_str(), "l");
    }

    legend->Draw();
    plateCanvas->RedrawAxis();
    plateCanvas->SaveAs("Fig/current/electron_energy_by_plate.png");

    TCanvas* countCanvas = new TCanvas(
      "electronEnergyCountByPlateCanvas",
      "Detected electrons by MCP plate",
      1000,
      700);
    countCanvas->SetLeftMargin(0.12);
    countCanvas->SetBottomMargin(0.16);
    countCanvas->SetGridy();

    TH1D* countHist = new TH1D(
      "electronCountByPlate",
      "Detected channel electrons by MCP plate;MCP plate;Electrons",
      static_cast<Int_t>(energiesByPlate.size()),
      0.0,
      static_cast<Double_t>(energiesByPlate.size()));
    countHist->SetDirectory(nullptr);
    countHist->SetFillColor(kAzure - 9);
    countHist->SetLineColor(kAzure + 2);
    countHist->SetLineWidth(2);

    Int_t bin = 1;
    for (std::map<PlateKey, std::vector<Double_t> >::const_iterator it =
           energiesByPlate.begin();
         it != energiesByPlate.end();
         ++it, ++bin) {
      countHist->SetBinContent(bin, it->second.size());
      countHist->GetXaxis()->SetBinLabel(
        bin, PlateLabel(it->first).c_str());
    }
    countHist->Draw("BAR");
    countCanvas->SaveAs("Fig/current/electron_count_by_plate_energy_macro.png");
  }

  std::cout << "\nSaved:" << std::endl;
  std::cout << "  electron_energy_all_mcp.png" << std::endl;
  if (!energiesByPlate.empty()) {
    std::cout << "  electron_energy_by_plate.png" << std::endl;
    std::cout << "  electron_count_by_plate_energy_macro.png"
              << std::endl;
  }

  file->Close();
  delete file;
}
