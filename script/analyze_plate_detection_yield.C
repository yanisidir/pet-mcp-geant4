#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{
  typedef std::pair<Int_t, Int_t> PlateKey;

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

void analyze_plate_detection_yield(
  const char* fileName = "build/mcp_output.root")
{
  gSystem->mkdir("Fig/current", kTRUE);

  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* gammaEntryTree = nullptr;
  TTree* plateTree = nullptr;
  file->GetObject("GammaMcpEntryTree", gammaEntryTree);
  file->GetObject("McpPlateStatsTree", plateTree);

  if (!gammaEntryTree || !plateTree) {
    std::cerr << "Missing GammaMcpEntryTree or McpPlateStatsTree in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t gammaSide = 0;
  Int_t gammaMcpIndex = -1;
  gammaEntryTree->SetBranchAddress("side", &gammaSide);
  gammaEntryTree->SetBranchAddress("mcpIndex", &gammaMcpIndex);

  std::map<PlateKey, Long64_t> gammaEntriesByPlate;
  std::set<PlateKey> allPlates;

  const Long64_t gammaEntryCount = gammaEntryTree->GetEntries();
  for (Long64_t entry = 0; entry < gammaEntryCount; ++entry) {
    gammaEntryTree->GetEntry(entry);

    const PlateKey plate(gammaSide, gammaMcpIndex);
    ++gammaEntriesByPlate[plate];
    allPlates.insert(plate);
  }

  Int_t plateSide = 0;
  Int_t plateMcpIndex = -1;
  Int_t electronChannelCount = 0;
  plateTree->SetBranchAddress("side", &plateSide);
  plateTree->SetBranchAddress("mcpIndex", &plateMcpIndex);
  plateTree->SetBranchAddress("electronChannelCount",
                              &electronChannelCount);

  std::map<PlateKey, Long64_t> eventsWithChannelElectronByPlate;
  std::map<PlateKey, Long64_t> channelElectronsByPlate;

  const Long64_t plateEntryCount = plateTree->GetEntries();
  for (Long64_t entry = 0; entry < plateEntryCount; ++entry) {
    plateTree->GetEntry(entry);

    const PlateKey plate(plateSide, plateMcpIndex);
    allPlates.insert(plate);
    channelElectronsByPlate[plate] += electronChannelCount;
    if (electronChannelCount > 0) {
      ++eventsWithChannelElectronByPlate[plate];
    }
  }

  const std::vector<PlateKey> plates(allPlates.begin(),
                                     allPlates.end());

  TH1D* eventYieldHistogram =
    new TH1D("hPlateDetectionEventYield",
             "Plate detection event yield;"
             "MCP plate;Event yield (%)",
             static_cast<Int_t>(plates.size()),
             0.0,
             static_cast<Double_t>(plates.size()));

  TH1D* electronYieldHistogram =
    new TH1D("hPlateElectronYield",
             "Plate electron yield;"
             "MCP plate;Electron yield (%)",
             static_cast<Int_t>(plates.size()),
             0.0,
             static_cast<Double_t>(plates.size()));

  std::cout << "File: " << fileName << std::endl;
  std::cout << "\n=== MCP plate detection yields ===" << std::endl;
  std::cout << std::left
            << std::setw(8) << "side"
            << std::setw(10) << "mcpIndex"
            << std::right
            << std::setw(16) << "gamma_entries"
            << std::setw(22) << "events_with_electron"
            << std::setw(20) << "channel_electrons"
            << std::setw(16) << "event_yield"
            << std::setw(16) << "event_yield_%"
            << std::setw(18) << "electron_yield"
            << std::setw(18) << "electron_yield_%"
            << std::endl;

  std::cout << std::fixed << std::setprecision(6);

  for (std::size_t index = 0; index < plates.size(); ++index) {
    const PlateKey plate = plates[index];
    const Long64_t gammaEntries =
      GetCount(gammaEntriesByPlate, plate);
    const Long64_t eventsWithElectron =
      GetCount(eventsWithChannelElectronByPlate, plate);
    const Long64_t channelElectrons =
      GetCount(channelElectronsByPlate, plate);

    const Double_t eventYield =
      gammaEntries > 0
        ? static_cast<Double_t>(eventsWithElectron)/gammaEntries
        : 0.0;
    const Double_t electronYield =
      gammaEntries > 0
        ? static_cast<Double_t>(channelElectrons)/gammaEntries
        : 0.0;

    const Double_t eventYieldPercent = 100.0*eventYield;
    const Double_t electronYieldPercent = 100.0*electronYield;
    const std::string label = PlateLabel(plate);

    eventYieldHistogram->GetXaxis()->SetBinLabel(
      static_cast<Int_t>(index + 1), label.c_str());
    eventYieldHistogram->SetBinContent(
      static_cast<Int_t>(index + 1), eventYieldPercent);

    electronYieldHistogram->GetXaxis()->SetBinLabel(
      static_cast<Int_t>(index + 1), label.c_str());
    electronYieldHistogram->SetBinContent(
      static_cast<Int_t>(index + 1), electronYieldPercent);

    std::cout << std::left
              << std::setw(8) << plate.first
              << std::setw(10) << plate.second
              << std::right
              << std::setw(16) << gammaEntries
              << std::setw(22) << eventsWithElectron
              << std::setw(20) << channelElectrons
              << std::setw(16) << eventYield
              << std::setw(16) << eventYieldPercent
              << std::setw(18) << electronYield
              << std::setw(18) << electronYieldPercent
              << std::endl;
  }

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.2f");

  eventYieldHistogram->SetFillColor(kBlue + 1);
  eventYieldHistogram->SetLineColor(kBlue + 2);
  eventYieldHistogram->SetLineWidth(2);
  eventYieldHistogram->SetMinimum(0.0);
  eventYieldHistogram->SetMaximum(
    1.20*eventYieldHistogram->GetMaximum());

  TCanvas* eventYieldCanvas =
    new TCanvas("cPlateDetectionEventYield",
                "Plate detection event yield",
                900,
                700);
  eventYieldCanvas->SetBottomMargin(0.16);
  eventYieldHistogram->Draw("HIST TEXT0");
  eventYieldCanvas->SaveAs("Fig/current/plate_detection_yield.png");

  electronYieldHistogram->SetFillColor(kGreen + 2);
  electronYieldHistogram->SetLineColor(kGreen + 3);
  electronYieldHistogram->SetLineWidth(2);
  electronYieldHistogram->SetMinimum(0.0);
  electronYieldHistogram->SetMaximum(
    1.20*electronYieldHistogram->GetMaximum());

  TCanvas* electronYieldCanvas =
    new TCanvas("cPlateElectronYield",
                "Plate electron yield",
                900,
                700);
  electronYieldCanvas->SetBottomMargin(0.16);
  electronYieldHistogram->Draw("HIST TEXT0");
  electronYieldCanvas->SaveAs("Fig/current/plate_electron_yield.png");

  std::cout << "\nSaved:" << std::endl;
  std::cout << "  plate_detection_yield.png" << std::endl;
  std::cout << "  plate_electron_yield.png" << std::endl;

  file->Close();
  delete file;
}
