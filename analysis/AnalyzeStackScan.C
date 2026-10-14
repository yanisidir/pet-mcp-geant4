#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Fichier cumulatif produit par macros/scan.mac.
void AnalyzeStackScan(const char* fileName = "scan_stacks.root")
{
  TFile file(fileName, "READ");
  if (file.IsZombie()) {
    std::cerr << "Erreur : impossible d'ouvrir " << fileName << '\n';
    return;
  }

  auto* scan = file.Get<TTree>("ScanTree");
  if (!scan) {
    std::cerr << "Erreur : arbre ScanTree absent dans " << fileName << '\n';
    return;
  }

  Int_t stacksPlus = 0;
  Int_t stacksMinus = 0;
  Int_t events = 0;
  Int_t electronHits = 0;
  Int_t electronHitsPlus = 0;
  Int_t electronHitsMinus = 0;
  Int_t coincidences = 0;

  scan->SetBranchAddress("stacksPlus", &stacksPlus);
  scan->SetBranchAddress("stacksMinus", &stacksMinus);
  scan->SetBranchAddress("events", &events);
  scan->SetBranchAddress("electronHits", &electronHits);
  scan->SetBranchAddress("electronHitsPlus", &electronHitsPlus);
  scan->SetBranchAddress("electronHitsMinus", &electronHitsMinus);
  scan->SetBranchAddress("coincidences", &coincidences);

  const Long64_t points = scan->GetEntries();
  if (points == 0) {
    std::cerr << "Erreur : ScanTree vide dans " << fileName << '\n';
    return;
  }

  std::vector<double> x, y, ey;

  std::cout << std::left
            << std::setw(10) << "plaques"
            << std::setw(12) << "evenements"
            << std::setw(12) << "electrons"
            << std::setw(10) << "+z"
            << std::setw(10) << "-z"
            << std::setw(16) << "e- / evenement"
            << "coincidences\n"
            << std::string(78, '-') << '\n';

  for (Long64_t i = 0; i < points; ++i) {
    scan->GetEntry(i);

    if (stacksPlus != stacksMinus) {
      std::cout << "(point asymetrique " << stacksPlus << "/" << stacksMinus
                << " ignore dans le graphe)\n";
    }

    const double perEvent =
      events > 0 ? double(electronHits)/events : 0.0;

    std::cout << std::left
              << std::setw(10) << stacksPlus
              << std::setw(12) << events
              << std::setw(12) << electronHits
              << std::setw(10) << electronHitsPlus
              << std::setw(10) << electronHitsMinus
              << std::setw(16) << std::fixed << std::setprecision(4) << perEvent
              << coincidences << '\n'
              << std::defaultfloat;

    if (stacksPlus != stacksMinus) continue;

    x.push_back(stacksPlus);
    y.push_back(electronHits);
    ey.push_back(std::sqrt(double(electronHits)));
  }

  if (x.empty()) {
    std::cerr << "Aucun point symetrique a tracer.\n";
    return;
  }

  auto* graph = new TGraphErrors(
    x.size(), x.data(), y.data(), nullptr, ey.data());
  graph->SetTitle(
    "Electrons atteignant un canal selon la profondeur de l'empilement;"
    "Nombre de plaques MCP par cote;Electrons atteignant un canal");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.2);

  auto* canvas = new TCanvas(
    "stackScanCanvas", "Balayage du nombre de plaques MCP", 900, 600);
  canvas->SetGrid();
  graph->Draw("APL");
  canvas->Update();
  gSystem->mkdir("figures", kTRUE);
  canvas->SaveAs("figures/scan_stacks.png");
}
