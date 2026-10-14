#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Taux de coincidence en fonction du nombre de plaques MCP par cote,
// lu dans le ScanTree cumulatif produit par un balayage /mcp/stacks.
//
// Definition du comptage (celle de la simulation, cf. EventAction.cc) :
// un evenement compte au maximum UNE coincidence, quel que soit le
// nombre d'electrons ou de plaques touches de chaque cote. Une
// coincidence = au moins un hit canal cote +z ET au moins un hit canal
// cote -z dans le meme evenement.
//
// Intervalle de confiance de Wilson plutot que l'approximation
// gaussienne p +/- z*sqrt(p(1-p)/n) : a taux faible l'intervalle de
// Wald sous-couvre et peut sortir de [0,1], alors que Wilson inverse le
// test du score et reste borne (cf. Wilson 1927, Agresti & Coull 1998).
void AnalyzeCoincidenceVsStacks(
  const char* fileName = "scan_stacks.root",
  const char* plotName = "coincidence_vs_stacks.png")
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
  Int_t eventsWithHitPlus = 0;
  Int_t eventsWithHitMinus = 0;
  Int_t coincidences = 0;

  scan->SetBranchAddress("stacksPlus", &stacksPlus);
  scan->SetBranchAddress("stacksMinus", &stacksMinus);
  scan->SetBranchAddress("events", &events);
  scan->SetBranchAddress("eventsWithHitPlus", &eventsWithHitPlus);
  scan->SetBranchAddress("eventsWithHitMinus", &eventsWithHitMinus);
  scan->SetBranchAddress("coincidences", &coincidences);

  const Long64_t points = scan->GetEntries();
  if (points == 0) {
    std::cerr << "Erreur : ScanTree vide dans " << fileName << '\n';
    return;
  }

  std::cout << "Fichier : " << fileName << '\n';
  std::cout << std::right
            << std::setw(10) << "N+N"
            << std::setw(12) << "evenements"
            << std::setw(12) << "eff +z %"
            << std::setw(12) << "eff -z %"
            << std::setw(10) << "k coinc"
            << std::setw(14) << "eff coinc %"
            << std::setw(26) << "IC 95% Wilson (%)"
            << std::setw(10) << "C"
            << '\n'
            << std::string(106, '-') << '\n';

  std::vector<double> x, y, eyLow, eyHigh;

  for (Long64_t i = 0; i < points; ++i) {
    scan->GetEntry(i);

    if (events <= 0) continue;

    const double effPlus = double(eventsWithHitPlus)/events;
    const double effMinus = double(eventsWithHitMinus)/events;
    const double effCoincidence = double(coincidences)/events;

    const double lo = TEfficiency::Wilson(events, coincidences, 0.95, false);
    const double hi = TEfficiency::Wilson(events, coincidences, 0.95, true);

    // C = eff_coinc / (eff_+z * eff_-z). C > 1 signale que les deux
    // cotes sont correles : les deux photons d'un evenement partagent
    // une seule direction d'emission (dos-a-dos), donc une direction
    // favorable a l'un l'est aussi a l'autre.
    const double independent = effPlus*effMinus;
    const double correlation =
      independent > 0.0 ? effCoincidence/independent : 0.0;

    std::ostringstream config;
    config << stacksPlus << "+" << stacksMinus;
    std::ostringstream ci;
    ci << "[" << std::fixed << std::setprecision(4)
       << 100.0*lo << ", " << 100.0*hi << "]";

    std::cout << std::right << std::fixed << std::setprecision(4)
              << std::setw(10) << config.str()
              << std::setw(12) << events
              << std::setw(12) << 100.0*effPlus
              << std::setw(12) << 100.0*effMinus
              << std::setw(10) << coincidences
              << std::setw(14) << 100.0*effCoincidence
              << std::setw(26) << ci.str()
              << std::setw(10) << std::setprecision(3) << correlation
              << '\n';

    if (stacksPlus != stacksMinus) {
      std::cout << "  (point asymetrique, exclu du graphe)\n";
      continue;
    }

    x.push_back(stacksPlus);
    y.push_back(100.0*effCoincidence);
    eyLow.push_back(100.0*(effCoincidence - lo));
    eyHigh.push_back(100.0*(hi - effCoincidence));
  }

  if (x.empty()) {
    std::cerr << "Aucun point symetrique a tracer.\n";
    return;
  }

  auto* graph = new TGraphAsymmErrors(
    x.size(), x.data(), y.data(), nullptr, nullptr,
    eyLow.data(), eyHigh.data());
  graph->SetTitle(
    "Taux de coincidence selon la profondeur des empilements;"
    "Nombre de plaques MCP par cote (N+N);Taux de coincidence (%)");
  graph->SetMarkerStyle(21);
  graph->SetMarkerSize(1.3);
  graph->SetMarkerColor(kRed + 1);
  graph->SetLineColor(kRed + 1);
  graph->SetLineWidth(2);
  graph->SetMinimum(0.0);

  auto* canvas = new TCanvas(
    "coincidenceVsStacks", "Coincidences selon le nombre de plaques",
    900, 650);
  canvas->SetGrid();
  graph->Draw("APL");
  canvas->Update();
  canvas->SaveAs(plotName);

  std::cout << "\nFigure : " << plotName << '\n';
}
