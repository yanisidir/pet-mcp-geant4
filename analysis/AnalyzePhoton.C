#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>

void AnalyzePhoton(const char* fileName = "mcp_output.root")
{
  TFile file(fileName, "READ");

  auto* hits = file.Get<TTree>("ElectronChannelHitTree");
  auto* events = file.Get<TTree>("EventSummaryTree");

  if (!hits || !events) {
    std::cerr << "Erreur : arbres ROOT absents dans " << fileName << '\n';
    return;
  }

  const Long64_t numberOfEvents = events->GetEntries();
  const Long64_t numberOfPhotons = numberOfEvents;
  const Long64_t numberOfElectrons = hits->GetEntries();

  const Long64_t photonsWithElectron =
    events->GetEntries("electronChannelPlusCount > 0") +
    events->GetEntries("electronChannelMinusCount > 0");

  const double rate = numberOfPhotons > 0
    ? 100.0 * photonsWithElectron / numberOfPhotons
    : 0.0;

  std::cout << "Photons de 511 keV envoyes : " << numberOfPhotons << '\n'
            << "Electrons atteignant un canal : " << numberOfElectrons << '\n'
            << "Photons produisant au moins un electron atteignant un canal : "
            << photonsWithElectron << '\n'
            << std::fixed << std::setprecision(4)
            << "Taux : " << rate << " %\n";

  double energy = 0.0;
  hits->SetBranchAddress("kineticEnergy_keV", &energy);

  auto* spectrum = new TH1D(
    "electronEnergy",
    "Spectre des electrons atteignant un canal;"
    "Energie cinetique [keV];Nombre d'electrons",
    100, 0.0, 511.0);
  spectrum->SetDirectory(nullptr);

  for (Long64_t i = 0; i < numberOfElectrons; ++i) {
    hits->GetEntry(i);
    spectrum->Fill(energy);
  }

  auto* canvas = new TCanvas(
    "electronEnergyCanvas", "Spectre en energie des electrons", 900, 600);
  spectrum->Draw();
  canvas->Update();
  gSystem->mkdir("figures", kTRUE);
  canvas->SaveAs("figures/spectre.png");

}
