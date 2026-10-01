#include "RootOutput.hh"
#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"
#include "TFile.h"
#include "TTree.h"
#include <cstring>
#include <string>
#include <utility>

RootOutput::~RootOutput()
{
  Close();
}

void RootOutput::Open(const char* fileName)
{
  Close();
  fFile = new TFile(fileName, "RECREATE");
  if (fFile->IsZombie()) {
    delete fFile;
    fFile = nullptr;
    G4Exception("RootOutput::Open", "RootOutput001", FatalException,
                "Cannot create the ROOT output file.");
    return;
  }
  fHitTree = new TTree("ElectronChannelHitTree", "Electrons entering MCP_channel");
  fEventTree = new TTree("EventSummaryTree", "One row per event");
  fHitTree->Branch("eventID", &fHit.eventID, "eventID/I");
  fHitTree->Branch("trackID", &fHit.trackID, "trackID/I");
  fHitTree->Branch("side", &fHit.side, "side/I");
  fHitTree->Branch("stackIndex", &fHit.stackIndex, "stackIndex/I");
  fHitTree->Branch("channelCopyNumber", &fHit.channelCopyNumber, "channelCopyNumber/I");
  fHitTree->Branch("kineticEnergy_keV", &fHit.kineticEnergy, "kineticEnergy_keV/D");
  fHitTree->Branch("globalTime_ns", &fHit.globalTime, "globalTime_ns/D");
  fHitTree->Branch("x_mm", &fPosition[0], "x_mm/D");
  fHitTree->Branch("y_mm", &fPosition[1], "y_mm/D");
  fHitTree->Branch("z_mm", &fPosition[2], "z_mm/D");
  fHitTree->Branch("dirX", &fDirection[0], "dirX/D");
  fHitTree->Branch("dirY", &fDirection[1], "dirY/D");
  fHitTree->Branch("dirZ", &fDirection[2], "dirZ/D");
  fHitTree->Branch("creatorProcessName", fCreatorProcess, "creatorProcessName/C");
  fEventTree->Branch("eventID", &fSummary.eventID, "eventID/I");
  fEventTree->Branch("electronChannelPlusCount", &fSummary.electronChannelPlusCount,
                     "electronChannelPlusCount/I");
  fEventTree->Branch("electronChannelMinusCount", &fSummary.electronChannelMinusCount,
                     "electronChannelMinusCount/I");
  fEventTree->Branch("isCoincidence", &fSummary.isCoincidence, "isCoincidence/O");
}

namespace
{
  void BindScanBranches(TTree* tree, ScanSummaryInfo& row, G4bool create)
  {
    const std::pair<const char*, G4int*> fields[] = {
      {"stacksPlus", &row.stacksPlus},
      {"stacksMinus", &row.stacksMinus},
      {"events", &row.events},
      {"electronHits", &row.electronHits},
      {"electronHitsPlus", &row.electronHitsPlus},
      {"electronHitsMinus", &row.electronHitsMinus},
      {"eventsWithHitPlus", &row.eventsWithHitPlus},
      {"eventsWithHitMinus", &row.eventsWithHitMinus},
      {"coincidences", &row.coincidences},
    };
    for (const auto& [name, address] : fields) {
      if (create) {
        tree->Branch(name, address, (std::string(name) + "/I").c_str());
      } else {
        tree->SetBranchAddress(name, address);
      }
    }
  }
}

// Recreated, so a new scan never inherits the rows of the previous one.
void RootOutput::ResetScanFile(const G4String& fileName)
{
  fScanFileName = fileName;

  TFile file(fileName.c_str(), "RECREATE");
  if (file.IsZombie()) {
    fScanFileName.clear();
    G4Exception("RootOutput::ResetScanFile", "RootOutput002", JustWarning,
                "Cannot create the scan ROOT file; scan recording disabled.");
    return;
  }
  TTree tree("ScanTree", "One row per run of a parameter scan");
  ScanSummaryInfo row;
  BindScanBranches(&tree, row, true);
  tree.Write();
}

void RootOutput::AppendScanRow(const ScanSummaryInfo& row)
{
  if (fScanFileName.empty()) return;

  TFile file(fScanFileName.c_str(), "UPDATE");
  if (file.IsZombie()) {
    G4Exception("RootOutput::AppendScanRow", "RootOutput003", JustWarning,
                "Cannot reopen the scan ROOT file; row dropped.");
    return;
  }

  auto* tree = file.Get<TTree>("ScanTree");
  ScanSummaryInfo buffer = row;
  if (!tree) {
    tree = new TTree("ScanTree", "One row per run of a parameter scan");
    BindScanBranches(tree, buffer, true);
  } else {
    BindScanBranches(tree, buffer, false);
  }
  tree->Fill();
  tree->Write("", TObject::kOverwrite);
}

void RootOutput::FillElectronChannelHit(const ElectronChannelHitInfo& hit)
{
  fHit = hit;
  fHit.kineticEnergy /= keV;
  fHit.globalTime /= ns;
  for (int i = 0; i < 3; ++i) {
    fPosition[i] = hit.position[i]/mm;
    fDirection[i] = hit.momentumDirection[i];
  }
  std::strncpy(fCreatorProcess, hit.creatorProcessName.c_str(), sizeof(fCreatorProcess) - 1);
  fCreatorProcess[sizeof(fCreatorProcess) - 1] = '\0';
  fHitTree->Fill();
}

void RootOutput::FillEventSummary(const EventSummaryInfo& summary)
{
  fSummary = summary;
  fEventTree->Fill();
}

void RootOutput::Close()
{
  if (!fFile) return;
  fFile->cd();
  fHitTree->Write();
  fEventTree->Write();
  delete fFile;  // The file owns its trees and closes on destruction.
  fFile = nullptr;
  fHitTree = fEventTree = nullptr;
}
