#ifndef ROOT_OUTPUT_HH
#define ROOT_OUTPUT_HH
#include "ElectronChannelHitInfo.hh"
#include "EventSummaryInfo.hh"
#include "ScanSummaryInfo.hh"
#include "globals.hh"
class TFile;
class TTree;

// Two ROOT trees, owned by RunAction. Conversion to keV/mm/ns occurs here.
class RootOutput
{
public:
  RootOutput() = default;
  ~RootOutput();
  RootOutput(const RootOutput&) = delete;
  RootOutput& operator=(const RootOutput&) = delete;
  void Open(const char* fileName);

  void ResetScanFile(const G4String& fileName);
  void AppendScanRow(const ScanSummaryInfo& row);
  G4bool HasScanFile() const { return !fScanFileName.empty(); }
  void FillElectronChannelHit(const ElectronChannelHitInfo& hit);
  void FillEventSummary(const EventSummaryInfo& summary);
  void Close();
private:
  TFile* fFile = nullptr;
  TTree* fHitTree = nullptr;
  TTree* fEventTree = nullptr;
  ElectronChannelHitInfo fHit;
  EventSummaryInfo fSummary;
  G4String fScanFileName;
  G4double fPosition[3]{};
  G4double fDirection[3]{};
  char fCreatorProcess[64]{};
};
#endif
