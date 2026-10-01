#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH
#include "RootOutput.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class OutputMessenger;

// Own the output and the run counters (serial execution).
class RunAction : public G4UserRunAction
{
public:
  explicit RunAction(const DetectorConstruction* detector);
  ~RunAction() override;
  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;
  void RecordEvent(const EventSummaryInfo& summary);
  RootOutput& GetOutput() { return fOutput; }

  void SetFileName(const G4String& name) { fFileName = name; }
  void ResetScanFile(const G4String& name) { fOutput.ResetScanFile(name); }

private:
  const DetectorConstruction* fDetector;
  OutputMessenger* fMessenger;
  RootOutput fOutput;
  G4String fFileName = "mcp_output.root";
  G4int fPlusCount = 0;
  G4int fMinusCount = 0;
  G4int fCoincidenceCount = 0;
  G4int fElectronHitsPlus = 0;
  G4int fElectronHitsMinus = 0;
};
#endif
