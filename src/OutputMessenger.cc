#include "OutputMessenger.hh"
#include "RunAction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

OutputMessenger::OutputMessenger(RunAction* runAction)
  : fRunAction(runAction)
{
  fDirectory = std::make_unique<G4UIdirectory>("/output/");
  fDirectory->SetGuidance("ROOT output configuration.");

  fFileCmd = std::make_unique<G4UIcmdWithAString>("/output/file", this);
  fFileCmd->SetGuidance(
    "Name of the per-run ROOT file; it is recreated at every beamOn.");
  fFileCmd->SetParameterName("name", false);
  fFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fScanFileCmd = std::make_unique<G4UIcmdWithAString>("/output/scanFile", this);
  fScanFileCmd->SetGuidance(
    "Start a scan: create this file and append one row per run to its "
    "ScanTree. Call it once, before the first beamOn of the scan.");
  fScanFileCmd->SetParameterName("name", false);
  fScanFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

OutputMessenger::~OutputMessenger() = default;

void OutputMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  if (command == fFileCmd.get()) {
    fRunAction->SetFileName(value);
  } else if (command == fScanFileCmd.get()) {
    fRunAction->ResetScanFile(value);
  }
}
