#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

namespace
{
  void DefineElectricFieldUnits()
  {
    static G4bool unitsDefined = false;
    if (unitsDefined) {
      return;
    }

    new G4UnitDefinition("kilovolt/mm",
                         "kV/mm",
                         "Electric field",
                         kilovolt/mm);
    new G4UnitDefinition("volt/mm",
                         "V/mm",
                         "Electric field",
                         volt/mm);
    unitsDefined = true;
  }
}

DetectorMessenger::DetectorMessenger(DetectorConstruction* detector)
  : G4UImessenger(),
    fDetector(detector),
    fMcpDirectory(new G4UIdirectory("/mcp/")),
    fMaterialCommand(nullptr),
    fEnableStackFieldsCommand(nullptr),
    fStackFieldMagnitudeCommand(nullptr),
    fFieldRegionRadialMarginCommand(nullptr),
    fFieldRegionZMarginCommand(nullptr)
{
  DefineElectricFieldUnits();

  fMcpDirectory->SetGuidance("MCP detector controls.");

  fMaterialCommand =
    new G4UIcmdWithAString("/mcp/setMaterial", this);
  fMaterialCommand->SetGuidance("Select the material used by all MCP bodies.");
  fMaterialCommand->SetParameterName("materialName", false);
  fMaterialCommand->SetCandidates(
    "GlassLead25Perc GlassLead30Perc GlassLead35Perc "
    "GlassLead40Perc GlassLead45Perc GlassLead50Perc "
    "GlassLead60Perc GlassLead75Perc");
  fMaterialCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEnableStackFieldsCommand =
    new G4UIcmdWithABool("/mcp/enableStackElectricFields", this);
  fEnableStackFieldsCommand->SetGuidance(
    "Enable or disable the local electric fields around MCP stacks.");
  fEnableStackFieldsCommand->SetParameterName("enabled", false);
  fEnableStackFieldsCommand->AvailableForStates(G4State_PreInit,
                                                G4State_Idle);

  fStackFieldMagnitudeCommand =
    new G4UIcmdWithADoubleAndUnit("/mcp/stackFieldMagnitude", this);
  fStackFieldMagnitudeCommand->SetGuidance(
    "Set the magnitude of the local MCP stack electric fields.");
  fStackFieldMagnitudeCommand->SetParameterName("magnitude", false);
  fStackFieldMagnitudeCommand->SetDefaultUnit("kV/mm");
  fStackFieldMagnitudeCommand->AvailableForStates(G4State_PreInit,
                                                  G4State_Idle);

  fFieldRegionRadialMarginCommand =
    new G4UIcmdWithADoubleAndUnit("/mcp/fieldRegionRadialMargin", this);
  fFieldRegionRadialMarginCommand->SetGuidance(
    "Set the radial margin around each MCP stack field region.");
  fFieldRegionRadialMarginCommand->SetParameterName("margin", false);
  fFieldRegionRadialMarginCommand->SetDefaultUnit("mm");
  fFieldRegionRadialMarginCommand->SetUnitCategory("Length");
  fFieldRegionRadialMarginCommand->AvailableForStates(G4State_PreInit,
                                                       G4State_Idle);

  fFieldRegionZMarginCommand =
    new G4UIcmdWithADoubleAndUnit("/mcp/fieldRegionZMargin", this);
  fFieldRegionZMarginCommand->SetGuidance(
    "Set the z margin before and after each MCP stack field region.");
  fFieldRegionZMarginCommand->SetParameterName("margin", false);
  fFieldRegionZMarginCommand->SetDefaultUnit("mm");
  fFieldRegionZMarginCommand->SetUnitCategory("Length");
  fFieldRegionZMarginCommand->AvailableForStates(G4State_PreInit,
                                                  G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fFieldRegionZMarginCommand;
  delete fFieldRegionRadialMarginCommand;
  delete fStackFieldMagnitudeCommand;
  delete fEnableStackFieldsCommand;
  delete fMaterialCommand;
  delete fMcpDirectory;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command,
                                    G4String newValue)
{
  if (!fDetector) {
    return;
  }

  if (command == fMaterialCommand) {
    fDetector->SetDetectorMaterial(newValue);
  } else if (command == fEnableStackFieldsCommand) {
    fDetector->SetStackElectricFieldsEnabled(
      fEnableStackFieldsCommand->GetNewBoolValue(newValue));
  } else if (command == fStackFieldMagnitudeCommand) {
    const G4double magnitude =
      fStackFieldMagnitudeCommand->GetNewDoubleValue(newValue);
    if (magnitude <= 0.0) {
      G4cout << "WARNING: /mcp/stackFieldMagnitude must be positive. "
             << "Keeping previous value." << G4endl;
      return;
    }
    fDetector->SetStackFieldMagnitude(magnitude);
  } else if (command == fFieldRegionRadialMarginCommand) {
    const G4double margin =
      fFieldRegionRadialMarginCommand->GetNewDoubleValue(newValue);
    if (margin < 0.0) {
      G4cout << "WARNING: /mcp/fieldRegionRadialMargin cannot be negative. "
             << "Keeping previous value." << G4endl;
      return;
    }
    fDetector->SetFieldRegionRadialMargin(margin);
  } else if (command == fFieldRegionZMarginCommand) {
    const G4double margin =
      fFieldRegionZMarginCommand->GetNewDoubleValue(newValue);
    if (margin < 0.0) {
      G4cout << "WARNING: /mcp/fieldRegionZMargin cannot be negative. "
             << "Keeping previous value." << G4endl;
      return;
    }
    fDetector->SetFieldRegionZMargin(margin);
  }
}
