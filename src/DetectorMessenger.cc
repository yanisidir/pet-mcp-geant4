#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

namespace
{
  template <typename Command>
  std::unique_ptr<Command> MakeCommand(const char* path,
                                       G4UImessenger* messenger,
                                       const char* guidance)
  {
    auto command = std::make_unique<Command>(path, messenger);
    command->SetGuidance(guidance);
    command->AvailableForStates(G4State_PreInit, G4State_Idle);
    return command;
  }

  template <typename Command>
  std::unique_ptr<Command> MakeLengthCommand(const char* path,
                                             G4UImessenger* messenger,
                                             const char* guidance,
                                             const char* parameterName,
                                             const char* defaultUnit)
  {
    auto command = MakeCommand<Command>(path, messenger, guidance);
    command->SetParameterName(parameterName, false);
    command->SetDefaultUnit(defaultUnit);
    return command;
  }
}

DetectorMessenger::DetectorMessenger(DetectorConstruction* detector)
  : fDetector(detector)
{
  fDirectory = std::make_unique<G4UIdirectory>("/mcp/");
  fDirectory->SetGuidance("MCP geometry configuration.");

  fStacksPlusCmd = MakeCommand<G4UIcmdWithAnInteger>(
    "/mcp/stacksPlus", this,
    "Number of MCP plates stacked on the +z side (0 removes the side).");
  fStacksPlusCmd->SetParameterName("count", false);
  fStacksPlusCmd->SetRange("count >= 0");

  fStacksMinusCmd = MakeCommand<G4UIcmdWithAnInteger>(
    "/mcp/stacksMinus", this,
    "Number of MCP plates stacked on the -z side (0 removes the side).");
  fStacksMinusCmd->SetParameterName("count", false);
  fStacksMinusCmd->SetRange("count >= 0");

  fStacksCmd = MakeCommand<G4UIcmdWithAnInteger>(
    "/mcp/stacks", this,
    "Shortcut: same number of plates on both sides.");
  fStacksCmd->SetParameterName("count", false);
  fStacksCmd->SetRange("count >= 0");

  fFirstPlateZCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/firstPlateZ", this,
    "Distance from the origin to the centre of the first plate of a stack.",
    "z", "cm");

  fStackGapCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/stackGap", this,
    "Vacuum gap between two consecutive plates of the same stack.",
    "gap", "um");

  fThicknessCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/plateThickness", this, "Thickness of one MCP plate.",
    "thickness", "um");

  fRadiusCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/plateRadius", this, "Radius of one MCP plate.", "radius", "mm");

  fChannelDiameterCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/channelDiameter", this, "Diameter of a vacuum channel.",
    "diameter", "um");

  fChannelPitchCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/channelPitch", this, "Hexagonal pitch between channels.",
    "pitch", "um");

  fChannelAngleCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/channelAngle", this,
    "Channel bias angle in the x-z plane; the sign follows the side.",
    "angle", "deg");

  fChevronCmd = MakeCommand<G4UIcmdWithABool>(
    "/mcp/chevron", this,
    "Alternate the channel angle sign from one plate to the next "
    "(chevron stack).");
  fChevronCmd->SetParameterName("chevron", true);
  fChevronCmd->SetDefaultValue(true);

  fEdgeMarginCmd = MakeCommand<G4UIcmdWithADouble>(
    "/mcp/edgeMarginFraction", this,
    "Fraction of the plate radius kept free of channels at the rim.");
  fEdgeMarginCmd->SetParameterName("fraction", false);
  fEdgeMarginCmd->SetRange("fraction >= 0. && fraction < 1.");

  fWorldSizeCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/worldSize", this,
    "Full side of the cubic world; enlarged automatically if too small.",
    "size", "cm");

  fProductionCutCmd = MakeLengthCommand<G4UIcmdWithADoubleAndUnit>(
    "/mcp/productionCut", this,
    "Secondary production cut inside the MCP region.", "cut", "um");

  fKillGammaCmd = MakeCommand<G4UIcmdWithABool>(
    "/mcp/killGammaOnExit", this,
    "Kill a gamma when it leaves the plate it is in. Set it to false for "
    "stacks of more than one plate so gammas can reach the next plate.");
  fKillGammaCmd->SetParameterName("kill", true);
  fKillGammaCmd->SetDefaultValue(true);

  fPrintCmd = MakeCommand<G4UIcmdWithoutParameter>(
    "/mcp/print", this, "Print the current MCP configuration.");
}

DetectorMessenger::~DetectorMessenger() = default;

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  G4bool geometryChanged = true;

  if (command == fStacksPlusCmd.get()) {
    fDetector->SetStackCount(1, G4UIcmdWithAnInteger::GetNewIntValue(value));
  } else if (command == fStacksMinusCmd.get()) {
    fDetector->SetStackCount(-1, G4UIcmdWithAnInteger::GetNewIntValue(value));
  } else if (command == fStacksCmd.get()) {
    const G4int count = G4UIcmdWithAnInteger::GetNewIntValue(value);
    fDetector->SetStackCount(1, count);
    fDetector->SetStackCount(-1, count);
  } else if (command == fFirstPlateZCmd.get()) {
    fDetector->SetFirstPlateZ(fFirstPlateZCmd->GetNewDoubleValue(value));
  } else if (command == fStackGapCmd.get()) {
    fDetector->SetStackGap(fStackGapCmd->GetNewDoubleValue(value));
  } else if (command == fThicknessCmd.get()) {
    fDetector->SetPlateThickness(fThicknessCmd->GetNewDoubleValue(value));
  } else if (command == fRadiusCmd.get()) {
    fDetector->SetPlateRadius(fRadiusCmd->GetNewDoubleValue(value));
  } else if (command == fChannelDiameterCmd.get()) {
    fDetector->SetChannelDiameter(
      fChannelDiameterCmd->GetNewDoubleValue(value));
  } else if (command == fChannelPitchCmd.get()) {
    fDetector->SetChannelPitch(fChannelPitchCmd->GetNewDoubleValue(value));
  } else if (command == fChannelAngleCmd.get()) {
    fDetector->SetChannelAngle(fChannelAngleCmd->GetNewDoubleValue(value));
  } else if (command == fChevronCmd.get()) {
    fDetector->SetChevron(G4UIcmdWithABool::GetNewBoolValue(value));
  } else if (command == fEdgeMarginCmd.get()) {
    fDetector->SetEdgeMarginFraction(
      G4UIcmdWithADouble::GetNewDoubleValue(value));
  } else if (command == fWorldSizeCmd.get()) {
    fDetector->SetWorldSize(fWorldSizeCmd->GetNewDoubleValue(value));
  } else if (command == fProductionCutCmd.get()) {
    fDetector->SetProductionCut(fProductionCutCmd->GetNewDoubleValue(value));
  } else if (command == fKillGammaCmd.get()) {
    fDetector->SetKillGammaOnExit(G4UIcmdWithABool::GetNewBoolValue(value));
    geometryChanged = false;  // Stepping rule only, the volumes are unchanged.
  } else if (command == fPrintCmd.get()) {
    fDetector->PrintConfiguration();
    geometryChanged = false;
  }

  // A geometry change after /run/initialize needs an explicit rebuild.
  const G4ApplicationState state =
    G4StateManager::GetStateManager()->GetCurrentState();
  if (geometryChanged && state != G4State_PreInit) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}
