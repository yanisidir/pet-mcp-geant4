#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
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
}

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
  PrimaryGeneratorAction* generator)
  : fGenerator(generator)
{
  fDirectory = std::make_unique<G4UIdirectory>("/source/");
  fDirectory->SetGuidance("Annihilation source configuration.");

  fParticleCmd = MakeCommand<G4UIcmdWithAString>(
    "/source/particle", this, "Name of the emitted particle.");
  fParticleCmd->SetParameterName("name", false);
  fParticleCmd->SetDefaultValue("gamma");

  fEnergyCmd = MakeCommand<G4UIcmdWithADoubleAndUnit>(
    "/source/energy", this, "Kinetic energy of each emitted particle.");
  fEnergyCmd->SetParameterName("energy", false);
  fEnergyCmd->SetDefaultUnit("keV");
  fEnergyCmd->SetRange("energy > 0.");

  fPositionCmd = MakeCommand<G4UIcmdWith3VectorAndUnit>(
    "/source/position", this, "Centre of the emission point.");
  fPositionCmd->SetParameterName("x", "y", "z", false);
  fPositionCmd->SetDefaultUnit("mm");

  fPositionSigmaCmd = MakeCommand<G4UIcmdWithADoubleAndUnit>(
    "/source/positionSigma", this,
    "Gaussian spread of the emission point on each axis; 0 means a point "
    "source.");
  fPositionSigmaCmd->SetParameterName("sigma", false);
  fPositionSigmaCmd->SetDefaultUnit("mm");
  fPositionSigmaCmd->SetRange("sigma >= 0.");

  fConeAxisCmd = MakeCommand<G4UIcmdWith3Vector>(
    "/source/coneAxis", this, "Axis of the emission cone (need not be unit).");
  fConeAxisCmd->SetParameterName("ux", "uy", "uz", false);

  fConeHalfAngleCmd = MakeCommand<G4UIcmdWithADoubleAndUnit>(
    "/source/coneHalfAngle", this,
    "Half angle of the emission cone; 180 deg gives an isotropic source.");
  fConeHalfAngleCmd->SetParameterName("angle", false);
  fConeHalfAngleCmd->SetDefaultUnit("deg");
  fConeHalfAngleCmd->SetRange("angle > 0.");

  fBackToBackCmd = MakeCommand<G4UIcmdWithABool>(
    "/source/backToBack", this,
    "Emit a second particle in the opposite direction from the same vertex.");
  fBackToBackCmd->SetParameterName("backToBack", true);
  fBackToBackCmd->SetDefaultValue(true);

  fMultiplicityCmd = MakeCommand<G4UIcmdWithAnInteger>(
    "/source/perEvent", this,
    "Number of annihilations (vertices, or pairs with backToBack) per event.");
  fMultiplicityCmd->SetParameterName("count", false);
  fMultiplicityCmd->SetRange("count >= 1");

  fPrintCmd = MakeCommand<G4UIcmdWithoutParameter>(
    "/source/print", this, "Print the current source configuration.");
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() = default;

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String value)
{
  if (command == fParticleCmd.get()) {
    fGenerator->SetParticle(value);
  } else if (command == fEnergyCmd.get()) {
    fGenerator->SetEnergy(fEnergyCmd->GetNewDoubleValue(value));
  } else if (command == fPositionCmd.get()) {
    fGenerator->SetPosition(fPositionCmd->GetNew3VectorValue(value));
  } else if (command == fPositionSigmaCmd.get()) {
    fGenerator->SetPositionSigma(fPositionSigmaCmd->GetNewDoubleValue(value));
  } else if (command == fConeAxisCmd.get()) {
    fGenerator->SetConeAxis(fConeAxisCmd->GetNew3VectorValue(value));
  } else if (command == fConeHalfAngleCmd.get()) {
    fGenerator->SetConeHalfAngle(fConeHalfAngleCmd->GetNewDoubleValue(value));
  } else if (command == fBackToBackCmd.get()) {
    fGenerator->SetBackToBack(G4UIcmdWithABool::GetNewBoolValue(value));
  } else if (command == fMultiplicityCmd.get()) {
    fGenerator->SetParticlesPerEvent(
      G4UIcmdWithAnInteger::GetNewIntValue(value));
  } else if (command == fPrintCmd.get()) {
    fGenerator->PrintConfiguration();
  }
}
