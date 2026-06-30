#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
  PrimaryGeneratorAction* generator)
  : G4UImessenger(),
    fGenerator(generator),
    fGunDirectory(new G4UIdirectory("/gun/")),
    fParticleCommand(nullptr),
    fEnergyCommand(nullptr),
    fPositionCommand(nullptr),
    fDirectionCommand(nullptr)
{
  fGunDirectory->SetGuidance("Primary particle gun control.");

  fParticleCommand = new G4UIcmdWithAString("/gun/particle", this);
  fParticleCommand->SetGuidance("Set the primary particle name.");
  fParticleCommand->SetGuidance("Examples: gamma, e-, mu-, proton.");
  fParticleCommand->SetParameterName("particleName", false); 
  fParticleCommand->AvailableForStates(G4State_PreInit, G4State_Idle); // Available in PreInit and Idle states

  fEnergyCommand =
    new G4UIcmdWithADoubleAndUnit("/gun/energy", this);
  fEnergyCommand->SetGuidance("Set the primary kinetic energy.");
  fEnergyCommand->SetParameterName("energy", false);
  fEnergyCommand->SetDefaultUnit("keV");
  fEnergyCommand->SetUnitCategory("Energy");
  fEnergyCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPositionCommand =
    new G4UIcmdWith3VectorAndUnit("/gun/position", this);
  fPositionCommand->SetGuidance("Set the primary source position.");
  fPositionCommand->SetParameterName("x", "y", "z", false);
  fPositionCommand->SetDefaultUnit("cm");
  fPositionCommand->SetUnitCategory("Length");
  fPositionCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDirectionCommand =
    new G4UIcmdWith3Vector("/gun/direction", this);
  fDirectionCommand->SetGuidance("Set the primary momentum direction.");
  fDirectionCommand->SetGuidance("The vector is normalized internally.");
  fDirectionCommand->SetParameterName("dx", "dy", "dz", false);
  fDirectionCommand->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fDirectionCommand;
  delete fPositionCommand;
  delete fEnergyCommand;
  delete fParticleCommand;
  delete fGunDirectory;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  if (!fGenerator) {
    return;
  }

  if (command == fParticleCommand) {
    fGenerator->SetParticleName(newValue);
  } else if (command == fEnergyCommand) {
    fGenerator->SetEnergy(
      fEnergyCommand->GetNewDoubleValue(newValue));
  } else if (command == fPositionCommand) {
    fGenerator->SetPosition(
      fPositionCommand->GetNew3VectorValue(newValue));
  } else if (command == fDirectionCommand) {
    fGenerator->SetDirection(
      fDirectionCommand->GetNew3VectorValue(newValue));
  }
}
