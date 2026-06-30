#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(new G4ParticleGun(1)),
    fMessenger(nullptr),
    fEnergy(511.0*keV),
    fPosition(0.0, 0.0, 0.0),
    fDirection(0.0, 0.0, 1.0)
{
  SetParticleName("gamma");
  SetEnergy(fEnergy);
  SetPosition(fPosition);
  SetDirection(fDirection);

  fMessenger = new PrimaryGeneratorMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fMessenger;
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if (!event) {
    return;
  }

  // fParticleGun->SetParticleEnergy(fEnergy);
  // fParticleGun->SetParticlePosition(fPosition);

  fParticleGun->SetParticleMomentumDirection(fDirection);
  fParticleGun->GeneratePrimaryVertex(event);

  fParticleGun->SetParticleMomentumDirection(-fDirection);

  fParticleGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::SetParticleName(
  const G4String& particleName)
{
  G4ParticleDefinition* particle =
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);

  if (!particle) {
    G4ExceptionDescription message;
    message << "Unknown particle name: " << particleName
            << ". The previous particle definition is kept.";
    G4Exception("PrimaryGeneratorAction::SetParticleName",
                "PrimaryGenerator001",
                JustWarning,
                message);
    return;
  }

  fParticleGun->SetParticleDefinition(particle);
}

void PrimaryGeneratorAction::SetEnergy(G4double energy)
{
  if (energy < 0.0) {
    G4ExceptionDescription message;
    message << "Invalid kinetic energy: " << energy/keV
            << " keV. The previous energy is kept.";
    G4Exception("PrimaryGeneratorAction::SetEnergy",
                "PrimaryGenerator002",
                JustWarning,
                message);
    return;
  }

  fEnergy = energy;
  fParticleGun->SetParticleEnergy(fEnergy);
}

void PrimaryGeneratorAction::SetPosition(
  const G4ThreeVector& position)
{
  fPosition = position;
  fParticleGun->SetParticlePosition(fPosition);
}

void PrimaryGeneratorAction::SetDirection(
  const G4ThreeVector& direction)
{
  if (direction.mag2() == 0.0) {
    G4Exception("PrimaryGeneratorAction::SetDirection",
                "PrimaryGenerator003",
                JustWarning,
                "Invalid zero momentum direction. "
                "The previous direction is kept.");
    return;
  }

  fDirection = direction.unit();
  fParticleGun->SetParticleMomentumDirection(fDirection);
}
