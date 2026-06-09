#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(new G4ParticleGun(1)),
    fMessenger(nullptr),
    fEnergy(511.0*keV),
    fPosition(0.0, 0.0, 0.0)
{
  fMessenger = new G4GenericMessenger(this,
                                      "/minimal/gun/",
                                      "Simplified PET source settings");

  fMessenger->DeclarePropertyWithUnit("energy", "MeV", fEnergy,
                                      "Energy of each annihilation photon");
  fMessenger->DeclarePropertyWithUnit("position", "cm", fPosition,
                                      "Common annihilation vertex");

  G4ParticleDefinition* gamma =
    G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  if (!gamma) {
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction",
                "PetSource001",
                FatalException,
                "The Geant4 gamma particle definition is unavailable.");
  }

  fParticleGun->SetParticleDefinition(gamma);
  fParticleGun->SetParticleEnergy(fEnergy);
  fParticleGun->SetParticlePosition(fPosition);
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

  fParticleGun->SetParticleEnergy(fEnergy);
  fParticleGun->SetParticlePosition(fPosition);

  // Premier photon: vers la pile MCP située du côté +z.
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  fParticleGun->GeneratePrimaryVertex(event);

  // Second photon: même événement et même vertex, direction opposée.
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));
  fParticleGun->GeneratePrimaryVertex(event);
}
