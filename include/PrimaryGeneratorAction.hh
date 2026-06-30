#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4ParticleGun;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  // Constructeur.
  PrimaryGeneratorAction();

  virtual ~PrimaryGeneratorAction();

  // Fonction appelée automatiquement à chaque événement.
  virtual void GeneratePrimaries(G4Event* event);

  void SetParticleName(const G4String& particleName);
  void SetEnergy(G4double energy);
  void SetPosition(const G4ThreeVector& position);
  void SetDirection(const G4ThreeVector& direction);

private:
  // =====================================================
  // OUTILS DE GENERATION
  // =====================================================

  G4ParticleGun* fParticleGun;

  PrimaryGeneratorMessenger* fMessenger;

  G4double fEnergy;

  G4ThreeVector fPosition;

  G4ThreeVector fDirection;
};

#endif
