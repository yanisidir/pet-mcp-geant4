#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4GenericMessenger;
class G4ParticleGun;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  // Constructeur.
  PrimaryGeneratorAction();

  virtual ~PrimaryGeneratorAction();

  // Fonction appelée automatiquement à chaque événement.
  //
  // Elle génère un photon gamma dirigé selon +z.
  virtual void GeneratePrimaries(G4Event* event);

private:
  // =====================================================
  // OUTILS DE GENERATION
  // =====================================================

  // Générateur du photon primaire.
  G4ParticleGun* fParticleGun;

  // Messenger permettant de modifier les paramètres
  // depuis une macro Geant4.
  G4GenericMessenger* fMessenger;

  // =====================================================
  // PARAMETRES DE LA SOURCE
  // =====================================================

  // Energie du photon.
  G4double fEnergy;

  // Position initiale du photon.
  G4ThreeVector fPosition;
};

#endif
