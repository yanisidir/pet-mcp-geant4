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
  // Elle génère deux photons gamma de même énergie,
  // au même point, dans les directions +z et -z.
  virtual void GeneratePrimaries(G4Event* event);

private:
  // =====================================================
  // OUTILS DE GENERATION
  // =====================================================

  // Générateur utilisé deux fois par événement pour créer
  // les deux photons de l'annihilation simplifiée.
  G4ParticleGun* fParticleGun;

  // Messenger permettant de modifier les paramètres
  // depuis une macro Geant4.
  G4GenericMessenger* fMessenger;

  // =====================================================
  // PARAMETRES DE LA SOURCE
  // =====================================================

  // Energie de chacun des deux photons.
  G4double fEnergy;

  // Vertex commun aux deux photons.
  G4ThreeVector fPosition;
};

#endif
