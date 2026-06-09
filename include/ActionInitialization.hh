#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

// ActionInitialization est la classe chargée d'enregistrer
// toutes les actions utilisateur Geant4.
//
// Son rôle est de dire au RunManager quelles classes utiliser
// pour gérer :
//
// - la génération des particules ;
// - le début et la fin des runs ;
// - le début et la fin des événements ;
// - le stepping ;
// - le tracking.
//
// Cette classe est appelée lors de l'initialisation de la simulation
// et constitue le point de connexion entre tous les modules.
class ActionInitialization : public G4VUserActionInitialization
{
public:

  // Constructeur.
  //
  // detector :
  //   pointeur vers la géométrie du détecteur.
  //
  // Il sera transmis à certaines actions
  // (par exemple SteppingAction) afin qu'elles puissent
  // accéder aux volumes du MCP.
  explicit ActionInitialization(
    const DetectorConstruction* detector);

  virtual ~ActionInitialization();

  // Fonction appelée par Geant4 pour construire
  // l'ensemble des actions utilisateur.
  //
  // Dans cette simulation elle enregistre :
  //
  // PrimaryGeneratorAction
  //     ↓
  // génération des particules primaires
  //
  // RunAction
  //     ↓
  // gestion du run et du fichier ROOT
  //
  // EventAction
  //     ↓
  // stockage des données événementielles
  //
  // SteppingAction
  //     ↓
  // analyse des interactions gamma
  //
  // TrackingAction
  //     ↓
  // enregistrement des électrons créés
  virtual void Build() const;
  virtual void BuildForMaster() const;

private:

  // Pointeur vers la géométrie du détecteur.
  //
  // Permet de transmettre l'accès au MCP
  // aux actions qui en ont besoin.
  const DetectorConstruction* fDetector;
};

#endif
