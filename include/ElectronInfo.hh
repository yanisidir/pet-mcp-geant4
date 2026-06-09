#ifndef ELECTRON_INFO_HH
#define ELECTRON_INFO_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// Structure décrivant les propriétés d'un électron observé pendant un événement.
struct ElectronInfo
{
  // Numéro de l'événement Geant4.
  G4int eventID;

  // Identifiant unique de la track électron.
  G4int trackID;

  // Identifiant de la particule mère.
  // Pour un électron créé par un gamma primaire,
  // parentID correspond généralement au trackID du gamma.
  G4int parentID;

  // Nom du volume dans lequel l'électron a été créé.
  // Exemples :
  // "MCP_body"
  // "MCP_channel"
  G4String volumeName;

  // Processus associé au step observé.
  // Actuellement souvent fixé à "track_start".
  G4String stepProcessName;

  // Nom du processus ayant créé l'électron.
  // Exemples :
  // "phot"
  // "compt"
  // "conv"
  G4String creatorProcessName;

  // Identifiant numérique du processus créateur.
  // 1 = photoélectrique
  // 2 = Compton
  // 3 = conversion de paire
  // 4 = Rayleigh
  // 0 = autre
  G4int creatorProcessID;

  // Energie cinétique initiale de l'électron au moment
  // de sa création.
  G4double kineticEnergy;

  // Temps global de création de l'électron.
  G4double globalTime;

  // Position de naissance de l'électron.
  // Correspond généralement à la position de l'interaction
  // gamma ayant produit cet électron.
  G4ThreeVector position;

  // Direction initiale de l'électron à sa création.
  // Cette information est particulièrement importante
  // pour un éventuel transport ultérieur dans CST.
  G4ThreeVector direction;
};

#endif