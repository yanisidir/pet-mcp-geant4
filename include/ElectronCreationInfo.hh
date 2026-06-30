#ifndef ELECTRON_CREATION_INFO_HH
#define ELECTRON_CREATION_INFO_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// Informations exactes à la naissance d'un électron secondaire
// produit par une interaction gamma dans le verre MCP.
struct ElectronCreationInfo
{
  G4int eventID;
  G4int electronTrackID;
  G4int parentGammaTrackID;
  G4int primaryGammaTrackID;
  G4String creatorProcessName;
  G4ThreeVector creationPosition;
  G4double creationTime;
  G4int creationSide;
  G4int creationMcpIndex;
};

#endif
