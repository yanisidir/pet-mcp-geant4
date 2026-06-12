#ifndef GAMMA_INTERACTION_INFO_HH
#define GAMMA_INTERACTION_INFO_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// Une ligne ROOT pour une interaction physique gamma dans un corps MCP.
struct GammaInteractionInfo
{
  G4int eventID;
  G4int trackID;
  G4int parentID;
  G4int side;
  G4int mcpIndex;
  G4String processName;
  G4double kineticEnergy;
  G4double globalTime;
  G4ThreeVector position;
  G4String volumeName;
};

#endif
