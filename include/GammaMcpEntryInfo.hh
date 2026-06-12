#ifndef GAMMA_MCP_ENTRY_INFO_HH
#define GAMMA_MCP_ENTRY_INFO_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// Informations enregistrées lorsqu'un gamma entre dans un corps MCP.
struct GammaMcpEntryInfo
{
  G4int eventID;
  G4int trackID;
  G4int parentID;
  G4int side;
  G4int mcpIndex;
  G4double kineticEnergy;
  G4double globalTime;
  G4ThreeVector position;
  G4String volumeName;
};

#endif
