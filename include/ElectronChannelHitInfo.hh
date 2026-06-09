#ifndef ELECTRON_CHANNEL_HIT_INFO_HH
#define ELECTRON_CHANNEL_HIT_INFO_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// Informations enregistrées lorsqu'un électron entre dans MCP_channel.
struct ElectronChannelHitInfo
{
  G4int eventID;
  G4int trackID;
  G4int parentID;
  G4int side;
  G4int mcpIndex;
  G4double kineticEnergy;
  G4double globalTime;
  G4ThreeVector position;
  G4ThreeVector momentumDirection;
  G4String volumeName;
  G4String creatorProcessName;
};

#endif
