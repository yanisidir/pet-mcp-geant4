#ifndef PHOTON_EXIT_INFO_HH
#define PHOTON_EXIT_INFO_HH

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// Informations enregistrées lorsqu'un gamma quitte le MCP.
struct PhotonExitInfo
{
  G4int eventID;
  G4int trackID;
  G4int parentID;
  G4int side;
  G4double kineticEnergy;
  G4double globalTime;
  G4ThreeVector position;
  G4ThreeVector momentumDirection;
  G4String volumeName;
  G4String stepProcessName;
};

#endif
