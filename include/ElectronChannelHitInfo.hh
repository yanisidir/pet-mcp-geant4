#ifndef ELECTRON_CHANNEL_HIT_INFO_HH
#define ELECTRON_CHANNEL_HIT_INFO_HH
#include "G4ThreeVector.hh"
#include "globals.hh"

// State at the first channel entry, in Geant4 internal units.
struct ElectronChannelHitInfo
{
  G4int eventID = 0;
  G4int trackID = 0;
  G4int side = 0;
  G4int stackIndex = 0;
  G4int channelCopyNumber = 0;
  G4double kineticEnergy = 0.;
  G4double globalTime = 0.;
  G4ThreeVector position;
  G4ThreeVector momentumDirection;
  G4String creatorProcessName;
};
#endif
