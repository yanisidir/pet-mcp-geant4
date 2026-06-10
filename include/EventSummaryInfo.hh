#ifndef EVENT_SUMMARY_INFO_HH
#define EVENT_SUMMARY_INFO_HH

#include "globals.hh"

// Une ligne minimale par événement.
struct EventSummaryInfo
{
  G4int eventID;
  G4int electronProducedCount;
  G4int electronChannelCount;
  G4int electronChannelPlusCount;
  G4int electronChannelMinusCount;
  G4bool hasElectronChannelPlus;
  G4bool hasElectronChannelMinus;
  G4bool isCoincidence;
  G4int photonExitCount;
  G4int photonExitPlusCount;
  G4int photonExitMinusCount;
};

#endif
