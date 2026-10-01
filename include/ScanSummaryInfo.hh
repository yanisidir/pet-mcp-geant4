#ifndef SCAN_SUMMARY_INFO_HH
#define SCAN_SUMMARY_INFO_HH
#include "globals.hh"

// One row per run of a parameter scan.
struct ScanSummaryInfo
{
  G4int stacksPlus = 0;
  G4int stacksMinus = 0;
  G4int events = 0;
  G4int electronHits = 0;        // electrons entering a channel, both sides
  G4int electronHitsPlus = 0;
  G4int electronHitsMinus = 0;
  G4int eventsWithHitPlus = 0;
  G4int eventsWithHitMinus = 0;
  G4int coincidences = 0;
};
#endif
