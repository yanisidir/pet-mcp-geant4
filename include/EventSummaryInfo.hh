#ifndef EVENT_SUMMARY_INFO_HH
#define EVENT_SUMMARY_INFO_HH
#include "globals.hh"

// One row per event, including events with no hit (efficiency denominator).
struct EventSummaryInfo
{
  G4int eventID = 0;
  G4int electronChannelPlusCount = 0;
  G4int electronChannelMinusCount = 0;
  G4bool isCoincidence = false;
};
#endif
