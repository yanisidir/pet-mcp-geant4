#include "EventAction.hh"
#include "ElectronChannelHitInfo.hh"
#include "RunAction.hh"
#include "G4Event.hh"

void EventAction::BeginOfEventAction(const G4Event* event)
{
  fSummary = {};
  fSummary.eventID = event->GetEventID();
}

void EventAction::RecordHit(const ElectronChannelHitInfo& hit)
{
  if (hit.side > 0) ++fSummary.electronChannelPlusCount;
  else ++fSummary.electronChannelMinusCount;
  fRunAction->GetOutput().FillElectronChannelHit(hit);
}

void EventAction::EndOfEventAction(const G4Event*)
{
  fSummary.isCoincidence = fSummary.electronChannelPlusCount > 0 &&
                           fSummary.electronChannelMinusCount > 0;
  fRunAction->RecordEvent(fSummary);
}
