#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH
#include "EventSummaryInfo.hh"
#include "G4UserEventAction.hh"
class RunAction;
struct ElectronChannelHitInfo;

class EventAction : public G4UserEventAction
{
public:
  explicit EventAction(RunAction* runAction) : fRunAction(runAction) {}
  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;
  void RecordHit(const ElectronChannelHitInfo& hit);
private:
  RunAction* fRunAction;  // Owned by Geant4.
  EventSummaryInfo fSummary;
};
#endif
