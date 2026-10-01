#ifndef STEPPING_ACTION_HH
#define STEPPING_ACTION_HH
#include "G4UserSteppingAction.hh"
class DetectorConstruction;
class EventAction;

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(const DetectorConstruction* detector, EventAction* eventAction)
    : fDetector(detector), fEventAction(eventAction) {}
  void UserSteppingAction(const G4Step* step) override;
private:
  const DetectorConstruction* fDetector;
  EventAction* fEventAction;
};
#endif
