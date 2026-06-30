#ifndef TRACKING_ACTION_HH
#define TRACKING_ACTION_HH

#include "G4UserTrackingAction.hh"

class DetectorConstruction;
class EventAction;
class G4Track;

// Compte chaque nouvelle track électron, sans modifier son transport.
class TrackingAction : public G4UserTrackingAction
{
public:
  TrackingAction(const DetectorConstruction* detector,
                 EventAction* eventAction);
  virtual ~TrackingAction();

  virtual void PreUserTrackingAction(const G4Track* track);

private:
  const DetectorConstruction* fDetector;
  EventAction* fEventAction;
};

#endif
