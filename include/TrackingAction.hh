#ifndef TRACKING_ACTION_HH
#define TRACKING_ACTION_HH

#include "G4UserTrackingAction.hh"

class EventAction;
class G4Track;

// Compte chaque nouvelle track électron, sans modifier son transport.
class TrackingAction : public G4UserTrackingAction
{
public:
  explicit TrackingAction(EventAction* eventAction);
  virtual ~TrackingAction();

  virtual void PreUserTrackingAction(const G4Track* track);

private:
  EventAction* fEventAction;
};

#endif
