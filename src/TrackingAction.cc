#include "TrackingAction.hh"

#include "EventAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"

TrackingAction::TrackingAction(EventAction* eventAction)
  : G4UserTrackingAction(),
    fEventAction(eventAction)
{
}

TrackingAction::~TrackingAction()
{
}

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  if (!track || !fEventAction || !track->GetParticleDefinition()) {
    return;
  }

  if (track->GetParticleDefinition()->GetParticleName() == "e-") {
    fEventAction->RecordProducedElectron(track->GetTrackID(),
                                         track->GetGlobalTime());
  }
}
