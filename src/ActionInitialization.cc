#include "ActionInitialization.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const DetectorConstruction* detector)
  : fDetector(detector),
    fGenerator(new PrimaryGeneratorAction),
    fRunAction(new RunAction(detector))
{
}

void ActionInitialization::Build() const
{
  // Geant4 takes ownership of every action registered here.
  SetUserAction(fGenerator);
  SetUserAction(fRunAction);
  auto* eventAction = new EventAction(fRunAction);
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(fDetector, eventAction));
}
