#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

ActionInitialization::ActionInitialization(const DetectorConstruction* detector)
  : G4VUserActionInitialization(),
    fDetector(detector)
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction);
  RunAction* runAction = new RunAction;
  SetUserAction(runAction);

  EventAction* eventAction = new EventAction(runAction);
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(fDetector, eventAction));
  SetUserAction(new TrackingAction(eventAction));
}

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}
