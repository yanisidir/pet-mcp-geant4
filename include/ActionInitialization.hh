#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH
#include "G4VUserActionInitialization.hh"
class DetectorConstruction;
class PrimaryGeneratorAction;
class RunAction;

// Register the actions of this serial application. The generator and the run
// action are built by the constructor, not by Build(), so that their
// /source/ and /output/ commands exist before /run/initialize.
class ActionInitialization : public G4VUserActionInitialization
{
public:
  explicit ActionInitialization(const DetectorConstruction* detector);
  void Build() const override;
private:
  const DetectorConstruction* fDetector;
  PrimaryGeneratorAction* fGenerator;
  RunAction* fRunAction;
};
#endif
