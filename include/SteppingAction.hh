#ifndef STEPPING_ACTION_HH
#define STEPPING_ACTION_HH

#include "ElectronChannelHitInfo.hh"
#include "PhotonExitInfo.hh"

#include "G4UserSteppingAction.hh"

class DetectorConstruction;
class EventAction;
class G4Step;
class G4VPhysicalVolume;

// SteppingAction observe uniquement deux passages de frontière utiles.
class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(const DetectorConstruction* detector,
                 EventAction* eventAction);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step* step);

private:
  G4bool GetChannelInfo(const G4VPhysicalVolume* volume,
                        G4int& side,
                        G4int& mcpIndex) const;
  G4int GetLastMcpSide(const G4VPhysicalVolume* volume) const;

  ElectronChannelHitInfo BuildElectronChannelHit(
    const G4Step* step,
    G4int side,
    G4int mcpIndex) const;

  PhotonExitInfo BuildPhotonExit(const G4Step* step,
                                 G4int side) const;

  const DetectorConstruction* fDetector;
  EventAction* fEventAction;
};

#endif
