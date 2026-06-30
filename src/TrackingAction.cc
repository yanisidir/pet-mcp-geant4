#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "ElectronCreationInfo.hh"
#include "EventAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

TrackingAction::TrackingAction(const DetectorConstruction* detector,
                               EventAction* eventAction)
  : G4UserTrackingAction(),
    fDetector(detector),
    fEventAction(eventAction)
{
}

TrackingAction::~TrackingAction()
{
}

// On enregistre les tracks gamma (trackID et parentID).
// On enregistre également tous les électrons produits (trackID et temps de création).
// Pour les électrons issus directement d'une interaction gamma
// (photoélectrique, Compton ou conversion) dans une plaque MCP,
// on enregistre des informations supplémentaires : parent gamma,
// gamma primaire, position et temps de création, processus créateur
// et plaque MCP de naissance.
void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  if (!track || !fEventAction || !track->GetParticleDefinition()) {
    return;
  }

  const G4String particleName =
    track->GetParticleDefinition()->GetParticleName();

  if (particleName == "gamma") {
    fEventAction->RecordGammaTrack(track->GetTrackID(),
                                   track->GetParentID());
  }

  if (particleName == "e-") {
    fEventAction->RecordProducedElectron(track->GetTrackID(),
                                         track->GetGlobalTime());

    const G4VProcess* creatorProcess = track->GetCreatorProcess();
    if (!creatorProcess || !fDetector) {
      return;
    }

    const G4String creatorProcessName =
      creatorProcess->GetProcessName();
    const G4bool comesFromGammaInteraction =
      creatorProcessName == "phot" ||
      creatorProcessName == "compt" ||
      creatorProcessName == "conv";

    if (!comesFromGammaInteraction) {
      return;
    }

    const G4VPhysicalVolume* creationVolume = track->GetVolume();
    if (!creationVolume) {
      return;
    }

    G4int creationSide = 0;
    G4int creationMcpIndex = -1;
    if (!fDetector->GetMcpInfo(creationVolume->GetLogicalVolume(),
                               creationSide,
                               creationMcpIndex)) {
      return;
    }

    const G4Event* event =
      G4EventManager::GetEventManager()->GetConstCurrentEvent();

    ElectronCreationInfo info;
    info.eventID = event ? event->GetEventID() : -1;
    info.electronTrackID = track->GetTrackID();
    info.parentGammaTrackID = track->GetParentID();
    info.primaryGammaTrackID =
      fEventAction->GetPrimaryGammaTrackID(track->GetParentID());
    info.creatorProcessName = creatorProcessName;
    info.creationPosition = track->GetVertexPosition();
    info.creationTime = track->GetGlobalTime();
    info.creationSide = creationSide;
    info.creationMcpIndex = creationMcpIndex;
    fEventAction->RecordElectronCreation(info);
  }
}
