#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "ElectronChannelHitInfo.hh"
#include "EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  const auto* post = step->GetPostStepPoint();
  if (post->GetStepStatus() != fGeomBoundary) return;

  auto* track = step->GetTrack();
  const auto* preVolume = step->GetPreStepPoint()->GetPhysicalVolume();
  const auto* postVolume = post->GetPhysicalVolume();
  const auto* preLogical = preVolume ? preVolume->GetLogicalVolume() : nullptr;
  const auto* postLogical = postVolume ? postVolume->GetLogicalVolume() : nullptr;
  const auto& particle = track->GetParticleDefinition()->GetParticleName();
  const auto side = fDetector->GetChannelSide(postLogical);

  if (particle == "e-" && side != 0) {
    ElectronChannelHitInfo hit;
    hit.eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
    hit.trackID = track->GetTrackID();
    hit.side = side;
    hit.stackIndex = fDetector->GetStackIndex(postLogical);
    hit.channelCopyNumber = postVolume->GetCopyNo();
    hit.kineticEnergy = post->GetKineticEnergy();
    hit.globalTime = track->GetGlobalTime();
    hit.position = post->GetPosition();
    hit.momentumDirection = post->GetMomentumDirection();
    const auto* creator = track->GetCreatorProcess();
    hit.creatorProcessName = creator ? creator->GetProcessName() : "primary";
    fEventAction->RecordHit(hit);

    // This simulation ends at channel entry; it does not model MCP multiplication.
    // A killed track cannot contribute a second hit, so no deduplication is needed.
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  // Preserve the existing transport boundary: a gamma leaving its MCP is stopped.
  const auto preSide = fDetector->GetMcpSide(preLogical);
  if (fDetector->KillGammaOnExit() && particle == "gamma" && preSide != 0 &&
      fDetector->GetMcpSide(postLogical) != preSide) {
    track->SetTrackStatus(fStopAndKill);
  }
}
