#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

SteppingAction::SteppingAction(const DetectorConstruction* detector,
                               EventAction* eventAction)
  : G4UserSteppingAction(),
    fDetector(detector),
    fEventAction(eventAction)
{
}

SteppingAction::~SteppingAction()
{
}

// On enregistre les hits d'électrons dans les canaux MCP et les sorties de photons
// du dernier MCP de la pile +z ou -z. On enregistre également les interactions
// gamma dans les plaques MCP. On enregistre aussi les entrées de photons gamma dans les plaques MCP.

// This Geant4 stepping callback processes each step to record gamma interactions inside MCP plates, 
// gamma entries into MCP bodies, electron hits when an electron enters an MCP channel, and photon 
// exits when a gamma leaves the last MCP on the +z or -z stack. It also stops and kills
// the track after recording an electron-channel hit or a photon-exit event, 
// while safely returning early if the step or geometry information is invalid.
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!step || !fEventAction) {
    return;
  }

  G4Track* track = step->GetTrack();
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();

  if (!track || !preStepPoint || !postStepPoint) {
    return;
  }

  const G4ParticleDefinition* particle =
    track->GetParticleDefinition();
  if (!particle) {
    return;
  }

  const G4VPhysicalVolume* preVolume =
    preStepPoint->GetPhysicalVolume();
  const G4VPhysicalVolume* postVolume =
    postStepPoint->GetPhysicalVolume();

  const G4String particleName = particle->GetParticleName();

  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.0) {
    G4int edepSide = 0;
    G4int edepMcpIndex = -1;
    if (GetMcpInfo(preVolume, edepSide, edepMcpIndex)) {
      fEventAction->AddMcpEnergyDeposit(edepSide,
                                        edepMcpIndex,
                                        edep);
    }
  }

  // Une interaction gamma peut avoir lieu au milieu d'un volume.
  if (particleName == "gamma") {
    const G4VProcess* process =
      postStepPoint->GetProcessDefinedStep(); // Le processus qui a défini le pas de l'étape (peut être nul si le pas est dû à une limite géométrique)
    if (process) {
      const G4String processName = process->GetProcessName();
      const G4bool isPhysicalGammaInteraction =
        processName == "phot" ||
        processName == "compt" ||
        processName == "Rayl" ||
        processName == "conv";

      G4int mcpSide = 0;
      G4int mcpIndex = -1;
      if (isPhysicalGammaInteraction &&
          GetMcpInfo(preVolume, mcpSide, mcpIndex)) {
        fEventAction->SaveGammaInteraction(
          BuildGammaInteraction(step,
                                mcpSide,
                                mcpIndex,
                                processName));
      }
    }
  }

  if (postStepPoint->GetStepStatus() != fGeomBoundary) {
    return;
  }

  G4int enteredMcpSide = 0;
  G4int enteredMcpIndex = -1;
  const G4bool entersMcpBody =
    particleName == "gamma" &&
    GetMcpInfo(postVolume, enteredMcpSide, enteredMcpIndex) &&
    (!preVolume ||
     preVolume->GetLogicalVolume() !=
       postVolume->GetLogicalVolume());

  if (entersMcpBody) {
    fEventAction->SaveGammaMcpEntry(
      BuildGammaMcpEntry(step,
                         enteredMcpSide,
                         enteredMcpIndex));
  }

  G4int channelSide = 0;
  G4int channelMcpIndex = -1;
  const G4bool entersChannel =
    GetChannelInfo(postVolume, channelSide, channelMcpIndex);

  // L'électron est enregistré lorsqu'il pénètre dans un canal MCP
  if (particleName == "e-" &&
      entersChannel) {
    fEventAction->SaveElectronChannelHit(BuildElectronChannelHit(step, channelSide, channelMcpIndex));
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  const G4int lastMcpSide = GetLastMcpSide(preVolume);

  // Le gamma est enregistré lorsqu'il quitte le dernier MCP
  // de la pile +z ou -z.
  if (particleName == "gamma" && lastMcpSide != 0 && GetLastMcpSide(postVolume) != lastMcpSide) {
    fEventAction->SavePhotonExit(BuildPhotonExit(step, lastMcpSide));
    track->SetTrackStatus(fStopAndKill);
  }
}

// Identifie le côté et la plaque MCP contenant le volume physique.
G4bool SteppingAction::GetMcpInfo(
  const G4VPhysicalVolume* volume,
  G4int& side,
  G4int& mcpIndex) const
{
  if (!fDetector || !volume) {
    side = 0;
    mcpIndex = -1;
    return false;
  }

  return fDetector->GetMcpInfo(volume->GetLogicalVolume(),
                               side,
                               mcpIndex);
}

// Identifie le côté et la plaque MCP contenant le canal.
G4bool SteppingAction::GetChannelInfo(
  const G4VPhysicalVolume* volume,
  G4int& side,
  G4int& mcpIndex) const
{
  if (!fDetector || !volume) {
    side = 0;
    mcpIndex = -1;
    return false;
  }

  return fDetector->GetChannelInfo(volume->GetLogicalVolume(), side, mcpIndex);
}

// Retourne +1 pour un volume du dernier MCP côté +z, -1 pour un volume du dernier MCP côté -z, 0 sinon
G4int SteppingAction::GetLastMcpSide(
  const G4VPhysicalVolume* volume) const
{
  if (!fDetector || !volume) {
    return 0;
  }

  return fDetector->GetLastMcpSide(volume->GetLogicalVolume());
}

// Construit les informations d'un hit d'électron dans un canal MCP
ElectronChannelHitInfo SteppingAction::BuildElectronChannelHit(
  const G4Step* step,
  G4int side,
  G4int mcpIndex) const
{
  ElectronChannelHitInfo hit;
  const G4Track* track = step->GetTrack();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4Event* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();

  hit.eventID = event ? event->GetEventID() : -1;
  hit.trackID = track->GetTrackID();
  hit.parentID = track->GetParentID();
  hit.parentGammaTrackID = -1;
  hit.primaryGammaTrackID = -1;
  hit.hasValidCreationInfo = false;
  hit.side = side;
  hit.mcpIndex = mcpIndex;
  hit.kineticEnergy = postStepPoint->GetKineticEnergy();
  hit.globalTime = track->GetGlobalTime();
  if (hit.globalTime == 0.0) {
    // Le transport dans le champ peut remettre le temps courant à zéro.
    // Le temps mémorisé à la création fournit alors un repli cohérent.
    hit.globalTime =
      fEventAction->GetProducedElectronGlobalTime(hit.trackID);
  }
  hit.position = postStepPoint->GetPosition();
  hit.momentumDirection = postStepPoint->GetMomentumDirection();
  hit.creationPosition = track->GetVertexPosition();
  hit.creationTime = -1.0;
  hit.creationSide = 0;
  hit.creationMcpIndex = -1;

  const G4VPhysicalVolume* volume =
    postStepPoint->GetPhysicalVolume();
  hit.volumeName = volume ? volume->GetName() : "unknown";

  const G4VProcess* creatorProcess = track->GetCreatorProcess();
  hit.creatorProcessName =
    creatorProcess ? creatorProcess->GetProcessName() : "primary";

  return hit;
}

// Construit les informations d'un photon quittant le dernier MCP
PhotonExitInfo SteppingAction::BuildPhotonExit(
  const G4Step* step,
  G4int side) const
{
  PhotonExitInfo exitInfo;
  const G4Track* track = step->GetTrack();
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4Event* event =
    G4EventManager::GetEventManager()->GetConstCurrentEvent();

  exitInfo.eventID = event ? event->GetEventID() : -1;
  exitInfo.trackID = track->GetTrackID();
  exitInfo.parentID = track->GetParentID();
  exitInfo.side = side;
  exitInfo.kineticEnergy = postStepPoint->GetKineticEnergy();
  exitInfo.globalTime = postStepPoint->GetGlobalTime();
  exitInfo.position = postStepPoint->GetPosition();
  exitInfo.momentumDirection = postStepPoint->GetMomentumDirection();

  const G4VPhysicalVolume* volume =
    preStepPoint->GetPhysicalVolume();
  exitInfo.volumeName = volume ? volume->GetName() : "unknown";

  const G4VProcess* stepProcess =
    postStepPoint->GetProcessDefinedStep();
  exitInfo.stepProcessName =
    stepProcess ? stepProcess->GetProcessName() : "unknown";

  return exitInfo;
}

GammaInteractionInfo SteppingAction::BuildGammaInteraction(
  const G4Step* step,
  G4int side,
  G4int mcpIndex,
  const G4String& processName) const
{
  GammaInteractionInfo info;
  const G4Track* track = step->GetTrack();
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4Event* event =
    G4EventManager::GetEventManager()->GetConstCurrentEvent();

  info.eventID = event ? event->GetEventID() : -1;
  info.trackID = track->GetTrackID();
  info.parentID = track->GetParentID();
  info.side = side;
  info.mcpIndex = mcpIndex;
  info.processName = processName;

  // Energie incidente juste avant que le processus physique soit appliqué.
  info.kineticEnergy = preStepPoint->GetKineticEnergy();
  info.globalTime = postStepPoint->GetGlobalTime();
  info.position = postStepPoint->GetPosition();

  const G4VPhysicalVolume* volume =
    preStepPoint->GetPhysicalVolume();
  info.volumeName = volume ? volume->GetName() : "unknown";

  return info;
}

GammaMcpEntryInfo SteppingAction::BuildGammaMcpEntry(
  const G4Step* step,
  G4int side,
  G4int mcpIndex) const
{
  GammaMcpEntryInfo entry;
  const G4Track* track = step->GetTrack();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();
  const G4Event* event =
    G4EventManager::GetEventManager()->GetConstCurrentEvent();

  entry.eventID = event ? event->GetEventID() : -1;
  entry.trackID = track->GetTrackID();
  entry.parentID = track->GetParentID();
  entry.side = side;
  entry.mcpIndex = mcpIndex;
  entry.kineticEnergy = postStepPoint->GetKineticEnergy();
  entry.globalTime = postStepPoint->GetGlobalTime();
  entry.position = postStepPoint->GetPosition();

  const G4VPhysicalVolume* volume =
    postStepPoint->GetPhysicalVolume();
  entry.volumeName = volume ? volume->GetName() : "unknown";

  return entry;
}
