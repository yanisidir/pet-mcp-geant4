#include "EventAction.hh"

#include "EventSummaryInfo.hh"
#include "RootOutput.hh"
#include "RunAction.hh"

#include "G4Event.hh"

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(),
    fRunAction(runAction),
    fElectronChannelPlusCount(0),
    fElectronChannelMinusCount(0),
    fPhotonExitPlusCount(0),
    fPhotonExitMinusCount(0)
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fProducedElectronGlobalTimes.clear();
  fElectronChannelTrackIDs.clear();
  fPhotonExitTrackIDs.clear();
  fElectronChannelCountsByMcp.clear();
  fElectronChannelPlusCount = 0;
  fElectronChannelMinusCount = 0;
  fPhotonExitPlusCount = 0;
  fPhotonExitMinusCount = 0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  EventSummaryInfo summary;
  summary.eventID = event ? event->GetEventID() : -1;
  summary.electronProducedCount =
    static_cast<G4int>(fProducedElectronGlobalTimes.size());
  summary.electronChannelCount =
    static_cast<G4int>(fElectronChannelTrackIDs.size());
  summary.electronChannelPlusCount = fElectronChannelPlusCount;
  summary.electronChannelMinusCount = fElectronChannelMinusCount;
  summary.electronChannelPlusMcp0Count =
    fElectronChannelCountsByMcp[std::make_pair(1, 0)];
  summary.electronChannelPlusMcp1Count =
    fElectronChannelCountsByMcp[std::make_pair(1, 1)];
  summary.electronChannelPlusMcp2Count =
    fElectronChannelCountsByMcp[std::make_pair(1, 2)];
  summary.electronChannelMinusMcp0Count =
    fElectronChannelCountsByMcp[std::make_pair(-1, 0)];
  summary.electronChannelMinusMcp1Count =
    fElectronChannelCountsByMcp[std::make_pair(-1, 1)];
  summary.electronChannelMinusMcp2Count =
    fElectronChannelCountsByMcp[std::make_pair(-1, 2)];
  summary.hasElectronChannelPlus = fElectronChannelPlusCount > 0;
  summary.hasElectronChannelMinus = fElectronChannelMinusCount > 0;
  summary.isCoincidence =
    summary.hasElectronChannelPlus &&
    summary.hasElectronChannelMinus;
  summary.photonExitCount =
    static_cast<G4int>(fPhotonExitTrackIDs.size());
  summary.photonExitPlusCount = fPhotonExitPlusCount;
  summary.photonExitMinusCount = fPhotonExitMinusCount;

  RootOutput::Instance()->FillEventSummary(summary);
  if (fRunAction) {
    fRunAction->CountEvent(summary.hasElectronChannelPlus,
                           summary.hasElectronChannelMinus);
  }

  G4cout << "-- Event " << summary.eventID
         << " ends with " << summary.electronProducedCount
         << " produced electrons, "
         << summary.electronChannelCount
         << " channel hits (+z=" << summary.electronChannelPlusCount
         << ", -z=" << summary.electronChannelMinusCount << ") and "
         << summary.photonExitCount
         << " photon exits (+z=" << summary.photonExitPlusCount
         << ", -z=" << summary.photonExitMinusCount << ")"
         << " | coincidence="
         << (summary.isCoincidence ? "yes" : "no")
         << "." << G4endl;
}

void EventAction::RecordProducedElectron(G4int trackID,
                                         G4double globalTime)
{
  fProducedElectronGlobalTimes.insert(
    std::make_pair(trackID, globalTime));
}

G4double EventAction::GetProducedElectronGlobalTime(G4int trackID) const
{
  const std::map<G4int, G4double>::const_iterator found =
    fProducedElectronGlobalTimes.find(trackID);

  return found != fProducedElectronGlobalTimes.end()
    ? found->second
    : 0.0;
}

G4bool EventAction::SaveElectronChannelHit(
  const ElectronChannelHitInfo& hit)
{
  const std::pair<std::set<G4int>::iterator, bool> inserted =
    fElectronChannelTrackIDs.insert(hit.trackID);

  if (!inserted.second) {
    return false;
  }

  if (hit.side > 0) {
    ++fElectronChannelPlusCount;
  } else if (hit.side < 0) {
    ++fElectronChannelMinusCount;
  }
  ++fElectronChannelCountsByMcp[
    std::make_pair(hit.side, hit.mcpIndex)];

  RootOutput::Instance()->FillElectronChannelHit(hit);
  return true;
}

G4bool EventAction::SavePhotonExit(const PhotonExitInfo& exitInfo)
{
  const std::pair<std::set<G4int>::iterator, bool> inserted =
    fPhotonExitTrackIDs.insert(exitInfo.trackID);

  if (!inserted.second) {
    return false;
  }

  if (exitInfo.side > 0) {
    ++fPhotonExitPlusCount;
  } else if (exitInfo.side < 0) {
    ++fPhotonExitMinusCount;
  }

  RootOutput::Instance()->FillPhotonExit(exitInfo);
  return true;
}
