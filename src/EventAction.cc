#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "EventSummaryInfo.hh"
#include "McpPlateStatsInfo.hh"
#include "RootOutput.hh"
#include "RunAction.hh"

#include "G4Event.hh"

EventAction::EventAction(const DetectorConstruction* detector,
                         RunAction* runAction)
  : G4UserEventAction(),
    fDetector(detector),
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

  const G4int numberOfMCPs =
    fDetector ? fDetector->GetNumberOfMCPs() : 0;
  const G4int sides[2] = {1, -1};
  for (G4int sideIndex = 0; sideIndex < 2; ++sideIndex) {
    const G4int side = sides[sideIndex];
    for (G4int mcpIndex = 0;
         mcpIndex < numberOfMCPs;
         ++mcpIndex) {
      const std::pair<G4int, G4int> key(side, mcpIndex); // pair of (side, mcpIndex) to look up the electron channel count for this plate
      const std::map<std::pair<G4int, G4int>, G4int>::const_iterator // map from (side, mcpIndex) to electron channel count
        found = fElectronChannelCountsByMcp.find(key);

      McpPlateStatsInfo plateStats;
      plateStats.eventID = summary.eventID;
      plateStats.side = side;
      plateStats.mcpIndex = mcpIndex;
      plateStats.electronChannelCount =
        found != fElectronChannelCountsByMcp.end()
          ? found->second
          : 0;
      RootOutput::Instance()->FillMcpPlateStats(plateStats);
    }
  }

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
