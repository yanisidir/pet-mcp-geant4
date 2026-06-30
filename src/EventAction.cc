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
    fPhotonExitMinusCount(0),
    fGammaInteractionCount(0),
    fGammaPhotCount(0),
    fGammaComptCount(0),
    fGammaRaylCount(0),
    fGammaConvCount(0),
    fEdepTotal(0.0),
    fEdepPlusTotal(0.0),
    fEdepMinusTotal(0.0)
{
}

EventAction::~EventAction()
{
}

// On commence un nouvel événement. On réinitialise les compteurs et les structures de données.
void EventAction::BeginOfEventAction(const G4Event*)
{
  fProducedElectronGlobalTimes.clear(); // Réinitialise la map des temps globaux des électrons produits
  fPrimaryGammaTrackIDs.clear(); // Réinitialise la map des trackIDs des gammas primaires
  fElectronCreations.clear(); // Réinitialise la map des informations de création des électrons
  fElectronChannelTrackIDs.clear(); // Réinitialise l'ensemble des trackIDs des hits d'électrons dans les canaux MCP
  fPhotonExitTrackIDs.clear(); // Réinitialise l'ensemble des trackIDs des sorties de photons du dernier MCP
  fGammaMcpEntryKeys.clear(); // Réinitialise l'ensemble des clés (trackID, side, mcpIndex) pour les entrées de photons gamma dans les plaques MCP
  fElectronChannelCountsByMcp.clear(); // Réinitialise la map des compteurs de hits d'électrons par plaque MCP
  fElectronChannelPlusCount = 0; // Réinitialise le compteur de hits d'électrons dans les canaux MCP côté +z
  fElectronChannelMinusCount = 0; // Réinitialise le compteur de hits d'électrons dans les canaux MCP côté -z
  fPhotonExitPlusCount = 0; // Réinitialise le compteur de sorties de photons du dernier MCP côté +z
  fPhotonExitMinusCount = 0; // Réinitialise le compteur de sorties de photons du dernier MCP côté -z
  fGammaInteractionCount = 0; // Réinitialise le compteur d'interactions gamma dans les plaques MCP
  fGammaPhotCount = 0; // Réinitialise le compteur d'interactions photoélectriques gamma
  fGammaComptCount = 0; // Réinitialise le compteur d'interactions Compton gamma
  fGammaRaylCount = 0; // Réinitialise le compteur d'interactions Rayleigh gamma
  fGammaConvCount = 0; // Réinitialise le compteur d'interactions de conversion gamma
  fEdepTotal = 0.0; // Réinitialise l'énergie déposée totale dans les MCP
  fEdepPlusTotal = 0.0; // Réinitialise l'énergie déposée côté +z
  fEdepMinusTotal = 0.0; // Réinitialise l'énergie déposée côté -z
  fEdepByMcp.clear(); // Réinitialise l'énergie déposée par plaque MCP
}

// On termine l'événement. 
// On enregistre les informations de résumé de l'événement et les statistiques des plaques MCP dans le fichier ROOT.
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
  summary.gammaInteractionCount = fGammaInteractionCount;
  summary.gammaPhotCount = fGammaPhotCount;
  summary.gammaComptCount = fGammaComptCount;
  summary.gammaRaylCount = fGammaRaylCount;
  summary.gammaConvCount = fGammaConvCount;
  summary.edepTotal = fEdepTotal;
  summary.edepPlusTotal = fEdepPlusTotal;
  summary.edepMinusTotal = fEdepMinusTotal;

  for (G4int i = 0; i < kMaxMcpPlatesInEventSummary; ++i) {
    summary.edepPlusByMcp[i] = 0.0;
    summary.edepMinusByMcp[i] = 0.0;
  }

  const G4int numberOfMCPs =
    fDetector ? fDetector->GetNumberOfMCPs() : 0;
  const G4int plateLimit =
    numberOfMCPs < kMaxMcpPlatesInEventSummary
      ? numberOfMCPs
      : kMaxMcpPlatesInEventSummary;
  for (G4int mcpIndex = 0; mcpIndex < plateLimit; ++mcpIndex) {
    const std::map<std::pair<G4int, G4int>, G4double>::const_iterator
      plusFound = fEdepByMcp.find(std::make_pair(1, mcpIndex));
    const std::map<std::pair<G4int, G4int>, G4double>::const_iterator
      minusFound = fEdepByMcp.find(std::make_pair(-1, mcpIndex));

    if (plusFound != fEdepByMcp.end()) {
      summary.edepPlusByMcp[mcpIndex] = plusFound->second;
    }
    if (minusFound != fEdepByMcp.end()) {
      summary.edepMinusByMcp[mcpIndex] = minusFound->second;
    }
  }

  RootOutput::Instance()->FillEventSummary(summary);

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
         << " and " << summary.gammaInteractionCount
         << " gamma interactions"
         << " | coincidence="
         << (summary.isCoincidence ? "yes" : "no")
         << "." << G4endl;
}

// Enregistre le temps global d'un électron produit, identifié par son trackID.
void EventAction::RecordProducedElectron(G4int trackID,
                                         G4double globalTime)
{
  fProducedElectronGlobalTimes.insert(
    std::make_pair(trackID, globalTime));
}

// Enregistre le trackID du gamma primaire associé à un trackID de gamma donné.
void EventAction::RecordGammaTrack(G4int trackID, G4int parentID)
{
  G4int primaryGammaTrackID = -1;
  if (parentID == 0) {
    primaryGammaTrackID = trackID;
  } else {
    const std::map<G4int, G4int>::const_iterator parentFound =
      fPrimaryGammaTrackIDs.find(parentID);
    if (parentFound != fPrimaryGammaTrackIDs.end()) {
      primaryGammaTrackID = parentFound->second;
    }
  }

  fPrimaryGammaTrackIDs[trackID] = primaryGammaTrackID;
}

// Enregistre les informations de création d'un électron produit, identifié par son trackID.
void EventAction::RecordElectronCreation(
  const ElectronCreationInfo& info)
{
  fElectronCreations[info.electronTrackID] = info;
}

// Retourne le trackID du gamma primaire associé à un trackID de gamma donné, ou -1 si non trouvé.
G4int EventAction::GetPrimaryGammaTrackID(G4int gammaTrackID) const
{
  const std::map<G4int, G4int>::const_iterator found =
    fPrimaryGammaTrackIDs.find(gammaTrackID);

  return found != fPrimaryGammaTrackIDs.end()
    ? found->second
    : -1;
}

// Retourne le temps global observé à la création de l'électron identifié par son trackID, ou 0.0 si non trouvé.
G4double EventAction::GetProducedElectronGlobalTime(G4int trackID) const
{
  const std::map<G4int, G4double>::const_iterator found =
    fProducedElectronGlobalTimes.find(trackID);

  return found != fProducedElectronGlobalTimes.end()
    ? found->second
    : 0.0;
}

// Retourne true seulement lors du premier enregistrement d'un hit d'électron dans un canal MCP identifié par son trackID.
G4bool EventAction::SaveElectronChannelHit(
  const ElectronChannelHitInfo& hit)
{
  ElectronChannelHitInfo enrichedHit(hit);
  const std::map<G4int, ElectronCreationInfo>::const_iterator
    creationFound = fElectronCreations.find(hit.trackID);

  // Si l'électron a été créé par une interaction gamma dans une plaque MCP,
  // on enrichit les informations du hit avec les informations de création.
  if (creationFound != fElectronCreations.end()) {
    enrichedHit.hasValidCreationInfo = true;
    enrichedHit.parentGammaTrackID = creationFound->second.parentGammaTrackID;
    enrichedHit.primaryGammaTrackID = creationFound->second.primaryGammaTrackID;
    enrichedHit.creationPosition = creationFound->second.creationPosition;
    enrichedHit.creationTime = creationFound->second.creationTime;
    enrichedHit.creationSide = creationFound->second.creationSide;
    enrichedHit.creationMcpIndex = creationFound->second.creationMcpIndex;
    enrichedHit.creatorProcessName = creationFound->second.creatorProcessName;
  }

  // On ne sauvegarde qu'une seule fois chaque trackID d'électron dans un canal MCP.
  const std::pair<std::set<G4int>::iterator, bool> inserted =
    fElectronChannelTrackIDs.insert(enrichedHit.trackID);

  // Si le trackID a déjà été enregistré, on ne fait rien et on retourne false.
  if (!inserted.second) {
    return false;
  }

  // On met à jour les compteurs de hits d'électrons par côté et par plaque MCP.
  if (enrichedHit.side > 0) {
    ++fElectronChannelPlusCount;
  } else if (enrichedHit.side < 0) {
    ++fElectronChannelMinusCount;
  }
  ++fElectronChannelCountsByMcp[std::make_pair(enrichedHit.side, enrichedHit.mcpIndex)];

  RootOutput::Instance()->FillElectronChannelHit(enrichedHit);
  return true;
}
// This function records a photon exit event for the current EventAction. It adds exitInfo.trackID to fPhotonExitTrackIDs 
// and returns false if that ID was already recorded; otherwise it updates the 
// plus/minus exit counters based on exitInfo.side and writes the event to RootOutput::Instance()->FillPhotonExit(exitInfo), then returns true.
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

void EventAction::SaveGammaInteraction(
  const GammaInteractionInfo& info)
{
  ++fGammaInteractionCount;

  if (info.processName == "phot") {
    ++fGammaPhotCount;
  } else if (info.processName == "compt") {
    ++fGammaComptCount;
  } else if (info.processName == "Rayl") {
    ++fGammaRaylCount;
  } else if (info.processName == "conv") {
    ++fGammaConvCount;
  }

  RootOutput::Instance()->FillGammaInteraction(info);
}

G4bool EventAction::SaveGammaMcpEntry(
  const GammaMcpEntryInfo& entry)
{
  const std::tuple<G4int, G4int, G4int> key(
    entry.trackID,
    entry.side,
    entry.mcpIndex);
  const std::pair<
    std::set<std::tuple<G4int, G4int, G4int> >::iterator,
    bool> inserted = fGammaMcpEntryKeys.insert(key);

  if (!inserted.second) {
    return false;
  }

  RootOutput::Instance()->FillGammaMcpEntry(entry);
  return true;
}

void EventAction::AddMcpEnergyDeposit(G4int side,
                                      G4int mcpIndex,
                                      G4double edep)
{
  if (edep <= 0.0 || mcpIndex < 0) {
    return;
  }

  fEdepTotal += edep;

  if (side > 0) {
    fEdepPlusTotal += edep;
  } else if (side < 0) {
    fEdepMinusTotal += edep;
  }

  fEdepByMcp[std::make_pair(side, mcpIndex)] += edep;
}
