#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include "ElectronChannelHitInfo.hh"
#include "GammaInteractionInfo.hh"
#include "GammaMcpEntryInfo.hh"
#include "PhotonExitInfo.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <map>
#include <set>
#include <tuple>

class G4Event;
class DetectorConstruction;
class RunAction;

// EventAction déduplique les tracks et écrit le résumé de l'événement.
class EventAction : public G4UserEventAction
{
public:
  EventAction(const DetectorConstruction* detector,
              RunAction* runAction);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  // Enregistre une seule fois chaque track électron créée.
  void RecordProducedElectron(G4int trackID, G4double globalTime);

  // Retourne le temps global observé à la création de l'électron.
  G4double GetProducedElectronGlobalTime(G4int trackID) const;

  // Retourne true seulement lors du premier enregistrement de ce trackID.
  G4bool SaveElectronChannelHit(const ElectronChannelHitInfo& hit);
  G4bool SavePhotonExit(const PhotonExitInfo& exitInfo);
  void SaveGammaInteraction(const GammaInteractionInfo& info);
  G4bool SaveGammaMcpEntry(const GammaMcpEntryInfo& entry);

private:
  const DetectorConstruction* fDetector;
  RunAction* fRunAction;
  std::map<G4int, G4double> fProducedElectronGlobalTimes;
  std::set<G4int> fElectronChannelTrackIDs;
  std::set<G4int> fPhotonExitTrackIDs;
  // EventAction est recréé/réinitialisé par événement. La clé contient
  // trackID, side et mcpIndex; eventID est donc implicite.
  std::set<std::tuple<G4int, G4int, G4int> > fGammaMcpEntryKeys;
  std::map<std::pair<G4int, G4int>, G4int>
    fElectronChannelCountsByMcp;
  G4int fElectronChannelPlusCount;
  G4int fElectronChannelMinusCount;
  G4int fPhotonExitPlusCount;
  G4int fPhotonExitMinusCount;
  G4int fGammaInteractionCount;
  G4int fGammaPhotCount;
  G4int fGammaComptCount;
  G4int fGammaRaylCount;
  G4int fGammaConvCount;
};

#endif
