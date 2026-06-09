#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include "ElectronChannelHitInfo.hh"
#include "PhotonExitInfo.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <map>
#include <set>

class G4Event;
class RunAction;

// EventAction déduplique les tracks et écrit le résumé de l'événement.
class EventAction : public G4UserEventAction
{
public:
  explicit EventAction(RunAction* runAction);
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

private:
  RunAction* fRunAction;
  std::map<G4int, G4double> fProducedElectronGlobalTimes;
  std::set<G4int> fElectronChannelTrackIDs;
  std::set<G4int> fPhotonExitTrackIDs;
  std::map<std::pair<G4int, G4int>, G4int>
    fElectronChannelCountsByMcp;
  G4int fElectronChannelPlusCount;
  G4int fElectronChannelMinusCount;
  G4int fPhotonExitPlusCount;
  G4int fPhotonExitMinusCount;
};

#endif
