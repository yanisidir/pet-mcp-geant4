#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH

#include "G4Run.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"

// G4Run spécialisé qui fusionne proprement les compteurs des workers.
class CoincidenceRun : public G4Run
{
public:
  CoincidenceRun();
  virtual ~CoincidenceRun();

  void CountEvent(G4bool hasPlusHit, G4bool hasMinusHit);
  virtual void Merge(const G4Run* run);

  G4int GetProcessedEventCount() const;
  G4int GetPlusDetectedEventCount() const;
  G4int GetMinusDetectedEventCount() const;
  G4int GetCoincidenceEventCount() const;

private:
  G4int fProcessedEventCount;
  G4int fPlusDetectedEventCount;
  G4int fMinusDetectedEventCount;
  G4int fCoincidenceEventCount;
};

class RunAction : public G4UserRunAction
{
public:

  // Constructeur.
  RunAction();

  virtual ~RunAction();

  // Fonction appelée automatiquement au début du run.
  //
  // Utilisée notamment pour :
  // - créer le fichier ROOT ;
  // - créer les TTrees ;
  // - initialiser les sorties.
  virtual void BeginOfRunAction(const G4Run* run);

  // Fonction appelée automatiquement à la fin du run.
  //
  // Utilisée notamment pour :
  // - écrire les données ROOT ;
  // - fermer le fichier ROOT ;
  // - afficher un résumé du run.
  virtual void EndOfRunAction(const G4Run* run);
  virtual G4Run* GenerateRun();

  // Transmet le résultat de l'événement au CoincidenceRun du worker.
  void CountEvent(G4bool hasPlusHit, G4bool hasMinusHit);
};

#endif
