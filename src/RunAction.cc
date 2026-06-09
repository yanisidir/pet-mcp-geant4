#include "RunAction.hh"

#include "RootOutput.hh"

#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4ios.hh"

CoincidenceRun::CoincidenceRun()
  : G4Run(),
    fProcessedEventCount(0),
    fPlusDetectedEventCount(0),
    fMinusDetectedEventCount(0),
    fCoincidenceEventCount(0)
{
}

CoincidenceRun::~CoincidenceRun()
{
}

void CoincidenceRun::CountEvent(G4bool hasPlusHit,
                                G4bool hasMinusHit)
{
  ++fProcessedEventCount;
  if (hasPlusHit) {
    ++fPlusDetectedEventCount;
  }
  if (hasMinusHit) {
    ++fMinusDetectedEventCount;
  }
  if (hasPlusHit && hasMinusHit) {
    ++fCoincidenceEventCount;
  }
}

void CoincidenceRun::Merge(const G4Run* run)
{
  const CoincidenceRun* localRun =
    dynamic_cast<const CoincidenceRun*>(run);
  if (localRun) {
    fProcessedEventCount += localRun->fProcessedEventCount;
    fPlusDetectedEventCount +=
      localRun->fPlusDetectedEventCount;
    fMinusDetectedEventCount +=
      localRun->fMinusDetectedEventCount;
    fCoincidenceEventCount +=
      localRun->fCoincidenceEventCount;
  }

  G4Run::Merge(run);
}

G4int CoincidenceRun::GetProcessedEventCount() const
{
  return fProcessedEventCount;
}

G4int CoincidenceRun::GetPlusDetectedEventCount() const
{
  return fPlusDetectedEventCount;
}

G4int CoincidenceRun::GetMinusDetectedEventCount() const
{
  return fMinusDetectedEventCount;
}

G4int CoincidenceRun::GetCoincidenceEventCount() const
{
  return fCoincidenceEventCount;
}

RunAction::RunAction()
  : G4UserRunAction()
{
}

RunAction::~RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  const G4bool isMultithreadedMaster =
    G4Threading::IsMultithreadedApplication() && IsMaster();
  if (!isMultithreadedMaster) {
    RootOutput::Instance()->Open("mcp_output.root");
  }

  G4cout << "### Run " << run->GetRunID() << " starts." << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  const G4bool isMultithreadedMaster =
    G4Threading::IsMultithreadedApplication() && IsMaster();

  if (!isMultithreadedMaster) {
    RootOutput::Instance()->Close();
  }

  if (IsMaster()) {
    const CoincidenceRun* coincidenceRun =
      dynamic_cast<const CoincidenceRun*>(run);
    const G4int eventCount =
      coincidenceRun
        ? coincidenceRun->GetProcessedEventCount()
        : 0;
    const G4int coincidenceCount =
      coincidenceRun
        ? coincidenceRun->GetCoincidenceEventCount()
        : 0;
    const G4double coincidenceEfficiency =
      eventCount > 0
        ? static_cast<G4double>(coincidenceCount)/eventCount
        : 0.0;

    G4cout << G4endl
           << "=== MCP/PET detection summary ===" << G4endl
           << "Processed events: " << eventCount << G4endl
           << "Events detected on +z: "
           << (coincidenceRun
                 ? coincidenceRun->GetPlusDetectedEventCount()
                 : 0)
           << G4endl
           << "Events detected on -z: "
           << (coincidenceRun
                 ? coincidenceRun->GetMinusDetectedEventCount()
                 : 0)
           << G4endl
           << "Coincidence events: " << coincidenceCount << G4endl
           << "Coincidence efficiency: "
           << coincidenceEfficiency
           << " (" << 100.0*coincidenceEfficiency << " %)"
           << G4endl
           << "=================================" << G4endl;
  }

  G4cout << "### Run " << run->GetRunID()
         << " ends after " << run->GetNumberOfEvent()
         << " events." << G4endl;
}

G4Run* RunAction::GenerateRun()
{
  return new CoincidenceRun;
}

void RunAction::CountEvent(G4bool hasPlusHit,
                           G4bool hasMinusHit)
{
  CoincidenceRun* run = dynamic_cast<CoincidenceRun*>(
    G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  if (run) {
    run->CountEvent(hasPlusHit, hasMinusHit);
  }
}
