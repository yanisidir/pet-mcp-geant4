#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "OutputMessenger.hh"
#include "ScanSummaryInfo.hh"

#include "G4Run.hh"
#include "G4ios.hh"

RunAction::RunAction(const DetectorConstruction* detector)
  : fDetector(detector),
    fMessenger(new OutputMessenger(this))
{
}

RunAction::~RunAction()
{
  delete fMessenger;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  fPlusCount = fMinusCount = fCoincidenceCount = 0;
  fElectronHitsPlus = fElectronHitsMinus = 0;
  fOutput.Open(fFileName.c_str());
}

void RunAction::RecordEvent(const EventSummaryInfo& summary)
{
  if (summary.electronChannelPlusCount > 0) ++fPlusCount;
  if (summary.electronChannelMinusCount > 0) ++fMinusCount;
  if (summary.isCoincidence) ++fCoincidenceCount;
  fElectronHitsPlus += summary.electronChannelPlusCount;
  fElectronHitsMinus += summary.electronChannelMinusCount;
  fOutput.FillEventSummary(summary);
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  fOutput.Close();
  const auto events = run->GetNumberOfEvent();
  const auto electronHits = fElectronHitsPlus + fElectronHitsMinus;

  G4cout << "Processed events: " << events << G4endl
         << "Electrons reaching a channel: " << electronHits << G4endl
         << "Events detected on +z: " << fPlusCount << G4endl
         << "Events detected on -z: " << fMinusCount << G4endl
         << "Coincidence events: " << fCoincidenceCount << G4endl
         << "Coincidence efficiency: "
         << (events > 0 ? G4double(fCoincidenceCount)/events : 0.) << G4endl;

  if (!fOutput.HasScanFile()) return;

  ScanSummaryInfo row;
  row.stacksPlus = fDetector ? fDetector->GetStackCount(1) : 0;
  row.stacksMinus = fDetector ? fDetector->GetStackCount(-1) : 0;
  row.events = events;
  row.electronHits = electronHits;
  row.electronHitsPlus = fElectronHitsPlus;
  row.electronHitsMinus = fElectronHitsMinus;
  row.eventsWithHitPlus = fPlusCount;
  row.eventsWithHitMinus = fMinusCount;
  row.coincidences = fCoincidenceCount;
  fOutput.AppendScanRow(row);

  G4cout << "Scan row appended: stacks " << row.stacksPlus << "/"
         << row.stacksMinus << ", electrons " << row.electronHits << G4endl;
}
