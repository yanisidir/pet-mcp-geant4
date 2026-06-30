#include "ElectricFieldSetup.hh"

#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4EqMagElectricField.hh"
#include "G4FieldManager.hh"
#include "G4IntegrationDriver.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField.hh"

ElectricFieldSetup::ElectricFieldSetup(
  const G4ThreeVector& fieldVector)
  : fElectricField(new G4UniformElectricField(fieldVector)),
    fEquation(new G4EqMagElectricField(fElectricField)),
    fStepper(new G4ClassicalRK4(fEquation, 8)),
    fChordFinder(nullptr),
    fFieldManager(nullptr)
{
  const G4double minimumStep = 0.001*mm;

  G4IntegrationDriver<G4ClassicalRK4>* integrationDriver =
    new G4IntegrationDriver<G4ClassicalRK4>(
      minimumStep,
      fStepper,
      fStepper->GetNumberOfVariables());

  fChordFinder = new G4ChordFinder(integrationDriver);
  fFieldManager = new G4FieldManager();
  fFieldManager->SetDetectorField(fElectricField);
  fFieldManager->SetChordFinder(fChordFinder);
}

ElectricFieldSetup::~ElectricFieldSetup()
{
  // G4FieldManager::SetChordFinder() stores the pointer but does not take
  // ownership when the chord finder was created manually. Detach the pointers
  // first so the manager cannot keep references to objects deleted below.
  if (fFieldManager) {
    fFieldManager->SetChordFinder(nullptr);
    fFieldManager->SetDetectorField(nullptr);
  }

  delete fFieldManager;
  delete fChordFinder;
  delete fStepper;
  delete fEquation;
  delete fElectricField;
}

G4FieldManager* ElectricFieldSetup::GetFieldManager() const
{
  return fFieldManager;
}
