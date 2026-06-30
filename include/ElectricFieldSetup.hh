#ifndef ELECTRIC_FIELD_SETUP_HH
#define ELECTRIC_FIELD_SETUP_HH

#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ChordFinder;
class G4ClassicalRK4;
class G4EqMagElectricField;
class G4FieldManager;
class G4UniformElectricField;

// Small owner class for one local uniform electric field.
//
// Ownership policy:
// - ElectricFieldSetup owns the electric field, equation, stepper,
//   chord finder and field manager below.
// - The integration driver is passed to G4ChordFinder; in Geant4 11.x
//   the chord finder owns that driver.
// - DetectorConstruction only attaches GetFieldManager() to a logical
//   volume. The logical volume uses this field manager but does not own it.
//
// Lifetime policy:
// An ElectricFieldSetup object must live at least as long as the logical
// volumes that use its G4FieldManager. DetectorConstruction keeps one setup
// per MCP stack and deletes it when rebuilding or destroying the detector.
class ElectricFieldSetup
{
public:
  explicit ElectricFieldSetup(const G4ThreeVector& fieldVector);
  ~ElectricFieldSetup();

  G4FieldManager* GetFieldManager() const;

private:
  ElectricFieldSetup(const ElectricFieldSetup&);
  ElectricFieldSetup& operator=(const ElectricFieldSetup&);

  G4UniformElectricField* fElectricField;
  G4EqMagElectricField* fEquation;
  G4ClassicalRK4* fStepper;
  G4ChordFinder* fChordFinder;
  G4FieldManager* fFieldManager;
};

#endif
