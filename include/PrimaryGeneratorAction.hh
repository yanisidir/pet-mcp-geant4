#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
class G4Event;
class G4ParticleGun;
class PrimaryGeneratorMessenger;

// Back-to-back 511 keV gammas by default, emitted in a cone;
// configured by PrimaryGeneratorMessenger.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event* event) override;

  void SetParticle(const G4String& name);
  void SetEnergy(G4double value) { fEnergy = value; }
  void SetPosition(const G4ThreeVector& value) { fPosition = value; }
  void SetPositionSigma(G4double value) { fPositionSigma = value; }
  void SetConeAxis(const G4ThreeVector& value);
  void SetConeHalfAngle(G4double value) { fConeHalfAngle = value; }
  void SetBackToBack(G4bool value) { fBackToBack = value; }
  void SetParticlesPerEvent(G4int value);
  void PrintConfiguration() const;

private:
  G4ThreeVector SampleDirectionInCone() const;
  G4ThreeVector SamplePosition() const;

  G4ParticleGun* fParticleGun;
  PrimaryGeneratorMessenger* fMessenger;

  G4double fEnergy = 511.0*keV;
  G4ThreeVector fPosition;
  G4double fPositionSigma = 0.0;
  G4ThreeVector fConeAxis{0.0, 0.0, 1.0};
  G4double fConeHalfAngle = 75.0*deg;
  G4int fParticlesPerEvent = 1;
  G4bool fBackToBack = true;
};
#endif
