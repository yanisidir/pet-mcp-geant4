#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : fParticleGun(new G4ParticleGun(1)),
    fMessenger(new PrimaryGeneratorMessenger(this))
{
  fParticleGun->SetParticleDefinition(G4Gamma::GammaDefinition());
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fMessenger;
  delete fParticleGun;
}

void PrimaryGeneratorAction::SetParticle(const G4String& name)
{
  auto* particle = G4ParticleTable::GetParticleTable()->FindParticle(name);
  if (!particle) {
    G4cerr << "Source: unknown particle '" << name << "', unchanged." << G4endl;
    return;
  }
  fParticleGun->SetParticleDefinition(particle);
}

void PrimaryGeneratorAction::SetConeAxis(const G4ThreeVector& value)
{
  if (value.mag2() <= 0.0) {
    G4cerr << "Source: null cone axis ignored." << G4endl;
    return;
  }
  fConeAxis = value.unit();
}

void PrimaryGeneratorAction::SetParticlesPerEvent(G4int value)
{
  if (value < 1) return;
  fParticlesPerEvent = value;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  fParticleGun->SetParticleEnergy(fEnergy);

  for (G4int i = 0; i < fParticlesPerEvent; ++i) {
    // Both gammas of a pair share one vertex.
    const G4ThreeVector position = SamplePosition();
    const G4ThreeVector direction = SampleDirectionInCone();

    fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->GeneratePrimaryVertex(event);

    if (fBackToBack) {
      fParticleGun->SetParticleMomentumDirection(-direction);
      fParticleGun->GeneratePrimaryVertex(event);
    }
  }
}

G4ThreeVector PrimaryGeneratorAction::SamplePosition() const
{
  if (fPositionSigma <= 0.0) return fPosition;

  return fPosition + G4ThreeVector(G4RandGauss::shoot(0.0, fPositionSigma),
                                   G4RandGauss::shoot(0.0, fPositionSigma),
                                   G4RandGauss::shoot(0.0, fPositionSigma));
}

G4ThreeVector PrimaryGeneratorAction::SampleDirectionInCone() const
{
  // Uniforme en angle solide, et non en angle polaire qui biaiserait vers
  // l'axe. Un demi-angle de 180 deg donne l'isotropie.
  const G4double cosThetaMin = std::cos(fConeHalfAngle);
  const G4double cosTheta =
    cosThetaMin + (1.0 - cosThetaMin)*G4UniformRand();
  const G4double sinTheta =
    std::sqrt(std::max(0.0, 1.0 - cosTheta*cosTheta));
  const G4double phi = twopi*G4UniformRand();

  const G4ThreeVector axis = fConeAxis.unit();
  const G4ThreeVector u = axis.orthogonal().unit();
  const G4ThreeVector v = axis.cross(u);

  return (sinTheta*(std::cos(phi)*u + std::sin(phi)*v)
          + cosTheta*axis).unit();
}

void PrimaryGeneratorAction::PrintConfiguration() const
{
  G4cout << "--- Source configuration ---" << G4endl
         << "  particle            : "
         << fParticleGun->GetParticleDefinition()->GetParticleName() << G4endl
         << "  energy              : " << fEnergy/keV << " keV" << G4endl
         << "  position            : (" << fPosition.x()/mm << ", "
         << fPosition.y()/mm << ", " << fPosition.z()/mm << ") mm" << G4endl
         << "  position sigma      : " << fPositionSigma/mm << " mm" << G4endl
         << "  cone axis           : (" << fConeAxis.x() << ", "
         << fConeAxis.y() << ", " << fConeAxis.z() << ")" << G4endl
         << "  cone half angle     : " << fConeHalfAngle/deg << " deg" << G4endl
         << "  back to back        : " << (fBackToBack ? "yes" : "no") << G4endl
         << "  vertices per event  : " << fParticlesPerEvent
         << (fBackToBack ? " pair(s)" : " particle(s)") << G4endl
         << "----------------------------" << G4endl;
}
