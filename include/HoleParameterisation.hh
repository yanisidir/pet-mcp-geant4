#ifndef HOLE_PARAMETERISATION_HH
#define HOLE_PARAMETERISATION_HH
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPVParameterisation.hh"
#include <vector>
class G4VPhysicalVolume;

// Shared channel shape; a position for each copy and one rotation per plate.
class HoleParameterisation : public G4VPVParameterisation
{
public:
  HoleParameterisation(const std::vector<G4ThreeVector>& positions, G4double angle);
  void ComputeTransformation(G4int copyNo, G4VPhysicalVolume* volume) const override;
private:
  std::vector<G4ThreeVector> fPositions;
  mutable G4RotationMatrix fRotation;  // Geant4 expects a non-const pointer.
};
#endif
