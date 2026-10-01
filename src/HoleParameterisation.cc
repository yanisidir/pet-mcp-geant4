#include "HoleParameterisation.hh"

#include "G4VPhysicalVolume.hh"

HoleParameterisation::HoleParameterisation(
  const std::vector<G4ThreeVector>& positions,
  G4double angle)
  : G4VPVParameterisation(),
    fPositions(positions),
    fRotation()
{
  // Incline l'axe local z des cylindres dans le plan x-z.
  fRotation.rotateY(angle);
}

void HoleParameterisation::ComputeTransformation(
  G4int copyNo,
  G4VPhysicalVolume* physicalVolume) const
{
  // Geant4 calls this once per copy to know where the channel goes.
  physicalVolume->SetTranslation(fPositions[copyNo]);
  physicalVolume->SetRotation(&fRotation);
}
