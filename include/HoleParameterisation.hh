#ifndef HOLE_PARAMETERISATION_HH
#define HOLE_PARAMETERISATION_HH

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VPVParameterisation.hh"
#include "globals.hh"

#include <vector>

class G4Tubs;
class G4VPhysicalVolume;

// HoleParameterisation est utilisée par Geant4 pour placer
// automatiquement un grand nombre de canaux MCP identiques.
//
// Au lieu de créer manuellement plusieurs centaines ou milliers
// de G4PVPlacement, cette classe fournit à Geant4 la position
// de chaque canal à partir d'une liste de coordonnées.
//
// Dans cette simulation :
// - tous les canaux ont la même géométrie ;
// - tous les canaux ont le même matériau ;
// - seule leur position change.
class HoleParameterisation : public G4VPVParameterisation
{
public:

  // Constructeur.
  //
  // positions :
  //   liste des centres des canaux MCP calculés
  //   dans DetectorConstruction.
  //
  // angle :
  //   inclinaison commune des canaux autour de l'axe y.
  explicit HoleParameterisation(
    const std::vector<G4ThreeVector>& positions,
    G4double angle);

  virtual ~HoleParameterisation();

  // Fonction appelée par Geant4 pour chaque copie du canal.
  //
  // copyNo :
  //   numéro du canal à placer.
  //
  // physicalVolume :
  //   volume physique correspondant au canal.
  //
  // Cette fonction définit :
  // - la position du canal ;
  // - son orientation éventuelle.
  virtual void ComputeTransformation(
    G4int copyNo,
    G4VPhysicalVolume* physicalVolume) const;

  // Fonction appelée par Geant4 pour définir
  // les dimensions d'une copie.
  //
  // Dans cette simulation, tous les canaux ont
  // exactement les mêmes dimensions.
  //
  // La fonction reste donc vide dans le .cc.
  virtual void ComputeDimensions(
    G4Tubs& channel,
    G4int copyNo,
    const G4VPhysicalVolume* physicalVolume) const;

private:

  // Liste des positions des centres des canaux MCP.
  //
  // fPositions[i]
  // correspond à la position du canal numéro i.
  //
  // Ces positions sont généralement organisées
  // selon un réseau hexagonal.
  std::vector<G4ThreeVector> fPositions;

  // Rotation commune à tous les canaux.
  // "mutable" permet de transmettre son adresse à Geant4 depuis
  // ComputeTransformation(), qui est une fonction const.
  mutable G4RotationMatrix fRotation;
};

#endif
