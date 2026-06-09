#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <vector>

class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

// DetectorConstruction décrit entièrement la géométrie
// et les matériaux de la simulation.
//
// Son rôle est de construire :
//
// - le monde Geant4 (World) ;
// - une pile MCP côté +z ;
// - une pile MCP symétrique côté -z ;
// - les canaux paramétrés de chaque MCP ;
// - les matériaux associés.
//
// C'est la classe qui définit l'environnement physique
// dans lequel les photons gamma vont se propager
// et interagir.
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  // Constructeur.
  DetectorConstruction();

  virtual ~DetectorConstruction();

  // Fonction principale appelée par Geant4
  // pour construire toute la géométrie.
  //
  // Elle construit successivement :
  //
  // World
  //   ↓
  // MCP stacks (+z et -z)
  //   ↓
  // MCP channels paramétrés
  virtual G4VPhysicalVolume* Construct();

  // Appelée par Geant4 après la construction de la géométrie.
  // Elle installe le champ électrique uniforme dans tout le World.
  virtual void ConstructSDandField();

  // Indique si le volume logique appartient à un canal de l'un des MCP.
  G4bool IsChannelLogicalVolume(
    const G4LogicalVolume* logicalVolume) const;

  // Retourne +1 pour un canal côté +z, -1 pour un canal côté -z,
  // et 0 si le volume n'est pas un canal MCP.
  G4int GetChannelSide(const G4LogicalVolume* logicalVolume) const;

  // Retourne le côté et l'index de la plaque contenant le canal.
  // Retourne false si le volume logique n'est pas un canal MCP.
  G4bool GetChannelInfo(const G4LogicalVolume* logicalVolume,
                        G4int& side,
                        G4int& mcpIndex) const;

  // Indique si le volume logique appartient au dernier MCP
  // de l'une des deux piles.
  G4bool IsLastMcpLogicalVolume(
    const G4LogicalVolume* logicalVolume) const;

  // Retourne le côté du dernier MCP: +1, -1 ou 0.
  G4int GetLastMcpSide(const G4LogicalVolume* logicalVolume) const;

private:

  // Définit les matériaux utilisés dans la simulation.
  //
  // Actuellement :
  // - air pour le monde ;
  // - verre au plomb pour le MCP ;
  // - vide (Galactic) pour les canaux.
  void DefineMaterials();

  // Construit le volume World.
  //
  // Le World contient l'ensemble de la simulation.
  G4VPhysicalVolume* BuildWorld();

  // Construit une pile MCP complète du côté demandé (+1 ou -1).
  void BuildMcpStack(G4int side);

  // Construit un étage MCP complet à la position et à l'angle demandés.
  void BuildMCP(G4int side,
                G4int mcpIndex,
                G4double positionZ,
                G4double channelAngle);

  // Construit les canaux dans le volume logique de l'étage concerné.
  G4int BuildChannels(G4LogicalVolume* mcpLogicalVolume,
                      G4int side,
                      G4int mcpIndex,
                      G4double channelAngle);

  // Calcule les positions des centres des canaux.
  //
  // Les positions sont organisées selon
  // un réseau hexagonal.
  //
  // detectorRadius :
  //   rayon du MCP.
  //
  // channelRadius :
  //   rayon d'un canal.
  //
  // channelPitch :
  //   distance entre deux canaux voisins.
  //
  // safeZone :
  //   marge de sécurité près du bord du MCP.
  std::vector<G4ThreeVector> ComputeChannelPositions(
    G4double detectorRadius,
    G4double channelRadius,
    G4double channelPitch,
    G4double safeZone) const;

  // =====================================================
  // VOLUMES LOGIQUES
  // =====================================================

  // Volume logique du monde.
  G4LogicalVolume* fWorldLogicalVolume;

  // Volumes logiques des canaux et côté correspondant.
  std::vector<G4LogicalVolume*> fChannelLogicalVolumes;
  std::vector<G4int> fChannelSides;
  std::vector<G4int> fChannelMcpIndices;

  // Corps et canaux du dernier étage de chaque pile.
  std::vector<G4LogicalVolume*> fLastMcpLogicalVolumes;
  std::vector<G4LogicalVolume*> fLastChannelLogicalVolumes;
  std::vector<G4int> fLastMcpSides;

  // =====================================================
  // MATERIAUX
  // =====================================================

  // Matériau du monde.
  //
  // Actuellement : G4_AIR.
  G4Material* fWorldMaterial;

  // Matériau du corps du MCP.
  //
  // Actuellement : G4_GLASS_LEAD.
  G4Material* fDetectorMaterial;

  // Matériau des canaux.
  //
  // Actuellement : G4_Galactic
  G4Material* fChannelMaterial;
};

#endif
