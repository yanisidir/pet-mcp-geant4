#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <map>
#include <vector>

class DetectorMessenger;
class G4LogicalVolume;
class G4Material;
class G4ProductionCuts;
class G4Region;
class G4VPhysicalVolume;
class HoleParameterisation;

// Stacks of lead-glass MCP plates with parameterised vacuum channels,
// configured by DetectorMessenger.
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override;
  G4VPhysicalVolume* Construct() override;

  // +1 / -1 identifies the side; 0 means outside the requested volumes.
  G4int GetChannelSide(const G4LogicalVolume* volume) const;
  G4int GetMcpSide(const G4LogicalVolume* volume) const;
  // Plate rank in its stack, counted from the origin; -1 outside.
  G4int GetStackIndex(const G4LogicalVolume* volume) const;

  void SetStackCount(G4int side, G4int count);
  G4int GetStackCount(G4int side) const;
  void SetPlateThickness(G4double value) { fPlateThickness = value; }
  void SetPlateRadius(G4double value) { fPlateRadius = value; }
  void SetFirstPlateZ(G4double value) { fFirstPlateZ = value; }
  void SetStackGap(G4double value) { fStackGap = value; }
  void SetChannelDiameter(G4double value) { fChannelDiameter = value; }
  void SetChannelPitch(G4double value) { fChannelPitch = value; }
  void SetChannelAngle(G4double value) { fChannelAngle = value; }
  void SetChevron(G4bool value) { fChevron = value; }
  void SetEdgeMarginFraction(G4double value) { fEdgeMarginFraction = value; }
  void SetWorldSize(G4double value) { fWorldSize = value; }
  void SetProductionCut(G4double value);
  void SetKillGammaOnExit(G4bool value) { fKillGammaOnExit = value; }
  G4bool KillGammaOnExit() const { return fKillGammaOnExit; }
  void PrintConfiguration() const;

private:
  void DefineMaterials();
  void CleanGeometry();
  G4VPhysicalVolume* BuildWorld();
  void BuildStack(G4int side);
  void BuildPlate(G4int side, G4int stackIndex);
  G4int BuildChannels(G4LogicalVolume* mother,
                      G4int side,
                      G4int stackIndex,
                      G4double angle);
  const std::vector<G4ThreeVector>& ChannelPositions();
  G4double PlateCentreZ(G4int side, G4int stackIndex) const;
  G4double RequiredWorldSize() const;

  // Volumes and materials belong to Geant4 stores; this class owns the cuts,
  // the messenger and the parameterisations.
  G4LogicalVolume* fWorldLogicalVolume = nullptr;
  std::vector<G4LogicalVolume*> fBodyVolumes;
  std::vector<HoleParameterisation*> fParameterisations;
  std::map<const G4LogicalVolume*, std::pair<G4int, G4int>> fBodyIds;
  std::map<const G4LogicalVolume*, std::pair<G4int, G4int>> fChannelIds;
  std::vector<G4ThreeVector> fChannelPositions;
  G4Material* fVacuum = nullptr;
  G4Material* fDetectorMaterial = nullptr;
  G4ProductionCuts* fMcpProductionCuts;
  G4Region* fMcpRegion = nullptr;
  DetectorMessenger* fMessenger;

  G4int fStackCountPlus = 1;
  G4int fStackCountMinus = 1;
  G4double fPlateThickness = 1500.0*micrometer;
  G4double fPlateRadius = 0.5*50000.0*micrometer;
  G4double fFirstPlateZ = 1.0*cm;
  G4double fStackGap = 100.0*micrometer;
  G4double fChannelDiameter = 25.0*micrometer;
  G4double fChannelPitch = 31.0*micrometer;
  G4double fChannelAngle = 8.0*deg;
  G4double fEdgeMarginFraction = 0.2;
  G4double fWorldSize = 10.0*cm;
  G4bool fChevron = true;
  G4bool fKillGammaOnExit = false;
};
#endif
