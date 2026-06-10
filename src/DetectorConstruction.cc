#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4FieldBuilder.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UniformElectricField.hh"
#include "G4VisAttributes.hh"
#include "HoleParameterisation.hh"

#include <cmath>
#include <sstream>
#include <vector>

namespace
{
  // MCP geometry
  const G4int kNumberOfMCPs = 3;
  const G4double kMcpLength = 3000.0*micrometer;
  const G4double kMcpDiameter = 50000.0*micrometer;
  const G4double kMcpRadius = 0.5*kMcpDiameter;
  const G4double kMcpGap = 50.0*micrometer;
  const G4double kFirstMcpZ = 20.0*cm;

  // MCP voltage / field
  const G4double kMcpVoltage = 3.0*kilovolt;
  const G4bool kEnableElectricField = false;

  // Channel geometry
  const G4double kChannelDiameter = 25.0*micrometer;
  const G4double kChannelPitch = 31.0*micrometer;
  const G4double kChannelAngle = 8.0*deg;
}

// namespace
// {
//   // =========================================================================
//   // DEBUG VISUALIZATION GEOMETRY
//   // =========================================================================

//   // MCP geometry
//   const G4int kNumberOfMCPs = 2;

//   const G4double kMcpLength   = 50.0*mm;
//   const G4double kMcpDiameter = 100.0*mm;
//   const G4double kMcpRadius   = 0.5*kMcpDiameter;

//   const G4double kMcpGap      = 20.0*mm;

//   // Bring stacks closer to the origin
//   const G4double kFirstMcpZ   = 5.0*cm;

//   // MCP voltage / field
//   const G4double kMcpVoltage = 3.0*kilovolt;
//   const G4bool kEnableElectricField = false;

//   // Channel geometry
//   const G4double kChannelDiameter = 5.0*mm;

//   // Very large pitch => only a few visible channels
//   const G4double kChannelPitch = 20.0*mm;

//   const G4double kChannelAngle = 8.0*deg;
// }

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fWorldLogicalVolume(nullptr),
    fWorldMaterial(nullptr),
    fDetectorMaterial(nullptr),
    fChannelMaterial(nullptr)
{
  G4FieldBuilder::Instance();
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();

  fChannelLogicalVolumes.clear();
  fChannelSides.clear();
  fChannelMcpIndices.clear();
  fLastMcpLogicalVolumes.clear();
  fLastChannelLogicalVolumes.clear();
  fLastMcpSides.clear();

  G4VPhysicalVolume* worldPhysical = BuildWorld();

  BuildMcpStack(1);
  BuildMcpStack(-1);

  return worldPhysical;
}

G4bool DetectorConstruction::IsChannelLogicalVolume(
  const G4LogicalVolume* logicalVolume) const
{
  return GetChannelSide(logicalVolume) != 0;
}

G4int DetectorConstruction::GetChannelSide( // returns +1 for plus side, -1 for minus side, 0 for not a channel
  const G4LogicalVolume* logicalVolume) const
{
  G4int side = 0;
  G4int mcpIndex = -1;
  GetChannelInfo(logicalVolume, side, mcpIndex);
  return side;
}

G4bool DetectorConstruction::GetChannelInfo(
  const G4LogicalVolume* logicalVolume,
  G4int& side,
  G4int& mcpIndex) const
{
  side = 0;
  mcpIndex = -1;
  if (!logicalVolume) {
    return false;
  }

  for (std::size_t i = 0; i < fChannelLogicalVolumes.size(); ++i) {
    if (fChannelLogicalVolumes[i] == logicalVolume) {
      side = fChannelSides[i];
      mcpIndex = fChannelMcpIndices[i];
      return true;
    }
  }

  return false;
}

G4int DetectorConstruction::GetNumberOfMCPs() const
{
  return kNumberOfMCPs;
}

G4bool DetectorConstruction::IsLastMcpLogicalVolume(
  const G4LogicalVolume* logicalVolume) const
{
  return GetLastMcpSide(logicalVolume) != 0;
}

G4int DetectorConstruction::GetLastMcpSide( // returns +1 for plus side, -1 for minus side, 0 for not the last MCP or its channels
  const G4LogicalVolume* logicalVolume) const
{
  if (!logicalVolume) {
    return 0;
  }

  for (std::size_t i = 0; i < fLastMcpLogicalVolumes.size(); ++i) {
    if (fLastMcpLogicalVolumes[i] == logicalVolume ||
        fLastChannelLogicalVolumes[i] == logicalVolume) {
      return fLastMcpSides[i];
    }
  }

  return 0;
}

void DetectorConstruction::ConstructSDandField()
{
  const G4double fieldMagnitude = kMcpVoltage/kMcpLength;

  G4ThreeVector fieldVector(0.0, 0.0, 0.0);

  if (kEnableElectricField) {
    fieldVector = G4ThreeVector(0.0, 0.0, -fieldMagnitude);
  }

  G4UniformElectricField* electricField =
    new G4UniformElectricField(fieldVector);

  G4FieldBuilder* fieldBuilder = G4FieldBuilder::Instance();
  fieldBuilder->SetGlobalField(electricField);
  fieldBuilder->SetFieldType(kElectroMagnetic);
  fieldBuilder->ConstructFieldSetup();

  if (kEnableElectricField) {
    G4cout << "Uniform electric field: "
           << fieldMagnitude/(kilovolt/cm)
           << " kV/cm along -z." << G4endl;
    G4cout << "Field limitation: this global direction is not symmetric; "
           << "the -z MCP stack would need an opposite local field."
           << G4endl;
  } else {
    G4cout << "Electric field disabled." << G4endl;
  }
}

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();

  fWorldMaterial = nist->FindOrBuildMaterial("G4_AIR");
  fDetectorMaterial = nist->FindOrBuildMaterial("G4_GLASS_LEAD");
  fChannelMaterial = nist->FindOrBuildMaterial("G4_Galactic");
}

G4VPhysicalVolume* DetectorConstruction::BuildWorld()
{
  const G4double worldSize = 1.0*m;

  G4Box* worldSolid = new G4Box("WorldSolid",
                                0.5*worldSize,
                                0.5*worldSize,
                                0.5*worldSize);

  fWorldLogicalVolume =
    new G4LogicalVolume(worldSolid, fWorldMaterial, "WorldLogical");

  G4VPhysicalVolume* worldPhysical =
    new G4PVPlacement(nullptr,
                      G4ThreeVector(),
                      fWorldLogicalVolume,
                      "World",
                      nullptr,
                      false,
                      0,
                      true);

  fWorldLogicalVolume->SetVisAttributes(G4VisAttributes::GetInvisible());

  return worldPhysical;
}

void DetectorConstruction::BuildMcpStack(G4int side)
{
  for (G4int mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    const G4double distanceFromSource =
      kFirstMcpZ + mcpIndex*(kMcpLength + kMcpGap);
    const G4double chevronAngle =
      (mcpIndex % 2 == 0) ? kChannelAngle : -kChannelAngle;

    BuildMCP(side,
             mcpIndex,
             side*distanceFromSource,
             side*chevronAngle);
  }
}

void DetectorConstruction::BuildMCP(
  G4int side,
  G4int mcpIndex,
  G4double positionZ,
  G4double channelAngle)
{
  const G4String sideName = side > 0 ? "plus" : "minus";

  std::ostringstream suffixStream;
  suffixStream << mcpIndex; // suffixStream recoit le numéro de l'étage MCP pour le rendre en string
  const G4String suffix = suffixStream.str();

  const G4String nameBase = "MCP_" + sideName + "_body_" + suffix;
  const G4String solidName = nameBase + "_solid";
  const G4String logicalName = nameBase + "_log";
  const G4String physicalName = nameBase;

  G4Tubs* mcpSolid = new G4Tubs(solidName,
                                0.0,
                                kMcpRadius,
                                0.5*kMcpLength,
                                0.0,
                                360.0*deg);

  G4LogicalVolume* mcpLogical =
    new G4LogicalVolume(mcpSolid, fDetectorMaterial, logicalName);

  const G4ThreeVector mcpPosition(0.0, 0.0, positionZ);
  const G4int copyNumber =
    side > 0 ? mcpIndex : kNumberOfMCPs + mcpIndex;

  new G4PVPlacement(nullptr,
                    mcpPosition,
                    mcpLogical,
                    physicalName,
                    fWorldLogicalVolume,
                    false,
                    copyNumber,
                    true);

  // Couleurs distinctes uniquement pour la visualisation.
  // La transparence permet de voir les canaux à l'intérieur du verre.
  const G4Colour mcpColour =
    side > 0
      ? G4Colour(0.10, 0.35, 1.00, 0.45) 
      : G4Colour(0.10, 0.80, 0.35, 0.45);
  G4VisAttributes* mcpVis = new G4VisAttributes(mcpColour);
  mcpVis->SetForceSolid(true);
  mcpLogical->SetVisAttributes(mcpVis);

  const G4int channelCount =
    BuildChannels(mcpLogical, side, mcpIndex, channelAngle);

  G4cout << "MCP volume=" << physicalName
         << " | side=" << (side > 0 ? "+z" : "-z")
         << " | index=" << mcpIndex
         << " | z=" << positionZ/cm << " cm"
         << " | channelAngle=" << channelAngle/deg << " deg"
         << " | channels=" << channelCount
         << G4endl;
}

G4int DetectorConstruction::BuildChannels(
  G4LogicalVolume* mcpLogicalVolume,
  G4int side,
  G4int mcpIndex,
  G4double channelAngle)
{
  const G4double channelRadius = 0.5*kChannelDiameter;
  const G4double safeZone = 0.2*kMcpRadius;

  std::vector<G4ThreeVector> channelPositions =
    ComputeChannelPositions(kMcpRadius,
                            channelRadius,
                            kChannelPitch,
                            safeZone);

  std::ostringstream suffixStream;
  suffixStream << mcpIndex;
  const G4String suffix = suffixStream.str();
  const G4String sideName = side > 0 ? "plus" : "minus";

  const G4String nameBase = "MCP_" + sideName + "_channel_" + suffix;
  const G4String solidName = nameBase + "_solid";
  const G4String logicalName = nameBase + "_log";
  const G4String physicalName = nameBase;

  G4Tubs* channelSolid = new G4Tubs(solidName,
                                    0.0,
                                    channelRadius,
                                    0.5*kMcpLength,
                                    0.0,
                                    360.0*deg);

  G4LogicalVolume* channelLogical =
    new G4LogicalVolume(channelSolid, fChannelMaterial, logicalName);

  G4VisAttributes* channelVis =
    new G4VisAttributes(G4Colour(1.00, 0.25, 0.05, 1.00));
  channelVis->SetForceSolid(true);
  channelLogical->SetVisAttributes(channelVis);

  fChannelLogicalVolumes.push_back(channelLogical);
  fChannelSides.push_back(side);
  fChannelMcpIndices.push_back(mcpIndex);

  if (mcpIndex == kNumberOfMCPs - 1) {
    fLastMcpLogicalVolumes.push_back(mcpLogicalVolume);
    fLastChannelLogicalVolumes.push_back(channelLogical);
    fLastMcpSides.push_back(side);
  }

  HoleParameterisation* holeParameterisation =
    new HoleParameterisation(channelPositions, channelAngle);

  new G4PVParameterised(physicalName,
                        channelLogical,
                        mcpLogicalVolume,
                        kUndefined,
                        static_cast<G4int>(channelPositions.size()),
                        holeParameterisation,
                        false);

  return static_cast<G4int>(channelPositions.size());
}

std::vector<G4ThreeVector> DetectorConstruction::ComputeChannelPositions(
  G4double detectorRadius,
  G4double channelRadius,
  G4double channelPitch,
  G4double safeZone) const
{
  std::vector<G4ThreeVector> channelPositions;

  const G4double usableRadius = detectorRadius - channelRadius - safeZone;
  const G4double rowSpacing = channelPitch*std::sqrt(3.0)/2.0;

  G4int row = 0;

  for (G4double y = -usableRadius; y <= usableRadius; y += rowSpacing) {
    G4double xMax = std::sqrt(usableRadius*usableRadius - y*y);
    G4double xOffset = (row % 2 == 0) ? 0.0 : 0.5*channelPitch;

    for (G4double x = -xMax; x <= xMax; x += channelPitch) {
      G4double xpos = x + xOffset;

      if (xpos*xpos + y*y <= usableRadius*usableRadius) {
        channelPositions.push_back(G4ThreeVector(xpos, y, 0.0));
      }
    }

    row++;
  }

  return channelPositions;
}
