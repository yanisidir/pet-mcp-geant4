#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "ElectricFieldSetup.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "HoleParameterisation.hh"

#include <cmath>
#include <sstream>
#include <vector>

namespace
{
  // MCP geometry
  const G4int kNumberOfMCPs = 1;
  const G4double kMcpLength = 3000.0*micrometer;
  const G4double kMcpDiameter = 50000.0*micrometer;
  const G4double kMcpRadius = 0.5*kMcpDiameter;
  const G4double kMcpGap = 1000.0*micrometer;
  const G4double kFirstMcpZ = 1.0*cm;

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
    fDetectorMaterialName("GlassLead60Perc"),
    fChannelMaterial(nullptr),
    fMessenger(nullptr),
    fPlusStackFieldSetup(nullptr),
    fMinusStackFieldSetup(nullptr),
    fEnableStackElectricFields(false),
    fStackFieldMagnitude(1.0*kilovolt/mm),
    fFieldRegionRadialMargin(1.0*mm),
    fFieldRegionZMargin(1.0*mm)
{
  fMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete fPlusStackFieldSetup;
  delete fMinusStackFieldSetup;
  delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();

  delete fPlusStackFieldSetup;
  delete fMinusStackFieldSetup;
  fPlusStackFieldSetup = nullptr;
  fMinusStackFieldSetup = nullptr;

  fChannelLogicalVolumes.clear();
  fChannelSides.clear();
  fChannelMcpIndices.clear();
  fDetectorLogicalVolumes.clear();
  fDetectorSides.clear();
  fDetectorMcpIndices.clear();
  fLastMcpLogicalVolumes.clear();
  fLastChannelLogicalVolumes.clear();
  fLastMcpSides.clear();

  G4VPhysicalVolume* worldPhysical = BuildWorld();

  G4LogicalVolume* plusMother = fWorldLogicalVolume;
  G4LogicalVolume* minusMother = fWorldLogicalVolume;
  G4double plusMotherCenterZ = 0.0;
  G4double minusMotherCenterZ = 0.0;

  if (fEnableStackElectricFields) {
    const G4double firstCenter = kFirstMcpZ;
    const G4double lastCenter = kFirstMcpZ + (kNumberOfMCPs - 1)*(kMcpLength + kMcpGap);
    const G4double stackMinZ = firstCenter - 0.5*kMcpLength;
    const G4double stackMaxZ = lastCenter + 0.5*kMcpLength;
    const G4double fieldCenterDistance = 0.5*(stackMinZ + stackMaxZ);

    plusMother = BuildElectricFieldRegion(1);
    minusMother = BuildElectricFieldRegion(-1);
    plusMotherCenterZ = fieldCenterDistance;
    minusMotherCenterZ = -fieldCenterDistance;
  }

  BuildMcpStack(1, plusMother, plusMotherCenterZ);
  BuildMcpStack(-1, minusMother, minusMotherCenterZ);

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

G4bool DetectorConstruction::GetMcpInfo(
  const G4LogicalVolume* logicalVolume,
  G4int& side,
  G4int& mcpIndex) const
{
  side = 0;
  mcpIndex = -1;
  if (!logicalVolume) {
    return false;
  }

  for (std::size_t i = 0; i < fDetectorLogicalVolumes.size(); ++i) {
    if (fDetectorLogicalVolumes[i] == logicalVolume) {
      side = fDetectorSides[i];
      mcpIndex = fDetectorMcpIndices[i];
      return true;
    }
  }

  return false;
}

void DetectorConstruction::SetDetectorMaterial(
  const G4String& materialName)
{
  fDetectorMaterialName = materialName;

  G4Material* material =
    G4Material::GetMaterial(materialName, false);
  if (!material) {
    // Permet d'utiliser /mcp/setMaterial avant /run/initialize.
    DefineMaterials();
    material = G4Material::GetMaterial(materialName, false);
  }

  if (!material) {
    G4ExceptionDescription message;
    message << "Unknown MCP material: " << materialName;
    G4Exception("DetectorConstruction::SetDetectorMaterial",
                "DetectorMaterial002",
                FatalException,
                message);
    return;
  }

  fDetectorMaterial = material;
  for (std::size_t i = 0;
       i < fDetectorLogicalVolumes.size();
       ++i) {
    fDetectorLogicalVolumes[i]->SetMaterial(fDetectorMaterial);
  }

  G4RunManager* runManager = G4RunManager::GetRunManager();
  if (runManager) {
    runManager->PhysicsHasBeenModified();
  }

  G4cout << "MCP material selected: "
         << fDetectorMaterial->GetName()
         << " | density="
         << fDetectorMaterial->GetDensity()/(g/cm3)
         << " g/cm3"
         << " | updated volumes="
         << fDetectorLogicalVolumes.size()
         << G4endl;
}

void DetectorConstruction::SetStackElectricFieldsEnabled(G4bool enabled)
{
  fEnableStackElectricFields = enabled;
  if (fWorldLogicalVolume) {
    G4cout << "WARNING: /mcp/enableStackElectricFields changed after "
           << "geometry construction. Run /run/initialize again or "
           << "reinitialize geometry before starting a new run."
           << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    if (runManager) {
      runManager->GeometryHasBeenModified();
    }
  }

  G4cout << "MCP stack electric fields "
         << (fEnableStackElectricFields ? "enabled" : "disabled")
         << "." << G4endl;
}

void DetectorConstruction::SetStackFieldMagnitude(G4double magnitude)
{
  if (magnitude <= 0.0) {
    G4cout << "WARNING: stack field magnitude must be positive. "
           << "Keeping previous value." << G4endl;
    return;
  }

  fStackFieldMagnitude = magnitude;
  if (fWorldLogicalVolume) {
    G4cout << "WARNING: /mcp/stackFieldMagnitude changed after "
           << "geometry construction. Run /run/initialize again or "
           << "reinitialize geometry before starting a new run."
           << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    if (runManager) {
      runManager->GeometryHasBeenModified();
    }
  }

  G4cout << "MCP stack field magnitude set to "
         << fStackFieldMagnitude/(kilovolt/mm)
         << " kV/mm." << G4endl;
}

void DetectorConstruction::SetFieldRegionRadialMargin(G4double margin)
{
  if (margin < 0.0) {
    G4cout << "WARNING: field region radial margin cannot be negative. "
           << "Keeping previous value." << G4endl;
    return;
  }

  fFieldRegionRadialMargin = margin;
  if (fWorldLogicalVolume) {
    G4cout << "WARNING: /mcp/fieldRegionRadialMargin changed after "
           << "geometry construction. Run /run/initialize again or "
           << "reinitialize geometry before starting a new run."
           << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    if (runManager) {
      runManager->GeometryHasBeenModified();
    }
  }

  G4cout << "MCP field region radial margin set to "
         << fFieldRegionRadialMargin/mm
         << " mm." << G4endl;
}

void DetectorConstruction::SetFieldRegionZMargin(G4double margin)
{
  if (margin < 0.0) {
    G4cout << "WARNING: field region z margin cannot be negative. "
           << "Keeping previous value." << G4endl;
    return;
  }

  fFieldRegionZMargin = margin;
  if (fWorldLogicalVolume) {
    G4cout << "WARNING: /mcp/fieldRegionZMargin changed after "
           << "geometry construction. Run /run/initialize again or "
           << "reinitialize geometry before starting a new run."
           << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    if (runManager) {
      runManager->GeometryHasBeenModified();
    }
  }

  G4cout << "MCP field region z margin set to "
         << fFieldRegionZMargin/mm
         << " mm." << G4endl;
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
  if (fEnableStackElectricFields) {
    G4cout << "Global electric field disabled. "
           << "MCP fields are attached locally to stack field regions."
           << G4endl;
  } else {
    G4cout << "MCP stack electric fields disabled." << G4endl;
  }
}


void DetectorConstruction::DefineMaterials()
{
  G4NistManager* nist = G4NistManager::Instance();

  G4Element* lead = nist->FindOrBuildElement("Pb");
  G4Element* silicon = nist->FindOrBuildElement("Si");
  G4Element* oxygen = nist->FindOrBuildElement("O");

  // Densities are approximate real lead-glass values.
  // Intermediate values are linearly interpolated.

  if (!G4Material::GetMaterial("GlassLead25Perc", false)) {
    G4Material* glassLead25 =
      new G4Material("GlassLead25Perc", 3.32*g/cm3, 3);
    glassLead25->AddElement(lead, 0.25);
    glassLead25->AddElement(silicon, 0.30);
    glassLead25->AddElement(oxygen, 0.45);
  }

  if (!G4Material::GetMaterial("GlassLead30Perc", false)) {
    G4Material* glassLead30 =
      new G4Material("GlassLead30Perc", 3.51*g/cm3, 3);
    glassLead30->AddElement(lead, 0.30);
    glassLead30->AddElement(silicon, 0.28);
    glassLead30->AddElement(oxygen, 0.42);
  }

  if (!G4Material::GetMaterial("GlassLead35Perc", false)) {
    G4Material* glassLead35 =
      new G4Material("GlassLead35Perc", 3.70*g/cm3, 3);
    glassLead35->AddElement(lead, 0.35);
    glassLead35->AddElement(silicon, 0.26);
    glassLead35->AddElement(oxygen, 0.39);
  }

  if (!G4Material::GetMaterial("GlassLead40Perc", false)) {
    G4Material* glassLead40 =
      new G4Material("GlassLead40Perc", 3.89*g/cm3, 3);
    glassLead40->AddElement(lead, 0.40);
    glassLead40->AddElement(silicon, 0.24);
    glassLead40->AddElement(oxygen, 0.36);
  }

  if (!G4Material::GetMaterial("GlassLead45Perc", false)) {
    G4Material* glassLead45 =
      new G4Material("GlassLead45Perc", 4.08*g/cm3, 3);
    glassLead45->AddElement(lead, 0.45);
    glassLead45->AddElement(silicon, 0.22);
    glassLead45->AddElement(oxygen, 0.33);
  }

  if (!G4Material::GetMaterial("GlassLead50Perc", false)) {
    G4Material* glassLead50 =
      new G4Material("GlassLead50Perc", 4.36*g/cm3, 3);
    glassLead50->AddElement(lead, 0.50);
    glassLead50->AddElement(silicon, 0.20);
    glassLead50->AddElement(oxygen, 0.30);
  }

  if (!G4Material::GetMaterial("GlassLead60Perc", false)) {
    G4Material* glassLead60 =
      new G4Material("GlassLead60Perc", 4.93*g/cm3, 3);
    glassLead60->AddElement(lead, 0.60);
    glassLead60->AddElement(silicon, 0.16);
    glassLead60->AddElement(oxygen, 0.24);
  }

  if (!G4Material::GetMaterial("GlassLead75Perc", false)) {
    G4Material* glassLead75 =
      new G4Material("GlassLead75Perc", 6.22*g/cm3, 3);
    glassLead75->AddElement(lead, 0.75);
    glassLead75->AddElement(silicon, 0.10);
    glassLead75->AddElement(oxygen, 0.15);
  }

  fWorldMaterial = nist->FindOrBuildMaterial("G4_Galactic");
  fChannelMaterial = nist->FindOrBuildMaterial("G4_Galactic");

  fDetectorMaterial =
    G4Material::GetMaterial(fDetectorMaterialName, false);

  if (!fDetectorMaterial) {
    G4ExceptionDescription message;
    message << "Unknown MCP material: " << fDetectorMaterialName
            << ". Available choices are GlassLead25Perc, "
            << "GlassLead30Perc, GlassLead35Perc, "
            << "GlassLead40Perc, GlassLead45Perc, "
            << "GlassLead50Perc, GlassLead60Perc and GlassLead75Perc.";
    G4Exception("DetectorConstruction::DefineMaterials",
                "DetectorMaterial001",
                FatalException,
                message);
    return;
  }

  const G4ElementVector* elements =
    fDetectorMaterial->GetElementVector();
  const G4double* fractions =
    fDetectorMaterial->GetFractionVector();

  G4cout << "MCP material: " << fDetectorMaterial->GetName()
         << " | density="
         << fDetectorMaterial->GetDensity()/(g/cm3)
         << " g/cm3 | mass fractions:";

  for (std::size_t i = 0;
       i < fDetectorMaterial->GetNumberOfElements();
       ++i) {
    G4cout << " " << (*elements)[i]->GetSymbol()
           << "=" << fractions[i];
  }

  G4cout << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::BuildWorld()
{
  const G4double worldSize = 10*cm;

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

G4LogicalVolume* DetectorConstruction::BuildElectricFieldRegion(G4int side)
{
  const G4double firstCenter = kFirstMcpZ;
  const G4double lastCenter =
    kFirstMcpZ + (kNumberOfMCPs - 1)*(kMcpLength + kMcpGap);
  const G4double stackMinZ = firstCenter - 0.5*kMcpLength;
  const G4double stackMaxZ = lastCenter + 0.5*kMcpLength;
  const G4double fieldLength =
    (stackMaxZ - stackMinZ) + 2.0*fFieldRegionZMargin;
  const G4double fieldCenterDistance = 0.5*(stackMinZ + stackMaxZ);
  const G4double fieldRadius =
    kMcpRadius + fFieldRegionRadialMargin;

  const G4String sideName = side > 0 ? "plus" : "minus";
  const G4String nameBase = "MCP_" + sideName + "_field_region";

  G4Tubs* fieldSolid = new G4Tubs(nameBase + "_solid",
                                  0.0,
                                  fieldRadius,
                                  0.5*fieldLength,
                                  0.0,
                                  360.0*deg);

  G4LogicalVolume* fieldLogical =
    new G4LogicalVolume(fieldSolid,
                        fChannelMaterial,
                        nameBase + "_log");

  const G4double fieldCenterZ = side*fieldCenterDistance;
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0.0, 0.0, fieldCenterZ),
                    fieldLogical,
                    nameBase,
                    fWorldLogicalVolume,
                    false,
                    side > 0 ? 1 : 2, // copy number: 1 for +z side, 2 for -z side
                    true);

  // The field region is invisible in the visualization.                  
  fieldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

  const G4ThreeVector fieldVector =
    side > 0
      ? G4ThreeVector(0.0, 0.0, -fStackFieldMagnitude)
      : G4ThreeVector(0.0, 0.0,  fStackFieldMagnitude);

  ElectricFieldSetup* fieldSetup =
    new ElectricFieldSetup(fieldVector);
  fieldLogical->SetFieldManager(fieldSetup->GetFieldManager(), true);

  if (side > 0) {
    fPlusStackFieldSetup = fieldSetup;
  } else {
    fMinusStackFieldSetup = fieldSetup;
  }

  G4cout << "Local MCP electric field enabled"
         << " | side=" << (side > 0 ? "+z" : "-z")
         << " | magnitude=" << fStackFieldMagnitude/(kilovolt/mm)
         << " kV/mm"
         << " | direction=("
         << fieldVector.x()/(kilovolt/mm) << ", "
         << fieldVector.y()/(kilovolt/mm) << ", "
         << fieldVector.z()/(kilovolt/mm) << ") kV/mm"
         << " | radius=" << fieldRadius/mm << " mm"
         << " | length=" << fieldLength/mm << " mm"
         << " | z=" << fieldCenterZ/cm << " cm"
         << G4endl;

  return fieldLogical;
}

void DetectorConstruction::BuildMcpStack(
  G4int side,
  G4LogicalVolume* motherLogicalVolume,
  G4double motherCenterZ)
{
  for (G4int mcpIndex = 0; mcpIndex < kNumberOfMCPs; ++mcpIndex) {
    const G4double distanceFromSource =
      kFirstMcpZ + mcpIndex*(kMcpLength + kMcpGap);
    const G4double chevronAngle =
      (mcpIndex % 2 == 0) ? kChannelAngle : -kChannelAngle;

    BuildMCP(side,
             mcpIndex,
             side*distanceFromSource,
             side*chevronAngle,
             motherLogicalVolume,
             motherCenterZ);
  }
}

void DetectorConstruction::BuildMCP(
  G4int side,
  G4int mcpIndex,
  G4double positionZ,
  G4double channelAngle,
  G4LogicalVolume* motherLogicalVolume,
  G4double motherCenterZ)
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
  fDetectorLogicalVolumes.push_back(mcpLogical);
  fDetectorSides.push_back(side);
  fDetectorMcpIndices.push_back(mcpIndex);

  const G4ThreeVector mcpPosition(
    0.0,
    0.0,
    positionZ - motherCenterZ);
  const G4int copyNumber =
    side > 0 ? mcpIndex : kNumberOfMCPs + mcpIndex;

  new G4PVPlacement(nullptr,
                    mcpPosition,
                    mcpLogical,
                    physicalName,
                    motherLogicalVolume,
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
