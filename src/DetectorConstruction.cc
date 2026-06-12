#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4FieldBuilder.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
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
  const G4int kNumberOfMCPs = 1;
  const G4double kMcpLength = 3000.0*micrometer;
  const G4double kMcpDiameter = 5000.0*micrometer;
  const G4double kMcpRadius = 0.5*kMcpDiameter;
  const G4double kMcpGap = 1000.0*micrometer;
  const G4double kFirstMcpZ = 1.0*cm;

  // MCP voltage / field
  const G4double kMcpVoltage = 3.0*kilovolt;
  const G4bool kEnableElectricField = false;

  // Channel geometry
  const G4double kChannelDiameter = 25.0*micrometer;
  const G4double kChannelPitch = 310.0*micrometer;
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
    fDetectorMaterialName("GlassLead25Perc"),
    fChannelMaterial(nullptr),
    fMessenger(nullptr)
{
  G4FieldBuilder::Instance();

  fMessenger = new G4GenericMessenger(
    this,
    "/mcp/",
    "MCP detector controls");

  G4GenericMessenger::Command& materialCommand =
    fMessenger->DeclareMethod(
      "setMaterial",
      &DetectorConstruction::SetDetectorMaterial,
      "Select the material used by all MCP bodies.");
  materialCommand.SetParameterName("materialName", false);

  materialCommand.SetCandidates(
    "GlassLead25Perc GlassLead30Perc GlassLead35Perc "
    "GlassLead40Perc GlassLead45Perc GlassLead50Perc "
    "GlassLead60Perc GlassLead75Perc");
}

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();

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

  BuildMcpStack(1);
  // BuildMcpStack(-1);

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
  fDetectorLogicalVolumes.push_back(mcpLogical);
  fDetectorSides.push_back(side);
  fDetectorMcpIndices.push_back(mcpIndex);

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
