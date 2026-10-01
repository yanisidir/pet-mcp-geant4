#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "HoleParameterisation.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <cmath>

DetectorConstruction::DetectorConstruction()
  : fMcpProductionCuts(new G4ProductionCuts),
    fMessenger(new DetectorMessenger(this))
{
  // Resolve secondary production at the scale of the glass walls.
  fMcpProductionCuts->SetProductionCut(1.0*um);
}

DetectorConstruction::~DetectorConstruction()
{
  for (auto* parameterisation : fParameterisations) delete parameterisation;
  delete fMessenger;
  delete fMcpProductionCuts;
}

void DetectorConstruction::SetStackCount(G4int side, G4int count)
{
  if (count < 0) return;
  if (side > 0) fStackCountPlus = count;
  else fStackCountMinus = count;
}

G4int DetectorConstruction::GetStackCount(G4int side) const
{
  return side > 0 ? fStackCountPlus : fStackCountMinus;
}

void DetectorConstruction::SetProductionCut(G4double value)
{
  fMcpProductionCuts->SetProductionCut(value);
  if (fMcpRegion) {
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  CleanGeometry();
  DefineMaterials();

  auto* world = BuildWorld();
  BuildStack(1);
  BuildStack(-1);

  // The region survives a rebuild, only its root volumes change.
  fMcpRegion = G4RegionStore::GetInstance()->GetRegion("MCP_Region", false);
  if (!fMcpRegion) {
    fMcpRegion = new G4Region("MCP_Region");
    fMcpRegion->SetProductionCuts(fMcpProductionCuts);
  }
  for (auto* body : fBodyVolumes) fMcpRegion->AddRootLogicalVolume(body);

  PrintConfiguration();

  return world;
}

// Construct() runs again whenever a /mcp/ command changes the geometry after
// initialisation, so the previous volumes must be released first.
void DetectorConstruction::CleanGeometry()
{
  if (fMcpRegion) {
    for (auto* body : fBodyVolumes) fMcpRegion->RemoveRootLogicalVolume(body);
  }
  for (auto* parameterisation : fParameterisations) delete parameterisation;
  fParameterisations.clear();
  fBodyVolumes.clear();
  fBodyIds.clear();
  fChannelIds.clear();
  fChannelPositions.clear();
  fWorldLogicalVolume = nullptr;

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
}

void DetectorConstruction::DefineMaterials()
{
  auto* nist = G4NistManager::Instance();
  fVacuum = nist->FindOrBuildMaterial("G4_Galactic");
  fDetectorMaterial = G4Material::GetMaterial("GlassLead60Perc", false);
  if (!fDetectorMaterial) {
    fDetectorMaterial = new G4Material("GlassLead60Perc", 4.93*g/cm3, 3);
    fDetectorMaterial->AddElement(nist->FindOrBuildElement("Pb"), 0.60);
    fDetectorMaterial->AddElement(nist->FindOrBuildElement("Si"), 0.16);
    fDetectorMaterial->AddElement(nist->FindOrBuildElement("O"), 0.24);
  }
}

G4int DetectorConstruction::GetChannelSide(const G4LogicalVolume* volume) const
{
  if (!volume) return 0;
  const auto entry = fChannelIds.find(volume);
  return entry == fChannelIds.end() ? 0 : entry->second.first;
}

G4int DetectorConstruction::GetMcpSide(const G4LogicalVolume* volume) const
{
  if (!volume) return 0;
  const auto entry = fBodyIds.find(volume);
  if (entry != fBodyIds.end()) return entry->second.first;
  return GetChannelSide(volume);
}

G4int DetectorConstruction::GetStackIndex(const G4LogicalVolume* volume) const
{
  if (!volume) return -1;
  const auto channel = fChannelIds.find(volume);
  if (channel != fChannelIds.end()) return channel->second.second;
  const auto body = fBodyIds.find(volume);
  if (body != fBodyIds.end()) return body->second.second;
  return -1;
}

G4double DetectorConstruction::PlateCentreZ(G4int side, G4int stackIndex) const
{
  return side*(fFirstPlateZ + stackIndex*(fPlateThickness + fStackGap));
}

G4double DetectorConstruction::RequiredWorldSize() const
{
  G4double reach = fPlateRadius;
  for (const G4int side : {1, -1}) {
    const G4int count = GetStackCount(side);
    if (count <= 0) continue;
    reach = std::max(reach,
                     std::abs(PlateCentreZ(side, count - 1))
                       + 0.5*fPlateThickness);
  }
  return 2.2*reach;
}

G4VPhysicalVolume* DetectorConstruction::BuildWorld()
{
  const G4double required = RequiredWorldSize();
  G4double worldSize = fWorldSize;
  if (worldSize < required) {
    G4cout << "MCP geometry: world enlarged from " << worldSize/cm
           << " cm to " << required/cm << " cm to contain the stacks."
           << G4endl;
    worldSize = required;
  }

  G4Box* worldSolid = new G4Box("WorldSolid",
                                0.5*worldSize,
                                0.5*worldSize,
                                0.5*worldSize);

  fWorldLogicalVolume =
    new G4LogicalVolume(worldSolid, fVacuum, "WorldLogical");

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

void DetectorConstruction::BuildStack(G4int side)
{
  const G4int count = GetStackCount(side);
  for (G4int stackIndex = 0; stackIndex < count; ++stackIndex) {
    BuildPlate(side, stackIndex);
  }
}

void DetectorConstruction::BuildPlate(G4int side, G4int stackIndex)
{
  const G4double angleSign =
    (fChevron && stackIndex % 2 == 1) ? -1.0 : 1.0;
  const G4double channelAngle = side*angleSign*fChannelAngle;
  const G4String sideName = side > 0 ? "plus" : "minus";

  const G4String nameBase =
    "MCP_" + sideName + "_body_" + std::to_string(stackIndex);

  G4Tubs* mcpSolid = new G4Tubs(nameBase + "_solid",
                                0.0,
                                fPlateRadius,
                                0.5*fPlateThickness,
                                0.0,
                                360.0*deg);

  G4LogicalVolume* mcpLogical =
    new G4LogicalVolume(mcpSolid, fDetectorMaterial, nameBase + "_log");
  fBodyVolumes.push_back(mcpLogical);
  fBodyIds[mcpLogical] = {side, stackIndex};

  const G4double positionZ = PlateCentreZ(side, stackIndex);

  new G4PVPlacement(nullptr,
                    G4ThreeVector(0.0, 0.0, positionZ),
                    mcpLogical,
                    nameBase,
                    fWorldLogicalVolume,
                    false,
                    stackIndex,
                    true);

  // Couleurs distinctes uniquement pour la visualisation.
  // La transparence permet de voir les canaux à l'intérieur du verre.
  const G4Colour mcpColour =
    side > 0
      ? G4Colour(0.10, 0.35, 1.00, 0.45)
      : G4Colour(0.10, 0.80, 0.35, 0.45);
  G4VisAttributes mcpVis(mcpColour);
  mcpVis.SetForceSolid(true);
  mcpLogical->SetVisAttributes(mcpVis);

  const G4int channelCount =
    BuildChannels(mcpLogical, side, stackIndex, channelAngle);

  G4cout << "MCP volume=" << nameBase
         << " | side=" << (side > 0 ? "+z" : "-z")
         << " | stack=" << stackIndex
         << " | z=" << positionZ/cm << " cm"
         << " | channelAngle=" << channelAngle/deg << " deg"
         << " | channels=" << channelCount
         << G4endl;
}

G4int DetectorConstruction::BuildChannels(
  G4LogicalVolume* mcpLogicalVolume,
  G4int side,
  G4int stackIndex,
  G4double channelAngle)
{
  const G4double channelRadius = 0.5*fChannelDiameter;
  const std::vector<G4ThreeVector>& channelPositions = ChannelPositions();

  const G4String sideName = side > 0 ? "plus" : "minus";
  const G4String nameBase =
    "MCP_" + sideName + "_channel_" + std::to_string(stackIndex);

  G4Tubs* channelSolid = new G4Tubs(nameBase + "_solid",
                                    0.0,
                                    channelRadius,
                                    0.5*fPlateThickness,
                                    0.0,
                                    360.0*deg);

  G4LogicalVolume* channelLogical =
    new G4LogicalVolume(channelSolid, fVacuum, nameBase + "_log");

  G4VisAttributes channelVis(G4Colour(1.00, 0.25, 0.05, 1.00));
  channelVis.SetForceSolid(true);
  channelLogical->SetVisAttributes(channelVis);

  fChannelIds[channelLogical] = {side, stackIndex};

  HoleParameterisation* holeParameterisation =
    new HoleParameterisation(channelPositions, channelAngle);
  fParameterisations.push_back(holeParameterisation);

  new G4PVParameterised(nameBase,
                        channelLogical,
                        mcpLogicalVolume,
                        kUndefined,
                        static_cast<G4int>(channelPositions.size()),
                        holeParameterisation,
                        false);

  return static_cast<G4int>(channelPositions.size());
}

// Shared by every plate, so computed once per geometry build.
const std::vector<G4ThreeVector>& DetectorConstruction::ChannelPositions()
{
  if (!fChannelPositions.empty()) return fChannelPositions;

  const G4double channelRadius = 0.5*fChannelDiameter;
  const G4double safeZone = fEdgeMarginFraction*fPlateRadius;
  const G4double usableRadius = fPlateRadius - channelRadius - safeZone;
  const G4double rowSpacing = fChannelPitch*std::sqrt(3.0)/2.0;

  if (usableRadius <= 0.0) return fChannelPositions;

  G4int row = 0;

  for (G4double y = -usableRadius; y <= usableRadius; y += rowSpacing) {
    G4double xMax = std::sqrt(usableRadius*usableRadius - y*y);
    G4double xOffset = (row % 2 == 0) ? 0.0 : 0.5*fChannelPitch;

    for (G4double x = -xMax; x <= xMax; x += fChannelPitch) {
      G4double xpos = x + xOffset;

      if (xpos*xpos + y*y <= usableRadius*usableRadius) {
        fChannelPositions.push_back(G4ThreeVector(xpos, y, 0.0));
      }
    }

    row++;
  }

  return fChannelPositions;
}

void DetectorConstruction::PrintConfiguration() const
{
  G4cout << "--- MCP configuration ---" << G4endl
         << "  stacks (+z / -z)    : " << fStackCountPlus << " / "
         << fStackCountMinus << G4endl
         << "  first plate centre  : " << fFirstPlateZ/cm << " cm" << G4endl
         << "  plate thickness     : " << fPlateThickness/um << " um" << G4endl
         << "  plate radius        : " << fPlateRadius/mm << " mm" << G4endl
         << "  gap between plates  : " << fStackGap/um << " um" << G4endl
         << "  channel diameter    : " << fChannelDiameter/um << " um" << G4endl
         << "  channel pitch       : " << fChannelPitch/um << " um" << G4endl
         << "  channel angle       : " << fChannelAngle/deg << " deg"
         << (fChevron ? " (chevron)" : "") << G4endl
         << "  edge margin         : " << fEdgeMarginFraction << " x radius"
         << G4endl
         << "  production cut      : "
         << fMcpProductionCuts->GetProductionCut(0)/um << " um" << G4endl
         << "  kill gamma on exit  : " << (fKillGammaOnExit ? "yes" : "no")
         << G4endl
         << "-------------------------" << G4endl;
}
