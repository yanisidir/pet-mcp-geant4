#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
const Double_t kMcpLengthMm = 3.0;
const Double_t kMcpGapMm = 1.0;
const Double_t kFirstMcpZMm = 10.0;
const Double_t kChannelAngleDeg = 8.0;
const Double_t kChannelDiameterUm = 25.0;
const Double_t kChannelRadiusUm = 12.5;
const Double_t kChannelPitchUm = 31.0;
const Double_t kCstMinimalMcpLengthUm = 3000.0;
const Double_t kCstMinimalChannelXUm = 0.14 * 3000.0 * 0.5;
const Double_t kCstMinimalChannelYUm = 0.0;

struct ElectronSeed
{
  Int_t eventID;
  Int_t trackID;
  Int_t side;
  Int_t mcpIndex;
  Double_t x_mm;
  Double_t y_mm;
  Double_t z_mm;
  Double_t dirX;
  Double_t dirY;
  Double_t dirZ;
  Double_t kineticEnergy_keV;
  Double_t time_ns;
  Double_t creationX_mm;
  Double_t creationY_mm;
  Double_t creationZ_mm;
  Double_t creationTime_ns;
  std::string creatorProcessName;
};

Double_t NaN()
{
  return std::numeric_limits<Double_t>::quiet_NaN();
}

// Returns the center z position of the MCP in mm, given the side (+1 or -1)
// and the MCP index (0, 1, 2, ...).
Double_t McpCenterZMm(Int_t side, Int_t mcpIndex)
{
  const Double_t distanceFromSource =
    kFirstMcpZMm + mcpIndex*(kMcpLengthMm + kMcpGapMm);
  return side*distanceFromSource;
}

// Returns the local z position in um for CST, given the electron seed.
Double_t CstFoldedLocalZUm(const ElectronSeed& seed, Bool_t warn = true)
{
  const Double_t centerZ_mm = McpCenterZMm(seed.side, seed.mcpIndex);
  const Double_t localOutward_mm = seed.side*(seed.z_mm - centerZ_mm);
  const Double_t cstZ_um =
    (0.5 - localOutward_mm/kMcpLengthMm)*kCstMinimalMcpLengthUm;

  if (warn && (cstZ_um < 0.0 || cstZ_um > kCstMinimalMcpLengthUm)) {
    std::cerr << "Warning: CST local z is outside [0, "
              << kCstMinimalMcpLengthUm << "] um"
              << " | eventID=" << seed.eventID
              << " trackID=" << seed.trackID
              << " side=" << seed.side
              << " z_mm=" << seed.z_mm
              << " cstZ_um=" << cstZ_um
              << std::endl;
  }

  return cstZ_um;
}

// Returns the local x position in um for CST, given the local z position in um.
Double_t CstMinimalChannelAxisXUm(Double_t cstZ_um)
{
  const Double_t angleRad = kChannelAngleDeg*3.14159265358979323846/180.0;
  return kCstMinimalChannelXUm +
    (cstZ_um - 0.5*kCstMinimalMcpLengthUm)*std::tan(angleRad);
}

// Returns the normalized direction vector for CST, given the electron seed.
Bool_t CstFoldedDirection(const ElectronSeed& seed,
                          Double_t& dx,
                          Double_t& dy,
                          Double_t& dz,
                          Bool_t warn = true)
{
  dx = seed.dirX;
  dy = seed.dirY;
  dz = -seed.side*seed.dirZ;

  Double_t norm = std::sqrt(dx*dx + dy*dy + dz*dz);
  if (norm <= 0.0) {
    if (warn) {
      std::cerr << "Warning: zero direction before CST normalization"
                << " | eventID=" << seed.eventID
                << " trackID=" << seed.trackID
                << " side=" << seed.side
                << std::endl;
    }
    return false;
  }

  dx /= norm;
  dy /= norm;
  dz /= norm;

  if (dz > 0.0) {
    if (warn) {
      std::cerr << "Warning: CST direction had dz > 0; flipping direction"
                << " | eventID=" << seed.eventID
                << " trackID=" << seed.trackID
                << " side=" << seed.side
                << std::endl;
    }
    dx = -dx;
    dy = -dy;
    dz = -dz;
  }

  return true;
}

// Validates the electron seed for CST export, printing warnings for any issues.
Bool_t ValidateSeedForExport(const ElectronSeed& seed)
{
  Bool_t valid = true;

  if (seed.side != 1 && seed.side != -1) {
    std::cerr << "Warning: invalid side for CST export"
              << " | eventID=" << seed.eventID
              << " trackID=" << seed.trackID
              << " side=" << seed.side
              << std::endl;
    valid = false;
  }

  if (seed.mcpIndex < 0) {
    std::cerr << "Warning: invalid MCP index for CST export"
              << " | eventID=" << seed.eventID
              << " trackID=" << seed.trackID
              << " mcpIndex=" << seed.mcpIndex
              << std::endl;
    valid = false;
  }

  if (seed.kineticEnergy_keV <= 0.0) {
    std::cerr << "Warning: non-positive electron kinetic energy"
              << " | eventID=" << seed.eventID
              << " trackID=" << seed.trackID
              << " E_keV=" << seed.kineticEnergy_keV
              << std::endl;
    valid = false;
  }

  const Double_t directionNorm =
    std::sqrt(seed.dirX*seed.dirX +
              seed.dirY*seed.dirY +
              seed.dirZ*seed.dirZ);
  if (directionNorm <= 0.0) {
    std::cerr << "Warning: zero Geant4 direction for CST export"
              << " | eventID=" << seed.eventID
              << " trackID=" << seed.trackID
              << std::endl;
    valid = false;
  }

  CstFoldedLocalZUm(seed, true);

  Double_t cstDirX = 0.0;
  Double_t cstDirY = 0.0;
  Double_t cstDirZ = 0.0;
  if (!CstFoldedDirection(seed, cstDirX, cstDirY, cstDirZ, true)) {
    valid = false;
  }

  return valid;
}

// Finds the first available leaf from a list of candidate names in the given tree.
TLeaf* FindLeaf(TTree* tree,
                const std::vector<std::string>& candidates,
                std::string& selectedName)
{
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    if (tree && tree->GetBranch(candidates[i].c_str())) {
      selectedName = candidates[i];
      return tree->GetLeaf(candidates[i].c_str());
    }
  }

  selectedName = "";
  return nullptr;
}

// Requires that the specified branch exists in the tree and returns its leaf.
TLeaf* RequireLeaf(TTree* tree, const std::string& branchName)
{
  if (!tree || !tree->GetBranch(branchName.c_str())) {
    std::cerr << "Missing required branch: " << branchName << std::endl;
    return nullptr;
  }
  return tree->GetLeaf(branchName.c_str());
}

// Finds the first event in the event tree that has isCoincidence == true.
bool FindFirstCoincidenceEvent(TTree* eventTree, Int_t& coincidenceEventID)
{
  if (!eventTree ||
      !eventTree->GetBranch("eventID") ||
      !eventTree->GetBranch("isCoincidence")) {
    std::cerr << "EventSummaryTree must contain eventID and isCoincidence."
              << std::endl;
    return false;
  }

  Int_t eventID = -1;
  Bool_t isCoincidence = false;
  eventTree->SetBranchAddress("eventID", &eventID);
  eventTree->SetBranchAddress("isCoincidence", &isCoincidence);

  const Long64_t entries = eventTree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    eventTree->GetEntry(entry);
    if (isCoincidence) {
      coincidenceEventID = eventID;
      eventTree->ResetBranchAddresses();
      return true;
    }
  }

  eventTree->ResetBranchAddresses();
  return false;
}

// Prints the details of the electron seed to standard output.
void PrintSeed(const ElectronSeed& seed)
{
  Double_t cstDirX = 0.0;
  Double_t cstDirY = 0.0;
  Double_t cstDirZ = 0.0;
  CstFoldedDirection(seed, cstDirX, cstDirY, cstDirZ, false);
  const Double_t cstZ_um = CstFoldedLocalZUm(seed, false);
  const Double_t cstX_um = CstMinimalChannelAxisXUm(cstZ_um);

  std::cout << "  side=" << (seed.side > 0 ? "+z" : "-z")
            << " trackID=" << seed.trackID
            << " mcpIndex=" << seed.mcpIndex
            << "\n    Geant4 position [mm] = ("
            << seed.x_mm << ", " << seed.y_mm << ", " << seed.z_mm << ")"
            << "\n    CST position [um]   = ("
            << cstX_um << ", "
            << kCstMinimalChannelYUm << ", "
            << cstZ_um << ")"
            << "\n    Geant4 direction    = ("
            << seed.dirX << ", " << seed.dirY << ", " << seed.dirZ << ")"
            << "\n    CST direction       = ("
            << cstDirX << ", " << cstDirY << ", " << cstDirZ << ")"
            << "\n    kinetic energy      = "
            << 1000.0*seed.kineticEnergy_keV << " eV"
            << "\n    creator process     = "
            << seed.creatorProcessName
            << std::endl;
}

// Writes the CSV header line to the output stream.
void WriteCsvHeader(std::ofstream& output)
{
  output
    << "eventID,trackID,side,mcpIndex,"
    << "cstLocalX_um,cstLocalY_um,cstLocalZ_um,"
    << "dirX,dirY,dirZ,"
    << "kineticEnergy_eV,time_ns,"
    << "channelAngle_deg,channelDiameter_um,channelRadius_um,channelPitch_um,"
    << "geant4X_mm,geant4Y_mm,geant4Z_mm,"
    << "geant4DirX,geant4DirY,geant4DirZ,"
    << "creationX_mm,creationY_mm,creationZ_mm,creationTime_ns,"
    << "creatorProcessName\n";
}

// Writes a single electron seed row to the CSV output steam.
void WriteCsvRow(std::ofstream& output, const ElectronSeed& seed)
{
  Double_t cstDirX = 0.0;
  Double_t cstDirY = 0.0;
  Double_t cstDirZ = 0.0;
  CstFoldedDirection(seed, cstDirX, cstDirY, cstDirZ, false);

  const Double_t cstZ_um = CstFoldedLocalZUm(seed, false);
  const Double_t cstX_um = CstMinimalChannelAxisXUm(cstZ_um);

  // CST script must use CHANNEL_RADIUS_UM = 12.5 to match this export.
  output << seed.eventID << ","
         << seed.trackID << ","
         << seed.side << ","
         << seed.mcpIndex << ","
         << cstX_um << ","
         << kCstMinimalChannelYUm << ","
         << cstZ_um << ","
         << cstDirX << ","
         << cstDirY << ","
         << cstDirZ << ","
         << 1000.0*seed.kineticEnergy_keV << ","
         << seed.time_ns << ","
         << kChannelAngleDeg << ","
         << kChannelDiameterUm << ","
         << kChannelRadiusUm << ","
         << kChannelPitchUm << ","
         << seed.x_mm << ","
         << seed.y_mm << ","
         << seed.z_mm << ","
         << seed.dirX << ","
         << seed.dirY << ","
         << seed.dirZ << ","
         << seed.creationX_mm << ","
         << seed.creationY_mm << ","
         << seed.creationZ_mm << ","
         << seed.creationTime_ns << ","
         << seed.creatorProcessName
         << "\n";
}

} // namespace


void ExportCstSeedElectrons(
  const char* fileName = "mcp_output.root",
  const char* outputName = "geant4_seed_electrons.csv")
{
  TFile* file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open ROOT file: " << fileName << std::endl;
    delete file;
    return;
  }

  TTree* eventTree = nullptr;
  TTree* electronTree = nullptr;
  file->GetObject("EventSummaryTree", eventTree);
  file->GetObject("ElectronChannelHitTree", electronTree);

  if (!eventTree || !electronTree) {
    std::cerr << "Missing EventSummaryTree or ElectronChannelHitTree in "
              << fileName << std::endl;
    file->Close();
    delete file;
    return;
  }

  Int_t selectedEventID = -1;
  if (!FindFirstCoincidenceEvent(eventTree, selectedEventID)) {
    std::cerr << "No event with isCoincidence == true was found."
              << std::endl;
    file->Close();
    delete file;
    return;
  }

  TLeaf* eventIDLeaf = RequireLeaf(electronTree, "eventID");
  TLeaf* trackIDLeaf = RequireLeaf(electronTree, "trackID");
  TLeaf* sideLeaf = RequireLeaf(electronTree, "side");
  TLeaf* mcpIndexLeaf = RequireLeaf(electronTree, "mcpIndex");
  TLeaf* energyLeaf = RequireLeaf(electronTree, "kineticEnergy_keV");
  TLeaf* timeLeaf = RequireLeaf(electronTree, "globalTime_ns");

  std::string xName;
  std::string yName;
  std::string zName;
  std::string dirXName;
  std::string dirYName;
  std::string dirZName;
  std::string creationXName;
  std::string creationYName;
  std::string creationZName;

  TLeaf* xLeaf = FindLeaf(electronTree, {"x_mm", "x", "posX"}, xName);
  TLeaf* yLeaf = FindLeaf(electronTree, {"y_mm", "y", "posY"}, yName);
  TLeaf* zLeaf = FindLeaf(electronTree, {"z_mm", "z", "posZ"}, zName);
  TLeaf* dirXLeaf = FindLeaf(electronTree, {"dirX", "px", "momX"}, dirXName);
  TLeaf* dirYLeaf = FindLeaf(electronTree, {"dirY", "py", "momY"}, dirYName);
  TLeaf* dirZLeaf = FindLeaf(electronTree, {"dirZ", "pz", "momZ"}, dirZName);

  TLeaf* creationXLeaf =
    FindLeaf(electronTree, {"creationX_mm", "creationX", "vertexX"}, creationXName);
  TLeaf* creationYLeaf =
    FindLeaf(electronTree, {"creationY_mm", "creationY", "vertexY"}, creationYName);
  TLeaf* creationZLeaf =
    FindLeaf(electronTree, {"creationZ_mm", "creationZ", "vertexZ"}, creationZName);
  TLeaf* creationTimeLeaf =
    electronTree->GetBranch("creationTime_ns")
      ? electronTree->GetLeaf("creationTime_ns")
      : nullptr;

  if (!eventIDLeaf || !trackIDLeaf || !sideLeaf || !mcpIndexLeaf ||
      !energyLeaf || !timeLeaf ||
      !xLeaf || !yLeaf || !zLeaf ||
      !dirXLeaf || !dirYLeaf || !dirZLeaf) {
    std::cerr << "ElectronChannelHitTree does not contain all required "
              << "branches for CST export." << std::endl;
    file->Close();
    delete file;
    return;
  }

  char creatorProcessName[128] = "unknown";
  const Bool_t hasCreatorProcessBranch =
    electronTree->GetBranch("creatorProcessName") != nullptr;
  if (hasCreatorProcessBranch) {
    electronTree->SetBranchAddress("creatorProcessName", creatorProcessName);
  }

  ElectronSeed plusSeed;
  ElectronSeed minusSeed;
  bool hasPlusSeed = false;
  bool hasMinusSeed = false;
  Int_t plusCount = 0;
  Int_t minusCount = 0;

  const Long64_t electronEntries = electronTree->GetEntries();
  // Loop over all electron entries to find the first valid seed for each side.
  for (Long64_t entry = 0; entry < electronEntries; ++entry) {
    std::snprintf(creatorProcessName,
                  sizeof(creatorProcessName),
                  "unknown");

    electronTree->GetEntry(entry);

    const Int_t eventID = static_cast<Int_t>(eventIDLeaf->GetValue());
    if (eventID != selectedEventID) {
      continue;
    }

    const Int_t side = static_cast<Int_t>(sideLeaf->GetValue());
    if (side > 0) {
      ++plusCount;
    } else if (side < 0) {
      ++minusCount;
    }

    if ((side > 0 && hasPlusSeed) ||
        (side < 0 && hasMinusSeed) ||
        side == 0) {
      continue;
    }

    ElectronSeed seed;
    seed.eventID = eventID;
    seed.trackID = static_cast<Int_t>(trackIDLeaf->GetValue());
    seed.side = side;
    seed.mcpIndex = static_cast<Int_t>(mcpIndexLeaf->GetValue());
    seed.x_mm = xLeaf->GetValue();
    seed.y_mm = yLeaf->GetValue();
    seed.z_mm = zLeaf->GetValue();
    seed.dirX = dirXLeaf->GetValue();
    seed.dirY = dirYLeaf->GetValue();
    seed.dirZ = dirZLeaf->GetValue();
    seed.kineticEnergy_keV = energyLeaf->GetValue();
    seed.time_ns = timeLeaf->GetValue();
    seed.creationX_mm = creationXLeaf ? creationXLeaf->GetValue() : NaN();
    seed.creationY_mm = creationYLeaf ? creationYLeaf->GetValue() : NaN();
    seed.creationZ_mm = creationZLeaf ? creationZLeaf->GetValue() : NaN();
    seed.creationTime_ns =
      creationTimeLeaf ? creationTimeLeaf->GetValue() : NaN();
    seed.creatorProcessName =
      hasCreatorProcessBranch && creatorProcessName[0] != '\0'
        ? creatorProcessName
        : "unknown";

    if (!ValidateSeedForExport(seed)) {
      continue;
    }

    if (side > 0) {
      plusSeed = seed;
      hasPlusSeed = true;
    } else if (side < 0) {
      minusSeed = seed;
      hasMinusSeed = true;
    }
  }

  electronTree->ResetBranchAddresses();

  if (!hasPlusSeed || !hasMinusSeed) {
    std::cerr << "Selected event " << selectedEventID
              << " is marked as coincidence, but one side could not be "
              << "retrieved from ElectronChannelHitTree." << std::endl;
    file->Close();
    delete file;
    return;
  }

  std::ofstream output(outputName);
  if (!output) {
    std::cerr << "Cannot create CSV file: " << outputName << std::endl;
    file->Close();
    delete file;
    return;
  }

  output << std::setprecision(12);
  WriteCsvHeader(output);
  WriteCsvRow(output, plusSeed);
  WriteCsvRow(output, minusSeed);
  output.close();

  std::cout << "\n=== CST seed electron export ===" << std::endl;
  std::cout << "ROOT file: " << fileName << std::endl;
  std::cout << "CSV file : " << outputName << std::endl;
  std::cout << "Selected coincidence eventID: "
            << selectedEventID << std::endl;
  std::cout << "Detected channel electrons in selected event:"
            << " +z=" << plusCount
            << " -z=" << minusCount << std::endl;
  std::cout << "Number of exported seeds: 2" << std::endl;
  std::cout << "\nExported electrons:" << std::endl;
  PrintSeed(plusSeed);
  PrintSeed(minusSeed);

  file->Close();
  delete file;
}
