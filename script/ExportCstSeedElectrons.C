#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const Int_t kNumberOfMCPs = 3;
const Double_t kMcpLengthMm = 3.0;
const Double_t kMcpGapMm = 1.0;
const Double_t kFirstMcpZMm = 10.0;
const Double_t kChannelAngleDeg = 8.0;
const Double_t kChannelDiameterUm = 25.0;
const Double_t kChannelPitchUm = 31.0;
const Double_t kCstReducedRadiusMm = 1.0;
const Double_t kCstMinimalMcpLengthUm = 500.0;
const Double_t kCstMinimalChannelXUm = 35.0;
const Double_t kCstMinimalChannelYUm = 0.0;

struct ElectronSeed
{
  Int_t eventID;
  Int_t trackID;
  Int_t parentGammaTrackID;
  Int_t primaryGammaTrackID;
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
  Bool_t hasValidCreationInfo;
};

std::string StackName(Int_t side)
{
  return side > 0 ? "+z" : "-z";
}

Double_t McpCenterZMm(Int_t side, Int_t mcpIndex)
{
  const Double_t distanceFromSource =
    kFirstMcpZMm + mcpIndex*(kMcpLengthMm + kMcpGapMm);
  return side*distanceFromSource;
}

Double_t ChannelAngleDeg(Int_t side, Int_t mcpIndex)
{
  const Double_t chevronAngle =
    (mcpIndex % 2 == 0) ? kChannelAngleDeg : -kChannelAngleDeg;
  return side*chevronAngle;
}

Double_t CstLocalZMm(const ElectronSeed& seed)
{
  return seed.z_mm - McpCenterZMm(seed.side, seed.mcpIndex);
}

Double_t CstMinimalLocalZUm(const ElectronSeed& seed)
{
  const Double_t normalizedZ =
    (CstLocalZMm(seed) + 0.5*kMcpLengthMm) / kMcpLengthMm;
  return normalizedZ*kCstMinimalMcpLengthUm;
}

Double_t CstMinimalChannelAxisXUm(const ElectronSeed& seed)
{
  const Double_t angleRad = kChannelAngleDeg*3.14159265358979323846/180.0;
  return kCstMinimalChannelXUm +
    (CstMinimalLocalZUm(seed) - 0.5*kCstMinimalMcpLengthUm)*std::sin(angleRad);
}

std::string PlateListAfterSeed(const ElectronSeed& seed)
{
  std::ostringstream plates;
  bool first = true;
  for (Int_t mcpIndex = seed.mcpIndex + 1;
       mcpIndex < kNumberOfMCPs;
       ++mcpIndex) {
    if (!first) {
      plates << ", ";
    }
    plates << "MCP " << mcpIndex;
    first = false;
  }

  return first ? "none" : plates.str();
}

void WriteGeometryBlock(std::ofstream& output,
                        const ElectronSeed& seed)
{
  const Double_t channelAngle = ChannelAngleDeg(seed.side, seed.mcpIndex);
  const Double_t mcpCenterZ = McpCenterZMm(seed.side, seed.mcpIndex);
  const Double_t localZ = CstLocalZMm(seed);

  output << "Stack: " << StackName(seed.side) << "\n";
  output << "Selected electron trackID: " << seed.trackID << "\n";
  output << "MCP where electron enters channel: MCP "
         << seed.mcpIndex << "\n";
  output << "MCP center z [mm]: " << mcpCenterZ << "\n";
  output << "Geant4 hit position:\n";
  output << "  global x,y,z [mm] = "
         << seed.x_mm << ", "
         << seed.y_mm << ", "
         << seed.z_mm << "\n";
  output << "CST local coordinate convention:\n";
  output << "  touched channel center is taken as local x = 0, y = 0\n";
  output << "  electron entry point is injected in the local channel frame\n";
  output << "  local z is kept relative to the MCP center\n";
  output << "  cst local x,y,z [mm] = 0, 0, "
         << localZ << "\n";
  output << "  note: do not place the reduced CST channel at Geant4 x,y; "
         << "Geant4 global coordinates are kept only for traceability\n";
  output << "Channel orientation:\n";
  output << "  rotation axis convention: channels are tilted by rotateY(angle)\n";
  output << "  channel angle [deg]: " << channelAngle << "\n";
  output << "  approximate channel axis direction: "
         << "(sin(angle), 0, cos(angle)); use both signs if CST treats "
         << "the cylinder axis as unoriented\n";
  output << "Electron direction at channel entry:\n";
  output << "  dirX,dirY,dirZ: "
         << seed.dirX << ", "
         << seed.dirY << ", "
         << seed.dirZ << "\n";
  output << "Plates to draw in CST:\n";
  output << "  local detailed geometry: MCP " << seed.mcpIndex
         << ", only the touched channel plus a small local glass environment\n";
  output << "  downstream plates in stack direction: "
         << PlateListAfterSeed(seed) << "\n";
  output << "  recommendation: keep a fuller or corridor-like channel geometry "
         << "only for these downstream plates, not the full 3+3 detector\n";
  output << "\n";
}

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

TLeaf* RequireLeaf(TTree* tree, const std::string& branchName)
{
  if (!tree || !tree->GetBranch(branchName.c_str())) {
    std::cerr << "Missing required branch: " << branchName << std::endl;
    return nullptr;
  }
  return tree->GetLeaf(branchName.c_str());
}

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

void PrintSeed(const ElectronSeed& seed)
{
  std::cout << "  side=" << (seed.side > 0 ? "+z" : "-z")
            << " trackID=" << seed.trackID
            << " mcpIndex=" << seed.mcpIndex
            << " E=" << seed.kineticEnergy_keV << " keV"
            << " t=" << seed.time_ns << " ns"
            << " pos=(" << seed.x_mm << ", "
            << seed.y_mm << ", "
            << seed.z_mm << ") mm"
            << " dir=(" << seed.dirX << ", "
            << seed.dirY << ", "
            << seed.dirZ << ")"
            << " creator=" << seed.creatorProcessName
            << " validCreation="
            << (seed.hasValidCreationInfo ? "true" : "false")
            << std::endl;
}

void WriteCsvRow(std::ofstream& output, const ElectronSeed& seed)
{
  output << seed.eventID << ","
         << seed.trackID << ","
         << CstMinimalChannelAxisXUm(seed) << ","
         << kCstMinimalChannelYUm << ","
         << CstMinimalLocalZUm(seed) << ","
         << seed.dirX << ","
         << seed.dirY << ","
         << seed.dirZ << ","
         << 1000.0*seed.kineticEnergy_keV << ","
         << seed.time_ns
         << "\n";
}

bool WriteGeometryRequest(const char* geometryRequestName,
                          const ElectronSeed& plusSeed,
                          const ElectronSeed& minusSeed,
                          Int_t plusCount,
                          Int_t minusCount)
{
  std::ofstream output(geometryRequestName);
  if (!output) {
    std::cerr << "Cannot create geometry request file: "
              << geometryRequestName << std::endl;
    return false;
  }

  output << std::setprecision(12);
  output << "CST geometry request for selected Geant4 coincidence event\n";
  output << "=========================================================\n\n";
  output << "Purpose\n";
  output << "-------\n";
  output << "Prepare a reduced CST geometry seeded by two Geant4 electrons: "
         << "one entering an MCP_channel on +z and one entering an MCP_channel "
         << "on -z.\n\n";

  output << "1. Selected Geant4 event information\n";
  output << "------------------------------------\n";
  output << "eventID: " << plusSeed.eventID << "\n";
  output << "detected channel electrons in event: +z="
         << plusCount << ", -z=" << minusCount << "\n";
  output << "Geant4 full detector note: the full Geant4 geometry contains "
         << "about 1.5 million parameterised channels per MCP. It must not "
         << "be reproduced directly in CST.\n\n";

  output << "Geant4 reference geometry used to interpret the hit\n";
  output << "--------------------------------------------------\n";
  output << "Number of MCPs per stack: " << kNumberOfMCPs << "\n";
  output << "MCP length [mm]: " << kMcpLengthMm << "\n";
  output << "MCP gap [mm]: " << kMcpGapMm << "\n";
  output << "First MCP center distance from source [mm]: "
         << kFirstMcpZMm << "\n";
  output << "Channel diameter [um]: " << kChannelDiameterUm << "\n";
  output << "Channel pitch [um]: " << kChannelPitchUm << "\n";
  output << "Chevron channel angle magnitude [deg]: "
         << kChannelAngleDeg << "\n";
  output << "Electric field in this Geant4 export configuration: disabled\n";
  output << "MCP material in this Geant4 export configuration: GlassLead60Perc\n";
  output << "Primary source: two back-to-back 511 keV photons\n\n";

  output << "2. Reduced CST geometry parameters\n";
  output << "----------------------------------\n";
  output << "CST geometry type: local/intermediate reduced MCP model\n";
  output << "Number of MCPs per stack to draw: " << kNumberOfMCPs << "\n";
  output << "MCP length [mm]: " << kMcpLengthMm << "\n";
  output << "MCP gap [mm]: " << kMcpGapMm << "\n";
  output << "Channel diameter [um]: " << kChannelDiameterUm << "\n";
  output << "Channel pitch [um]: " << kChannelPitchUm << "\n";
  output << "Channel angle [deg]: +/-" << kChannelAngleDeg << "\n";
  output << "Reduced transverse radius around touched channel [mm]: "
         << kCstReducedRadiusMm << "\n";
  output << "Channel selection rule: draw a limited set of channels around "
         << "local x=0, y=0 using the pitch above; do not draw the full "
         << "Geant4 channel population.\n\n";

  output << "Important limitation\n";
  output << "--------------------\n";
  output << "ElectronChannelHitTree stores the electron entry point into a "
         << "channel, but not the Geant4 parameterised channel copy number. "
         << "Therefore CST local x,y are set to 0,0 at the touched channel "
         << "center. Geant4 global x,y,z are kept separately for "
         << "traceability only.\n\n";

  output << "+z seed geometry\n";
  output << "----------------\n";
  WriteGeometryBlock(output, plusSeed);

  output << "-z seed geometry\n";
  output << "----------------\n";
  WriteGeometryBlock(output, minusSeed);

  output << "Reduced CST drawing strategy\n";
  output << "----------------------------\n";
  output << "1. Do not draw the full 3+3 MCP detector for this seed study.\n";
  output << "2. For the MCP where the selected electron enters, draw only:\n";
  output << "   - the touched channel centered near the listed entry point;\n";
  output << "   - a small local lead-glass environment around that channel;\n";
  output << "   - the channel orientation listed above.\n";
  output << "3. For MCPs downstream of that entry point in the stack direction, "
         << "draw either:\n";
  output << "   - a fuller local corridor of channels that the particle may cross; "
         << "or\n";
  output << "   - a simplified volume sufficient for transport studies.\n";
  output << "4. Avoid drawing unrelated plates from the opposite stack unless "
         << "the selected seed electron belongs to that stack.\n";
  output << "5. Use geant4_seed_electrons.csv as the particle initial condition file.\n";

  return true;
}
}

void ExportCstSeedElectrons(
  const char* fileName = "build/mcp_output.root",
  const char* outputName = "geant4_seed_electrons.csv",
  const char* geometryRequestName = "cst_geometry_request.txt")
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
  TLeaf* parentGammaLeaf = RequireLeaf(electronTree, "parentGammaTrackID");
  TLeaf* primaryGammaLeaf = RequireLeaf(electronTree, "primaryGammaTrackID");
  TLeaf* energyLeaf = RequireLeaf(electronTree, "kineticEnergy_keV");
  TLeaf* timeLeaf = RequireLeaf(electronTree, "globalTime_ns");
  TLeaf* validCreationLeaf =
    RequireLeaf(electronTree, "hasValidCreationInfo");
  TLeaf* creationTimeLeaf = RequireLeaf(electronTree, "creationTime_ns");

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

  if (!eventIDLeaf || !trackIDLeaf || !sideLeaf || !mcpIndexLeaf ||
      !parentGammaLeaf || !primaryGammaLeaf || !energyLeaf ||
      !timeLeaf || !validCreationLeaf || !creationTimeLeaf ||
      !xLeaf || !yLeaf || !zLeaf ||
      !dirXLeaf || !dirYLeaf || !dirZLeaf ||
      !creationXLeaf || !creationYLeaf || !creationZLeaf ||
      !electronTree->GetBranch("creatorProcessName")) {
    std::cerr << "ElectronChannelHitTree does not contain all branches "
              << "needed for CST export." << std::endl;
    file->Close();
    delete file;
    return;
  }

  char creatorProcessName[128] = "";
  electronTree->SetBranchAddress("creatorProcessName", creatorProcessName);

  ElectronSeed plusSeed;
  ElectronSeed minusSeed;
  bool hasPlusSeed = false;
  bool hasMinusSeed = false;
  Int_t plusCount = 0;
  Int_t minusCount = 0;

  const Long64_t electronEntries = electronTree->GetEntries();
  for (Long64_t entry = 0; entry < electronEntries; ++entry) {
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
    seed.parentGammaTrackID =
      static_cast<Int_t>(parentGammaLeaf->GetValue());
    seed.primaryGammaTrackID =
      static_cast<Int_t>(primaryGammaLeaf->GetValue());
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
    seed.creationX_mm = creationXLeaf->GetValue();
    seed.creationY_mm = creationYLeaf->GetValue();
    seed.creationZ_mm = creationZLeaf->GetValue();
    seed.creationTime_ns = creationTimeLeaf->GetValue();
    seed.creatorProcessName = creatorProcessName;
    seed.hasValidCreationInfo =
      static_cast<Bool_t>(validCreationLeaf->GetValue());

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
  output
    << "eventID,trackID,"
    << "cstLocalX_um,cstLocalY_um,cstLocalZ_um,"
    << "dirX,dirY,dirZ,"
    << "kineticEnergy_eV,time_ns\n";
  WriteCsvRow(output, plusSeed);
  WriteCsvRow(output, minusSeed);
  output.close();

  const bool wroteGeometryRequest =
    WriteGeometryRequest(geometryRequestName,
                         plusSeed,
                         minusSeed,
                         plusCount,
                         minusCount);

  std::cout << "\n=== CST seed electron export ===" << std::endl;
  std::cout << "ROOT file: " << fileName << std::endl;
  std::cout << "CSV file : " << outputName << std::endl;
  if (wroteGeometryRequest) {
    std::cout << "Geometry request: "
              << geometryRequestName << std::endl;
  }
  std::cout << "Selected coincidence eventID: "
            << selectedEventID << std::endl;
  std::cout << "Detected channel electrons in selected event:"
            << " +z=" << plusCount
            << " -z=" << minusCount << std::endl;
  std::cout << "\nExported electrons:" << std::endl;
  PrintSeed(plusSeed);
  PrintSeed(minusSeed);

  file->Close();
  delete file;
}
