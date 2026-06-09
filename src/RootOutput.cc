#include "RootOutput.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "TFile.h"
#include "TTree.h"

#include <cstring>
#include <sstream>
#include <string>

namespace // Copie sécurisée d'une G4String dans un tableau de caractères C.
{
  void CopyText(char* destination,
                std::size_t destinationSize,
                const G4String& source)
  {
    std::strncpy(destination, source.c_str(), destinationSize);
    destination[destinationSize - 1] = '\0';
  }
}

G4ThreadLocal RootOutput* RootOutput::fInstance = nullptr;

RootOutput* RootOutput::Instance()
{
  if (!fInstance) {
    fInstance = new RootOutput;
  }
  return fInstance;
}

RootOutput::RootOutput()
  : fFile(nullptr),
    fElectronChannelHitTree(nullptr),
    fPhotonExitTree(nullptr),
    fEventSummaryTree(nullptr),
    fElectronEventID(0),
    fElectronTrackID(0),
    fElectronParentID(0),
    fElectronSide(0),
    fElectronMcpIndex(-1),
    fElectronKineticEnergy(0.0),
    fElectronGlobalTime(0.0),
    fElectronX(0.0),
    fElectronY(0.0),
    fElectronZ(0.0),
    fElectronDirX(0.0),
    fElectronDirY(0.0),
    fElectronDirZ(0.0),
    fPhotonEventID(0),
    fPhotonTrackID(0),
    fPhotonParentID(0),
    fPhotonSide(0),
    fPhotonKineticEnergy(0.0),
    fPhotonGlobalTime(0.0),
    fPhotonX(0.0),
    fPhotonY(0.0),
    fPhotonZ(0.0),
    fPhotonDirX(0.0),
    fPhotonDirY(0.0),
    fPhotonDirZ(0.0),
    fSummaryEventID(0),
    fSummaryElectronProducedCount(0),
    fSummaryElectronChannelCount(0),
    fSummaryElectronChannelPlusCount(0),
    fSummaryElectronChannelMinusCount(0),
    fSummaryElectronChannelPlusMcp0Count(0),
    fSummaryElectronChannelPlusMcp1Count(0),
    fSummaryElectronChannelPlusMcp2Count(0),
    fSummaryElectronChannelMinusMcp0Count(0),
    fSummaryElectronChannelMinusMcp1Count(0),
    fSummaryElectronChannelMinusMcp2Count(0),
    fSummaryHasElectronChannelPlus(false),
    fSummaryHasElectronChannelMinus(false),
    fSummaryIsCoincidence(false),
    fSummaryPhotonExitCount(0),
    fSummaryPhotonExitPlusCount(0),
    fSummaryPhotonExitMinusCount(0)
{
  fElectronVolumeName[0] = '\0';
  fElectronCreatorProcessName[0] = '\0';
  fPhotonVolumeName[0] = '\0';
  fPhotonStepProcessName[0] = '\0';
}

RootOutput::~RootOutput()
{
  Close();
}

void RootOutput::Open(const char* fileName)
{
  Close();

  const std::string threadFileName = BuildThreadFileName(fileName);
  fFile = new TFile(threadFileName.c_str(), "RECREATE");

  if (!fFile || fFile->IsZombie()) {
    G4ExceptionDescription message;
    message << "Cannot create ROOT file: " << threadFileName;
    delete fFile;
    fFile = nullptr;
    G4Exception("RootOutput::Open",
                "RootOutput001",
                FatalException,
                message);
    return;
  }

  fElectronChannelHitTree =
    new TTree("ElectronChannelHitTree",
              "Electrons entering MCP_channel");
  fPhotonExitTree =
    new TTree("PhotonExitTree",
              "Gamma photons leaving the MCP");
  fEventSummaryTree =
    new TTree("EventSummaryTree",
              "One minimal row per event");

  fElectronChannelHitTree->Branch(
    "eventID", &fElectronEventID, "eventID/I");
  fElectronChannelHitTree->Branch(
    "trackID", &fElectronTrackID, "trackID/I");
  fElectronChannelHitTree->Branch(
    "parentID", &fElectronParentID, "parentID/I");
  fElectronChannelHitTree->Branch(
    "side", &fElectronSide, "side/I");
  fElectronChannelHitTree->Branch(
    "mcpIndex", &fElectronMcpIndex, "mcpIndex/I");
  fElectronChannelHitTree->Branch(
    "kineticEnergy_keV",
    &fElectronKineticEnergy,
    "kineticEnergy_keV/D");
  fElectronChannelHitTree->Branch(
    "globalTime_ns", &fElectronGlobalTime, "globalTime_ns/D");
  fElectronChannelHitTree->Branch(
    "x_cm", &fElectronX, "x_cm/D");
  fElectronChannelHitTree->Branch(
    "y_cm", &fElectronY, "y_cm/D");
  fElectronChannelHitTree->Branch(
    "z_cm", &fElectronZ, "z_cm/D");
  fElectronChannelHitTree->Branch(
    "dirX", &fElectronDirX, "dirX/D");
  fElectronChannelHitTree->Branch(
    "dirY", &fElectronDirY, "dirY/D");
  fElectronChannelHitTree->Branch(
    "dirZ", &fElectronDirZ, "dirZ/D");
  fElectronChannelHitTree->Branch(
    "volumeName", fElectronVolumeName, "volumeName/C");
  fElectronChannelHitTree->Branch(
    "creatorProcessName",
    fElectronCreatorProcessName,
    "creatorProcessName/C");

  fPhotonExitTree->Branch(
    "eventID", &fPhotonEventID, "eventID/I");
  fPhotonExitTree->Branch(
    "trackID", &fPhotonTrackID, "trackID/I");
  fPhotonExitTree->Branch(
    "parentID", &fPhotonParentID, "parentID/I");
  fPhotonExitTree->Branch(
    "side", &fPhotonSide, "side/I");
  fPhotonExitTree->Branch(
    "kineticEnergy_keV",
    &fPhotonKineticEnergy,
    "kineticEnergy_keV/D");
  fPhotonExitTree->Branch(
    "globalTime_ns", &fPhotonGlobalTime, "globalTime_ns/D");
  fPhotonExitTree->Branch(
    "x_cm", &fPhotonX, "x_cm/D");
  fPhotonExitTree->Branch(
    "y_cm", &fPhotonY, "y_cm/D");
  fPhotonExitTree->Branch(
    "z_cm", &fPhotonZ, "z_cm/D");
  fPhotonExitTree->Branch(
    "dirX", &fPhotonDirX, "dirX/D");
  fPhotonExitTree->Branch(
    "dirY", &fPhotonDirY, "dirY/D");
  fPhotonExitTree->Branch(
    "dirZ", &fPhotonDirZ, "dirZ/D");
  fPhotonExitTree->Branch(
    "volumeName", fPhotonVolumeName, "volumeName/C");
  fPhotonExitTree->Branch(
    "stepProcessName",
    fPhotonStepProcessName,
    "stepProcessName/C");

  fEventSummaryTree->Branch(
    "eventID", &fSummaryEventID, "eventID/I");
  fEventSummaryTree->Branch(
    "electronProducedCount",
    &fSummaryElectronProducedCount,
    "electronProducedCount/I");
  fEventSummaryTree->Branch(
    "electronChannelCount",
    &fSummaryElectronChannelCount,
    "electronChannelCount/I");
  fEventSummaryTree->Branch(
    "electronChannelPlusCount",
    &fSummaryElectronChannelPlusCount,
    "electronChannelPlusCount/I");
  fEventSummaryTree->Branch(
    "electronChannelMinusCount",
    &fSummaryElectronChannelMinusCount,
    "electronChannelMinusCount/I");
  fEventSummaryTree->Branch(
    "electronChannelPlusMcp0Count",
    &fSummaryElectronChannelPlusMcp0Count,
    "electronChannelPlusMcp0Count/I");
  fEventSummaryTree->Branch(
    "electronChannelPlusMcp1Count",
    &fSummaryElectronChannelPlusMcp1Count,
    "electronChannelPlusMcp1Count/I");
  fEventSummaryTree->Branch(
    "electronChannelPlusMcp2Count",
    &fSummaryElectronChannelPlusMcp2Count,
    "electronChannelPlusMcp2Count/I");
  fEventSummaryTree->Branch(
    "electronChannelMinusMcp0Count",
    &fSummaryElectronChannelMinusMcp0Count,
    "electronChannelMinusMcp0Count/I");
  fEventSummaryTree->Branch(
    "electronChannelMinusMcp1Count",
    &fSummaryElectronChannelMinusMcp1Count,
    "electronChannelMinusMcp1Count/I");
  fEventSummaryTree->Branch(
    "electronChannelMinusMcp2Count",
    &fSummaryElectronChannelMinusMcp2Count,
    "electronChannelMinusMcp2Count/I");
  fEventSummaryTree->Branch(
    "hasElectronChannelPlus",
    &fSummaryHasElectronChannelPlus,
    "hasElectronChannelPlus/O");
  fEventSummaryTree->Branch(
    "hasElectronChannelMinus",
    &fSummaryHasElectronChannelMinus,
    "hasElectronChannelMinus/O");
  fEventSummaryTree->Branch(
    "isCoincidence",
    &fSummaryIsCoincidence,
    "isCoincidence/O");
  fEventSummaryTree->Branch(
    "photonExitCount",
    &fSummaryPhotonExitCount,
    "photonExitCount/I");
  fEventSummaryTree->Branch(
    "photonExitPlusCount",
    &fSummaryPhotonExitPlusCount,
    "photonExitPlusCount/I");
  fEventSummaryTree->Branch(
    "photonExitMinusCount",
    &fSummaryPhotonExitMinusCount,
    "photonExitMinusCount/I");

  G4cout << "ROOT output opened: " << threadFileName << G4endl;
}

void RootOutput::FillElectronChannelHit(
  const ElectronChannelHitInfo& hit)
{
  if (!fElectronChannelHitTree) {
    return;
  }

  fElectronEventID = hit.eventID;
  fElectronTrackID = hit.trackID;
  fElectronParentID = hit.parentID;
  fElectronSide = hit.side;
  fElectronMcpIndex = hit.mcpIndex;
  fElectronKineticEnergy = hit.kineticEnergy/keV;
  fElectronGlobalTime = hit.globalTime/ns;
  fElectronX = hit.position.x()/cm;
  fElectronY = hit.position.y()/cm;
  fElectronZ = hit.position.z()/cm;
  fElectronDirX = hit.momentumDirection.x();
  fElectronDirY = hit.momentumDirection.y();
  fElectronDirZ = hit.momentumDirection.z();

  CopyText(fElectronVolumeName,
           sizeof(fElectronVolumeName),
           hit.volumeName);
  CopyText(fElectronCreatorProcessName,
           sizeof(fElectronCreatorProcessName),
           hit.creatorProcessName);

  fElectronChannelHitTree->Fill();
}

void RootOutput::FillPhotonExit(const PhotonExitInfo& exitInfo)
{
  if (!fPhotonExitTree) {
    return;
  }

  fPhotonEventID = exitInfo.eventID;
  fPhotonTrackID = exitInfo.trackID;
  fPhotonParentID = exitInfo.parentID;
  fPhotonSide = exitInfo.side;
  fPhotonKineticEnergy = exitInfo.kineticEnergy/keV;
  fPhotonGlobalTime = exitInfo.globalTime/ns;
  fPhotonX = exitInfo.position.x()/cm;
  fPhotonY = exitInfo.position.y()/cm;
  fPhotonZ = exitInfo.position.z()/cm;
  fPhotonDirX = exitInfo.momentumDirection.x();
  fPhotonDirY = exitInfo.momentumDirection.y();
  fPhotonDirZ = exitInfo.momentumDirection.z();

  CopyText(fPhotonVolumeName,
           sizeof(fPhotonVolumeName),
           exitInfo.volumeName);
  CopyText(fPhotonStepProcessName,
           sizeof(fPhotonStepProcessName),
           exitInfo.stepProcessName);

  fPhotonExitTree->Fill();
}

void RootOutput::FillEventSummary(const EventSummaryInfo& summary)
{
  if (!fEventSummaryTree) {
    return;
  }

  fSummaryEventID = summary.eventID;
  fSummaryElectronProducedCount = summary.electronProducedCount;
  fSummaryElectronChannelCount = summary.electronChannelCount;
  fSummaryElectronChannelPlusCount =
    summary.electronChannelPlusCount;
  fSummaryElectronChannelMinusCount =
    summary.electronChannelMinusCount;
  fSummaryElectronChannelPlusMcp0Count =
    summary.electronChannelPlusMcp0Count;
  fSummaryElectronChannelPlusMcp1Count =
    summary.electronChannelPlusMcp1Count;
  fSummaryElectronChannelPlusMcp2Count =
    summary.electronChannelPlusMcp2Count;
  fSummaryElectronChannelMinusMcp0Count =
    summary.electronChannelMinusMcp0Count;
  fSummaryElectronChannelMinusMcp1Count =
    summary.electronChannelMinusMcp1Count;
  fSummaryElectronChannelMinusMcp2Count =
    summary.electronChannelMinusMcp2Count;
  fSummaryHasElectronChannelPlus =
    summary.hasElectronChannelPlus;
  fSummaryHasElectronChannelMinus =
    summary.hasElectronChannelMinus;
  fSummaryIsCoincidence = summary.isCoincidence;
  fSummaryPhotonExitCount = summary.photonExitCount;
  fSummaryPhotonExitPlusCount = summary.photonExitPlusCount;
  fSummaryPhotonExitMinusCount = summary.photonExitMinusCount;
  fEventSummaryTree->Fill();
}

void RootOutput::Close()
{
  if (!fFile) {
    fElectronChannelHitTree = nullptr;
    fPhotonExitTree = nullptr;
    fEventSummaryTree = nullptr;
    return;
  }

  fFile->cd();
  fElectronChannelHitTree->Write();
  fPhotonExitTree->Write();
  fEventSummaryTree->Write();
  fFile->Close();

  delete fFile;
  fFile = nullptr;
  fElectronChannelHitTree = nullptr;
  fPhotonExitTree = nullptr;
  fEventSummaryTree = nullptr;

  G4cout << "ROOT output closed." << G4endl;
}

std::string RootOutput::BuildThreadFileName(const char* fileName)
{
  const G4int threadID = G4Threading::G4GetThreadId();

  if (!G4Threading::IsMultithreadedApplication() ||
      threadID == G4Threading::SEQUENTIAL_ID) {
    return std::string(fileName);
  }

  std::string baseName(fileName);
  std::string extension;
  const std::string::size_type dotPosition = baseName.rfind('.');

  if (dotPosition != std::string::npos) {
    extension = baseName.substr(dotPosition);
    baseName = baseName.substr(0, dotPosition);
  }

  std::ostringstream outputName;
  outputName << baseName;
  if (threadID == G4Threading::MASTER_ID) {
    outputName << "_master";
  } else {
    outputName << "_t" << threadID;
  }
  outputName << extension;

  return outputName.str();
}
