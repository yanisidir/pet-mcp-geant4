#ifndef ROOT_OUTPUT_HH
#define ROOT_OUTPUT_HH

#include "ElectronChannelHitInfo.hh"
#include "EventSummaryInfo.hh"
#include "GammaInteractionInfo.hh"
#include "GammaMcpEntryInfo.hh"
#include "McpPlateStatsInfo.hh"
#include "PhotonExitInfo.hh"

#include "G4Threading.hh"

#include <string>

class TFile;
class TTree;

// Sortie ROOT : arbres indépendants pour les observables de simulation.
class RootOutput
{
public:
  static RootOutput* Instance();

  void Open(const char* fileName);
  void FillElectronChannelHit(const ElectronChannelHitInfo& hit);
  void FillPhotonExit(const PhotonExitInfo& exitInfo);
  void FillGammaInteraction(const GammaInteractionInfo& info);
  void FillGammaMcpEntry(const GammaMcpEntryInfo& entry);
  void FillEventSummary(const EventSummaryInfo& summary);
  void FillMcpPlateStats(const McpPlateStatsInfo& stats);
  void Close();

private:
  RootOutput();
  ~RootOutput();

  RootOutput(const RootOutput&);
  RootOutput& operator=(const RootOutput&);

  static std::string BuildThreadFileName(const char* fileName);
  static G4ThreadLocal RootOutput* fInstance;

  TFile* fFile;
  TTree* fElectronChannelHitTree;
  TTree* fPhotonExitTree;
  TTree* fGammaInteractionTree;
  TTree* fGammaMcpEntryTree;
  TTree* fEventSummaryTree;
  TTree* fMcpPlateStatsTree;

  // ElectronChannelHitTree
  int fElectronEventID;
  int fElectronTrackID;
  int fElectronParentID;
  int fElectronSide;
  int fElectronMcpIndex;
  double fElectronKineticEnergy;
  double fElectronGlobalTime;
  double fElectronX;
  double fElectronY;
  double fElectronZ;
  double fElectronDirX;
  double fElectronDirY;
  double fElectronDirZ;
  char fElectronVolumeName[64];
  char fElectronCreatorProcessName[64];

  // PhotonExitTree
  int fPhotonEventID;
  int fPhotonTrackID;
  int fPhotonParentID;
  int fPhotonSide;
  double fPhotonKineticEnergy;
  double fPhotonGlobalTime;
  double fPhotonX;
  double fPhotonY;
  double fPhotonZ;
  double fPhotonDirX;
  double fPhotonDirY;
  double fPhotonDirZ;
  char fPhotonVolumeName[64];
  char fPhotonStepProcessName[64];

  // GammaInteractionTree
  int fGammaEventID;
  int fGammaTrackID;
  int fGammaParentID;
  int fGammaSide;
  int fGammaMcpIndex;
  char fGammaProcessName[64];
  double fGammaKineticEnergy;
  double fGammaGlobalTime;
  double fGammaX;
  double fGammaY;
  double fGammaZ;
  char fGammaVolumeName[64];

  // GammaMcpEntryTree
  int fGammaEntryEventID;
  int fGammaEntryTrackID;
  int fGammaEntryParentID;
  int fGammaEntrySide;
  int fGammaEntryMcpIndex;
  double fGammaEntryKineticEnergy;
  double fGammaEntryGlobalTime;
  double fGammaEntryX;
  double fGammaEntryY;
  double fGammaEntryZ;
  char fGammaEntryVolumeName[64];

  // EventSummaryTree
  int fSummaryEventID;
  int fSummaryElectronProducedCount;
  int fSummaryElectronChannelCount;
  int fSummaryElectronChannelPlusCount;
  int fSummaryElectronChannelMinusCount;
  bool fSummaryHasElectronChannelPlus;
  bool fSummaryHasElectronChannelMinus;
  bool fSummaryIsCoincidence;
  int fSummaryPhotonExitCount;
  int fSummaryPhotonExitPlusCount;
  int fSummaryPhotonExitMinusCount;
  int fSummaryGammaInteractionCount;
  int fSummaryGammaPhotCount;
  int fSummaryGammaComptCount;
  int fSummaryGammaRaylCount;
  int fSummaryGammaConvCount;

  // McpPlateStatsTree
  int fPlateEventID;
  int fPlateSide;
  int fPlateMcpIndex;
  int fPlateElectronChannelCount;
};

#endif
