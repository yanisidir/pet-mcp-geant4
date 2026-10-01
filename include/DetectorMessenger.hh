#ifndef DETECTOR_MESSENGER_HH
#define DETECTOR_MESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"
#include <memory>

class DetectorConstruction;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIdirectory;

// /mcp/... macro commands.
class DetectorMessenger : public G4UImessenger
{
public:
  explicit DetectorMessenger(DetectorConstruction* detector);
  ~DetectorMessenger() override;
  void SetNewValue(G4UIcommand* command, G4String value) override;

private:
  DetectorConstruction* fDetector;

  std::unique_ptr<G4UIdirectory> fDirectory;
  std::unique_ptr<G4UIcmdWithAnInteger> fStacksPlusCmd;
  std::unique_ptr<G4UIcmdWithAnInteger> fStacksMinusCmd;
  std::unique_ptr<G4UIcmdWithAnInteger> fStacksCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fFirstPlateZCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fStackGapCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fThicknessCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fRadiusCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fChannelDiameterCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fChannelPitchCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fChannelAngleCmd;
  std::unique_ptr<G4UIcmdWithABool> fChevronCmd;
  std::unique_ptr<G4UIcmdWithADouble> fEdgeMarginCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fWorldSizeCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fProductionCutCmd;
  std::unique_ptr<G4UIcmdWithABool> fKillGammaCmd;
  std::unique_ptr<G4UIcmdWithoutParameter> fPrintCmd;
};
#endif
