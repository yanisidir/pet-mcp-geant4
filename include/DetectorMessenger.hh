#ifndef DETECTOR_MESSENGER_HH
#define DETECTOR_MESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIdirectory;

class DetectorMessenger : public G4UImessenger
{
public:
  explicit DetectorMessenger(DetectorConstruction* detector);
  virtual ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  DetectorMessenger(const DetectorMessenger&);
  DetectorMessenger& operator=(const DetectorMessenger&);

  DetectorConstruction* fDetector;
  G4UIdirectory* fMcpDirectory;
  G4UIcmdWithAString* fMaterialCommand;
  G4UIcmdWithABool* fEnableStackFieldsCommand;
  G4UIcmdWithADoubleAndUnit* fStackFieldMagnitudeCommand;
  G4UIcmdWithADoubleAndUnit* fFieldRegionRadialMarginCommand;
  G4UIcmdWithADoubleAndUnit* fFieldRegionZMarginCommand;
};

#endif
