#ifndef OUTPUT_MESSENGER_HH
#define OUTPUT_MESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"
#include <memory>

class G4UIcmdWithAString;
class G4UIdirectory;
class RunAction;

// /output/... macro commands.
class OutputMessenger : public G4UImessenger
{
public:
  explicit OutputMessenger(RunAction* runAction);
  ~OutputMessenger() override;
  void SetNewValue(G4UIcommand* command, G4String value) override;

private:
  RunAction* fRunAction;
  std::unique_ptr<G4UIdirectory> fDirectory;
  std::unique_ptr<G4UIcmdWithAString> fFileCmd;
  std::unique_ptr<G4UIcmdWithAString> fScanFileCmd;
};
#endif
