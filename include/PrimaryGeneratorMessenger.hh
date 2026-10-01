#ifndef PRIMARY_GENERATOR_MESSENGER_HH
#define PRIMARY_GENERATOR_MESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"
#include <memory>

class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIdirectory;
class PrimaryGeneratorAction;

// /source/... macro commands. The gun's own /gun/ commands stay available
// but bypass the cone sampling.
class PrimaryGeneratorMessenger : public G4UImessenger
{
public:
  explicit PrimaryGeneratorMessenger(PrimaryGeneratorAction* generator);
  ~PrimaryGeneratorMessenger() override;
  void SetNewValue(G4UIcommand* command, G4String value) override;

private:
  PrimaryGeneratorAction* fGenerator;

  std::unique_ptr<G4UIdirectory> fDirectory;
  std::unique_ptr<G4UIcmdWithAString> fParticleCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fEnergyCmd;
  std::unique_ptr<G4UIcmdWith3VectorAndUnit> fPositionCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fPositionSigmaCmd;
  std::unique_ptr<G4UIcmdWith3Vector> fConeAxisCmd;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fConeHalfAngleCmd;
  std::unique_ptr<G4UIcmdWithABool> fBackToBackCmd;
  std::unique_ptr<G4UIcmdWithAnInteger> fMultiplicityCmd;
  std::unique_ptr<G4UIcmdWithoutParameter> fPrintCmd;
};
#endif
