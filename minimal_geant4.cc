#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"

#include "FTFP_BERT.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"

int main(int argc, char** argv)
{
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  G4SteppingVerbose::UseBestUnit(4);

  // ROOT objects are created by Geant4 worker threads, so initialize ROOT
  // once on the main thread before the run starts.
  ROOT::EnableThreadSafety();
  TFile::Class();
  TTree::Class();

  //auto* runManager =
    //G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  auto* runManager = G4RunManagerFactory::CreateRunManager();

  G4int nThreads = 1;
  runManager->SetNumberOfThreads(nThreads);

  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  // Keep the complete FTFP_BERT reference list, but replace its standard
  // electromagnetic physics with the Livermore low-energy EM models.
  FTFP_BERT* physicsList = new FTFP_BERT(0);  // 0 means "don't print verbose output"
  physicsList->ReplacePhysics(new G4EmLivermorePhysics(0));
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization(detector));

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* uiManager = G4UImanager::GetUIpointer();

  if (ui) {
    uiManager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    uiManager->ApplyCommand(command + fileName);
  }

  delete visManager;
  delete runManager;

  return 0;
}
