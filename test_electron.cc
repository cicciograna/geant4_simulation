/// \file test_electrons.cc
/// \brief Modified version of B2a to fire electrons

#include "TEDetectorConstruction.hh"
#include "TEActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
// #include "QGSP_BERT.hh"
// #include "G4EmPenelopePhysics.hh"
// #include "G4EmLivermorePhysics.hh"
// #include "G4EmStandardPhysics_option4.hh"
// #include "G4StepLimiterPhysics.hh"

#include "TEPhysicsList.hh"
#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include<fstream>
#include<stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED  
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new TEDetectorConstruction());

//   G4VModularPhysicsList* physicsList = new QGSP_BERT;
//   physicsList->RegisterPhysics(new G4EmPenelopePhysics());
//   physicsList->RegisterPhysics(new G4EmLivermorePhysics());
//   physicsList->RegisterPhysics(new G4EmStandardPhysics_option4());
//   physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(new TEPhysicsList);
    
  // Set user action classes
  runManager->SetUserInitialization(new TEActionInitialization());
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }
  # ifdef GEANT4_ROOT
  system("tac /home/tug26830/geant4_fullkinematics/local_output/g4root/out_nparts.txt > /home/tug26830/geant4_fullkinematics/local_output/g4root/out_nparts_rev.txt");
  system("mv /home/tug26830/geant4_fullkinematics/local_output/g4root/out_nparts_rev.txt /home/tug26830/geant4_fullkinematics/local_output/g4root/out_nparts.txt");
  std::ofstream   *out_particles = new std::ofstream( "/home/tug26830/geant4_fullkinematics/local_output/g4root/out_particles.txt", std::ofstream::app );
  *out_particles << "1.1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1" << G4endl;
  out_particles->close();
  # endif

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  //
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
