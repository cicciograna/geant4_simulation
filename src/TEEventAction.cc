/// \file TEEventAction.cc
/// \brief Implementation of the TEEventAction class

#include "TEEventAction.hh"
#include "TETrackerHit.hh"

#include "TEAnalysis.hh"
#include "G4SDManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
// #include "G4TrajectoryContainer.hh"
// #include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found 
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl; 
      G4Exception("B5EventAction::EndOfEventAction()",
                  "B5Code001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl; 
    G4Exception("B5EventAction::EndOfEventAction()",
                "B5Code001", JustWarning, msg);
  }
  return hc;  
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TEEventAction::TEEventAction()
: G4UserEventAction()
{
  index = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TEEventAction::~TEEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEEventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEEventAction::EndOfEventAction(const G4Event* event)
{
# if defined(ISO_CURVES_GENERATOR) || defined(ENERGIES_RANGE_ISO_CURVES_GENERATOR) || defined(FIXED_PARTICLES)
  auto runManager =  G4RunManager::GetRunManager();
  G4int run_id = runManager->GetCurrentRun()->GetRunID();
//   delete runManager;
  
  if(index < 0){ index = G4SDManager::GetSDMpointer()->GetCollectionID("TrackerHitsCollection"); }
  
  auto analysisManager = G4AnalysisManager::Instance();

  auto hc = GetHC(event, index);
  for (int i = 0; i < hc->GetSize(); ++i) {
    auto hit = static_cast<TETrackerHit*>(hc->GetHit(i));
//     if(hit->GetCharge() < 0 && hit->GetInitKE() < 150){ // only electrons with standard energies for now
    if(hit->GetCharge() < 0          && 
      (hit->GetInitPos().x() >   -3. && hit->GetInitPos().x() <    3.) &&
      (hit->GetInitPos().y() >   -3. && hit->GetInitPos().y() <    3.) &&
      (hit->GetInitPos().z() > originZpos -3. && hit->GetInitPos().z() < originZpos + 3.) ){   // only electrons generated in the source region
      analysisManager->FillNtupleDColumn( 0, hit->GetTrackID()          );  // ionn
      analysisManager->FillNtupleDColumn( 1, hit->GetMass()             );  // mass
      analysisManager->FillNtupleDColumn( 2, hit->GetCharge()           );  // charge
      analysisManager->FillNtupleDColumn( 3, hit->GetTOF()              );  // tof
      analysisManager->FillNtupleDColumn( 4, hit->GetInitPos().x()      );  // init_x
      analysisManager->FillNtupleDColumn( 5, hit->GetInitPos().y()      );  // init_y
      analysisManager->FillNtupleDColumn( 6, hit->GetInitPos().z()      );  // init_z
      analysisManager->FillNtupleDColumn( 7, hit->GetInitAzm()          );  // init_azm
      analysisManager->FillNtupleDColumn( 8, hit->GetInitElev()         );  // init_elev
      analysisManager->FillNtupleDColumn( 9, hit->GetInitDircos()       );  // init_dircos
      analysisManager->FillNtupleDColumn(10, hit->GetInitVel().x()      );  // init_vx
      analysisManager->FillNtupleDColumn(11, hit->GetInitVel().y()      );  // init_vy
      analysisManager->FillNtupleDColumn(12, hit->GetInitVel().z()      );  // init_vz
      analysisManager->FillNtupleDColumn(13, hit->GetInitVel().mag()    );  // init_vtot
      analysisManager->FillNtupleDColumn(14, hit->GetInitVel().perp()   );  // init_vrad
      analysisManager->FillNtupleDColumn(15, hit->GetInitVel().z()      );  // init_vlong
      analysisManager->FillNtupleDColumn(16, hit->GetInitMom().x()      );  // init_px
      analysisManager->FillNtupleDColumn(17, hit->GetInitMom().y()      );  // init_py
      analysisManager->FillNtupleDColumn(18, hit->GetInitMom().z()      );  // init_pz
      analysisManager->FillNtupleDColumn(19, hit->GetInitMom().mag()    );  // init_ptot
      analysisManager->FillNtupleDColumn(20, hit->GetInitMom().perp()   );  // init_prad
      analysisManager->FillNtupleDColumn(21, hit->GetInitMom().z()      );  // init_plong
      analysisManager->FillNtupleDColumn(22, hit->GetInitKE()           );  // init_ke
      analysisManager->FillNtupleDColumn(23, hit->GetSplatPos().x()     );  // splat_x
      analysisManager->FillNtupleDColumn(24, hit->GetSplatPos().y()     );  // splat_y
      analysisManager->FillNtupleDColumn(25, hit->GetSplatPos().z()     );  // splat_z
      analysisManager->FillNtupleDColumn(26, hit->GetSplatAzm()         );  // splat_azm
      analysisManager->FillNtupleDColumn(27, hit->GetSplatElev()        );  // splat_elev
      analysisManager->FillNtupleDColumn(28, hit->GetSplatDircos()      );  // splat_dircos
      analysisManager->FillNtupleDColumn(29, hit->GetSplatVel().x()     );  // splat_vx
      analysisManager->FillNtupleDColumn(30, hit->GetSplatVel().y()     );  // splat_vy
      analysisManager->FillNtupleDColumn(31, hit->GetSplatVel().z()     );  // splat_vz
      analysisManager->FillNtupleDColumn(32, hit->GetSplatVel().mag()   );  // splat_vtot
      analysisManager->FillNtupleDColumn(33, hit->GetSplatVel().perp()  );  // splat_vrad
      analysisManager->FillNtupleDColumn(34, hit->GetSplatVel().z()     );  // splat_vlong
      analysisManager->FillNtupleDColumn(35, hit->GetSplatMom().x()     );  // splat_px
      analysisManager->FillNtupleDColumn(36, hit->GetSplatMom().y()     );  // splat_py
      analysisManager->FillNtupleDColumn(37, hit->GetSplatMom().z()     );  // splat_pz
      analysisManager->FillNtupleDColumn(38, hit->GetSplatMom().mag()   );  // splat_ptot
      analysisManager->FillNtupleDColumn(39, hit->GetSplatMom().perp()  );  // splat_prad
      analysisManager->FillNtupleDColumn(40, hit->GetSplatMom().z()     );  // splat_plong
      analysisManager->FillNtupleDColumn(41, hit->GetSplatKE()          );  // splat_ke
      analysisManager->FillNtupleDColumn(42, hit->GetSplatPos().perp()  );  // splat_radius
      analysisManager->FillNtupleDColumn(43, 64501111.                );  // run_id
      
      
      if((hit->GetSplatPos().x() > -60. && hit->GetSplatPos().x() < 60. ) &&
         (hit->GetSplatPos().y() > -60. && hit->GetSplatPos().y() < 60. ) && 
         (hit->GetSplatPos().z() >= (originZpos + targetZ - targetLength - 0.01) && hit->GetSplatPos().z() <= (originZpos + targetZ - targetLength + 0.01) )) {
           analysisManager->FillNtupleIColumn(44, 1);}
      else{analysisManager->FillNtupleIColumn(44, 0);}                        // ontarget
      
      if(hit->GetSplatPos().z() >= (originZpos + targetZ - targetLength - 0.01) && hit->GetSplatPos().z() <= (originZpos + targetZ - targetLength + 0.01) ) {
           analysisManager->FillNtupleIColumn(45, 1);}
      else{analysisManager->FillNtupleIColumn(45, 0);}                        // onmcpplane
      
      analysisManager->AddNtupleRow();
    }
  }
# endif
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
