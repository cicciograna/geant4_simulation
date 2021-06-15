/// \file TESteppingAction.cc
/// \brief Implementation of the TESteppingAction class

#include "TESteppingAction.hh"
#include "TEEventAction.hh"
#include "TEDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4VSteppingVerbose.hh"
#include "G4SteppingManager.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TESteppingAction::TESteppingAction(TEEventAction* eventAction)
: G4UserSteppingAction(), 
  fEventAction(eventAction),
  out_particles(NULL),
  oss_init(NULL),
  oss_splat(NULL)
//   fScoringVolume(0)
{
# ifdef GEANT4_ROOT
  out_particles = new std::ofstream( "/home/tug26830/geant4_fullkinematics/local_output/g4root/out_particles.txt", std::ofstream::out );
  oss_init = new std::ostringstream;
  oss_splat = new std::ostringstream;
# endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TESteppingAction::~TESteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TESteppingAction::UserSteppingAction(const G4Step* aStep)
{
  # ifdef GEANT4_ROOT
  std::cout << std::fixed;
  std::cout << std::setprecision(10);
  // energy deposit
  G4String  prePhysVol = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  G4String postPhysVol;
  if( aStep->GetPostStepPoint()->GetPhysicalVolume() ){
    postPhysVol = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
  }
  else{
    postPhysVol = "OutOfWorld";
  }
  
  if( aStep->GetTrack()->GetParentID() == 0 ){
    if ( aStep->GetTrack()->GetCurrentStepNumber() == 1 ){
      oss_init  ->str(""); oss_init ->clear();
      oss_splat ->str(""); oss_splat->clear();
      
      
      TrackID     = aStep->GetTrack()       ->GetTrackID(); 
      TOF         = 0.;
      Mass        = aStep->GetTrack()       ->GetDynamicParticle()->GetMass(); 
      Charge      = aStep->GetTrack()       ->GetDynamicParticle()->GetCharge();
      InitX       = aStep->GetPreStepPoint()->GetPosition().getX();
      InitY       = aStep->GetPreStepPoint()->GetPosition().getY();
      InitZ       = aStep->GetPreStepPoint()->GetPosition().getZ();
      InitAzm     = aStep->GetPreStepPoint()->GetMomentum().phi()/deg; if(InitAzm < 0.){ InitAzm += 360.; }
      InitElev    = aStep->GetPreStepPoint()->GetMomentum().theta()/deg;
      if(Mass > 0.) { InitVx = ( aStep->GetPreStepPoint()->GetMomentum().getX()*(keV/CLHEP::c_light) ) / ( Mass*(keV/CLHEP::c_squared) ); } else{ InitVx = CLHEP::c_light * sin(aStep->GetPreStepPoint()->GetMomentum().theta()) * cos(aStep->GetPreStepPoint()->GetMomentum().phi()); };
      if(Mass > 0.) { InitVy = ( aStep->GetPreStepPoint()->GetMomentum().getY()*(keV/CLHEP::c_light) ) / ( Mass*(keV/CLHEP::c_squared) ); } else{ InitVy = CLHEP::c_light * sin(aStep->GetPreStepPoint()->GetMomentum().theta()) * sin(aStep->GetPreStepPoint()->GetMomentum().phi()); };
      if(Mass > 0.) { InitVz = ( aStep->GetPreStepPoint()->GetMomentum().getZ()*(keV/CLHEP::c_light) ) / ( Mass*(keV/CLHEP::c_squared) ); } else{ InitVz = CLHEP::c_light * cos(aStep->GetPreStepPoint()->GetMomentum().theta())                                                     ; };
      InitKE      = aStep->GetPreStepPoint()->GetKineticEnergy();
      
//       G4cout << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << " " << TrackID << " " << TOF << " " << Mass << " " << Charge << " " << InitX << " " << InitY << " " << InitZ << " " << InitAzm << " " << InitElev << " " << InitDircos << " " << InitVx << " " << InitVy << " " << InitVz << " " << InitVtot << " " << InitVrad << " " << InitVlong << " " << InitPx << " " << InitPy << " " << InitPz << " " << InitPtot << " " << InitPrad << " " << InitPlong << " " << InitKE << G4endl;
      
      *oss_init << std::fixed << std::setprecision(10) << TrackID << " " << TOF/us << " " << Mass/CLHEP::amu_c2 << " " << Charge << " " << InitX << " " << InitY << " " << InitZ << " " << InitAzm << " " << InitElev << " " << InitVx/(mm/us) << " " << InitVy/(mm/us) << " " << InitVz/(mm/us) << " 1 1 1 1 1 1 1 1 1 " << InitKE/eV;
    }
    
    
    if( aStep->GetPostStepPoint()->GetKineticEnergy() == 0 || postPhysVol == "OutOfWorld" || aStep->GetTrack()->GetLocalTime() >= 5000.){
      TrackID     = aStep->GetTrack()       ->GetTrackID(); 
      TOF         = aStep->GetTrack()       ->GetLocalTime();
      Mass        = aStep->GetTrack()       ->GetDynamicParticle()->GetMass();
      Charge      = aStep->GetTrack()       ->GetDynamicParticle()->GetCharge();
      SplatX      = aStep->GetPostStepPoint()->GetPosition().getX();
      SplatY      = aStep->GetPostStepPoint()->GetPosition().getY();
      SplatZ      = aStep->GetPostStepPoint()->GetPosition().getZ();
      SplatAzm    = aStep->GetPostStepPoint()->GetMomentum().phi()/deg; if(InitAzm < 0.){ InitAzm += 360.; }
      SplatElev   = aStep->GetPostStepPoint()->GetMomentum().theta()/deg;
      if(Mass > 0.) { SplatVx = ( aStep->GetPostStepPoint()->GetMomentum().getX()*(keV/CLHEP::c_light) ) / ( Mass*(keV/CLHEP::c_squared) ); } else{ SplatVx = CLHEP::c_light * sin(aStep->GetPostStepPoint()->GetMomentum().theta()) * cos(aStep->GetPreStepPoint()->GetMomentum().phi()); };
      if(Mass > 0.) { SplatVy = ( aStep->GetPostStepPoint()->GetMomentum().getY()*(keV/CLHEP::c_light) ) / ( Mass*(keV/CLHEP::c_squared) ); } else{ SplatVy = CLHEP::c_light * sin(aStep->GetPostStepPoint()->GetMomentum().theta()) * sin(aStep->GetPreStepPoint()->GetMomentum().phi()); };
      if(Mass > 0.) { SplatVz = ( aStep->GetPostStepPoint()->GetMomentum().getZ()*(keV/CLHEP::c_light) ) / ( Mass*(keV/CLHEP::c_squared) ); } else{ SplatVz = CLHEP::c_light * cos(aStep->GetPostStepPoint()->GetMomentum().theta())                                                     ; };
      SplatKE     = aStep->GetPostStepPoint()->GetKineticEnergy();

//       *out_particles << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << " morta" << G4endl;
      *oss_splat << std::fixed << std::setprecision(10) << TrackID << " " << TOF/us << " " << Mass/CLHEP::amu_c2 << " " << Charge << " " << SplatX << " " << SplatY << " " << SplatZ << " " << SplatAzm << " " << SplatElev << " " << SplatVx/(mm/us) << " " << SplatVy/(mm/us) << " " << SplatVz/(mm/us) << " 1 1 1 1 1 1 1 1 1 " << SplatKE/eV;

//       G4cout << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << " " << TrackID << " " << TOF << " " << Mass << " " << Charge << " " << InitX << " " << InitY << " " << InitZ << " " << InitAzm << " " << InitElev << " " << InitDircos << " " << InitVx << " " << InitVy << " " << InitVz << " " << InitVtot << " " << InitVrad << " " << InitVlong << " " << InitPx << " " << InitPy << " " << InitPz << " " << InitPtot << " " << InitPrad << " " << InitPlong << " " << InitKE << G4endl;
      *out_particles << oss_init->str().c_str() << G4endl << oss_splat->str().c_str() << G4endl;
      if(TOF > 5000.){ aStep->GetTrack()->SetTrackStatus(fStopAndKill); }
    }
  }
  
  if(aStep->IsFirstStepInVolume()){
//     G4cout << "First step" << G4endl << aStep->GetPreStepPoint()->GetMomentum() << G4endl << G4endl;
  }
//   ++stepcounter;
  if(aStep->IsLastStepInVolume()){
//     G4cout << "Last step " << stepcounter << G4endl << aStep->GetPreStepPoint()->GetMomentum() << G4endl << aStep->GetPostStepPoint()->GetMomentum() << G4endl;
  }
   
# endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

