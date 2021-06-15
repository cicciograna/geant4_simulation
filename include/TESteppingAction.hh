/// \file TESteppingAction.hh
/// \brief Definition of the TESteppingAction class

#ifndef TESteppingAction_h
#define TESteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

# include <fstream>
# include <sstream>
#include "G4StepStatus.hh"

class TEEventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class TESteppingAction : public G4UserSteppingAction
{
  public:
    TESteppingAction(TEEventAction* eventAction);
    virtual ~TESteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    TEEventAction*  fEventAction;
//     G4LogicalVolume* fScoringVolume;
    G4int stepcounter = 0;
    
    std::ofstream   *out_particles;
    std::ostringstream *oss_init;
    std::ostringstream *oss_splat;
    G4int            nparticles;
    G4StepStatus    stepStatus;
    
      G4int         TrackID     ; 
      G4double      TOF         ;
      G4double      Mass        ;
      G4int         Charge      ;
      G4double      InitX       ;
      G4double      InitY       ;
      G4double      InitZ       ;
      G4double      InitAzm     ;
      G4double      InitElev    ;
      G4double      InitDircos  ;
      G4double      InitVx      ;
      G4double      InitVy      ;
      G4double      InitVz      ;
      G4double      InitVtot    ;
      G4double      InitVrad    ;
      G4double      InitVlong   ;
      G4double      InitPx      ;
      G4double      InitPy      ;
      G4double      InitPz      ;
      G4double      InitPtot    ;
      G4double      InitPrad    ;
      G4double      InitPlong   ;
      G4double      InitKE      ;
      G4double      SplatX      ;
      G4double      SplatY      ;
      G4double      SplatZ      ;
      G4double      SplatAzm    ;
      G4double      SplatElev   ;
      G4double      SplatDircos ;
      G4double      SplatVx     ;
      G4double      SplatVy     ;
      G4double      SplatVz     ;
      G4double      SplatVtot   ;
      G4double      SplatVrad   ;
      G4double      SplatVlong  ;
      G4double      SplatPx     ;
      G4double      SplatPy     ;
      G4double      SplatPz     ;
      G4double      SplatPtot   ;
      G4double      SplatPrad   ;
      G4double      SplatPlong  ;
      G4double      SplatKE     ;
      G4double      SplatRadius ;
      G4bool        OnTarget    ;
      G4bool        OnMCPPlane  ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
