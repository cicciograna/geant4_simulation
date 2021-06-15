/// \file TEDetectorConstruction.hh
/// \brief Definition of the TEDetectorConstruction class

#ifndef TEDetectorConstruction_h
#define TEDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

#include "TFile.h"
#include "TTree.h"

#include "G4AssemblyVolume.hh"

#include "G4Cache.hh"
#include "G4ElectroMagneticField.hh"

#include "TEDefinitions.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class TEDetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class TEDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    TEDetectorConstruction();
    virtual ~TEDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // Set methods
    void SetTargetMaterial (G4String );
    void SetMaxStep (G4double );
    void SetCheckOverlaps(G4bool );

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    
    // POINTERS TO LOGICAL VOLUMES
    // Container structure that will hold the various electrodes
    // and that will mirror aroudn the source
    G4AssemblyVolume *fSpectrometer;
    
    // Electrodes
    G4LogicalVolume  *fLVE01      ; G4int nE01      =  1; G4double dzE01      =  0. ;
    G4LogicalVolume **fLVE02_E13  ; G4int nE02_E13  = 12; G4double dzE02_E13  = 14. ;
    G4LogicalVolume **fLVE14_E41  ; G4int nE14_E41  = 28; G4double dzE14_E41  = 19. ;
    G4LogicalVolume  *fLVE42      ; G4int nE42      =  1; G4double dzE42      =  9. ;
    G4LogicalVolume **fLVE43_E46  ; G4int nE43_E46  =  4; G4double dzE43_E46  = 11. ;
    G4LogicalVolume **fLVE47_E48  ; G4int nE47_E48  =  2; G4double dzE47_E48  =  5. ;
    G4LogicalVolume  *fLVE49      ; G4int nE49      =  1; G4double dzE49      =  2. ;
    G4LogicalVolume **fLVE50_E52  ; G4int nE50_E52  =  3; G4double dzE50_E52  =  6. ;
    G4LogicalVolume  *fLVE53      ; G4int nE53      =  1; G4double dzE53      =  2. ;
    G4LogicalVolume **fLVE54_E57  ; G4int nE54_E57  =  4; G4double dzE54_E57  =  6. ;
    G4LogicalVolume  *fLVCHLeft   ; G4int nCHLeft   =  1; G4double dzCHLeft   =  1.5;
    G4LogicalVolume  *fLVCHExtraE ; G4int nCHExtraE =  1; G4double dzCHExtraE =  1. ;
    G4LogicalVolume  *fLVCHRight  ; G4int nCHRight  =  1; G4double dzCHRight  =  1. ;
    G4LogicalVolume  *fLVE58      ; G4int nE58      =  1; G4double dzE58      = 29. ;
                                                          G4double dzSource   = 31. ; // Between E57 and E58 there's the source
    // MCP
    G4LogicalVolume*   fLVMCP;     // pointer to the logical Target

    // Vessel
    G4Material*        fTargetMaterial;  // pointer to the target material
    G4Material*        fVesselMaterial;  // pointer to the vessel material

    G4UserLimits* fStepLimit;            // pointer to user step limits

//     TEDetectorMessenger*  fMessenger;   // messenger

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                         // magnetic field messenger
    G4Transform3D nulltrans;
    G4Transform3D reflect3D;
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 
    
    TTree *tfield;
    
    G4double                  zOffset = 0.;
    G4Cache<G4ElectroMagneticField*> fBField;  // pointer to the thread-local fields
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
