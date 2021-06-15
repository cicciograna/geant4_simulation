/// \file TEDetectorConstruction.hh
/// \brief Definition of the TEDetectorConstruction class

#ifndef TEDetectorConstruction_h
#define TEDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

#include "TFile.h"
#include "TTree.h"

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
    
    G4int nE27_E31 =  5;
    G4int nE32_E34 =  4;
    G4int nE45_E49 =  5;
    G4int nE51_E67 = 17;
    G4int nE68_E87 = 20;
    
    G4LogicalVolume **fLVE001_E013      ; G4int nE001_E013 = 13;
    G4LogicalVolume **fLVE014_E041      ; G4int nE014_E041 = 28;
    G4LogicalVolume  *fLVE042           ;
    G4LogicalVolume **fLVE043_E046      ; G4int nE043_E046 =  4;
    G4LogicalVolume **fLVE047_E048      ; G4int nE047_E048 =  2;
    G4LogicalVolume  *fLVE049           ;
    G4LogicalVolume **fLVE050_E052      ; G4int nE050_E052 =  3;
    G4LogicalVolume  *fLVE053           ;
    G4LogicalVolume **fLVE054_E057      ; G4int nE054_E057 =  4;
    G4LogicalVolume  *fLVIonCoilHolder  ;
    G4LogicalVolume  *fLVE058;
    // Between E58 and E59 there's the source
    G4LogicalVolume  *fLVE059           ;
    G4LogicalVolume  *fLVEleCoilHolder  ;
    G4LogicalVolume **fLVE060_E064      ; G4int nE060_E064 =  4;
    G4LogicalVolume  *fLVE065           ;
    G4LogicalVolume **fLVE066_E068      ; G4int nE066_E068 =  3;
    G4LogicalVolume  *fLVE069           ;
    G4LogicalVolume **fLVE070_E071      ; G4int nE070_E071 =  2;
    G4LogicalVolume **fLVE072_E076      ; G4int nE072_E076 =  4;
    G4LogicalVolume  *fLVE077           ;
    G4LogicalVolume **fLVE078_E105      ; G4int nE078_E105 = 28;
    G4LogicalVolume **fLVE106_E117      ; G4int nE106_E117 = 13;
  
    // data members
    G4LogicalVolume*   fLogicTarget;     // pointer to the logical Target

    G4Material*        fTargetMaterial;  // pointer to the target material
    G4Material*        fVesselMaterial;  // pointer to the vessel material

    G4UserLimits* fStepLimit;            // pointer to user step limits

//     TEDetectorMessenger*  fMessenger;   // messenger

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                         // magnetic field messenger
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 
    
    TTree *tfield;
    
    G4double                  zOffset = 0.;
    G4Cache<G4ElectroMagneticField*> fBField;  // pointer to the thread-local fields
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
