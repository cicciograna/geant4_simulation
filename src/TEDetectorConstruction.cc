/// \file TEDetectorConstruction.cc
/// \brief Implementation of the TEDetectorConstruction class

#include "TEDetectorConstruction.hh"
// #include "TEDetectorMessenger.hh"
#include "TETrackerSD.hh"

#include "TEMagTabulatedField3D.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
// #include "G4SimpleRunge.hh"
// #include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4PropagatorInField.hh"

#include "G4EquationOfMotion.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4GenericPolycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include <cstdio>

// #define MAGSTANDARD

// #define ELECTRODE_ZPOS

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* TEDetectorConstruction::fMagFieldMessenger = 0;

TEDetectorConstruction::TEDetectorConstruction() :G4VUserDetectorConstruction(),
 fLVMCP(NULL),
 fTargetMaterial(NULL),
 fVesselMaterial(NULL),
 fStepLimit(NULL),
 fSpectrometer(NULL),
 fLVE02_E13(NULL),
 fLVE14_E41(NULL),
 fLVE43_E46(NULL),
 fLVE47_E48(NULL),
 fLVE50_E52(NULL),
 fLVE54_E57(NULL),
 fCheckOverlaps(true)
{
  fBField.Put(0);
  fSpectrometer = new G4AssemblyVolume();
  fLVE02_E13 = new G4LogicalVolume*[nE02_E13];
  fLVE14_E41 = new G4LogicalVolume*[nE14_E41];
  fLVE43_E46 = new G4LogicalVolume*[nE43_E46];
  fLVE47_E48 = new G4LogicalVolume*[nE47_E48];
  fLVE50_E52 = new G4LogicalVolume*[nE50_E52];
  fLVE54_E57 = new G4LogicalVolume*[nE54_E57];
  nulltrans = HepGeom::Transform3D(); 
  reflect3D = HepGeom::ReflectZ3D();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TEDetectorConstruction::~TEDetectorConstruction()
{
//  delete [] fSpectrometer;
//  delete [] fLVE02_E13;
//  delete [] fLVE14_E41;
//  delete [] fLVE43_E46;
//  delete [] fLVE47_E48;
//  delete [] fLVE50_E52;
//  delete [] fLVE54_E57;
//   delete fStepLimit;
//   delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TEDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEDetectorConstruction::DefineMaterials()
{
  // Material definition

  G4NistManager* nistManager = G4NistManager::Instance();

  // Vacuum defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_Galactic");

  // Gaseous Helium, initialized from NIST Helium, but modified to change its pressure
  G4Material *NIST_He = nistManager->FindOrBuildMaterial("G4_He");
//   G4double atm = 1e-30*atmosphere;
  G4double atm = 1.25e-6;
  G4Material *gHe = new G4Material("gHe",                   // Name
//                                    4./25400.*atm,   // Density
                                   NIST_He->GetDensity()*atm,   // Density
//                                    NIST_He->GetDensity(),   // Density
                                   NIST_He,                 // Base material
                                   kStateUndefined,         // State
                                   CLHEP::STP_Temperature,  // Temperature
                                   atm*atmosphere);                    // Pressure

  // Stainless steel defined using NIST Manager
  fVesselMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

  // Copper defined using NIST Manager
  fTargetMaterial = nistManager->FindOrBuildMaterial("G4_Cu");


  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TEDetectorConstruction::DefineVolumes()
{
  G4Material* vacuum  = G4Material::GetMaterial("G4_Galactic");
  G4Material* helium  = G4Material::GetMaterial("G4_He");
  G4Material* gHe     = G4Material::GetMaterial("gHe");
  G4cout << "Vacuum pressure: " << vacuum->GetPressure()/atmosphere << " atm" << G4endl;
  G4cout << "Vacuum density: "  << vacuum->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
  G4cout << "G4_He pressure: " << helium->GetPressure()/atmosphere << " atm" << G4endl;
  G4cout << "G4_He density: "  << helium->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
  G4cout << "Gaseous Helium pressure: " << gHe->GetPressure()/atmosphere << " atm" << G4endl;
  G4cout << "Gaseous Helium density: "  << gHe->GetDensity()/(g/cm3) << " g/cm3" << G4endl;
//   G4Material* ssteel  = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  G4Material* copper  = G4Material::GetMaterial("G4_Cu");

  // Sizes of the principal geometrical components (solids)

  // Sizes included in the file TEDefinitions.hh

  // Visualization attributes

  G4VisAttributes* worldVisAtt  = new G4VisAttributes( G4Colour(1.0,1.0,1.0,0.0 ) );
  G4VisAttributes* boxVisAtt    = new G4VisAttributes( G4Colour(1.0,1.0,1.0,0.0 ) );
  G4VisAttributes* VesselVisAtt = new G4VisAttributes( G4Colour(1.0,1.0,1.0,0.0 ) );
  G4VisAttributes* InitVisAtt   = new G4VisAttributes( G4Colour(1.0,0.0,0.0,1) );
  G4VisAttributes* LargeVisAtt  = new G4VisAttributes( G4Colour(0.0,1.0,0.0,1) );
  G4VisAttributes* SmallVisAtt  = new G4VisAttributes( G4Colour(0.0,0.0,1.0,1) );
  G4VisAttributes* MCPVisAtt    = new G4VisAttributes( G4Colour(1.0,1.0,0.0,1.0) );
  G4VisAttributes* MCPBoxVisAtt = new G4VisAttributes( G4Colour(1.0,1.0,1.0,0.50 ) );

  // Definitions of Solids, Logical Volumes, Physical Volumes
  // World
  G4ThreeVector positionWorld = G4ThreeVector(0,0,0);   // with this vector, ther source will be in (0,0,0)

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = " << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm << " mm" << G4endl;

  G4Box* worldS = new G4Box("world",                                                       //its name
                            (worldWidth/2.)*mm, (worldWidth/2.)*mm, (worldLength/2.*mm));   //its size

  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS,    //its solid
                                                 vacuum,    //its material
                                                 "World");  //its name

  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume* worldPV = new G4PVPlacement(0,               // no rotation
                                                 positionWorld,   // at (0,0,0)
                                                 worldLV,         // its logical volume
                                                 "World",         // its name
                                                 0,               // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps

  // Vessel
  // Outside
//   G4ThreeVector positionVesselOut = G4ThreeVector(0,0,0);
//   G4ThreeVector positionVesselOut = G4ThreeVector(0,0,(vesselLength/2. + vesselThick)*mm);   // with this vector, ther source will be in (0,0,0)
  G4ThreeVector positionVesselOut = G4ThreeVector(0,0,vesselZoffset);   // with this vector, ther source will be in (0,0,0)

  G4Tubs* vesselOutS = new G4Tubs("vesselOut",        // name
                               0.*mm,                 // inner radius
                               vesselRadius*mm,       // outer radius
                               (vesselLength/2. + vesselThick)*mm,  // length of outer vessel
                               0.*deg,                // init phi
                               360.*deg);             // end phi

  G4LogicalVolume *vesselOutLV = new G4LogicalVolume(vesselOutS,    // its solid
                                                  fVesselMaterial,  // its material
                                                  "VesselOut",      // its name
                                                  0,                // its FieldManager (will have to change)
                                                  0,                // its sensitive detector
                                                  0);               // its user cuts

  new G4PVPlacement(0,                  // no rotation
                    positionVesselOut,  // at (0,0,0)
                    vesselOutLV,        // its logical volume
                    "VesselOut",        // its name
                    worldLV,            // its mother  volume
                    false,              // no boolean operations
                    0,                  // copy number
                    fCheckOverlaps);    // checking overlaps

  // Inside
  G4ThreeVector positionVesselIn = G4ThreeVector(0,0,0);

  G4Tubs* vesselInS = new G4Tubs("vesselIn",                              // name
                                 0.*mm,                                   // inner radius
                                 (vesselRadius-vesselThick)*mm,           // outer radius
                                 (vesselLength/2)*mm,    // length of inner vessel
                                 0.*deg,                                  // init phi
                                 360.*deg);                               // end phi

  G4LogicalVolume *vesselInLV = new G4LogicalVolume(vesselInS,    // its solid
//                                                     gHe,          // its material
                                                    vacuum,       // its material
                                                    "VesselIn",   // its name
                                                    0,            // its FieldManager (will have to change)
                                                    0,            // its sensitive detector
                                                    0);           // its user cuts

  new G4PVPlacement(0,                  // no rotation
                    positionVesselIn,   // at (0,0,0)
                    vesselInLV,         // its logical volume
                    "VesselIn",         // its name
                    vesselOutLV,        // its mother  volume
                    false,              // no boolean operations
                    0,                  // copy number
                    fCheckOverlaps);    // checking overlaps

  
  // POSITIONING OF ELECTRODES
  // Remember that the reference position from which coordinates are calculated is 
  // the center of VesselIn, which is in and by itself at (0., 0., 1120.)
  
  // E01 (Initial electrode)
  G4double E01Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall ) + 
                      nE49     * ( dzE49     + eThicksmall ) + 
                      nE47_E48 * ( dzE47_E48 + eThicksmall ) + 
                      nE43_E46 * ( dzE43_E46 + eThicksmall ) + 
                      nE42     * ( dzE42     + eThicklarge ) + 
                      nE14_E41 * ( dzE14_E41 + eThicklarge ) + 
                      nE02_E13 * ( dzE02_E13 + eThicklarge ) + 
                      nE01     * ( dzE01     + eThicklarge ); 
  G4double         E01_Zpos = -E01Zpos + dzE01;
  G4ThreeVector    posE01   = G4ThreeVector(0,0,E01_Zpos + eThicklarge/2.);
  # ifdef ELECTRODE_ZPOS
  G4cout << "E01 : " << E01_Zpos << G4endl;
  # endif
  
  G4Tubs  *E01S   = new G4Tubs("E01S", eRifinal*mm, eRolarge*mm, (eThicklarge/2)*mm, 0.*deg, 360.*deg);
           fLVE01 = new G4LogicalVolume(E01S, copper, "E01", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVE01, posE01, 0 );
  
  fLVE01->SetVisAttributes(LargeVisAtt);

  
  // E02 - E13 (Large electrodes with large ID)
  G4double E02Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall ) + 
                      nE49     * ( dzE49     + eThicksmall ) + 
                      nE47_E48 * ( dzE47_E48 + eThicksmall ) + 
                      nE43_E46 * ( dzE43_E46 + eThicksmall ) + 
                      nE42     * ( dzE42     + eThicklarge ) + 
                      nE14_E41 * ( dzE14_E41 + eThicklarge ) + 
                      nE02_E13 * ( dzE02_E13 + eThicklarge ); 
  for(G4int copyNo=0; copyNo < nE02_E13; ++copyNo) {
    G4double      E02_E13_Zpos  = -E02Zpos + dzE02_E13 + copyNo*(dzE02_E13 + eThicklarge);
    G4ThreeVector posE02_E13    = G4ThreeVector(0,0,E02_E13_Zpos + eThicklarge/2.);
    # ifdef ELECTRODE_ZPOS
    G4cout << "E" <<  2+copyNo << " : " << E02_E13_Zpos << G4endl;
    # endif
    
    G4Tubs* E02_E13S    = new G4Tubs("E02_E13S", eRilarge*mm, eRolarge*mm, (eThicklarge/2)*mm, 0.*deg, 360.*deg);
    fLVE02_E13[copyNo]  = new G4LogicalVolume(E02_E13S,copper,"E02_E13LV",0,0,0);
    fSpectrometer->AddPlacedVolume( fLVE02_E13[copyNo], posE02_E13, 0 );

    fLVE02_E13[copyNo]->SetVisAttributes(LargeVisAtt);
  }

  
  // E14 - E41 (Large electrodes with small ID)
  G4double E14Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall ) + 
                      nE49     * ( dzE49     + eThicksmall ) + 
                      nE47_E48 * ( dzE47_E48 + eThicksmall ) + 
                      nE43_E46 * ( dzE43_E46 + eThicksmall ) + 
                      nE42     * ( dzE42     + eThicklarge ) + 
                      nE14_E41 * ( dzE14_E41 + eThicklarge );
  for(G4int copyNo=0; copyNo < nE14_E41; ++copyNo) {
    G4double      E14_E41_Zpos  = -E14Zpos + dzE14_E41 + copyNo*(dzE14_E41 + eThicklarge);
    G4ThreeVector posE14_E41    = G4ThreeVector(0,0,E14_E41_Zpos + eThicklarge/2.);
    # ifdef ELECTRODE_ZPOS
    G4cout << "E" << 14+copyNo << " : " << E14_E41_Zpos << G4endl;
    # endif

    G4Tubs* E14_E41S    = new G4Tubs("E14_E41S", eRismall*mm, eRosmall*mm, (eThicklarge/2)*mm, 0.*deg, 360.*deg);
    fLVE14_E41[copyNo]  = new G4LogicalVolume(E14_E41S,copper,"E14_E41LV",0,0,0);
    fSpectrometer->AddPlacedVolume( fLVE14_E41[copyNo], posE14_E41, 0 );
    
    fLVE14_E41[copyNo]->SetVisAttributes(LargeVisAtt);
  }
  
  
  // E42 (Separator between source and ion/ele regions)
  G4double E42Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall ) + 
                      nE49     * ( dzE49     + eThicksmall ) + 
                      nE47_E48 * ( dzE47_E48 + eThicksmall ) + 
                      nE43_E46 * ( dzE43_E46 + eThicksmall ) + 
                      nE42     * ( dzE42     + eThicklarge );
  G4double         E42_Zpos = -E42Zpos + dzE42;
  G4ThreeVector    posE42  = G4ThreeVector(0,0,E42_Zpos + eThicklarge/2.);
//  G4cout << "E42 : " << E42_Zpos << G4endl;
  
  G4Tubs  *E42S   = new G4Tubs("E42S", eRiinit*mm, eRosmall*mm, (eThicklarge/2)*mm, 0.*deg, 360.*deg);
           fLVE42 = new G4LogicalVolume(E42S, copper, "E42", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVE42, posE42, 0 );
  
  fLVE42->SetVisAttributes(LargeVisAtt);

  
  // E43 - E46 (Source electrodes)
  G4double E43Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall ) + 
                      nE49     * ( dzE49     + eThicksmall ) + 
                      nE47_E48 * ( dzE47_E48 + eThicksmall ) + 
                      nE43_E46 * ( dzE43_E46 + eThicksmall );
  for(G4int copyNo=0; copyNo < nE43_E46; ++copyNo) {
    G4double      E43_E46_Zpos  = -E43Zpos + dzE43_E46 + copyNo*(dzE43_E46 + eThicklarge);
    G4ThreeVector posE43_E46    = G4ThreeVector(0,0,E43_E46_Zpos + eThicksmall/2.);
    # ifdef ELECTRODE_ZPOS
    G4cout << "E" << 43+copyNo << " : " << E43_E46_Zpos << G4endl;
    # endif

    G4Tubs* E43_E46S      = new G4Tubs("E43_E46S", eRiinit*mm, eRoinit*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
    fLVE43_E46[copyNo]    = new G4LogicalVolume(E43_E46S,copper,"E43_E46LV",0,0,0);
    fSpectrometer->AddPlacedVolume( fLVE43_E46[copyNo], posE43_E46, 0 );
    
    fLVE43_E46[copyNo]->SetVisAttributes(SmallVisAtt);
  }

  
  // E47 - E48 (Source electrodes)
  G4double E47Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall ) + 
                      nE49     * ( dzE49     + eThicksmall ) + 
                      nE47_E48 * ( dzE47_E48 + eThicksmall );
  for(G4int copyNo=0; copyNo < nE47_E48; ++copyNo) {
    G4double      E47_E48_Zpos  = -E47Zpos + dzE47_E48 + copyNo*(dzE47_E48 + eThicklarge);
    G4ThreeVector posE47_E48    = G4ThreeVector(0,0,E47_E48_Zpos + eThicksmall/2.);
    # ifdef ELECTRODE_ZPOS
    G4cout << "E" << 47+copyNo << " : " << E47_E48_Zpos << G4endl;
    # endif

    G4Tubs* E47_E48S    = new G4Tubs("E47_E48S", eRiinit*mm, eRoinit*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
    fLVE47_E48[copyNo]  = new G4LogicalVolume(E47_E48S,copper,"E47_E48LV",0,0,0);
    fSpectrometer->AddPlacedVolume( fLVE47_E48[copyNo], posE47_E48, 0 );
    
    fLVE47_E48[copyNo]->SetVisAttributes(SmallVisAtt);
  }
  
  
  // E49 (Electrostatic lens)
  G4double E49Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall ) + 
                      nE49     * ( dzE49     + eThicksmall );
  G4double         E49_Zpos = -E49Zpos + dzE49;
  G4ThreeVector    posE49  = G4ThreeVector(0,0,E49_Zpos + eThicksmall/2.);
  # ifdef ELECTRODE_ZPOS
  G4cout << "E49 : " << E49_Zpos << G4endl;
  # endif
  
  G4Tubs  *E49S   = new G4Tubs("E49S", eRiinit*mm, eRoinit*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
           fLVE49 = new G4LogicalVolume(E49S, copper, "E49", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVE49,posE49, 0 );
  
  fLVE49->SetVisAttributes(SmallVisAtt);

  
  // E50 - E52 (Source electrodes)
  G4double E50Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall ) + 
                      nE50_E52 * ( dzE50_E52 + eThicksmall );
  for(G4int copyNo=0; copyNo < nE50_E52; ++copyNo) {
    G4double      E50_E52_Zpos  = -E50Zpos + dzE50_E52 + copyNo*(dzE50_E52 + eThicklarge);
    G4ThreeVector posE50_E52    = G4ThreeVector(0,0,E50_E52_Zpos + eThicksmall/2.);
    # ifdef ELECTRODE_ZPOS
    G4cout << "E" << 50+copyNo << " : " << E50_E52_Zpos << G4endl;
    # endif

    G4Tubs* E50_E52S = new G4Tubs("E50_E52S", eRiinit*mm, eRoinit*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
    fLVE50_E52[copyNo] = new G4LogicalVolume(E50_E52S,copper,"E50_E52LV",0,0,0);
    fSpectrometer->AddPlacedVolume( fLVE50_E52[copyNo], posE50_E52, 0 );
    
    fLVE50_E52[copyNo]->SetVisAttributes(SmallVisAtt);
  }
  
  
  // E53 (Source electrode)
  G4double E53Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall ) + 
                      nE53     * ( dzE53     + eThicksmall );
  G4double         E53_Zpos = -E53Zpos + dzE53;
  G4ThreeVector    posE53 = G4ThreeVector(0,0,E53_Zpos + eThicksmall/2.);
  # ifdef ELECTRODE_ZPOS
  G4cout << "E53 : " << E53_Zpos << G4endl;
  # endif
  
  G4Tubs  *E53S   = new G4Tubs("E53S", eRiinit*mm, eRoinit*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
           fLVE53 = new G4LogicalVolume(E53S, copper, "E53", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVE53, posE53, 0 );
  
  fLVE53->SetVisAttributes(SmallVisAtt);

  
  // E54 - E57 (Source electrodes)
  G4double E54Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall ) + 
                      nE54_E57 * ( dzE54_E57 + eThicksmall );
  for(G4int copyNo=0; copyNo < nE54_E57; ++copyNo) {
    G4double      E54_E57_Zpos  = -E54Zpos + dzE54_E57 + copyNo*(dzE54_E57 + eThicklarge);
    G4ThreeVector posE54_E57    = G4ThreeVector(0,0,E54_E57_Zpos + eThicksmall/2.);
    # ifdef ELECTRODE_ZPOS
    G4cout << "E" << 54+copyNo << " : " << E54_E57_Zpos << G4endl;
    # endif

    G4Tubs* E54_E57S = new G4Tubs("E54_E57S", eRiinit*mm, eRoinit*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
    fLVE54_E57[copyNo] = new G4LogicalVolume(E54_E57S,copper,"E54_E57LV",0,0,0);
    fSpectrometer->AddPlacedVolume( fLVE54_E57[copyNo], posE54_E57, 0 );
    
    fLVE54_E57[copyNo]->SetVisAttributes(SmallVisAtt);
  }
  
  
  // COIL HOLDER - LEFT PIECE
  G4double ECHLeftZpos  = dzSource + 
                          nE58     * ( dzE58     + eThicksmall );
  G4double         ECHLeft_Zpos = -ECHLeftZpos + dzCHLeft;
  G4ThreeVector    posECHLeft = G4ThreeVector(0,0,ECHLeft_Zpos + eCHheight/2.);
  # ifdef ELECTRODE_ZPOS
  G4cout << "ECHLeft : " << ECHLeft_Zpos << G4endl;
  # endif
  
  G4Tubs    *ECHLeftFullS = new G4Tubs("ECHLeftFullS" ,  eCHRi     *mm,  eCHRo    *mm, (  eCHheight     /2)*mm, 0.*deg, 360.*deg);
  G4Tubs    *ECHsubLeftS  = new G4Tubs("ECHsubLeftSub", (eCHRi+1.) *mm, (eCHRo+1.)*mm, ( (eCHheight-2.) /2)*mm, 0.*deg, 360.*deg);
  G4VSolid  *ECHLeftS     = new G4SubtractionSolid("ECHLeftS", ECHLeftFullS, ECHsubLeftS, 0, G4ThreeVector(0.,0.,-2*mm));
             fLVCHLeft    = new G4LogicalVolume(ECHLeftS, copper, "ECHLeft", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVCHLeft, posECHLeft, 0 );
  
  fLVCHLeft->SetVisAttributes(MCPVisAtt);
  
  // COIL HOLDER - EXTRA ELECTRODE
  G4double ECHExtraESZpos = dzSource + 
                            nE58     * ( dzE58     + eThicksmall );
  G4double         ECHExtraES_Zpos = -ECHLeftZpos + dzCHLeft + eCHheight + dzCHExtraE;
  G4ThreeVector    posECHExtraE = G4ThreeVector(0,0,ECHExtraES_Zpos + eThicksmall/2.);
  # ifdef ELECTRODE_ZPOS
  G4cout << "ECHExtraE : " << ECHExtraES_Zpos << G4endl;
  # endif
  
  G4Tubs  *ECHExtraES  = new G4Tubs("ECHExtraES", eRiinit*mm, eRoinitSpecial*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
           fLVCHExtraE = new G4LogicalVolume(ECHExtraES, copper, "ECHExtraE", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVCHExtraE, posECHExtraE, 0 );
  
  fLVCHExtraE->SetVisAttributes(MCPVisAtt);
  
  // COIL HOLDER - RIGHT PIECE
  G4double ECHRightZpos  = dzSource + 
                           nE58     * ( dzE58     + eThicksmall );
  G4double         ECHRight_Zpos = -ECHLeftZpos + dzCHLeft + eCHheight + dzCHExtraE + eThicksmall + dzCHRight;
  G4ThreeVector    posECHRight = G4ThreeVector(0,0,ECHRight_Zpos + eCHheight/2.);
  # ifdef ELECTRODE_ZPOS
  G4cout << "ECHRight : " << ECHRight_Zpos << G4endl;
  # endif
  
  G4Tubs    *ECHRightFullS = new G4Tubs("ECHRightFullS" ,  eCHRi     *mm,  eCHRo    *mm, (  eCHheight     /2)*mm, 0.*deg, 360.*deg);
  G4Tubs    *ECHsubRightS  = new G4Tubs("ECHsubRightSub", (eCHRi+1.) *mm, (eCHRo+1.)*mm, ( (eCHheight-2.) /2)*mm, 0.*deg, 360.*deg);
  G4VSolid  *ECHRightS     = new G4SubtractionSolid("ECHRightS", ECHRightFullS, ECHsubRightS, 0, G4ThreeVector(0.,0., 2*mm));
             fLVCHRight    = new G4LogicalVolume(ECHRightS, copper, "ECHRight", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVCHRight, posECHRight, 0 );
  
  fLVCHRight->SetVisAttributes(MCPVisAtt);
  
  
  // E58 (Source electrode)
  G4double E58Zpos  = dzSource + 
                      nE58     * ( dzE58     + eThicksmall );
  G4double         E58_Zpos = -E58Zpos + dzE58;
  G4ThreeVector    posE58 = G4ThreeVector(0,0,E58_Zpos + eThicksmall/2.);
  # ifdef ELECTRODE_ZPOS
  G4cout << "E58 : " << E58_Zpos << G4endl;
  # endif
  
  G4Tubs  *E58S   = new G4Tubs("E58S", eRiinit*mm, eRoinitSpecial*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
           fLVE58 = new G4LogicalVolume(E58S, copper, "E58", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( fLVE58, posE58, 0 );
  
  fLVE58->SetVisAttributes(SmallVisAtt);
//   
//   fSpectrometer->MakeImprint( vesselInLV, nulltrans, 1);
//   fSpectrometer->MakeImprint( vesselInLV, reflect3D, 2);

//   // Source (TEMPORARY!!!)
//
//   G4ThreeVector positionSource = G4ThreeVector(0,0,0);
//
//
//   G4Sphere* sourceS = new G4Sphere("source",          // name
//                                    0.*mm,             // inner radius
//                                    1*mm,   // outer radius
//                                    0.*deg,            // init phi
//                                    360.*deg,          // end phi
//                                    0.*deg,            // init theta
//                                    180.*deg);         // end theta
//
//   G4LogicalVolume *fLogicSource = new G4LogicalVolume(sourceS,
//                                                       vacuum,
//                                                       "Source",
//                                                       0,0,0);
//
//   new G4PVPlacement(0,               // no rotation
//                     positionSource,  // at (x,y,z)
//                     fLogicSource,    // its logical volume
//                     "Source",      // its name
//                     vesselInLV,      // its mother volume
//                     false,           // no boolean operations
//                     0,               // copy number
//                     fCheckOverlaps); // checking overlaps
//
//   G4cout << "Target is " << 2*targetLength/cm << " cm of "
//          << fTargetMaterial->GetName() << G4endl;

  // MCP box
  // Front
  G4double BoxFrontZpos = E01Zpos + eThicklarge/2.;
  G4ThreeVector  posBoxFront    = G4ThreeVector(0,0,BoxFrontZpos);
  G4Tubs          *BoxFrontS    = new G4Tubs("BoxFrontS", targetRadius*mm, 150.*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
  G4LogicalVolume *BoxFrontLV   = new G4LogicalVolume(BoxFrontS, copper, "BoxFront", 0, 0, 0);
//   new G4PVPlacement(0, posBoxFront, BoxFrontLV, "BoxFront", vesselInLV, false, 0, fCheckOverlaps);
  fSpectrometer->AddPlacedVolume( BoxFrontLV, posBoxFront, 0 );
  BoxFrontLV->SetVisAttributes(MCPBoxVisAtt);

  // Body
  G4double BoxBodyLength = 124.;
  G4double BoxBodyZpos = E01Zpos + 1. + BoxBodyLength/2.;
  G4ThreeVector  posBoxBody    = G4ThreeVector(0,0,BoxBodyZpos);
  G4Tubs          *BoxBodyS    = new G4Tubs("BoxBodyS", 149*mm, 150.*mm, (BoxBodyLength/2)*mm, 0.*deg, 360.*deg);
  G4LogicalVolume *BoxBodyLV   = new G4LogicalVolume(BoxBodyS, copper, "BoxBody", 0, 0, 0);
  fSpectrometer->AddPlacedVolume( BoxBodyLV, posBoxBody, 0 );
//   new G4PVPlacement(0, posBoxBody, BoxBodyLV, "BoxBody", vesselInLV, false, 0, fCheckOverlaps);
  BoxBodyLV->SetVisAttributes(MCPBoxVisAtt);

  // Back
  G4double BoxBackZpos = E01Zpos + 1. + BoxBodyLength + eThicksmall/2.;
  G4ThreeVector  posBoxBack    = G4ThreeVector(0,0,BoxBackZpos);
  G4Tubs          *BoxBackS    = new G4Tubs("BoxBackS", 0.*mm, 150.*mm, (eThicksmall/2)*mm, 0.*deg, 360.*deg);
  G4LogicalVolume *BoxBackLV   = new G4LogicalVolume(BoxBackS, copper, "BoxBack", 0, 0, 0);
//   new G4PVPlacement(0, posBoxBack, BoxBackLV, "BoxBack", vesselInLV, false, 0, fCheckOverlaps);
  fSpectrometer->AddPlacedVolume( BoxBackLV, posBoxBack, 0 );
  BoxBackLV->SetVisAttributes(MCPBoxVisAtt);

  // Target

  G4ThreeVector posMCP = G4ThreeVector(0,0,targetZ*mm);
//   G4ThreeVector positionTarget = G4ThreeVector(0,0,0.);


  G4Tubs* fMCPS = new G4Tubs("MCPS",          // name
                             0.*mm,           // inner radius
                             targetRadius*mm, // outer radius
                             targetLength*mm, // length of target
                             0.*deg,          // init phi
                             360.*deg);       // end phi

  fLVMCP          = new G4LogicalVolume(fMCPS,
//                                         vacuum,
                                        fTargetMaterial,
                                        "MCP",
                                        0,0,0);
  fSpectrometer->AddPlacedVolume( fLVMCP, posMCP, 0 );

//   new G4PVPlacement(0,               // no rotation
//                     positionTarget,  // at (x,y,z)
//                     fLVMCP,    // its logical volume
//                     "MCP",      // its name
//                     vesselInLV,      // its mother volume
//                     false,           // no boolean operations
//                     0,               // copy number
//                     fCheckOverlaps); // checking overlaps
  
  
  fSpectrometer->MakeImprint( vesselInLV, nulltrans, 1);
  fSpectrometer->MakeImprint( vesselInLV, reflect3D, 2);










//
  // TEST User limits
  G4double maxStep = 100*mm;
//   G4double maxStep = 10*MeV;
//   G4double maxStep = 0.1* vesselLength;
  fStepLimit = new G4UserLimits(maxStep);
//   fStepLimit = new G4UserLimits();
//   fStepLimit->SetUserMinEkine(maxStep);
//   vesselInLV->SetUserLimits(fStepLimit);
//   vesselOutLV->SetUserLimits(fStepLimit);
//   fLogicTarget->SetUserLimits(fStepLimit);



  worldLV      ->SetVisAttributes(worldVisAtt);
  vesselOutLV  ->SetVisAttributes(VesselVisAtt);
  vesselInLV   ->SetVisAttributes(VesselVisAtt);
//   fLogicSource ->SetVisAttributes(boxVisAtt);
  fLVMCP ->SetVisAttributes(MCPVisAtt);

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  TETrackerSD* aTrackerSD = new TETrackerSD("TargetSD", "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector("VesselIn", aTrackerSD, true);
  SetSensitiveDetector("VesselOut", aTrackerSD, true);
  SetSensitiveDetector("World", aTrackerSD, true);
//   SetSensitiveDetector("Target", aTrackerSD, true);

#ifdef MAGSTANDARD
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
#else
  // STUFF FROM PURGING_MAGNET
  //Field grid in A9.TABLE. File must be in accessible from run urn directory.
  

//  TEMagTabulatedField3D* TEMagField= new TEMagTabulatedField3D("/home/tug26830/geant4_fullkinematics/xunzhen_map_reader/EMfieldMap_Src-360_FullSpectrometer_Radius220mm_Real0.400Vmm_COMSOL8G.root", zOffset);  // Default field map
//   TEMagTabulatedField3D* TEMagField= new TEMagTabulatedField3D("/home/tug26830/geant4_fullkinematics/xunzhen_map_reader/EMfieldMap_NoSleeves_WithShield_Src-1000_FullSpectrometer_Radius240mm_Real0.400Vmm_COMSOL8G.root", zOffset);  // NoSleeves_WithShield
  
  std::string source_dir = "/home/tug26830/geant4_fullkinematics/xunzhen_map_reader";
  
//   std::string bfield_map = "EMfieldMap_FeFe_NoSleeves_NoShield_Src-1000_1200_Radius240mm_Real0.400Vmm_COMSOL8G.root"                         ;    // 1 - Fe Fe - No Sleeves   - No PMT shield    - QuadInt
//   std::string bfield_map = "EMfieldMap_FeFe_NoSleeves_WithShield_Src-1000_1200_Radius240mm_Real0.400Vmm_COMSOL8G_LinearInterpolation.root"   ;    // 2 - Fe Fe - No Sleeves   - With PMT shield  - LinInt
//   std::string bfield_map = "EMfieldMap_MuFe_NoSleeves_NoShield_Src-1000_1200_Radius240mm_Real0.400Vmm_COMSOL8G.root"                         ;    // 3 - Mu Fe - No Sleeves   - No PMT shield    - Quadint
//   std::string bfield_map = "EMfieldMap_MuFe_NoSleeves_NoShield_Src-1000_1200_Radius240mm_Real0.400Vmm_COMSOL8G_LinearInterpolation.root"     ;    // 33 - Mu Fe - No Sleeves   - No PMT shield    - LinInt
//   std::string bfield_map = "EMfieldMap_MuFe_NoSleeves_WithShield_Src-1000_1200_Radius240mm_Real0.400Vmm_COMSOL8G.root"                       ;    // 4 - Mu Fe - No Sleeves   - With PMT shield  - QuadInt
//   std::string bfield_map = "EMfieldMap_MuFe_WithSleeves_NoShield_Src-1000_1200_Radius240mm_Real0.400Vmm_COMSOL8G_LinearInterpolation.root"   ;    // 5 - Mu Fe - With Sleeves - No PMT shield    - LinInt
//   std::string bfield_map = "EMfieldMap_MuFe_WithSleeves_WithShield_Src-1000_1200_Radius240mm_Real0.400Vmm_COMSOL8G_LinearInterpolation.root" ;    // 6 - Mu Fe - With Sleeves - With PMT shield  - LinInt

//   std::string bfield_map = "EMfieldMap_MuFe_WithSleeves_NoShield_Src-1000_1200_Radius240mm_Real0.400Vmm_pa20210111_COMSOL8G_LinearInterpolation.root"              ;    // 5400 - Mu Fe - With Sleeves - No PMT shield    - LinInt - PMT @ 400mm
//   std::string bfield_map = "EMfieldMap_MuFe_WithSleeves_WithShield_Src-1000_1200_Radius240mm_Real0.400Vmm_pa20210111_COMSOL8G_LinearInterpolation.root"            ;    // 6400 - Mu Fe - With Sleeves - With PMT shield  - LinInt - PMT @ 400mm
   std::string bfield_map = "EMfieldMap_MuFe_WithSleeves_WithShield_PMTat450mm_Src-1000_1200_Radius240mm_Real0.400Vmm_pa20210111_COMSOL8G_LinearInterpolation.root"      ;    // 6450 - Mu Fe - With Sleeves - With PMT shield  - LinInt - PMT @ 450mm
//   std::string bfield_map = "EMfieldMap_MuFe_WithSleeves_WithShield_PMTat500mm_Src-1000_1200_Radius240mm_Real0.400Vmm_pa20210111_COMSOL8G_LinearInterpolation.root" ;    // 6500 - Mu Fe - With Sleeves - With PMT shield  - LinInt - PMT @ 500mm
//   std::string bfield_map = "EMfieldMap_MuFe_WithSleeves_WithShield_PMTat600mm_Src-1000_1200_Radius240mm_Real0.400Vmm_pa20210111_COMSOL8G_LinearInterpolation.root" ;    // 6600 - Mu Fe - With Sleeves - With PMT shield  - LinInt - PMT @ 600mm

  
  char bfield_to_open[300];
  sprintf( bfield_to_open, "%s/%s", source_dir.c_str(), bfield_map.c_str() );
  
  G4cout << "Opening source file:" << G4endl << bfield_to_open << G4endl << G4endl;
  TEMagTabulatedField3D* TEMagField= new TEMagTabulatedField3D(bfield_to_open, zOffset);
  
  
  fBField.Put(TEMagField);

  G4EqMagElectricField* equation = new G4EqMagElectricField(TEMagField);

  //This is thread-local
  G4FieldManager* pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  G4cout<< "DeltaStep "<< pFieldMgr->GetDeltaOneStep()/mm <<"mm" <<G4endl;

//   G4MagIntegratorStepper* stepper = new G4SimpleRunge (equation,12);
//   G4MagIntegratorStepper* stepper = new G4SimpleHeum (equation,12);
  G4MagIntegratorStepper* stepper = new G4ClassicalRK4 (equation,12);
  G4double minStep           = 1e-7*mm;
  G4ChordFinder* chordFinder = new G4ChordFinder((G4MagneticField*)TEMagField,minStep,stepper);

  G4TransportationManager* transportManager = G4TransportationManager::GetTransportationManager();

//   G4PropagatorInField* fieldPropagator = transportManager->GetPropagatorInField();
//
//   G4double epsMin            = 2.5e-7*mm;
//   G4double epsMax            = 2.5e-6*mm;
//
//   fieldPropagator->SetMinimumEpsilonStep(epsMin);
//   fieldPropagator->SetMaximumEpsilonStep(epsMax);
//
  pFieldMgr->SetChordFinder(chordFinder);

//   G4double chordPrecision = 1e-5*mm;
//   pFieldMgr->GetChordFinder()->SetDeltaChord( chordPrecision );
  pFieldMgr->SetDetectorField(fBField.Get());
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

//   if (fTargetMaterial != pttoMaterial) {
//      if ( pttoMaterial ) {
//         fTargetMaterial = pttoMaterial;
//         if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
//         G4cout
//           << G4endl
//           << "----> The target is made of " << materialName << G4endl;
//      } else {
//         G4cout
//           << G4endl
//           << "-->  WARNING from SetTargetMaterial : "
//           << materialName << " not found" << G4endl;
//      }
//   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}
