/// \file TEPrimaryGeneratorAction.cc
/// \brief Implementation of the TEPrimaryGeneratorAction class

#include "TEPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "G4IonTable.hh"

#include "G4Tubs.hh"
#include <cstdlib>
#include <cstdio>

#include <set>

#include "Randomize.hh"

// # include "TROOT.h"
# include "TSystem.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TEPrimaryGeneratorAction::TEPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
 fin(NULL),
 te(NULL),
 uinit_ke(0),
 uinit_x(0),
 uinit_y(0),
 uinit_z(0),
 uinit_px_dir(0),
 uinit_py_dir(0),
 uinit_pz_dir(0),
 tsource_particles(NULL),
 lxenon(NULL),
 lxray(NULL),
 lnu(NULL),
 lauger0(NULL),
 lauger1(NULL),
 out_nparts(NULL),
 xe_plus(NULL),
 xe_plusplus(NULL),
 nparticles(0)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic

  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");

  fParticleGun->SetParticleDefinition(particleDefinition);
  // Basic Auger
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,originZpos*mm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.75,1.));
  fParticleGun->SetParticleEnergy(122.0*eV);
  
  #ifdef GEANT4_ROOT
  out_nparts = new std::ofstream( "/home/tug26830/geant4_fullkinematics/local_output/g4root/out_nparts.txt", std::ofstream::out );
  gSystem->Load( "/home/tug26830/geant4_fullkinematics/source_particle_class_hh.so" );
  
  std::string input_file = "/home/tug26830/geant4_fullkinematics/local_output/g4root/fg4root.root";
  
  fin = new TFile( input_file.c_str() ,"OPEN" );
  G4cout << "Input file: " << input_file.c_str() << G4endl;
  
  tsource_particles = (TTree*)fin->Get("tsource_particles");
  
  tsource_particles->SetBranchAddress( "xenon"  ,  &lxenon );
  tsource_particles->SetBranchAddress( "xray"   ,   &lxray );
  tsource_particles->SetBranchAddress( "nu"     ,     &lnu );
  tsource_particles->SetBranchAddress( "auger0" , &lauger0 );
  tsource_particles->SetBranchAddress( "auger1" , &lauger1 );

  G4cout << "The output file ./local_output/g4root/out_particles.txt will contain " << tsource_particles->GetEntries() << " events" << G4endl;
  
  xe_plusplus = new G4ParticleDefinition(
    "XePlusPlus",               //  const G4String &aName,
    121936.2909053164721*MeV,   //  G4double  mass,
    0.0*MeV,                    //  G4double  width,
    2*eplus,                    //  G4double  charge,
    1,                          //  G4int  iSpin,
    1,                          //  G4int  iParity,
    0,                          //  G4int  iConjugation,
    1,                          //  G4int  iIsospin,
    1,                          //  G4int  iIsospinZ,
    0,                          //  G4int  gParity,
    "nucleus",                  //  const G4String &pType,
    0,                          //  G4int  lepton,
    131,                        //  G4int  baryon,
    0,                          //  G4int  encoding,
    true,                       //  G4bool  stable,
    -1.0,                       //  G4double  lifetime,
    NULL,                       //  G4DecayTable *decaytable,
    false,                      //  G4bool  shortlived = false,
    "",                         //  const G4String &subType = ""
    0,                          //  G4int  anti_encoding = 0,
    0.0                         //  G4double  magneticMoment = 0.0
  );
  xe_plus = new G4ParticleDefinition(
    "XePlus",                   //  const G4String &aName,
    121936.8018921231925*MeV,   //  G4double  mass,
    0.0*MeV,                    //  G4double  width,
    eplus,                      //  G4double  charge,
    1,                          //  G4int  iSpin,
    1,                          //  G4int  iParity,
    0,                          //  G4int  iConjugation,
    1,                          //  G4int  iIsospin,
    1,                          //  G4int  iIsospinZ,
    0,                          //  G4int  gParity,
    "nucleus",                  //  const G4String &pType,
    0,                          //  G4int  lepton,
    131,                        //  G4int  baryon,
    0,                          //  G4int  encoding,
    true,                       //  G4bool  stable,
    -1.0,                       //  G4double  lifetime,
    NULL,                       //  G4DecayTable *decaytable,
    false,                      //  G4bool  shortlived = false,
    "",                         //  const G4String &subType = ""
    0,                          //  G4int  anti_encoding = 0,
    0.0                         //  G4double  magneticMoment = 0.0
  );

  
//  for(int i=0;i<tsource_particles->GetEntries();++i){
//  tsource_particles->GetEntry(i);
//  G4cout << " xenon:\t"  <<  lxenon->charge << "\t" <<  lxenon->init_ke << "\t" <<  lxenon->init_x << "\t" <<  lxenon->init_y << "\t" <<  lxenon->init_z << "\t" <<  lxenon->init_px_dir << "\t" <<  lxenon->init_py_dir << "\t" <<  lxenon->init_pz_dir << G4endl;
//  G4cout << "  xray:\t"  <<   lxray->charge << "\t" <<   lxray->init_ke << "\t" <<   lxray->init_x << "\t" <<   lxray->init_y << "\t" <<   lxray->init_z << "\t" <<   lxray->init_px_dir << "\t" <<   lxray->init_py_dir << "\t" <<   lxray->init_pz_dir << G4endl;
//  G4cout << "    nu:\t"  <<     lnu->charge << "\t" <<     lnu->init_ke << "\t" <<     lnu->init_x << "\t" <<     lnu->init_y << "\t" <<     lnu->init_z << "\t" <<     lnu->init_px_dir << "\t" <<     lnu->init_py_dir << "\t" <<     lnu->init_pz_dir << G4endl;
//  G4cout << "auger0:\t"  << lauger0->charge << "\t" << lauger0->init_ke << "\t" << lauger0->init_x << "\t" << lauger0->init_y << "\t" << lauger0->init_z << "\t" << lauger0->init_px_dir << "\t" << lauger0->init_py_dir << "\t" << lauger0->init_pz_dir << G4endl;
//  G4cout << "auger1:\t"  << lauger1->charge << "\t" << lauger1->init_ke << "\t" << lauger1->init_x << "\t" << lauger1->init_y << "\t" << lauger1->init_z << "\t" << lauger1->init_px_dir << "\t" << lauger1->init_py_dir << "\t" << lauger1->init_pz_dir << G4endl;
//  G4cout << G4endl;
//  }
  # elif FIXED_PARTICLES
//   std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_standardE_100k_felev_30.0_PointSrc_20201214.root";
   std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_realProbabilities_GaussSrc_20210106.root";    // DEFAULT TO GENERATE GAUSSIAN-DISTRIBUTED AUGERS
  
//   std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_minE90_felev_5-15-30-60_GaussSrc_20201029.root";
//   std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_GaussSrc_20201014.root";
//  std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_10k_felev_30.0_PointSrc_20201214.root";
//  std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_minE20-maxE130_100k_felev_30.0_PointSrc_20201214.root";

//  std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_standardE_100k_felev_15.0_PointSrc_20201214.root";
//  std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_standardE_100k_felev_30.0_PointSrc_20201214.root";
//  std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_standardE_100k_felev_60.0_PointSrc_20201214.root";

//   std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_standardE_100k_felev_24.13_PointSrc_20210223.root";
//   std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_standardE_100k_felev_31.31_PointSrc_20210226.root";
//   std::string input_file = "/home/tug26830/geant4_fullkinematics/generate_electrons/source_electrons_standardE_100k_felev_43.12_PointSrc_20210223.root";
  
  fin = new TFile( input_file.c_str(), "OPEN");

  G4cout << "Input file: " << input_file.c_str() << G4endl;

  tsource_particles = (TTree*)fin->Get("tsource_particles");
  tsource_particles->SetBranchAddress("init_ke"    , &uinit_ke    );
  tsource_particles->SetBranchAddress("init_x"     , &uinit_x     );
  tsource_particles->SetBranchAddress("init_y"     , &uinit_y     );
  tsource_particles->SetBranchAddress("init_z"     , &uinit_z     );
  tsource_particles->SetBranchAddress("init_px_dir", &uinit_px_dir);
  tsource_particles->SetBranchAddress("init_py_dir", &uinit_py_dir);
  tsource_particles->SetBranchAddress("init_pz_dir", &uinit_pz_dir);
  #endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TEPrimaryGeneratorAction::~TEPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4double vesselZHalfLength = 0;
  G4LogicalVolume* vesselLV = G4LogicalVolumeStore::GetInstance()->GetVolume("VesselIn");
  G4Tubs* vesselTube = NULL;
  if ( vesselLV ) vesselTube = dynamic_cast<G4Tubs*>(vesselLV->GetSolid());
  if ( vesselTube ) vesselZHalfLength = vesselTube->GetZHalfLength();
  else  {
    G4cerr << "Vessel volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }
  
#ifdef GEANT4_ROOT
G4double naug = 1;
  // Fires particles from ROOT file
  for(int i=0;i<tsource_particles->GetEntries();++i){
    naug = 1;
    tsource_particles->GetEntry(i);
    nparticles = 0;
        
    //////////////
    /// AUGERS ///
    //////////////
    
    
    // Definition and firing of lauger1
    if(lauger1->init_ke > 1e-10){  // checks if lauger1 exists
    naug = 2;
//     nparticles = 5;
    fParticleGun->SetParticleDefinition( G4ParticleTable::GetParticleTable()->FindParticle("e-") ) ;
    fParticleGun->SetParticleEnergy( lauger1->init_ke*keV ) ;
    fParticleGun->SetParticlePosition( G4ThreeVector( lauger1->init_x*mm, lauger1->init_y*mm, lauger1->init_z*mm ) );
    fParticleGun->SetParticleMomentumDirection( G4ThreeVector( lauger1->init_px_dir, lauger1->init_py_dir, lauger1->init_pz_dir ) );
    fParticleGun->GeneratePrimaryVertex(anEvent);
    ++nparticles;
    }
    
    // Definition and firing of lauger0
    fParticleGun->SetParticleDefinition( G4ParticleTable::GetParticleTable()->FindParticle("e-") ) ;
    fParticleGun->SetParticleEnergy( lauger0->init_ke*keV ) ;
    fParticleGun->SetParticlePosition( G4ThreeVector( lauger0->init_x*mm, lauger0->init_y*mm, lauger0->init_z*mm ) );
    fParticleGun->SetParticleMomentumDirection( G4ThreeVector( lauger0->init_px_dir, lauger0->init_py_dir, lauger0->init_pz_dir ) );
    fParticleGun->GeneratePrimaryVertex(anEvent);
    ++nparticles;
            
    ////////////////
    /// NEUTRINO ///
    ////////////////
    
    // Definition and firing
    fParticleGun->SetParticleDefinition( G4ParticleTable::GetParticleTable()->FindParticle("nu_e") ) ;
    fParticleGun->SetParticleEnergy( lnu->init_ke*keV ) ;
    fParticleGun->SetParticlePosition( G4ThreeVector( lnu->init_x*mm, lnu->init_y*mm, lnu->init_z*mm ) );
    fParticleGun->SetParticleMomentumDirection( G4ThreeVector( lnu->init_px_dir, lnu->init_py_dir, lnu->init_pz_dir ) );
    fParticleGun->GeneratePrimaryVertex(anEvent);
    ++nparticles;
            
    ////////////
    /// XRAY ///
    ////////////
    
    // Definition and firing
    fParticleGun->SetParticleDefinition( G4ParticleTable::GetParticleTable()->FindParticle("gamma") ) ;
    fParticleGun->SetParticleEnergy( lxray->init_ke*keV ) ;
    fParticleGun->SetParticlePosition( G4ThreeVector( lxray->init_x*mm, lxray->init_y*mm, lxray->init_z*mm ) );
    fParticleGun->SetParticleMomentumDirection( G4ThreeVector( lxray->init_px_dir, lxray->init_py_dir, lxray->init_pz_dir ) );
    fParticleGun->GeneratePrimaryVertex(anEvent);
    ++nparticles;
    
    /////////////
    /// XENON ///
    /////////////
    
    // Definition
    G4int xe_Z = 54, xe_A = 131;
//    G4double xe_charge = 1.*eplus; if( lauger1->init_ke > 1. ){ G4cout << "boggio" << G4endl; xe_charge = 2.*eplus; }
    G4double xe_charge = naug*eplus;
    G4double xe_excitEnergy = 0.*keV;
    

    if(lauger1->init_ke > 1e-10){
    fParticleGun->SetParticleDefinition(xe_plusplus);
    }
    else{
    fParticleGun->SetParticleDefinition(xe_plus);
    }
    
    // Firing
    fParticleGun->SetParticleEnergy( lxenon->init_ke*keV );
    fParticleGun->SetParticlePosition( G4ThreeVector( lxenon->init_x*mm, lxenon->init_y*mm, lxenon->init_z*mm ) );
    fParticleGun->SetParticleMomentumDirection( G4ThreeVector( lxenon->init_px_dir, lxenon->init_py_dir, lxenon->init_pz_dir ) );
    fParticleGun->GeneratePrimaryVertex(anEvent);
    ++nparticles;

    // Writing into out_nparts
    *out_nparts << nparticles << G4endl;
  }
#elif ISO_CURVES_GENERATOR
  // Fires multiple particles at increasing theta/phi and energy (or momentum)
  G4double   azimuth_step = 1; // degrees
  G4double elevation_step = 1; // degrees
  G4int totelevations = 180/elevation_step + 1;     // DOES NOT GENERATE UNIFORMLY DISTRIBUTED ELEVATIONS
  G4int   totazimuths = 360/azimuth_step   + 1;     // DEGREE-BASED METHOD

  for(int en=1;en<=150;en+=1){
    fParticleGun->SetParticleEnergy(en*eV);

    for(int j=0;j<totazimuths;++j){
      G4double phi = 2.*CLHEP::pi*( ((G4double)j/(G4double)(totazimuths-1)) );
//       G4double phi = j*azimuth_step;
      for(int i=0;i<totelevations;++i){
        G4double theta = CLHEP::pi*( ((G4double)i/(G4double)(totelevations-1)) );
//         G4double theta = acos(startacos + i*elevation_step);
        fParticleGun->SetParticleMomentumDirection( G4ThreeVector( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) ) );
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }
  }
#elif ENERGIES_RANGE_ISO_CURVES_GENERATOR
  // Fires multiple particles at increasing theta/phi and energy (or momentum)
  G4double   azimuth_step = 1; // degrees
  G4double elevation_step = 0.5; // degrees
  G4int totelevations = 180/elevation_step + 1;     // DOES NOT GENERATE UNIFORMLY DISTRIBUTED ELEVATIONS
  G4int   totazimuths = 360/azimuth_step   + 1;     // DEGREE-BASED METHOD

  int energy_range = 5; // amplitude of energy range centered on standard energies (e.g., if 5 around 23, it will be 21, 22, 23, 24, 25)
  G4double energies[] = {23, 34, 45, 54, 65, 100, 111, 122};
  std::set<double> *set_energies = new std::set<double>;
  
  int en_mult = 4;
  for(int i=0; i<sizeof(energies)/sizeof(double); ++i){
    for(int j=0; j<en_mult*energy_range; ++j){ set_energies->insert( energies[i] - ( energy_range ) / 2. + (double)j/en_mult ); }
  }

  for(std::set<double>::iterator it = set_energies->begin() ; it != set_energies->end() ; ++it){
    fParticleGun->SetParticleEnergy(*it*eV);

    for(int j=0;j<totazimuths;++j){
      G4double phi = 2.*CLHEP::pi*( ((G4double)j/(G4double)(totazimuths-1)) );
//       G4double phi = j*azimuth_step;
      for(int i=0;i<totelevations;++i){
        G4double theta = CLHEP::pi*( ((G4double)i/(G4double)(totelevations-1)) );
//         G4double theta = acos(startacos + i*elevation_step);
        fParticleGun->SetParticleMomentumDirection( G4ThreeVector( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) ) );
        fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }
  }
  # elif FIXED_PARTICLES
  G4cout << tsource_particles->GetEntries() << " electron inital positions" << G4endl;
  int aaa = 0;
  // Fires particles from ROOT file
  for(int i=tsource_particles->GetEntries()-1;i>=0;--i){
    tsource_particles->GetEntry(i);
    fParticleGun->SetParticleEnergy(uinit_ke*eV);
    fParticleGun->SetParticlePosition( G4ThreeVector(uinit_x*mm, uinit_y*mm, uinit_z*mm) );
    fParticleGun->SetParticleMomentumDirection( G4ThreeVector(uinit_px_dir, uinit_py_dir, uinit_pz_dir) );
    if(aaa < 5){G4cout << uinit_ke << " " << uinit_x << " " << uinit_y << " " << uinit_z << " " << uinit_px_dir << " " << uinit_py_dir << " " << uinit_pz_dir << G4endl;}
    ++aaa;
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  
#else
  fParticleGun->GeneratePrimaryVertex(anEvent);
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
