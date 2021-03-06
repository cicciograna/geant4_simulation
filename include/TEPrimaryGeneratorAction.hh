/// \file TEPrimaryGeneratorAction.hh
/// \brief Definition of the TEPrimaryGeneratorAction class

#ifndef TEPrimaryGeneratorAction_h
#define TEPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "TFile.h"
#include "TTree.h"
# include "/home/tug26830/geant4_fullkinematics/source_particle_class.hh"
# include <fstream>

#include "TEDefinitions.hh"

#include "G4ParticleDefinition.hh"

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class TEPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    TEPrimaryGeneratorAction();
    virtual ~TEPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* );

    G4ParticleGun* GetParticleGun() {return fParticleGun;}

    // Set methods
    void SetRandomFlag(G4bool );

  private:
    G4ParticleGun*          fParticleGun; // G4 particle gun

    TFile *fin;
    TTree *te;
    TTree *tsource_particles;
    
    double uinit_ke    ;
    double uinit_x     ;
    double uinit_y     ;
    double uinit_z     ;
    double uinit_px_dir;
    double uinit_py_dir;
    double uinit_pz_dir;
    
    source_particle *lxenon  = new source_particle;
    source_particle *lxray   = new source_particle;
    source_particle *lnu     = new source_particle;
    source_particle *lauger0 = new source_particle;
    source_particle *lauger1 = new source_particle;
    
    G4ParticleDefinition* xe_plus;
    G4ParticleDefinition* xe_plusplus;
    
    std::ofstream   *out_nparts;
    G4int            nparticles;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
