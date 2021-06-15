/// \file TETrackerHit.hh
/// \brief Definition of the TETrackerHit class

#ifndef TETrackerHit_h
#define TETrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
#include <cmath>
#include "G4SystemOfUnits.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class TETrackerHit : public G4VHit
{
  public:
    TETrackerHit();
    TETrackerHit(const TETrackerHit&);
    virtual ~TETrackerHit();

    // operators
    const TETrackerHit& operator=(const TETrackerHit&);
    G4bool operator==(const TETrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID     (G4int track)               { fTrackID      = track;                    };
    void SetTOF         (G4double tof)              { fTOF          = tof;                      };
    void SetMass        (G4double mass)             { fMass         = mass;                     };
    void SetCharge      (G4double charge)           { fCharge       = charge;                   };
    void SetInitPos     (G4ThreeVector init_xyz)    { fInitPos      = init_xyz;                 };
//     void SetInitAzm     (G4ThreeVector init_pxyz)   { fInitAzm      = atan2(init_pxyz.y(),init_pxyz.x())/deg;          };
//     void SetInitAzm     (G4ThreeVector init_pxyz)   { fInitAzm      = atan2(init_pxyz.x(),init_pxyz.z())/deg;          };
    void SetInitAzm     (G4ThreeVector init_pxyz)   { fInitAzm      = init_pxyz.phi()/deg; if(fInitAzm < 0.){fInitAzm+=360.;}     };
    void SetInitElev    (G4ThreeVector init_pxyz)   { fInitElev     = init_pxyz.theta()/deg;    };
    void SetInitDircos  (G4double init_dircos)      { fInitDircos   = init_dircos;              };
    void SetInitVel     (G4ThreeVector init_vxyz)   { fInitVel      = init_vxyz;                };
    void SetInitMom     (G4ThreeVector init_pxyz)   { fInitMom      = init_pxyz;                };
    void SetInitKE      (G4double init_ke)          { fInitKE       = init_ke;                  };
    void SetSplatPos    (G4ThreeVector splat_xyz)   { fSplatPos     = splat_xyz;                };
    void SetSplatAzm    (G4ThreeVector splat_xyz)  { fSplatAzm     = atan2(splat_xyz.y(),splat_xyz.x())/deg;         };
//     void SetSplatAzm    (G4ThreeVector splat_pxyz)  { fSplatAzm     = splat_pxyz.phi()/deg;         };
    void SetSplatElev   (G4ThreeVector splat_pxyz)  { fSplatElev    = splat_pxyz.theta()/deg;       };
    void SetSplatDircos (G4double splat_dircos)     { fSplatDircos  = splat_dircos;             };
    void SetSplatVel    (G4ThreeVector splat_vxyz)  { fSplatVel     = splat_vxyz;               };
    void SetSplatMom    (G4ThreeVector splat_pxyz)  { fSplatMom     = splat_pxyz;               };
    void SetSplatKE     (G4double splat_ke)         { fSplatKE      = splat_ke;                 };
    void SetEdep        (G4double de)               { fEdep         = de;                       };

    // Get methods
    G4int           GetTrackID()      const { return fTrackID     ;};
    G4double        GetTOF()          const { return fTOF         ;};
    G4double        GetMass()         const { return fMass        ;};
    G4double        GetCharge()       const { return fCharge      ;};
    G4ThreeVector   GetInitPos()      const { return fInitPos     ;};
    G4double        GetInitAzm()      const { return fInitAzm     ;};
    G4double        GetInitElev()     const { return fInitElev    ;};
    G4double        GetInitDircos()   const { return fInitDircos  ;};
    G4ThreeVector   GetInitVel()      const { return fInitVel     ;};
    G4ThreeVector   GetInitMom()      const { return fInitMom     ;};
    G4double        GetInitKE()       const { return fInitKE      ;};
    G4ThreeVector   GetSplatPos()     const { return fSplatPos    ;};
    G4double        GetSplatAzm()     const { return fSplatAzm    ;};
    G4double        GetSplatElev()    const { return fSplatElev   ;};
    G4double        GetSplatDircos()  const { return fSplatDircos ;};
    G4ThreeVector   GetSplatVel()     const { return fSplatVel    ;};
    G4ThreeVector   GetSplatMom()     const { return fSplatMom    ;};
    G4double        GetSplatKE()      const { return fSplatKE     ;};
    G4double        GetEdep()         const { return fEdep        ;};

  private:
    
    G4int           fTrackID      ;
    G4double        fTOF          ;
    G4double        fMass         ;
    G4int           fCharge       ;
    G4ThreeVector   fInitPos      ;
    G4double        fInitAzm      ;
    G4double        fInitElev     ;
    G4double        fInitDircos   ;
    G4ThreeVector   fInitVel      ;
    G4ThreeVector   fInitMom      ;
    G4double        fInitKE       ;
    G4ThreeVector   fSplatPos     ;
    G4double        fSplatAzm     ;
    G4double        fSplatElev    ;
    G4double        fSplatDircos  ;
    G4ThreeVector   fSplatVel     ;
    G4ThreeVector   fSplatMom     ;
    G4double        fSplatKE      ;
    G4double        fEdep         ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TETrackerHit> TETrackerHitsCollection;

extern G4ThreadLocal G4Allocator<TETrackerHit>* TETrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TETrackerHit::operator new(size_t)
{
  if(!TETrackerHitAllocator)
      TETrackerHitAllocator = new G4Allocator<TETrackerHit>;
  return (void *) TETrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TETrackerHit::operator delete(void *hit)
{
  TETrackerHitAllocator->FreeSingle((TETrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
