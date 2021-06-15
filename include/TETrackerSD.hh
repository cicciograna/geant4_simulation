/// \file TETrackerSD.hh
/// \brief Definition of the TETrackerSD class

#ifndef TETrackerSD_h
#define TETrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4StepStatus.hh"

#include "TETrackerHit.hh"

#include <vector>

# include <fstream>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// TETracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class TETrackerSD : public G4VSensitiveDetector
{
  public:
    TETrackerSD(const G4String& name, 
                const G4String& hitsCollectionName);
    virtual ~TETrackerSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    TETrackerHitsCollection* fHitsCollection;
    std::ofstream   *out_nparts;
    std::ofstream   *out_particles;
    G4int            nparticles;
    G4StepStatus    stepStatus;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
