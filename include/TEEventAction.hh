/// \file TEEventAction.hh
/// \brief Definition of the TEEventAction class

#ifndef TEEventAction_h
#define TEEventAction_h 1

#include "G4UserEventAction.hh"

#include "G4THitsMap.hh"

#include "globals.hh"

#include "TEDefinitions.hh"

/// Event action class

class TEEventAction : public G4UserEventAction
{
  public:
    TEEventAction();
    virtual ~TEEventAction();

    virtual void  BeginOfEventAction(const G4Event* );
    virtual void    EndOfEventAction(const G4Event* );
    
    G4int index;
    
  private:
    G4THitsMap<G4double>* GetHitsCollection(G4int hcID, const G4Event* event) const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
