/// \file TERunAction.hh
/// \brief Definition of the TERunAction class

#ifndef TERunAction_h
#define TERunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

/// Run action class

class TERunAction : public G4UserRunAction
{
  public:
    TERunAction();
    virtual ~TERunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void   EndOfRunAction(const G4Run* run);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
