/// \file TEActionInitialization.cc
/// \brief Implementation of the TEActionInitialization class

#include "TEActionInitialization.hh"
#include "TEPrimaryGeneratorAction.hh"
#include "TERunAction.hh"
#include "TEEventAction.hh"
#include "TESteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TEActionInitialization::TEActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TEActionInitialization::~TEActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEActionInitialization::BuildForMaster() const
{
  SetUserAction(new TERunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TEActionInitialization::Build() const
{
  SetUserAction(new TEPrimaryGeneratorAction);
  TERunAction* runAction = new TERunAction;
  SetUserAction(runAction);
  
  TEEventAction* eventAction = new TEEventAction;
  SetUserAction(eventAction);
  
  SetUserAction(new TESteppingAction(eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
