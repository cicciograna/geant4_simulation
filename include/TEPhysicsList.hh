//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Code developed by:
//  S.Larsson
//
//    ********************************
//    *                              *
//    *    TEPhysicsList.hh     *
//    *                              *
//    ********************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef TEPhysicsList_h
#define TEPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4GammaConversion;     
class G4ComptonScattering;   
class G4PhotoElectricEffect;       

class G4eIonisation;          
class G4eBremsstrahlung;     


class TEPhysicsList: public G4VUserPhysicsList
{
public:
  TEPhysicsList();
  ~TEPhysicsList();
  
protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  
  void SetCuts();
  
public: 
  // Set Cuts
  void SetGammaCut(G4double);
  void SetElectronCut(G4double);
  void SetPositronCut(G4double);
  void SetNeutrinoECut(G4double);
  void SetProtonCut(G4double);
  void SetGenericIonCut(G4double);
  
  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);
  void SetPositronLowLimit(G4double);
  void SetNeutrinoELowLimit(G4double);
  void SetProtonLowLimit(G4double);
  void SetGenericIonLowLimit(G4double);
  void SetGEPLowLimit(G4double);

  void SetGELowLimit(G4double);
  
private:
  
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForNeutrinoE;
  G4double cutForProton;
  G4double cutForGenericIon;
  
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructBarions();
  void ConstructIons();
  
protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  
private:
  
};
#endif



