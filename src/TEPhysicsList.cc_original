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
//    *    TEPhysicsList.cc     *
//    *                              *
//    ********************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TEPhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "G4IonConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TEPhysicsList::TEPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1*millimeter;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForProton    = defaultCutValue;

  SetVerboseLevel(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TEPhysicsList::~TEPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
  ConstructIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructBosons()
{

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructBarions()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructIons()
{
  //  generic ion
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
// Bosons
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
// Leptons
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
// GenericIon
// #include "G4MultipleScattering.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructEM()
{
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() )
    {
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();


      if (particleName == "gamma") {
  //gamma
  pmanager->AddDiscreteProcess(new G4GammaConversion);
  pmanager->AddDiscreteProcess(new G4ComptonScattering);
  pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);

      } else if (particleName == "e-") {
  //electron
  pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
  pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
  pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
//   pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
//   pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
//   pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);

      } else if (particleName == "e+") {
  //positron
  pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
  pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
  pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
  pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);

      } else if (particleName == "proton") {
  //proton
  pmanager->AddProcess(new G4hMultipleScattering,  -1,1,1);

      } else if (particleName == "anti_proton") {
  //antiproton
  pmanager->AddProcess(new G4hMultipleScattering,  -1,1,1);
      
      } else if (particleName == "GenericIon") {
  //generic ion
  pmanager->AddProcess(new G4ionIonisation,        -1, 2,2);
      }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::ConstructGeneral()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "TEPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  SetCutValue(cutForNeutrinoE, "nu_e");

  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton
  SetCutValue(cutForProton, "proton");
  SetCutValue(cutForProton, "anti_proton");

  SetCutValue(cutForGenericIon, "GenericIon");
  
    //  SetCutValueForOthers(defaultCutValue);

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "TEPhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  // G4Gamma::SetEnergyRange(lowcut,1e5);
  SetGELowLimit(lowcut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "TEPhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  // G4Electron::SetEnergyRange(lowcut,1e5);
  SetGELowLimit(lowcut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetPositronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "TEPhysicsList::SetCuts:";
    G4cout << "Positron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  G4cerr << "TEPhysicsList::SetPositronLowLimit: Not currently able to set Positron LowLimit." << G4endl;
  G4Exception("TEPhysicsList::SetPositronLowLimit()","PurMag001",
        FatalException,"Positron Low Limit: not implemented in TEPhysicsList");
  //
  // G4Positron::SetEnergyRange(lowcut,1e5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetNeutrinoELowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "TEPhysicsList::SetCuts:";
    G4cout << "NeutrinoE cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  // G4Electron::SetEnergyRange(lowcut,1e5);
  SetGELowLimit(lowcut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetProtonLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "TEPhysicsList::SetCuts:";
    G4cout << "Proton cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  G4cerr << "TEPhysicsList::SetProtonLowLimit: Not currently able to set Proton LowLimit." << G4endl;
  G4Exception("TEPhysicsList::SetProtonLowLimit()","PurMag002",
        FatalException,"Proton Low Limit: not implemented in TEPhysicsList");
  //
  // G4Proton::SetEnergyRange(lowcut,1e5);
  // G4AntiProton::SetEnergyRange(lowcut,1e5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetGEPLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "TEPhysicsList::SetGEPLowLimit:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  // G4Gamma::SetEnergyRange(lowcut,1e5);
  // G4Electron::SetEnergyRange(lowcut,1e5);
  // G4Positron::SetEnergyRange(lowcut,1e5);
  this->SetGELowLimit(lowcut);

  G4cerr << " SetGEPLowLimit : Uncertain whether setting Positron low limit " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetGenericIonLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "TEPhysicsList::SetCuts:";
    G4cout << "GenericIon cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  G4cerr << "TEPhysicsList::SetGenericIonLowLimit: Not currently able to set GenericIon LowLimit." << G4endl;
  G4Exception("TEPhysicsList::SetGenericIonLowLimit()","PurMag002",
        FatalException,"GenericIon Low Limit: not implemented in TEPhysicsList");
  //
  // G4Proton::SetEnergyRange(lowcut,1e5);
  // G4AntiProton::SetEnergyRange(lowcut,1e5);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void TEPhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "TEPhysicsList::SetGELowLimit:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowcut,1e5);
}
void TEPhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetElectronCut(G4double val)
{
  //  ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetPositronCut(G4double val)
{
  //  ResetCuts();
  cutForPositron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetNeutrinoECut(G4double val)
{
  //  ResetCuts();
  cutForNeutrinoE = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetProtonCut(G4double val)
{
  //ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TEPhysicsList::SetGenericIonCut(G4double val)
{
  //ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






