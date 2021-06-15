/// \file TETrackerSD.cc
/// \brief Implementation of the TETrackerSD class

#include "TETrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TETrackerSD::TETrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);    // collectionName is a member of G4VSensitiveDetector ( http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4VSensitiveDetector_8hh-source.html#l00099 )
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TETrackerSD::~TETrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TETrackerSD::Initialize(G4HCofThisEvent* hce)    // Initialize is a member of G4VSensitiveDetector ( http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4VSensitiveDetector_8hh-source.html#l00069 )
{
  // Create hits collection

  fHitsCollection = new TETrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TETrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)   // ProcessHits is a member of G4VSensitiveDetector ( http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4VSensitiveDetector_8hh-source.html#l00085 )
{
# if defined(ISO_CURVES_GENERATOR) || defined(ENERGIES_RANGE_ISO_CURVES_GENERATOR) || defined(FIXED_PARTICLES)
//   double mmus= mm/microsecond;
  G4double c = 299792.46;
  // energy deposit
  G4double edep = aStep->GetPreStepPoint()->GetKineticEnergy();
//   G4double edep = aStep->GetTotalEnergyDeposit();
//   G4cout << "TETrackerHit: " << edep << G4endl;

  if (edep==0.) return false;                                                                
  G4int         TrackID   = aStep->GetTrack()        ->GetTrackID();                     
  G4double      Mass      = aStep->GetTrack()        ->GetDynamicParticle()->GetMass();  
  G4int         Charge    = aStep->GetTrack()        ->GetDynamicParticle()->GetCharge();
  G4double      TOF       = aStep->GetTrack()        ->GetLocalTime();                   
  G4ThreeVector InitPos   = aStep->GetPreStepPoint() ->GetPosition();                    
  G4ThreeVector InitMom   = aStep->GetPreStepPoint() ->GetMomentum();                    
  G4double      InitKE    = aStep->GetPreStepPoint() ->GetKineticEnergy();               
  G4ThreeVector SplatPos  = aStep->GetPostStepPoint()->GetPosition();                    
  G4ThreeVector SplatMom  = aStep->GetPostStepPoint()->GetMomentum();                    
  G4double      SplatKE   = aStep->GetPostStepPoint()->GetKineticEnergy();               
                                                                                            
  G4double      InitDircos  = InitMom.z()/InitMom.mag();
  G4double      SplatDircos = SplatMom.z()/SplatMom.mag();

  TETrackerHit* newHit = new TETrackerHit();
  newHit->SetTrackID     ( TrackID          );
  newHit->SetMass        ( Mass     /keV    );
  newHit->SetCharge      ( Charge   /eplus  );
  newHit->SetTOF         ( TOF      /us     );  // Check https://scicomex.com/global-local-and-proper-geant4-time/ for difference between Global/Local/Proper Time
  newHit->SetInitPos     ( InitPos  /mm     );
  newHit->SetInitAzm     ( InitMom          );
  newHit->SetInitElev    ( InitMom          );
  newHit->SetInitDircos  ( InitDircos       );
  newHit->SetInitVel     ( InitMom/Mass*c   );
  newHit->SetInitMom     ( InitMom  /keV    );
  newHit->SetInitKE      ( InitKE   /eV     );
  newHit->SetSplatPos    ( SplatPos /mm     );
  newHit->SetSplatAzm    ( SplatPos         );
  newHit->SetSplatElev   ( SplatMom         );
  newHit->SetSplatDircos ( SplatDircos      );
  newHit->SetSplatVel    ( SplatMom/Mass*c  );
  newHit->SetSplatMom    ( SplatMom /keV    );
  newHit->SetSplatKE     ( SplatKE  /eV     );

  fHitsCollection->insert( newHit );

  //newHit->Print();
# endif
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TETrackerSD::EndOfEvent(G4HCofThisEvent*)    // EndOfEvent is a member of G4VSensitiveDetector ( http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4VSensitiveDetector_8hh-source.html#l00070 )
{
  if ( verboseLevel>-1 ) { 
//      G4int nofHits = fHitsCollection->entries();
//      G4cout << G4endl
//             << "-------->Hits Collection: in this event they are " << nofHits 
//             << " hits in the tracker chambers: " << G4endl;
//      for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
