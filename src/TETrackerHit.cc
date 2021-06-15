/// \file TETrackerHit.cc
/// \brief Implementation of the TETrackerHit class

#include "TETrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<TETrackerHit>* TETrackerHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TETrackerHit::TETrackerHit()
 : G4VHit(),
  fTrackID(-1),
  fMass(0.),
  fCharge(0),
  fTOF(0),
  fInitPos(G4ThreeVector(0.,0.,0.)),
  fInitAzm(0.),
  fInitElev(0.),
  fInitDircos(0.),
  fInitVel(G4ThreeVector(0.,0.,0.)),
  fInitMom(G4ThreeVector(0.,0.,0.)),
  fInitKE(0.),
  fSplatPos(G4ThreeVector(0.,0.,0.)),
  fSplatAzm(0.),
  fSplatElev(0.),
  fSplatDircos(0.),
  fSplatVel(G4ThreeVector(0.,0.,0.)),
  fSplatMom(G4ThreeVector(0.,0.,0.)),
  fSplatKE(0.),
  fEdep(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TETrackerHit::~TETrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TETrackerHit::TETrackerHit(const TETrackerHit& right)
  : G4VHit()
{
  fTrackID      = right.fTrackID      ;
  fMass         = right.fMass         ;
  fCharge       = right.fCharge       ;
  fTOF          = right.fTOF          ;
  fInitPos      = right.fInitPos      ;
  fInitAzm      = right.fInitAzm      ;
  fInitElev     = right.fInitElev     ;
  fInitDircos   = right.fInitDircos   ;
  fInitVel      = right.fInitVel      ;
  fInitMom      = right.fInitMom      ;
  fInitKE       = right.fInitKE       ;
  fSplatPos     = right.fSplatPos     ;
  fSplatAzm     = right.fSplatAzm     ;
  fSplatElev    = right.fSplatElev    ;
  fSplatDircos  = right.fSplatDircos  ;
  fSplatVel     = right.fSplatVel     ;
  fSplatMom     = right.fSplatMom     ;
  fSplatKE      = right.fSplatKE      ;
  fEdep         = right.fEdep         ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TETrackerHit& TETrackerHit::operator=(const TETrackerHit& right)
{
  fTrackID      = right.fTrackID      ;
  fMass         = right.fMass         ;
  fCharge       = right.fCharge       ;
  fTOF          = right.fTOF          ;
  fInitPos      = right.fInitPos      ;
  fInitAzm      = right.fInitAzm      ;
  fInitElev     = right.fInitElev     ;
  fInitDircos   = right.fInitDircos   ;
  fInitVel      = right.fInitVel      ;
  fInitMom      = right.fInitMom      ;
  fInitKE       = right.fInitKE       ;
  fSplatPos     = right.fSplatPos     ;
  fSplatAzm     = right.fSplatAzm     ;
  fSplatElev    = right.fSplatElev    ;
  fSplatDircos  = right.fSplatDircos  ;
  fSplatVel     = right.fSplatVel     ;
  fSplatMom     = right.fSplatMom     ;
  fSplatKE      = right.fSplatKE      ;
  fEdep         = right.fEdep         ;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TETrackerHit::operator==(const TETrackerHit& right) const
{
  return ( this == &right ) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TETrackerHit::Draw()
{
//   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
//   if(pVVisManager)
//   {
//     G4Circle circle(fSplatPos);
//     circle.SetScreenSize(4.);
//     circle.SetFillStyle(G4Circle::filled);
//     if(fInitKE > 110.9 && fInitKE < 111.1){
//       G4Colour colour(0.,1.,0.);
//       G4VisAttributes attribs(colour);
//       circle.SetVisAttributes(attribs);
//     }
//     else if(fInitKE > 121.9 && fInitKE < 122.1){
//       G4Colour colour(1.,0.,0.);
//       G4VisAttributes attribs(colour);
//       circle.SetVisAttributes(attribs);
//     }
//     else {
//       G4Colour colour(0.,0.,1.);
//       G4VisAttributes attribs(colour);
//       circle.SetVisAttributes(attribs);
//     }
//     pVVisManager->Draw(circle);
//   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TETrackerHit::Print()
{
//   G4cout
//      << "  trackID: "                       << fTrackID     << G4endl
//      << "Kinetic energy: "  << std::setw(7) << fInitKE      << G4endl
//      << " TOF: "            << std::setw(7) << fTOF/us      << G4endl
//      << " Init Position: "  << std::setw(7) << fInitPos/mm  << G4endl
//      << " Splat Position: " << std::setw(7) << fSplatPos/mm << G4endl
//      << " Init Momentum: "  << std::setw(7) << fInitMom/keV << G4endl
//      << G4endl
//      << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
