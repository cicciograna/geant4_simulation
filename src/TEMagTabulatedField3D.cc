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
//  S.Larsson and J. Generowicz.
//
//    *************************************
//    *                                   *
//    *    TEMagTabulatedField3D.cc     *
//    *                                   *
//    *************************************
//
//

#include "TEMagTabulatedField3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

# include "TFile.h"
# include "TTree.h"

namespace{
  G4Mutex myTEMagTabulatedField3DLock = G4MUTEX_INITIALIZER;
}

using namespace std;

TEMagTabulatedField3D::TEMagTabulatedField3D(const char* filename, double zOffset ) 
  :fZoffset(zOffset),invertX(false),invertY(false),invertZ(false)
{    
  double lenUnit= millimeter;
  double BfieldUnit= gauss; 
  double EfieldUnit= volt/millimeter; 
  G4cout  << "\n-----------------------------------------------------------"
          << "\n      Magnetic field"
          << "\n-----------------------------------------------------------";
    
  G4cout << "\n ---> " "Reading the field grid from " << filename << " ... " << G4endl; 

  //
  //This is a thread-local class and we have to avoid that all workers open the 
  //file at the same time
  G4AutoLock lock(&myTEMagTabulatedField3DLock);

  TFile *file = new TFile(filename);
  
//   if (!file.is_open())
  if (!file->IsOpen())
    {
      G4ExceptionDescription ed;
      ed << "Could not open input file " << filename << std::endl;
      G4Exception("TEMagTabulatedField3D::TEMagTabulatedField3D", "pugmag001",FatalException,ed);
    }
  
  // Read table dimensions 
  double txpos, typos, tzpos;
  double tbx, tby, tbz;
  double tex, tey, tez;
  
  TTree *tout = (TTree*)file->Get("tout");
  tout->SetBranchAddress("xpos", &txpos);
  tout->SetBranchAddress("ypos", &typos);
  tout->SetBranchAddress("zpos", &tzpos);
  tout->SetBranchAddress("bx"  , &tbx  );
  tout->SetBranchAddress("by"  , &tby  );
  tout->SetBranchAddress("bz"  , &tbz  );
  tout->SetBranchAddress("ex"  , &tex  );
  tout->SetBranchAddress("ey"  , &tey  );
  tout->SetBranchAddress("ez"  , &tez  );
  tout->GetEntry(0);
  nx = (int)txpos;
  ny = (int)typos;
  nz = (int)tzpos;
  
  G4cout << "  [ Number of values x,y,z: "  << nx << " " << ny << " " << nz << " ] " << G4endl;

  // Set up storage space for table
  xBField.resize( nx );
  yBField.resize( nx );
  zBField.resize( nx );
  xEField.resize( nx );
  yEField.resize( nx );
  zEField.resize( nx );
  int ix, iy, iz;
  for (ix=0; ix<nx; ix++) {
    xBField[ix].resize(ny);
    yBField[ix].resize(ny);
    zBField[ix].resize(ny);
    xEField[ix].resize(ny);
    yEField[ix].resize(ny);
    zEField[ix].resize(ny);
    for (iy=0; iy<ny; iy++) {
      xBField[ix][iy].resize(nz);
      yBField[ix][iy].resize(nz);
      zBField[ix][iy].resize(nz);
      xEField[ix][iy].resize(nz);
      yEField[ix][iy].resize(nz);
      zEField[ix][iy].resize(nz);
    }
  }
  
  // Read in the data
  G4int   origin_z = 1120;    // NOTE!!! THIS VALUE MUST BE THE SAME USED IN XUNZHEN FIELD MAP CREATOR!!!
  G4int ion_offset = 1000;    // NOTE!!! THIS VALUE MUST BE THE SAME USED IN XUNZHEN FIELD MAP CREATOR!!!
//  G4int ion_offset = 128;    // NOTE!!! THIS VALUE MUST BE THE SAME USED IN XUNZHEN FIELD MAP CREATOR!!!
  G4int delta_array = origin_z - ion_offset;
  
  G4cout << tout->GetEntries() << G4endl;
  for (long numb = 1;numb < tout->GetEntries();++numb){
//   for (long numb = 1;numb < 966;++numb){
    tout->GetEntry(numb);
    if((int)txpos > maxx){maxx = (int)txpos;} if((int)txpos < minx){minx = (int)txpos;} 
    if((int)typos > maxy){maxy = (int)typos;} if((int)typos < miny){miny = (int)typos;} 
    if((int)tzpos > maxz){maxz = (int)tzpos;} if((int)tzpos < minz){minz = (int)tzpos;}
//     G4cout << numb << "\t" << txpos << "\t" << (int)txpos+(nx)/2 << "\t" << typos << "\t" << (int)typos+(ny)/2 << "\t" << tzpos << "\t" << (int)tzpos-992 << "\t"/* << G4endl*/;
    xBField[(int)txpos+(nx)/2][(int)typos+(ny)/2][(int)tzpos-delta_array] = tbx * BfieldUnit;  // THAT 10 MUST BE EQUAL TO "origin_z - ion_offset"
    yBField[(int)txpos+(nx)/2][(int)typos+(ny)/2][(int)tzpos-delta_array] = tby * BfieldUnit;  // THAT 10 MUST BE EQUAL TO "origin_z - ion_offset"
    zBField[(int)txpos+(nx)/2][(int)typos+(ny)/2][(int)tzpos-delta_array] = tbz * BfieldUnit;  // THAT 10 MUST BE EQUAL TO "origin_z - ion_offset"
    xEField[(int)txpos+(nx)/2][(int)typos+(ny)/2][(int)tzpos-delta_array] = tex * EfieldUnit;  // THAT 10 MUST BE EQUAL TO "origin_z - ion_offset"
    yEField[(int)txpos+(nx)/2][(int)typos+(ny)/2][(int)tzpos-delta_array] = tey * EfieldUnit;  // THAT 10 MUST BE EQUAL TO "origin_z - ion_offset"
    zEField[(int)txpos+(nx)/2][(int)typos+(ny)/2][(int)tzpos-delta_array] = tez * EfieldUnit;  // THAT 10 MUST BE EQUAL TO "origin_z - ion_offset"
  }
  file->Close();

  lock.unlock();

  minx *= lenUnit; maxx *= lenUnit;
  miny *= lenUnit; maxy *= lenUnit;
  minz *= lenUnit; maxz *= lenUnit;

  G4cout << "\n ---> ... done reading " << G4endl;

  // G4cout << " Read values of field from file " << filename << G4endl; 
  G4cout  << " ---> assumed the order:  x, y, z, Bx, By, Bz "
          << "\n ---> Min values x,y,z: " << minx/mm << " " << miny/mm << " " << minz/mm << " mm "
          << "\n ---> Max values x,y,z: " << maxx/mm << " " << maxy/mm << " " << maxz/mm << " mm "
          << "\n ---> The field will be offset by " << zOffset/mm << " mm " << G4endl;

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;} 
  if (maxy < miny) {swap(maxy,miny); invertY = true;} 
  if (maxz < minz) {swap(maxz,minz); invertZ = true;} 
  G4cout  << "\nAfter reordering if neccesary"
          << "\n ---> Min values x,y,z: " << minx/mm << " " << miny/mm << " " << minz/mm << " mm "
          << "\n ---> Max values x,y,z: " << maxx/mm << " " << maxy/mm << " " << maxz/mm << " mm ";

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout  << "\n ---> Dif values x,y,z (range): " << dx/mm << " " << dy/mm << " " << dz/mm << " mm in z "
          << "\n-----------------------------------------------------------" << G4endl;
}

void TEMagTabulatedField3D::GetFieldValue(const double point[4], double *Bfield ) const
{
  double x = point[0];
  double y = point[1];
  double z = point[2] + fZoffset;

  // Check that the point is within the defined region 
  if ( x>=minx && x<=maxx &&
       y>=miny && y<=maxy && 
       z>=minz && z<=maxz ) {
    
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy; 
    double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    

    

// #ifdef DEBUG_INTERPOLATING_FIELD
//     G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << G4endl;
//     G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << G4endl;
//     double valx0z0, mulx0z0, valx1z0, mulx1z0;
//     double valx0z1, mulx0z1, valx1z1, mulx1z1;
//     valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
//     valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
//     valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
//     valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
// #endif

    // Full 3-dimensional version (see http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4ElectroMagneticField_8hh-source.html )
    // B field
    Bfield[0] =
      xBField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xBField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xBField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xBField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xBField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xBField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xBField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xBField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[1] =
      yBField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yBField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yBField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yBField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yBField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yBField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yBField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yBField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[2] =
      zBField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zBField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zBField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zBField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zBField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zBField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zBField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zBField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
      
    // E field  
    Bfield[3] =
      xEField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xEField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xEField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xEField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xEField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xEField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xEField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xEField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[4] =
      yEField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yEField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yEField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yEField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yEField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yEField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yEField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yEField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[5] =
      zEField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zEField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zEField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zEField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zEField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zEField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zEField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zEField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

#ifdef DEBUG_MAGNETIC_FIELD
    G4cout << "Basic x,y,z: " << x      << " " << y      << " " << z      << G4endl;
    G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << G4endl;
    G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << G4endl;
    
    G4cout << "xBField = " << Bfield[0] << " = " << G4endl;
    G4cout << "xBField[" << xindex   << "][" << yindex   << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xBField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "xBField[" << xindex   << "][" << yindex   << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << xBField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "xBField[" << xindex   << "][" << yindex+1 << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xBField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "xBField[" << xindex   << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << xBField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  << G4endl;
    G4cout << "xBField[" << xindex+1 << "][" << yindex   << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xBField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "xBField[" << xindex+1 << "][" << yindex   << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << xBField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "xBField[" << xindex+1 << "][" << yindex+1 << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xBField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "xBField[" << xindex+1 << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << xBField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal  << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
    
    G4cout << "yBField = " << Bfield[1] << " = " << G4endl;
    G4cout << "yBField[" << xindex   << "][" << yindex   << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yBField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "yBField[" << xindex   << "][" << yindex   << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << yBField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "yBField[" << xindex   << "][" << yindex+1 << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yBField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "yBField[" << xindex   << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << yBField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  << G4endl;
    G4cout << "yBField[" << xindex+1 << "][" << yindex   << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yBField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "yBField[" << xindex+1 << "][" << yindex   << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << yBField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "yBField[" << xindex+1 << "][" << yindex+1 << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yBField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "yBField[" << xindex+1 << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << yBField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal  << G4endl;
    
    G4cout << "zBField = " << Bfield[2] << " = " << G4endl;
    G4cout << "zBField[" << xindex   << "][" << yindex   << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zBField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "zBField[" << xindex   << "][" << yindex   << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << zBField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "zBField[" << xindex   << "][" << yindex+1 << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zBField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "zBField[" << xindex   << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << zBField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  << G4endl;
    G4cout << "zBField[" << xindex+1 << "][" << yindex   << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zBField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "zBField[" << xindex+1 << "][" << yindex   << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << zBField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "zBField[" << xindex+1 << "][" << yindex+1 << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zBField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "zBField[" << xindex+1 << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << zBField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal  << G4endl;
#endif
    
#ifdef DEBUG_ELECTRIC_FIELD
    G4cout << "xEField = " << Bfield[3] << " = " << G4endl;
    G4cout << "xEField[" << xindex   << "][" << yindex   << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xEField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "xEField[" << xindex   << "][" << yindex   << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << xEField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "xEField[" << xindex   << "][" << yindex+1 << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xEField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "xEField[" << xindex   << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << xEField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  << G4endl;
    G4cout << "xEField[" << xindex+1 << "][" << yindex   << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xEField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "xEField[" << xindex+1 << "][" << yindex   << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << xEField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "xEField[" << xindex+1 << "][" << yindex+1 << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << xEField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "xEField[" << xindex+1 << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << xEField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal  << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
    
    G4cout << "yEField = " << Bfield[4] << " = " << G4endl;
    G4cout << "yEField[" << xindex   << "][" << yindex   << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yEField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "yEField[" << xindex   << "][" << yindex   << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << yEField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "yEField[" << xindex   << "][" << yindex+1 << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yEField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "yEField[" << xindex   << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << yEField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  << G4endl;
    G4cout << "yEField[" << xindex+1 << "][" << yindex   << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yEField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "yEField[" << xindex+1 << "][" << yindex   << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << yEField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "yEField[" << xindex+1 << "][" << yindex+1 << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << yEField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "yEField[" << xindex+1 << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << yEField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal  << G4endl;
    
    G4cout << "zEField = " << Bfield[5] << " = " << G4endl;
    G4cout << "zEField[" << xindex   << "][" << yindex   << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zEField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "zEField[" << xindex   << "][" << yindex   << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << zEField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "zEField[" << xindex   << "][" << yindex+1 << "][" << zindex   << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zEField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "zEField[" << xindex   << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (1-xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << zEField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  << G4endl;
    G4cout << "zEField[" << xindex+1 << "][" << yindex   << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zEField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) << G4endl;
    G4cout << "zEField[" << xindex+1 << "][" << yindex   << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (1-ylocal) << " ) * ( " << (  zlocal) << " ) = " << zEField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  << G4endl;
    G4cout << "zEField[" << xindex+1 << "][" << yindex+1 << "][" << zindex   << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (1-zlocal) << " ) = " << zEField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) << G4endl;
    G4cout << "zEField[" << xindex+1 << "][" << yindex+1 << "][" << zindex+1 << "] * ( " << (  xlocal) << " ) * ( " << (  ylocal) << " ) * ( " << (  zlocal) << " ) = " << zEField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal  << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
    G4cout << "----------------------------------------------------------------------------------" << G4endl;
#endif
    
  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
    Bfield[3] = 0.0;
    Bfield[4] = 0.0;
    Bfield[5] = 0.0;
  }
}

