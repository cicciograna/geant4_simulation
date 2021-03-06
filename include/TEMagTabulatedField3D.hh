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
//    *    PurgMagTabulatedField3D.hh     *
//    *                                   *
//    *************************************
//
//

#include "globals.hh"
#include "G4ElectroMagneticField.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>

class TEMagTabulatedField3D : public G4ElectroMagneticField
{
  
  // Storage space for the table
  std::vector< std::vector< std::vector< double > > > xBField;
  std::vector< std::vector< std::vector< double > > > yBField;
  std::vector< std::vector< std::vector< double > > > zBField;
  std::vector< std::vector< std::vector< double > > > xEField;
  std::vector< std::vector< std::vector< double > > > yEField;
  std::vector< std::vector< std::vector< double > > > zEField;
  // The dimensions of the table
  int nx,ny,nz; 
  // The physical limits of the defined region
  int minx =  10000; int miny =  10000; int minz =  10000;
  int maxx = -10000; int maxy = -10000; int maxz = -10000;
  // The physical extent of the defined region
  double dx, dy, dz;
  double fZoffset;
  bool invertX, invertY, invertZ;

public:
  virtual G4bool DoesFieldChangeEnergy() const { return true; };
  
  TEMagTabulatedField3D(const char* filename, double zOffset );
  void  GetFieldValue( const  double Point[4], double *Bfield          ) const;
};

