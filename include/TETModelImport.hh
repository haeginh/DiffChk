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
// TETModelImport.hh
// \file   MRCP_GEANT4/External/include/TETModelImport.hh
// \author Haegin Han
//

#ifndef TETModelImport_h
#define TETModelImport_h 1

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4Tet.hh"

// *********************************************************************
// This class is to import the phantom data from *.ele, *.node, and
// *.material files.
// -- DataRead: Construct G4Tet by reading data from *.ele and *.node
//              files
// -- MaterialRead: Construct G4Material by reading material data from
//                  *.material file
// -- ColourRead: Construct std::map that contains G4Colour data
//                according to data in colour.dat file
// -- PrintMaterialInformation: Print a table that contains organ ID,
//                              number of tetrahedrons, volume, density,
//                              mass, and name for each organ
// *********************************************************************
using namespace std;

class TETModelImport
{
public:
	TETModelImport(G4String phantomName);
    virtual ~TETModelImport() {}

	// get methods
	int           GetNumTotTet()         {return tetVector.size();}
	int           GetNumTet(int idx)        { return numTetMap[idx];}
	double        GetVolume(int idx)       { return volumeMap[idx]; }
	int           GetMaterialIndex(int idx){ return organID[idx]; }
	G4Tet*        GetTetrahedron(int idx)  { return tetVector[idx]; }
	G4ThreeVector GetPhantomSize()           { return phantomSize; }
	G4ThreeVector GetPhantomBoxMin()         { return boundingBox_Min; }
	G4ThreeVector GetPhantomBoxMax()         { return boundingBox_Max; }
    std::map<G4int, G4double> GetRBMRatio()  {return rbmRatio;}

private:
	void DataRead(G4String, G4String);
    void ReadRBMnBS(G4String);
	void PrintInfomation();

	G4String phantomName;

	G4ThreeVector boundingBox_Min;
	G4ThreeVector boundingBox_Max;
	G4ThreeVector phantomSize;

	std::vector<G4ThreeVector> vertexVector;
	std::vector<G4Tet*>        tetVector;
	std::vector<int*>        eleVector;
	std::vector<int>         organID;
	std::map<int, int>     numTetMap;
	std::map<int, double>  volumeMap;

    std::map<G4int, G4double>  rbmRatio;
};

#endif
