/*
 * ImportVoxelPhantom.hh
 *
 *  Created on: Dec 04, 2019
 *      Author: haeginh
 */

#ifndef INCLUDE_IMPORTVOXELPHANTOM_HH_
#define INCLUDE_IMPORTVOXELPHANTOM_HH_

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <string>
#include <tuple>
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

using namespace std;
typedef tuple<int, int, int> INT3;

class VOXModelImport {
public:
	VOXModelImport(string _phantomFile);
	virtual ~VOXModelImport();

	void SetPhantomInfo(string infoFile);
	void ImportPhantomVoxelData(string voxFile);
	void PrintInfomation();

	INT3          GetVoxelDim()          {return voxelDim;}
	G4ThreeVector GetVoxelSize()         {return voxelSize;}
	double        GetVoxelVol()          {return voxelVol;}
	int           GetNumVoxel(G4int idx) {return organVoxels[idx];}
	int           GetVoxelIdx(G4int idx) {return voxelData[idx];}
	G4ThreeVector GetPhantomSize()       {return phantomSize;}

private:
	//voxel info.
	string         phantomName;
	G4ThreeVector  voxelSize;
	double         voxelVol;
	INT3           voxelDim;
	int *          voxelData;
	map<int, int>  organVoxels;
	G4ThreeVector  phantomSize;

	//air separator
	bool chkET, chkStomach;
	int etSepZ, stmSepZ;

	string Filename;

};

#endif /* INCLUDE_IMPORTVOXELPHANTOM_HH_ */
