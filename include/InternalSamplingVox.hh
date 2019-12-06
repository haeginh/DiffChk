/*
 * ExternalBeam.hh
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#ifndef SRC_INTERNALSAMPLINGVOX_HH_
#define SRC_INTERNALSAMPLINGVOX_HH_

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "VOXModelImport.hh"

#include <vector>
#include <set>

class InternalSourceVox
{
public:
	InternalSourceVox(VOXModelImport* voxData);
	~InternalSourceVox();

	void SetSource(int source);
	void GetAprimaryPos(G4ThreeVector &pos);

	int GetSource() 	        const{return sourceID;}

	bool IsInside(G4ThreeVector point);

private:
	G4ThreeVector RandomSamplingInAVox(INT3 ijk);

private:
    int                   sourceID;
    VOXModelImport*       voxData;
    std::vector<INT3>     sourceVoxels;
    std::set<INT3>        sourceVoxelsSet;
    G4ThreeVector         voxelSize;
    G4ThreeVector         voxelSizeInv;
    G4ThreeVector         trans;
};



#endif /* SRC_EXTERNALBEAM_HH_ */
