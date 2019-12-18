/*
 * ExternalBeam.hh
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#ifndef SRC_BEAMGENERATOR_HH_
#define SRC_BEAMGENERATOR_HH_

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Tet.hh"

#include <vector>

class    TETModelImport;
typedef  std::pair<G4double, G4int> VOLPICK;
class InternalSource
{
public:
	InternalSource(TETModelImport* tetData);
	~InternalSource();

	void SetSource(vector<int> source);
	void GetAprimaryPos(G4ThreeVector &pos);

private:
	G4ThreeVector RandomSamplingInTet(G4Tet* tet);

private:
    TETModelImport*       tetData;
    std::vector<VOLPICK>  tetPick;
};



#endif /* SRC_EXTERNALBEAM_HH_ */
