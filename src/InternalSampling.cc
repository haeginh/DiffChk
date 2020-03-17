/*
 * ExternalBeam.cc
 *
 *  Created on: Oct 18, 2019
 *      Author: hhg
 */

#include "TETModelImport.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include <set>
#include <algorithm>
#include <fstream>
#include "InternalSampling.hh"

InternalSource::InternalSource(TETModelImport* _tetData)
:tetData(_tetData)
{}

InternalSource::~InternalSource()
{}

void InternalSource::SetSource(std::vector<int> sources)
{
	//Cout
    tetPick.clear(); isRBM = false;
    if(sources[0]==0){isRBM = true;}

	//Extract source tet IDs
	std::set<int> sourceSet(sources.begin(), sources.end());
    if(isRBM) sourceSet.erase(0);

    if(!isRBM){
        for(G4int i=0;i<tetData->GetNumTotTet();i++){
            if(sourceSet.find(tetData->GetMaterialIndex(i)) != sourceSet.end())
                tetPick.push_back(VOLPICK(tetData->GetTetrahedron(i)->GetCubicVolume(), i));
        }
    }else{
        auto rbmRatio = tetData->GetRBMRatio();
	for(auto ss:sourceSet){
            if(rbmRatio.find(ss)==rbmRatio.end()){
                G4cerr<<ss<<" is not included in RBMnBS file!"<<G4endl; exit(0);
            }
        }
        for(G4int i=0;i<tetData->GetNumTotTet();i++){
            if(sourceSet.find(tetData->GetMaterialIndex(i)) != sourceSet.end())
                tetPick.push_back(VOLPICK(tetData->GetTetrahedron(i)->GetCubicVolume()*rbmRatio[tetData->GetMaterialIndex(i)], i));
        }
    }

	//Arrange volumes
	std::sort(tetPick.begin(), tetPick.end());
	std::reverse(tetPick.begin(), tetPick.end());

	G4double previousVol(0.);
	for(auto &tp:tetPick) {
		tp.first += previousVol;
		previousVol = tp.first;
	}

	for(auto &tp:tetPick) tp.first /= previousVol;

	//cout<<source<<" -> "<<tetPick.size()<<endl;
}

void InternalSource::GetAprimaryPos(G4ThreeVector &position)
{
	G4double rand = G4UniformRand();
	for(auto tp:tetPick){
		if(rand>tp.first) continue;
		position = RandomSamplingInTet(tetData->GetTetrahedron(tp.second)); break;
	}
}

G4ThreeVector InternalSource::RandomSamplingInTet(G4Tet* tet){

	G4double varS = G4UniformRand();
	G4double varT = G4UniformRand();
	G4double varU = G4UniformRand();

	if (varS+varT>1.0){

		varS = 1.0 - varS;
		varT = 1.0 - varT;

	}
	if (varT+varU>1.0){

		double tmp = varU;
		varU = 1.0 - varS - varT;
		varT = 1.0 -tmp;
	} else if (varS+varT+varU>1.0){

		double tmp = varU;
		varU = varS + varT + varU - 1.0;
		varS = 1 - varT - tmp;
	}

	double a = 1 - varS - varT - varU;

	G4ThreeVector SampledPosition = a*(tet->GetVertices()[0])+varS*(tet->GetVertices()[1])+varT*(tet->GetVertices()[2])+varU*(tet->GetVertices()[3]);
	return SampledPosition;
}


