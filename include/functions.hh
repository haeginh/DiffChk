/*
 * functions.hh
 *
 *  Created on: Dec 3, 2019
 *      Author: hhg
 */

#ifndef INCLUDE_FUNCTIONS_HH_
#define INCLUDE_FUNCTIONS_HH_

#include "VOXModelImport.hh"
#include "TETModelImport.hh"
#include "InternalSampling.hh"
#include "InternalSamplingVox.hh"
#include "G4Timer.hh"
#include <omp.h>
using namespace std;

vector<pair<string, pair<int, int>>> ReadOrganFile
  (string fileName, TETModelImport* tetModel, VOXModelImport* voxModel)
{
	ifstream ifs(fileName);
	if(!ifs.is_open()){
		cerr<<"There is no "<<fileName<<"!!"<<endl;
		exit(1);
	}

	string dump;
	vector<pair<string, pair<int, int>>> selected;
	while(getline(ifs, dump)){
		stringstream ss(dump);
		int tet, vox;
		ss>>dump>>tet>>vox;
		selected.push_back(make_pair(dump, make_pair(tet, vox)));
	}
	ifs.close();

	// Print the overall information for each organ
	//
	double voxelVol = voxModel->GetVoxelSize().getX()
			         *voxModel->GetVoxelSize().getY()
					 *voxModel->GetVoxelSize().getZ();
	cout << endl
		 << setw(15) << "Organ Name"
		 << setw(11) << "tet ID"
		 << setw(11) << "tet #"
	     << setw(11) << "tet vol"
		 << setw(11) << "vox ID"
		 << setw(11) << "vox #"
	     << setw(11) << "vox vol"<<endl;
	cout << "---------------------------------------------------------------------------------"<<endl;
	for(auto iter:selected){
		cout<<setw(15)<<iter.first
			<<setw(11)<<iter.second.first
			<<setw(11)<<tetModel->GetNumTet(iter.second.first)
			<<setw(11)<<tetModel->GetVolume(iter.second.first)/cm3
			<<setw(11)<<iter.second.second
			<<setw(11)<<voxModel->GetNumVoxel(iter.second.second)
			<<setw(11)<<voxModel->GetNumVoxel(iter.second.second)*voxelVol/cm3<<endl;
	}
	cout<<endl;
	return selected;
}

vector<double> CalculateCD(vector<pair<string, pair<int, int>>> selected,
		          TETModelImport* tetPhan, VOXModelImport* voxPhan, int samplingNum)
{
	vector<double> cdVec;
	for(auto iter:selected){
		InternalSource tetInternal(tetPhan);
		tetInternal.SetSource(iter.second.first);
		InternalSourceVox voxInternal(voxPhan);
		voxInternal.SetSource(iter.second.second);
		G4ThreeVector tetPoint,voxPoint;
		double tx(0.), ty(0.), tz(0.), vx(0.), vy(0.), vz(0.);
#pragma omp parallel for private(tetPoint, voxPoint) reduction(+:tx, ty, tz, vx, vy, vz)
		for(int i=0;i<samplingNum;i++){
			tetInternal.GetAprimaryPos(tetPoint);
			voxInternal.GetAprimaryPos(voxPoint);
			tx += tetPoint.getX(); ty += tetPoint.getY(); tz += tetPoint.getZ();
			vx += voxPoint.getX(); vy += voxPoint.getY(); vz += voxPoint.getZ();
		}
		double cd=(G4ThreeVector(tx, ty, tz)-G4ThreeVector(vx, vy, vz)).mag()/(double)samplingNum;
		cdVec.push_back(cd);
	}
	return cdVec;
}

vector<double> CalculateDI(vector<pair<string, pair<int, int>>> selected,
		          TETModelImport* tetPhan, VOXModelImport* voxPhan, int samplingNum)
{
	vector<double> diVec;
	for(auto iter:selected){
		InternalSource tetInternal(tetPhan);
		tetInternal.SetSource(iter.second.first);
		InternalSourceVox voxInternal(voxPhan);
		voxInternal.SetSource(iter.second.second);
		G4ThreeVector tetPoint;
		int count(0);
#pragma omp parallel for private(tetPoint) reduction(+:count)
		for(int i=0;i<samplingNum;i++){
			tetInternal.GetAprimaryPos(tetPoint);
			if(voxInternal.IsInside(tetPoint)) count++;
		}
		double tetVol = tetPhan->GetVolume(iter.second.first);
		double voxVol = (double)voxPhan->GetNumVoxel(iter.second.second)*voxPhan->GetVoxelVol();
		double di = ((double)count/(double)samplingNum)*tetVol*2./(tetVol+voxVol);
		diVec.push_back(di);
	}
	return diVec;
}
#endif /* INCLUDE_FUNCTIONS_HH_ */
