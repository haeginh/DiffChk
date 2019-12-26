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
#include <iomanip>
using namespace std;

typedef vector<pair<string, pair<vector<int>, vector<int>>>> SELECTED;
SELECTED ReadOrganFile
  (string fileName, TETModelImport* tetModel, VOXModelImport* voxModel)
{
	ifstream ifs(fileName);
	if(!ifs.is_open()){
		cerr<<"There is no "<<fileName<<"!!"<<endl;
		exit(1);
	}

	string dump;
	SELECTED selected;
	while(getline(ifs, dump)){
		stringstream ss(dump);
		string sTet; vector<int> tet, vox;
		ss>>dump;
		while(ss>>sTet){
			if(sTet=="/") break;
			tet.push_back(atoi(sTet.c_str()));
		}
		int intTemp;
		while(ss>>intTemp) vox.push_back(intTemp);
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
		 << setw(15) << "tet ID"
		 << setw(11) << "tet #"
	     << setw(11) << "tet vol"
		 << setw(15) << "vox ID"
		 << setw(11) << "vox #"
	     << setw(11) << "vox vol"<<endl;
	cout << "-----------------------------------------------------------------------------------------"<<endl;
	for(auto it = selected.begin();it!=selected.end();){
		int tetNum(0); double tetVol(0.); string tetString;
		for(auto idx:it->second.first){
			tetNum+=tetModel->GetNumTet(idx);
			tetVol+=tetModel->GetVolume(idx);
			tetString += to_string(idx) + "/";
		}
		int voxNum(0); string voxString;
		for(auto idx:it->second.second){
			voxNum+=voxModel->GetNumVoxel(idx);
			voxString += to_string(idx) + "/";
		}

		cout<<setw(15)<<it->first
			<<setw(15)<<tetString
			<<setw(11)<<tetNum
			<<setw(11)<<tetVol/cm3
			<<setw(15)<<voxString
			<<setw(11)<<voxNum
			<<setw(11)<<voxNum*voxelVol/cm3;

		if(tetNum==0 || voxNum==0){
			it = selected.erase(it);
			cout<<" > EXCLUDED";
		}
		else ++it;
		cout<<endl;
	}
	cout<<endl;
	return selected;
}

pair<vector<double>, vector<G4ThreeVector>> CalculateCD(SELECTED selected,
		          TETModelImport* tetPhan, VOXModelImport* voxPhan, int samplingNum)
{
	pair<vector<double>, vector<G4ThreeVector>> cdVec;
	for(auto iter:selected){
		cout<<"\t"<<iter.first<<"..."<<flush;
		G4Timer timer; timer.Start();
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
		G4ThreeVector trans = (G4ThreeVector(tx, ty, tz)-G4ThreeVector(vx, vy, vz))/(double)samplingNum;
		double cd=trans.mag();
		cdVec.first.push_back(cd);
		cdVec.second.push_back(trans);
		timer.Stop();
		cout<<timer.GetRealElapsed()<<" -> "<<cd<<" "<<trans<<endl;
	}
	return cdVec;
}

vector<double> CalculateDI(SELECTED selected,
		          TETModelImport* tetPhan, VOXModelImport* voxPhan, int samplingNum)
{
	vector<double> diVec;
	for(auto iter:selected){
		cout<<"\t"<<iter.first<<"..."<<flush;
		G4Timer timer; timer.Start();
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
		double tetVol(0.); double voxVol(0.);
		for(auto id:iter.second.first) tetVol += tetPhan->GetVolume(id);
		for(auto id:iter.second.second) voxVol += (double)voxPhan->GetNumVoxel(id)*voxPhan->GetVoxelVol();
		double di = ((double)count/(double)samplingNum)*tetVol*2./(tetVol+voxVol);
		diVec.push_back(di);
		timer.Stop();
		cout<<timer.GetRealElapsed()<<" -> "<<di<<endl;
	}
	return diVec;
}
#endif /* INCLUDE_FUNCTIONS_HH_ */
