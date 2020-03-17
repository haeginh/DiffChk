/*
 * functions.hh
 *
 *  Created on: Dec 3, 2019
 *      Author: hhg
 */

#ifndef INCLUDE_FUNCTIONS_HH_
#define INCLUDE_FUNCTIONS_HH_

#include "TETModelImport.hh"
#include "InternalSampling.hh"
#include "G4Timer.hh"
#include <omp.h>
#include <iomanip>
using namespace std;

vector<pair<string, pair<vector<int>, vector<int>>>> ReadOrganFile
  (string fileName)
{
	ifstream ifs(fileName);
	if(!ifs.is_open()){
		cerr<<"There is no "<<fileName<<"!!"<<endl;
		exit(1);
	}

	string dump;
	vector<pair<string, pair<vector<int>, vector<int>>>> selected;
	while(getline(ifs, dump)){
        if(dump.empty()) continue;
		stringstream ss(dump);
		string a, b;
		//ss>>quoted(dump)>>quoted(a)>>quoted(b);
		ss>>dump>>quoted(a)>>quoted(b);
		stringstream ssa(a), ssb(b);
		vector<int> aVec, bVec;
		int intTemp;
		while(ssa>>intTemp) aVec.push_back(intTemp);
		while(ssb>>intTemp) bVec.push_back(intTemp);
		selected.push_back(make_pair(dump, make_pair(aVec, bVec)));
	}
	ifs.close();

	return selected;
}

void CalculateCLD(string fileName, vector<pair<string, pair<vector<int>, vector<int>>>> selected,
		          TETModelImport* tetPhan, int samplingNum)
{
    G4Timer timer; timer.Start();
    int count(1);
    for(auto iter:selected){
        cout<<'\r'<<count++<<"/"<<selected.size()<<" : "<<iter.first<<"...source setting        "<<flush;
		InternalSource internalA(tetPhan);
		internalA.SetSource(iter.second.first);
		InternalSource internalB(tetPhan);
		internalB.SetSource(iter.second.second);
		G4ThreeVector a, b;
		vector<map<int, int>> distBin;
		map<int, int> initialMap;
        cout<<'\r'<<count++<<"/"<<selected.size()<<" : "<<iter.first<<"...calculating CLD         "<<flush;

#pragma omp parallel
		{
#pragma omp single
			{
				int numThread = omp_get_num_threads();
				for(int i=0;i<numThread;i++) distBin.push_back(initialMap);
			}
			int tNum = omp_get_thread_num();
#pragma omp for private(a, b)
			for(int i=0;i<samplingNum;i++){
				internalA.GetAprimaryPos(a);
				internalB.GetAprimaryPos(b);
				int bin = floor((a-b).mag());
				distBin[tNum][bin] = distBin[tNum][bin] +1;
			}
		}

		//merge
		map<int, int> finalBin;
		for(auto db:distBin){
			for(auto aBin:db)
				finalBin[aBin.first] = finalBin[aBin.first] + aBin.second;
		}

		string file = fileName + "_" + iter.first + ".cld";
		ofstream ofs(file);
		for(auto aa:iter.second.first)  ofs<<aa<<" "; ofs<<"/ ";
		for(auto bb:iter.second.second) ofs<<bb<<" "; ofs<<endl;
		int max = finalBin.crbegin()->first;
		for(int i=0;i<max+2;i++)
			ofs<<i<<"\t"<<finalBin[i]/(double)samplingNum<<"\n";
		ofs.close();
	}
    timer.Stop();
    cout<<endl<<"--> "<<timer.GetRealElapsed()<<" s"<<endl;
	return;
}

#endif /* INCLUDE_FUNCTIONS_HH_ */
