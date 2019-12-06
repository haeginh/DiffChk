/*
 * DiffChk.cc
 *
 *  Created on: Dec 5, 2019
 *      Author: hhg
 */

#include "functions.hh"

void PrintUsage(){

}
int main(int argc, char** argv){
	bool cdSwitch (false);
	bool diSwitch(false);
	string organFile, tetName, voxName;
	int samplingNum(-1);
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-cd")
			cdSwitch=true;
		else if(string(argv[i])=="-di")
			diSwitch=true;
		else if(string(argv[i])=="-f")
			organFile = string(argv[(i++)+1]);
		else if(string(argv[i])=="-n")
			samplingNum = atoi(argv[(i++)+1]);
		else if(string(argv[i])=="-t")
			omp_set_num_threads(atoi(argv[(i++)+1]));
		else if(i==argc-2)
			tetName = string(argv[i]);
		else if(i==argc-1)
			voxName = string(argv[i]);
		else{
			PrintUsage();
			return 1;
		}
	}

	TETModelImport* tetPhan = new TETModelImport(tetName);
	VOXModelImport* voxPhan = new VOXModelImport(voxName);
	auto selected = ReadOrganFile(organFile, tetPhan, voxPhan);

	vector<string> contents;
	vector<vector<double>> values;
	if(cdSwitch){
		cout<<"Calculate cd.."<<flush;
		G4Timer timer; timer.Start();
		contents.push_back("CD");
		if(samplingNum<0) samplingNum = 10000000;
		values.push_back(CalculateCD(selected, tetPhan, voxPhan, samplingNum));
		timer.Stop(); cout<<timer.GetRealElapsed()<<endl;
	}
	if(diSwitch){
		cout<<"Calculate di.."<<flush;
		G4Timer timer; timer.Start();
		contents.push_back("DI");
		if(samplingNum<0) samplingNum = 10000000;
		values.push_back(CalculateDI(selected, tetPhan, voxPhan, samplingNum));
		timer.Stop(); cout<<timer.GetRealElapsed()<<endl;
	}
}


