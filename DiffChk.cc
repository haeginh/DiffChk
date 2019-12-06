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
	string organFile, tetName, voxName;
	int samplingNum(-1);
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-f")
			organFile = string(argv[(i++)+1]);
		else if(string(argv[i])=="-n")
			samplingNum = atoi(argv[(i++)+1]);
		else if(string(argv[i])=="-t")
			omp_set_num_threads(atoi(argv[(i++)+1]));
		else if(i==argc-1)
			tetName = string(argv[i]);
		else{
			PrintUsage();
			return 1;
		}
	}

	TETModelImport* tetPhan = new TETModelImport(tetName);
	auto selected = ReadOrganFile(organFile);
	cout<<"CLD will be calculated for "<<selected.size()<<" pairs"<<endl;

	if(samplingNum<0) samplingNum = 10000000;
	CalculateCLD(organFile, selected, tetPhan, samplingNum);

}


