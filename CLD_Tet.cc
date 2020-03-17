/*
 * DiffChk.cc
 *
 *  Created on: Dec 5, 2019
 *      Author: hhg
 */

#include "functions.hh"

void PrintUsage(){
	cout<<"USAGE: ./CLD_TET [options] (tet. phan. w/o suffix)"<<endl;
	cout<<"OPTIONS:"<<endl;
	cout<<setw(4)<<"-f "<<"file name that contains (see organs.txt)"<<endl;
	cout<<setw(4)<<"-n "<<"sampling number (default:1E7)"<<endl;
	cout<<setw(4)<<"-t "<<"number of threads (default:1)"<<endl;
	cout<<setw(4)<<"-o "<<"prefix of output file name (optional)"<<endl;
}
int main(int argc, char** argv){
	string organFile, tetName, outName;
	int samplingNum(-1);
	omp_set_num_threads(1);
	if(argc==1){PrintUsage(); return 1;}
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-f")
			organFile = string(argv[++i]);
		else if(string(argv[i])=="-n")
			samplingNum = atoi(argv[++i]);
		else if(string(argv[i])=="-t")
			omp_set_num_threads(atoi(argv[++i]));
		else if(string(argv[i])=="-o")
			outName = string(argv[++i]);
		else if(i==argc-1)
			tetName = string(argv[i]);
		else{
			PrintUsage();
			return 1;
		}
	}
	if(organFile.empty()) {PrintUsage(); return 1;}

	TETModelImport* tetPhan = new TETModelImport(tetName);
	auto selected = ReadOrganFile(organFile);
	cout<<"CLD will be calculated for "<<selected.size()<<" pairs"<<endl;

	if(samplingNum<0) samplingNum = 10000000;
	CalculateCLD(outName, selected, tetPhan, samplingNum);

}

