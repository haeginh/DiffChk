/*
 * DiffChk.cc
 *
 *  Created on: Dec 5, 2019
 *      Author: hhg
 */

#include "functions.hh"
#include "HDCalculator.hh"

void PrintUsage(){
	cout<<"USAGE: ./DiffChk [options] (tet. phan. w/o suffix) (vox. phan.)"<<endl;
	cout<<"OPTIONS:"<<endl;
	cout<<setw(5)<<"-cd "<<"option for calculating CD"<<endl;
	cout<<setw(5)<<"-di "<<"option for calculating DI"<<endl;
	cout<<setw(5)<<"-hd "<<"option for calculating HD"<<endl;
	cout<<setw(5)<<"-f "<<"file name that contains organ info. (see organs.txt)"<<endl;
	cout<<setw(5)<<"-n "<<"sampling number (default: 1E7)"<<endl;
	cout<<setw(5)<<"-t "<<"the number of threads (default: 1)"<<endl;
	cout<<setw(5)<<"-o "<<"output file name (default: out.txt)"<<endl;
}
int main(int argc, char** argv){
	bool cdSwitch (false);
	bool diSwitch(false);
	bool hdSwitch(false);
	string organFile, outName("out.txt"), tetName, voxName;
	int samplingNum(1e7);
	omp_set_num_threads(1);
	if(argc==1) {PrintUsage(); return 1;}
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-cd")
			cdSwitch=true;
		else if(string(argv[i])=="-di")
			diSwitch=true;
		else if(string(argv[i])=="-hd")
			hdSwitch=true;
		else if(string(argv[i])=="-f")
			organFile = string(argv[(i++)+1]);
		else if(string(argv[i])=="-n")
			samplingNum = atoi(argv[(i++)+1]);
		else if(string(argv[i])=="-t")
			omp_set_num_threads(atoi(argv[(i++)+1]));
		else if(string(argv[i])=="-o")
			outName = string(argv[(i++)+1]);
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
	vector<G4ThreeVector> cdCoordVec;
	if(cdSwitch){
		cout<<"Calculate CD.."<<endl;
		contents.push_back("CD\tx\ty\tz");
		auto val = CalculateCD(selected, tetPhan, voxPhan, samplingNum);
		values.push_back(val.first);
		cdCoordVec = val.second;
	}
	if(diSwitch){
		cout<<"Calculate DI.."<<endl;
		contents.push_back("DI");
		values.push_back(CalculateDI(selected, tetPhan, voxPhan, samplingNum));
	}
	if(hdSwitch){
		cout<<"Calculate HD.."<<endl;
		contents.push_back("HD");
		auto hdCal = new HDCalculator(selected, tetPhan, voxPhan, samplingNum);
		values.push_back(hdCal->CalculateHDs());
	}

	ofstream ofs(outName);
	for(auto c:contents) ofs<<"\t"<<c;
	ofs<<endl;

	for(size_t i=0;i<selected.size();i++){
		ofs<<selected[i].first<<flush;
		for(size_t j=0;j<values.size();j++){
			ofs<<"\t"<<values[j][i]<<flush;
			if(contents[j]=="CD\tx\ty\tz"){
				G4ThreeVector p = cdCoordVec[i];
				ofs<<"\t"<<p.getX()<<"\t"<<p.getY()<<"\t"<<p.getZ();
			}
		}
		ofs<<endl;
	}
	ofs.close();
}


