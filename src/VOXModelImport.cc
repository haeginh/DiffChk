/*
 * ImportVoxelPhantom.cc
 *
 *  Created on: May 20, 2016
 *      Author: mchan
 */

#include "VOXModelImport.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Timer.hh"

VOXModelImport::VOXModelImport(string _phantomName)
:phantomName(_phantomName)
{
	G4Timer timer; timer.Start();
	cout << "Importing tet. phantom (" << phantomName << ")..."<< flush;

	string InfoFile = "Phantom_information";
	string voxelFile = _phantomName;

	SetPhantomInfo(InfoFile);
	ImportPhantomVoxelData(voxelFile);
//	PrintInfomation();
	timer.Stop();
	cout << timer.GetRealElapsed()<<endl;
}

VOXModelImport::~VOXModelImport() {
}


void VOXModelImport::SetPhantomInfo(string filename){
	std::ifstream ifp;
	ifp.open(filename);

	if(!ifp.is_open()) {
		cerr << filename << " not found!!" << endl;
		exit(1);
	}

	string dump;
	while(getline(ifp, dump))
	{
		if(dump.empty()) continue;
		stringstream ss(dump);
		ss >> dump;
		if(phantomName.find(dump)==string::npos) continue;
		int fNx, fNy, fNz;
		ss >> fNx >> fNy >> fNz;
		voxelDim = INT3(fNx, fNy, fNz);

		double xDim, yDim, zDim;
		ss>>xDim>>yDim>>zDim;
		voxelSize = G4ThreeVector(xDim, yDim, zDim) * mm;
		voxelVol = xDim*yDim*zDim*mm3;

		phantomSize.setX((double)fNx*voxelSize.x());
		phantomSize.setY((double)fNy*voxelSize.y());
		phantomSize.setZ((double)fNz*voxelSize.z());
		break;
	}

	ifp.close();
}

void VOXModelImport::ImportPhantomVoxelData(string voxelFile){
	std::ifstream is(voxelFile, std::ifstream::binary);
	if (!is) {
		cerr << "***error: file is not open." << endl;
		exit(1);
	}

	// get length of file:
	is.seekg(0, is.end);
	int64_t length = is.tellg();
	is.seekg(0, is.beg);

	voxelData = new char[length];
	int fNx, fNy, fNz;
	tie(fNx, fNy, fNz) = voxelDim;

	if (length != fNx * fNy * fNz) {
		cerr << "***error: wrong dimension." << endl;
		exit(1);
	}

//	std::cout << "Reading " << length << " characters... ";

	// read data as a block:
	is.read(voxelData, length);

//	if (is)
//		std::cout << "all characters read successfully." << endl;
	if(!is)
		std::cout << "error: only " << is.gcount() << " could be read";
	is.close();

	for(int i=0;i<length;i++){
		int idx = voxelData[i];
		if(idx<0) voxelData[i] += 256;
		organVoxels[voxelData[i]] = organVoxels[voxelData[i]] +1;
	}
}

void VOXModelImport::PrintInfomation(){
	// Print the overall information for each organ
	//
	cout<<setw(15)<<"Voxel size"<<voxelSize<<" mm3"<<endl;
	cout<<setw(15)<<"Voxel dim."<<get<0>(voxelDim)<<" * "
			                    <<get<1>(voxelDim)<<" * "
								<<get<2>(voxelDim)<<endl;
	cout << endl
		 << setw(9)  << "Organ ID"
		 << setw(11) << "# of Vox"
		 << setw(11) << "vol [cm3]"<<endl;
	cout << "----------------------------------------------------------------"<<endl;

	cout<<std::setiosflags(std::ios::fixed);
	cout.precision(3);
	double voxVol = voxelSize.getX()*voxelSize.getY()*voxelSize.getZ();
	for(auto ov:organVoxels)
	{
		cout << std::setw(9)  << ov.first                    // organ ID
			   << std::setw(11) << ov.second         // # of tetrahedrons
			   << std::setw(11) << ov.second*voxVol/cm3 << endl;    // organ volume
	}

}



