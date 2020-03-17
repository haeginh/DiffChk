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
#include <set>
#include <algorithm>
#include <fstream>
#include "InternalSamplingVox.hh"

InternalSourceVox::InternalSourceVox(VOXModelImport* _voxData)
:voxData(_voxData)
{
	voxelSize = voxData->GetVoxelSize();
	voxelSizeInv = G4ThreeVector(1./voxelSize.getX(), 1./voxelSize.getY(), 1./voxelSize.getZ());
	trans = -voxData->GetPhantomSize()*0.5;
}

InternalSourceVox::~InternalSourceVox()
{}

void InternalSourceVox::SetSource(vector<int> source)
{
	sourceVoxels.clear();
	INT3 voxDim = voxData->GetVoxelDim();
	set<int> sourceSet(source.begin(), source.end());
	for(int k=0;k<get<2>(voxDim);k++){
		for(int j=0;j<get<1>(voxDim);j++){
			for(int i=0;i<get<0>(voxDim);i++){
				int idx = i + j*get<0>(voxDim) + k*get<0>(voxDim)*get<1>(voxDim);
				if(sourceSet.find(voxData->GetVoxelIdx(idx))!=sourceSet.end())
					sourceVoxels.push_back(INT3(i, j, k));
			}
		}
	}
	std::set<INT3> toSet(sourceVoxels.begin(), sourceVoxels.end());
	sourceVoxelsSet.clear();
	sourceVoxelsSet = toSet;
}

void InternalSourceVox::GetAprimaryPos(G4ThreeVector &position)
{
	G4double rand = G4UniformRand();
	rand *= (double)sourceVoxels.size();
	INT3 selected = sourceVoxels[floor(rand)];
	position = RandomSamplingInAVox(selected);
}

G4ThreeVector InternalSourceVox::RandomSamplingInAVox(INT3 selected){
	double x, y, z;
	tie(x, y, z) = selected;
	G4ThreeVector sampledPos = G4ThreeVector((x+G4UniformRand())*voxelSize.getX(),
			                                 (y+G4UniformRand())*voxelSize.getY(),
											 (z+G4UniformRand())*voxelSize.getZ()) + trans;
	return sampledPos;
}

bool InternalSourceVox::IsInside(G4ThreeVector point){
	point -= trans;
	int fNx = floor(point.getX()*voxelSizeInv.getX());
	int fNy = floor(point.getY()*voxelSizeInv.getY());
	int fNz = floor(point.getZ()*voxelSizeInv.getZ());
	if(sourceVoxelsSet.find(INT3(fNx, fNy, fNz))==sourceVoxelsSet.end()) return false;
	return true;
}


