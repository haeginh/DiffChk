//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// TETModelImport.cc
// \file   MRCP_GEANT4/External/src/TETModelImport.cc
// \author Haegin Han
//

#include "TETModelImport.hh"
#include "G4Timer.hh"
TETModelImport::TETModelImport(G4String _phantomName)
{
	// set phantom name
	phantomName = _phantomName;

	G4Timer timer; timer.Start();
	cout << "Importing tet. phantom (" << phantomName << ")..."<< flush;

	G4String eleFile      =  phantomName + "ele";
	G4String nodeFile     =  phantomName + "node";

	// read phantom data files (*. ele, *.node)
	DataRead(eleFile, nodeFile);
	// print the summary of phantom information
	//PrintInfomation();
	timer.Stop();
	cout << timer.GetRealElapsed()<<endl;
}

void TETModelImport::DataRead(G4String eleFile, G4String nodeFile)
{
	G4String tempStr;
	int tempInt;

	// Read *.node file
	//
	std::ifstream ifpNode;

	ifpNode.open(nodeFile);
	if(!ifpNode.is_open()) {
		// exception for the case when there is no *.node file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + nodeFile ).c_str());
	}

	int numVertex;
	double xPos, yPos, zPos;
	double xMin(DBL_MAX), yMin(DBL_MAX), zMin(DBL_MAX);
	double xMax(DBL_MIN), yMax(DBL_MIN), zMax(DBL_MIN);

	ifpNode >> numVertex >> tempInt >> tempInt >> tempInt;

	for(int i=0; i<numVertex; i++)
	{
		ifpNode >> tempInt >> xPos >> yPos >> zPos;

		// set the unit
		xPos*=cm;
		yPos*=cm;
		zPos*=cm;

		// save the node data as the form of std::vector<G4ThreeVector>
		vertexVector.push_back(G4ThreeVector(xPos, yPos, zPos));

		// to get the information of the bounding box of phantom
		if (xPos < xMin) xMin = xPos;
		if (xPos > xMax) xMax = xPos;
		if (yPos < yMin) yMin = yPos;
		if (yPos > yMax) yMax = yPos;
		if (zPos < zMin) zMin = zPos;
		if (zPos > zMax) zMax = zPos;
	}

	// set the variables for the bounding box and phantom size
	boundingBox_Min = G4ThreeVector(xMin,yMin,zMin);
	boundingBox_Max = G4ThreeVector(xMax,yMax,zMax);
	G4ThreeVector center = (boundingBox_Max+boundingBox_Min)*0.5;
	phantomSize = G4ThreeVector(xMax-xMin,yMax-yMin,zMax-zMin);

	ifpNode.close();

	// Read *.ele file
	//
	std::ifstream ifpEle;

	ifpEle.open(eleFile);
	if(!ifpEle.is_open()) {
		// exception for the case when there is no *.ele file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + eleFile ).c_str());
	}

	int numEle;
	ifpEle >> numEle  >> tempInt >> tempInt;

	for(int i=0; i<numEle; i++)
	{
		ifpEle >> tempInt;
		int* ele = new int[4];
		for(int j=0;j<4;j++){
			ifpEle >> tempInt;
			ele[j]=tempInt;
		}
		eleVector.push_back(ele);
		ifpEle >> tempInt;
		organID.push_back(tempInt);

		// save the element (tetrahedron) data as the form of std::vector<G4Tet*>
		tetVector.push_back(new G4Tet("Tet_Solid",
							   		  vertexVector[ele[0]]-center,
									  vertexVector[ele[1]]-center,
									  vertexVector[ele[2]]-center,
									  vertexVector[ele[3]]-center));

		// calculate the total volume and the number of tetrahedrons for each organ
		std::map<int, double>::iterator FindIter = volumeMap.find(organID[i]);

		if(FindIter!=volumeMap.end()){
			FindIter->second += tetVector[i]->GetCubicVolume();
			numTetMap[organID[i]]++;
		}
		else {
			volumeMap[organID[i]] = tetVector[i]->GetCubicVolume();
			numTetMap[organID[i]] = 1;
		}
	}
	ifpEle.close();
}

void TETModelImport::ReadRBMnBS(G4String bonefile){
    std::ifstream ifs(bonefile);
    if(!ifs.is_open()) {
        // exception for the case when there is no *.material file
        G4Exception("TETModelImport::RBMBSRead","",JustWarning,
                G4String("      There is no " + bonefile ).c_str());
        return;
    }
    G4int idx;
    G4double rbm, bs;
    while(ifs>>idx>>rbm>>bs){
        rbmRatio[idx]=rbm;
    }
}

void TETModelImport::PrintInfomation()
{
	// Print the overall information for each organ
	//
	cout << endl
		 << setw(9)  << "Organ ID"
		 << setw(11) << "# of Tet"
		 << setw(11) << "vol [cm3]"<<endl;
	cout << "----------------------------------------------------------------"<<endl;

	cout<<std::setiosflags(std::ios::fixed);
	cout.precision(3);
	for(auto vm:volumeMap)
	{
		cout << std::setw(9)  << vm.first                    // organ ID
			   << std::setw(11) << numTetMap[vm.first]         // # of tetrahedrons
			   << std::setw(11) << vm.second/cm3 << endl;    // organ volume
	}
}

