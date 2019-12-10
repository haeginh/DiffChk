/*
 * HDCalculator.hh
 *
 *  Created on: Dec 8, 2019
 *      Author: hhg
 */

#ifndef SRC_HDCALCULATOR_HH_
#define SRC_HDCALCULATOR_HH_

#include "VOXModelImport.hh"
#include "TETModelImport.hh"
#include "G4TessellatedSolid.hh"

using namespace std;
typedef std::tuple<int, int, int> DOT_INDEX;
typedef map<int, vector<pair<pair<int, int>, pair<int, int>>>> COMPRESSED;
enum  Axis {xAxis, yAxis, zAxis};

class HDCalculator {
public:
	HDCalculator(vector<pair<string, pair<int, int>>> selected, TETModelImport* tetPhan, VOXModelImport* voxPhan, int samplingNum);
	virtual ~HDCalculator();

//	void   ExtractBoundaries();
	vector<double> CalculateHDs();
private:
	double CalculateHD(int idx);
	vector<pair<double, int>> ArrangeFacetArea(G4TessellatedSolid* tess);
	//tet
	G4TessellatedSolid* ExtractTetBoundary(G4int idx);
	G4ThreeVector       GetAPointOnTetSurface(G4TessellatedSolid* tess);

	//vox
	G4TessellatedSolid* ExtractVoxBoundary(G4int idx);
	G4ThreeVector       GetAPointOnVoxSurface(G4TessellatedSolid* tess);
	COMPRESSED Compress2D(map<int, vector<pair<int, int>>> boundaryPix);
	void GetVertices(COMPRESSED compressed, vector<DOT_INDEX> &vertices, Axis axis);

	//variables
	vector<pair<string, pair<int, int>>> selected;
	int samplingNum;
	//tet
	TETModelImport* tetPhan;
	G4TessellatedSolid* tetTess;
	vector<pair<double, int>>   tetFacetArea;

	//vox
	VOXModelImport* voxPhan;
	G4ThreeVector   voxTrans;
	G4TessellatedSolid* voxTess;
	vector<pair<double, int>>   voxFacetArea;

};

#endif /* SRC_HDCALCULATOR_HH_ */
