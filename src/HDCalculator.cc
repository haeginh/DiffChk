/*
 * HDCalculator.cc
 *
 *  Created on: Dec 8, 2019
 *      Author: hhg
 */

#include "HDCalculator.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include <iterator>

HDCalculator::HDCalculator(vector<pair<string, pair<int, int>>> _selected, TETModelImport* _tetPhan, VOXModelImport* _voxPhan, int _samplingNum)
: selected(_selected), samplingNum(_samplingNum), tetPhan(_tetPhan), voxPhan(_voxPhan)
{
	voxTrans = -voxPhan->GetPhantomSize()*0.5;
}

HDCalculator::~HDCalculator() {}

vector<double> HDCalculator::CalculateHDs()
{
	vector<double> values;
	for(size_t i=0; i<selected.size();i++){
		cout<<"\tCalculating for "<<selected[i].first<<"...I:"<<flush;
		values.push_back(CalculateHD(i));
	}
	return values;
}

double HDCalculator::CalculateHD(int idx)
{
	G4Timer timer; timer.Start();
	//Boundary extraction
	tetTess = ExtractTetBoundary(selected[idx].second.first);
	voxTess = ExtractVoxBoundary(selected[idx].second.second);
	//Arranging facet area (very fast)
	tetFacetArea.clear(); tetFacetArea = ArrangeFacetArea(tetTess);
	voxFacetArea.clear(); voxFacetArea = ArrangeFacetArea(voxTess);
	timer.Stop(); cout<<timer.GetRealElapsed()<<"/C:"<<flush;

	//Calculate HD (tet->vox)
	timer.Start();
	double hd(0.);
#pragma omp parallel for reduction(max:hd)
	for(int i=0;i<samplingNum;i++){
		G4ThreeVector tetP = GetAPointOnTetSurface(tetTess);
		double distTet = voxTess->DistanceToIn(tetP)+voxTess->DistanceToOut(tetP);
		G4ThreeVector voxP = GetAPointOnVoxSurface(voxTess);
		double distVox = tetTess->DistanceToIn(voxP)+tetTess->DistanceToOut(voxP);
		double dist = max(distTet, distVox);
		if(hd<dist) hd = dist;
	}
	timer.Stop(); cout<<timer.GetRealElapsed()<<" -> "<<hd<<endl;
	delete tetTess;
	delete voxTess;
	return hd;
}

vector<pair<double, int>> HDCalculator::ArrangeFacetArea(G4TessellatedSolid* aTess)
{
	vector<pair<double, int>> facetArea;
	for(int i=0;i<aTess->GetNumberOfFacets();i++)
		facetArea.push_back({aTess->GetFacet(i)->GetArea(), i});
	sort(facetArea.begin(), facetArea.end(), greater<pair<double, int>>());
	double previous(0);
	for(auto &aFacet:facetArea){
		aFacet.first += previous;
		previous = aFacet.first;
	}
	for(auto &aFacet:facetArea)
		aFacet.first /= previous;
	return facetArea;
}

G4TessellatedSolid* HDCalculator::ExtractTetBoundary(G4int idx)
{
	map<INT3, bool> facePool;
	for(int i=0;i<tetPhan->GetNumTotTet();i++){
		if(tetPhan->GetMaterialIndex(i)!=idx) continue;
		auto ele = tetPhan->GetAnElement(i);
		INT3 a = INT3(ele[1], ele[2], ele[3]);
		INT3 b = INT3(ele[0], ele[2], ele[3]);
		INT3 c = INT3(ele[0], ele[1], ele[3]);
		INT3 d = INT3(ele[0], ele[1], ele[2]);
		auto ret = facePool.insert({a, true}); //when it is not paired (will be selected as boundary)
		if(!ret.second) ret.first->second = false; // when it was paired
		ret = facePool.insert({b, true});
		if(!ret.second) ret.first->second = false;
		ret = facePool.insert({c, true});
		if(!ret.second) ret.first->second = false;
		ret = facePool.insert({d, true});
		if(!ret.second) ret.first->second = false;
	}
	G4TessellatedSolid* aTess = new G4TessellatedSolid();
	for(auto aFace:facePool){
		if(!aFace.second) continue;
		int a, b, c; tie(a, b, c) = aFace.first;
		aTess->AddFacet(new G4TriangularFacet(tetPhan->GetAVertiex(a),
				                              tetPhan->GetAVertiex(b),
											  tetPhan->GetAVertiex(c), ABSOLUTE));
	}
	aTess->SetSolidClosed(true);
	return aTess;
}

G4ThreeVector HDCalculator::GetAPointOnTetSurface(G4TessellatedSolid* tess){
	G4ThreeVector point;
	double rand1 = G4UniformRand();
	for(auto aFacet: tetFacetArea){
		if(rand1>aFacet.first) continue;
		auto   facet = tess->GetFacet(aFacet.second);
		G4ThreeVector vVec = facet->GetVertex(0) - facet->GetVertex(1);
		G4ThreeVector uVec = facet->GetVertex(2) - facet->GetVertex(1);
		double rand2(1.), rand3(1.);
		while(rand2 + rand3 > 1){
			rand2 = G4UniformRand();
			rand3 = G4UniformRand();
		}
		vVec = vVec * rand2;
		uVec = uVec * rand3;
		point = facet->GetVertex(1) + vVec + uVec;
		break;
	}
	return point;
}
G4TessellatedSolid* HDCalculator::ExtractVoxBoundary(G4int idx)
{
	int dimX, dimY, dimZ;
	tie(dimX, dimY, dimZ) = voxPhan->GetVoxelDim();
	map<int, vector<pair<int, int>>> xPos, xNeg, yPos, yNeg, zPos, zNeg;
	for (int k = 0; k < dimZ; k++) {
		for (int j = 0; j < dimY; j++) {
			for (int i = 0; i < dimX; i++) {
				int index = (dimX*dimY)*k + dimX * j + i;
				if (voxPhan->GetVoxelIdx(index) != idx) continue;
				//-x 비교
				if (i == 0)
					xNeg[i].push_back(make_pair(j, k));
				else if (voxPhan->GetVoxelIdx(index) != voxPhan->GetVoxelIdx(index - 1))
					xNeg[i].push_back(make_pair(j, k));
				//+x 비교
				if (i == dimX - 1)
					xPos[i+1].push_back(make_pair(j, k));
				else if (voxPhan->GetVoxelIdx(index) != voxPhan->GetVoxelIdx(index + 1))
					xPos[i+1].push_back(make_pair(j, k));

				//-y 확인
				if (j == 0)
					yNeg[j].push_back(make_pair(k, i));
				if (voxPhan->GetVoxelIdx(index) != voxPhan->GetVoxelIdx(index - dimX))
					yNeg[j].push_back(make_pair(k, i));

				//+y 비교
				if (j == dimY - 1)
					yPos[j+1].push_back(make_pair(k, i));
				else if (voxPhan->GetVoxelIdx(index) != voxPhan->GetVoxelIdx(index + dimX))
					yPos[j+1].push_back(make_pair(k, i));

				//-z 확인
				if (k == 0)
					zNeg[k].push_back(make_pair(i, j));
				else if (voxPhan->GetVoxelIdx(index) != voxPhan->GetVoxelIdx(index - dimX * dimY))
					zNeg[k].push_back(make_pair(i, j));

				//+z 비교
				if (k == dimZ - 1)
					zPos[k+1].push_back(make_pair(i, j));
				else if (voxPhan->GetVoxelIdx(index) != voxPhan->GetVoxelIdx(index + (dimX*dimY)))
					zPos[k+1].push_back(make_pair(i, j));
			}
		}
	}
	map<int, vector<pair<pair<int, int>, pair<int, int>>>> xPosCom, xNegCom, yPosCom, yNegCom, zPosCom, zNegCom;
	xPosCom = Compress2D(xPos); xNegCom = Compress2D(xNeg);
	yPosCom = Compress2D(yPos); yNegCom = Compress2D(yNeg);
	zPosCom = Compress2D(zPos); zNegCom = Compress2D(zNeg);
	int numFace(0);
	for (auto l : xPosCom) numFace += l.second.size();
	for (auto l : yPosCom) numFace += l.second.size();
	for (auto l : zPosCom) numFace += l.second.size();
	for (auto l : xNegCom) numFace += l.second.size();
	for (auto l : yNegCom) numFace += l.second.size();
	for (auto l : zNegCom) numFace += l.second.size();

	vector<DOT_INDEX> vertices;
	GetVertices(xPosCom, vertices, xAxis); GetVertices(xNegCom, vertices, xAxis);
	GetVertices(yPosCom, vertices, yAxis); GetVertices(yNegCom, vertices, yAxis);
	GetVertices(zPosCom, vertices, zAxis); GetVertices(zNegCom, vertices, zAxis);
	sort(vertices.begin(), vertices.end());
	vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());

	map<DOT_INDEX, int> vertexID;
	double vX = voxPhan->GetVoxelSize().getX();
	double vY = voxPhan->GetVoxelSize().getY();
	double vZ = voxPhan->GetVoxelSize().getZ();

	G4TessellatedSolid*   aTess = new G4TessellatedSolid();
	vector<G4ThreeVector> verCoord;
	for (size_t i = 0; i < vertices.size(); i++) {
		vertexID[vertices[i]] = i;
		G4ThreeVector thePoint = G4ThreeVector((double)get<0>(vertices[i])*vX, (double)get<1>(vertices[i])*vY, (double)get<2>(vertices[i])*vZ);
		thePoint += voxTrans;
		verCoord.push_back(thePoint);
	}
	for (auto l : xPosCom) {
		for (auto f : l.second) {
			int a = vertexID[DOT_INDEX(l.first, f.first.first, f.first.second)];
			int b = vertexID[DOT_INDEX(l.first, f.first.first + f.second.first, f.first.second)];
			int c = vertexID[DOT_INDEX(l.first, f.first.first + f.second.first, f.first.second + f.second.second)];
			int d = vertexID[DOT_INDEX(l.first, f.first.first, f.first.second + f.second.second)];
			aTess->AddFacet(new G4QuadrangularFacet(verCoord[a], verCoord[b], verCoord[c], verCoord[d], ABSOLUTE));
		}
	}
	for (auto l : xNegCom) {
		for (auto f : l.second) {
			int a = vertexID[DOT_INDEX(l.first, f.first.first, f.first.second)];
			int b = vertexID[DOT_INDEX(l.first, f.first.first, f.first.second + f.second.second)];
			int c = vertexID[DOT_INDEX(l.first, f.first.first + f.second.first, f.first.second + f.second.second)];
			int d = vertexID[DOT_INDEX(l.first, f.first.first + f.second.first, f.first.second)];
			aTess->AddFacet(new G4QuadrangularFacet(verCoord[a], verCoord[b], verCoord[c], verCoord[d], ABSOLUTE));
		}
	}
	for (auto l : yPosCom) {
		for (auto f : l.second) {
			int a = vertexID[DOT_INDEX(f.first.second, l.first, f.first.first)];
			int b = vertexID[DOT_INDEX(f.first.second, l.first, f.first.first + f.second.first)];
			int c = vertexID[DOT_INDEX(f.first.second + f.second.second, l.first, f.first.first + f.second.first)];
			int d = vertexID[DOT_INDEX(f.first.second + f.second.second, l.first, f.first.first)];
			aTess->AddFacet(new G4QuadrangularFacet(verCoord[a], verCoord[b], verCoord[c], verCoord[d], ABSOLUTE));
		}
	}
	for (auto l : yNegCom) {
		for (auto f : l.second) {
			int a = vertexID[DOT_INDEX(f.first.second, l.first, f.first.first)];
			int b = vertexID[DOT_INDEX(f.first.second + f.second.second, l.first, f.first.first)];
			int c = vertexID[DOT_INDEX(f.first.second + f.second.second, l.first, f.first.first + f.second.first)];
			int d = vertexID[DOT_INDEX(f.first.second, l.first, f.first.first + f.second.first)];
			aTess->AddFacet(new G4QuadrangularFacet(verCoord[a], verCoord[b], verCoord[c], verCoord[d], ABSOLUTE));
		}
	}
	for (auto l : zPosCom) {
		for (auto f : l.second) {
			int a = vertexID[DOT_INDEX(f.first.first, f.first.second, l.first)];
			int b = vertexID[DOT_INDEX(f.first.first + f.second.first, f.first.second, l.first)];
			int c = vertexID[DOT_INDEX(f.first.first + f.second.first, f.first.second + f.second.second, l.first)];
			int d = vertexID[DOT_INDEX(f.first.first, f.first.second + f.second.second, l.first)];
			aTess->AddFacet(new G4QuadrangularFacet(verCoord[a], verCoord[b], verCoord[c], verCoord[d], ABSOLUTE));
		}
	}
	for (auto l : zNegCom) {
		for (auto f : l.second) {
			int a = vertexID[DOT_INDEX(f.first.first, f.first.second, l.first)];
			int b = vertexID[DOT_INDEX(f.first.first, f.first.second + f.second.second, l.first)];
			int c = vertexID[DOT_INDEX(f.first.first + f.second.first, f.first.second + f.second.second, l.first)];
			int d = vertexID[DOT_INDEX(f.first.first + f.second.first, f.first.second, l.first)];
			aTess->AddFacet(new G4QuadrangularFacet(verCoord[a], verCoord[b], verCoord[c], verCoord[d], ABSOLUTE));
		}
	}
	aTess->SetSolidClosed(true);
	return aTess;
}

G4ThreeVector HDCalculator::GetAPointOnVoxSurface(G4TessellatedSolid* tess)
{
	G4ThreeVector point;
	double rand1 = G4UniformRand();
	for(auto aFacet: tetFacetArea){
		if(rand1>aFacet.first) continue;
		auto   facet = tess->GetFacet(aFacet.second);
		G4ThreeVector vVec = facet->GetVertex(1) - facet->GetVertex(0);
		G4ThreeVector uVec = facet->GetVertex(2) - facet->GetVertex(0);
		vVec = vVec * G4UniformRand();
		uVec = uVec * G4UniformRand();
		point = facet->GetVertex(0) + vVec + uVec;
		break;
	}
	return point;
}

COMPRESSED HDCalculator::Compress2D(map<int, vector<pair<int, int>>> boundaryPix) {
	COMPRESSED compressed;

	for (auto layer : boundaryPix) {
		map<int, vector<int>> pixels;
		map<int, vector<pair<int, int>>> compressedPix1;
		for (auto pix : layer.second) {
			if (pixels.find(pix.first) == pixels.end()) pixels[pix.first] = {pix.second};
			else pixels[pix.first].push_back(pix.second);
		}
		for (auto line : pixels) {
			vector<int> dots = line.second;
			sort(dots.begin(), dots.end());
			int start=dots[0];
			vector<pair<int, int>> comVec;
			for (size_t i = 1; i < dots.size(); i++) {
				if (dots[i] - dots[i-1] == 1) continue;
				comVec.push_back(make_pair(start, dots[i-1]-start+1));
				start = dots[i];
			}
			comVec.push_back(make_pair(start, dots.back()-start+1));
			compressedPix1[line.first] = comVec;
		}
		vector<pair<pair<int, int>, pair<int, int>>> compressedPix2;
		map<pair<int, int>, pair<int, bool>> compressedPix1Map;
		for (auto pix : compressedPix1) {
			for (auto p : pix.second) {
				compressedPix1Map[make_pair(pix.first, p.first)] = make_pair(p.second, false);
			}
		}
		for (auto &pix : compressedPix1Map) {
			if (pix.second.second) continue;

			for (int i = 1; 1; i++) {
				auto nextPix = make_pair(pix.first.first + i, pix.first.second);
				if (compressedPix1Map.find(nextPix)==compressedPix1Map.end() || compressedPix1Map[nextPix].first != pix.second.first) {
					compressedPix2.push_back(make_pair(pix.first, make_pair(i, pix.second.first)));
					break;
				}
				else compressedPix1Map[nextPix].second = true;
			}
		}
		compressed[layer.first] = compressedPix2;
	}
	return compressed;
}

void HDCalculator::GetVertices(COMPRESSED compressed, vector<DOT_INDEX> &vertices, Axis axis) {
	if(axis ==xAxis)
	for (auto l : compressed) {
		for (auto f : l.second) {
			vertices.push_back(DOT_INDEX(l.first, f.first.first, f.first.second));
			vertices.push_back(DOT_INDEX(l.first, f.first.first + f.second.first, f.first.second));
			vertices.push_back(DOT_INDEX(l.first, f.first.first, f.first.second + f.second.second));
			vertices.push_back(DOT_INDEX(l.first, f.first.first + f.second.first, f.first.second + f.second.second));
		}
	}
	else if(axis==yAxis)
		for (auto l : compressed) {
			for (auto f : l.second) {
				vertices.push_back(DOT_INDEX(f.first.second, l.first, f.first.first));
				vertices.push_back(DOT_INDEX(f.first.second, l.first, f.first.first + f.second.first));
				vertices.push_back(DOT_INDEX(f.first.second + f.second.second, l.first, f.first.first));
				vertices.push_back(DOT_INDEX(f.first.second + f.second.second, l.first, f.first.first + f.second.first));
			}
		}
	else if (axis == zAxis)
		for (auto l : compressed) {
			for (auto f : l.second) {
				vertices.push_back(DOT_INDEX(f.first.first, f.first.second, l.first));
				vertices.push_back(DOT_INDEX(f.first.first + f.second.first, f.first.second, l.first));
				vertices.push_back(DOT_INDEX(f.first.first, f.first.second + f.second.second, l.first));
				vertices.push_back(DOT_INDEX(f.first.first + f.second.first, f.first.second + f.second.second, l.first));
			}
		}
}
