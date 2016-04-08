
#ifndef DATASTRUCTURE_H_H
#define DATASTRUCTURE_H_H



#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
class CSR{
public:
	vector<double> A;
	vector<double> b;
	CSR(){};
	CSR(Mesh mesh){
		A.resize(mesh.JA.size());
		b.resize(mesh.cellNum);
	}

	static CSR plus(CSR E1, CSR E2){
		CSR E3(E1);
		//transform(E2.A.begin(), E2.A.end(), E2.A.begin(), coefficents, multiplies<double>());
		transform(E1.A.begin(), E1.A.end(), E2.A.begin(), E3.A.begin(), std::plus<double>());
		transform(E1.b.begin(), E1.b.end(), E2.b.begin(), E3.b.begin(), std::plus<double>());
		return E3;
	}

	static CSR minus(CSR E1, CSR E2){
		CSR E3(E1);
		transform(E2.A.begin(), E2.A.end(), E2.A.begin(), negate<double>());
		transform(E1.A.begin(), E1.A.end(), E2.A.begin(), E3.A.begin(), std::plus<double>());
		transform(E2.b.begin(), E2.b.end(), E2.b.begin(), negate<double>());
		transform(E1.b.begin(), E1.b.end(), E2.b.begin(), E3.b.begin(), std::plus<double>());
		return E3;
	}
};

class CellField{
public:
	vector<double> inner;
	vector<double> boundary;
	vector<vector<double>> boundaryCondition;

	CellField(){}
	CellField(Mesh mesh){
		inner.resize(mesh.cellNum);
		boundary.resize(mesh.boundaryFaceNum);
		//boundaryCondition.resize(mesh.boundaryEdge.size());
	}
	void setUniformValue(Mesh mesh, double value){
		inner.resize(mesh.cellNum);
		boundary.resize(mesh.boundaryFaceNum);
		inner.assign(inner.size(), value);
		boundary.assign(boundary.size(), value);

	}
	void assignBoundary(Mesh mesh) {
		for (int p = 0; p<mesh.bFbegin.size() - 1; p++){
			for (int f = mesh.bFbegin[p]; f<mesh.bFbegin[p + 1]; f++){
				int boundaryF = mesh.bF_1[f];
				double a = boundaryCondition[p][0];
				double b = boundaryCondition[p][1];
				double c = boundaryCondition[p][1];
				double d = mesh.F_d[boundaryF];
				int ownner = mesh.FC[boundaryF * 2 + 0];
				double phi = inner[ownner];
				boundary[f] = (a / d*phi + c) / (a / d + b);
			}
		}
	}
};

#endif