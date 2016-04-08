#ifndef PRINTER_H_H
#define PRINTER_H_H


#include "mesh.h"
#include <fstream>
#include <vector>
#include "dataStructure.h"
#include "serviceLocator.h"
using namespace std;
class  A_Printer:public interface_t{

public:
	string cases = "cases ";
	virtual void print(string fileName, Mesh mesh, vector<CellField> phi) = 0;
};

class  Printer :public A_Printer,public member_t<Printer>{
public:
	string cases = "cases ";
	void printFileHead(ofstream &fout, Mesh mesh, vector<CellField> phi){
		fout << "variables = \"x\", \"y\"," << " ";
		for (int vn = 0; vn < phi.size(); vn++){
			fout << " \"v" << vn << "\",";
		}
		fout << endl;
		fout << "ZONE I = " << mesh.Ny << " J = " << mesh.Nx << " F = POINT";
		fout << endl;
	}
	void printBoundary(ofstream &fout, Mesh mesh, vector<CellField> phi){
		for (int f = 0; f < mesh.boundaryFaceNum; f++){
			int boundaryF = mesh.bF_1[f];
			fout << mesh.F_centroid[boundaryF * 2 + 0] << " " << mesh.F_centroid[boundaryF * 2 + 1] << " ";
			for (int vn = 0; vn < phi.size(); vn++){
				fout << phi[vn].boundary[f] << " ";
			}
			fout << endl;
		}
	}
	void printInner(ofstream &fout, Mesh mesh, vector<CellField> phi){
		for (int c = 0; c < mesh.cellNum; ++c){
			fout << mesh.C_centroid[c * 2 + 0] << " " << mesh.C_centroid[c * 2 + 1] << " ";
			for (int vn = 0; vn < phi.size(); vn++){
				fout << phi[vn].inner[c] << " ";
			}
			fout << endl;
		}
	}
	void print(string fileName, Mesh mesh, vector<CellField> phi)
	{
		ofstream fout(cases + fileName);

		printFileHead(fout, mesh, phi);
		printInner(fout, mesh, phi);


	}
};


class  Printer1 :public Printer, public member_t<Printer>{
public:
	void print(string fileName, Mesh mesh, vector<CellField> phi)
	{
		ofstream fout(cases + fileName);
		printFileHead(fout, mesh, phi);
		printBoundary(fout, mesh, phi);
		printInner(fout, mesh, phi);


	}
};


#endif