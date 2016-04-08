#include <string>

template<class T> class  A_Printer{
public:
	string cases="cases ";
	virtual void print(string fileName, Mesh<T> mesh, vector<CellField<T>*> phi) = 0;
};

template<class T> class  Printer:public A_Printer<T>{
public:
	string cases = "cases ";
	void printFileHead(ofstream &fout, Mesh<T> mesh, vector<CellField<T>*> phi){
		fout << "variables = \"x\", \"y\"," << " ";
		for (int vn = 0; vn < phi.size(); vn++){
			fout << " \"v" << vn << "\",";
		}
		fout << endl;
		fout << "ZONE I = " << mesh.Nx << " J = " << mesh.Ny << " F = POINT";
		fout << endl;
	}
	void printBoundary(ofstream &fout, Mesh<T> mesh, vector<CellField<T>*> phi){
		for (int f = 0; f < mesh.boundaryFaceNum; f++){
			int boundaryF = mesh.bF_1[f];
			fout << mesh.F_centroid[boundaryF * 2 + 0] << " " << mesh.F_centroid[boundaryF * 2 + 1] << " ";
			for (int vn = 0; vn < phi.size(); vn++){
				fout << phi[vn]->boundary[f] << " ";
			}
			fout << endl;
		}
	}
	void printInner(ofstream &fout, Mesh<T> mesh, vector<CellField<T>*> phi){
		for (int c = 0; c < mesh.cellNum; ++c){
			fout << mesh.C_centroid[c * 2 + 0] << " " << mesh.C_centroid[c * 2 + 1] << " ";
			for (int vn = 0; vn < phi.size(); vn++){
				fout << phi[vn]->inner[c] << " ";
			}
			fout << endl;
		}
	}
	void print(string fileName, Mesh<T> mesh, vector<CellField<T>*> phi)
	{
		ofstream fout(cases + fileName);

		printFileHead(fout, mesh,phi);
		printInner(fout, mesh, phi);
		
		
	}
};


template<class T> class  Printer1 :public Printer<T>{
public:
	void print(string fileName, Mesh<T> mesh, vector<CellField<T>*> phi)
	{
		ofstream fout(cases + fileName);
		printFileHead(fout, mesh, phi);
		printBoundary(fout, mesh, phi);
		printInner(fout, mesh, phi);


	}
};