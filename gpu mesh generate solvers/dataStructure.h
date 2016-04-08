#include <vector>
template<class T> class CSR{
public :	
	vector<T> A;
	vector<T> b;
	CSR<T>(){};
	CSR(Mesh<T> mesh){
		A.resize(mesh.JA.size());
		b.resize(mesh.cellNum);
	}

	static CSR<T> plus(CSR<T> E1, CSR<T> E2){
		CSR<T> E3(E1);
		//transform(E2.A.begin(), E2.A.end(), E2.A.begin(), coefficents, multiplies<T>());
		transform(E1.A.begin(), E1.A.end(), E2.A.begin(),  E3.A.begin(), std::plus<T>());
		transform(E1.b.begin(), E1.b.end(), E2.b.begin(), E3.b.begin(), std::plus<T>());
		return E3;
	}

	static CSR<T> minus(CSR<T> E1, CSR<T> E2){
		CSR<T> E3 (E1);
		transform(E2.A.begin(), E2.A.end(), E2.A.begin(), negate<T>());
		transform(E1.A.begin(), E1.A.end(), E2.A.begin(), E3.A.begin(), std::plus<T>());
		transform(E2.b.begin(), E2.b.end(), E2.b.begin(), negate<T>());
		transform(E1.b.begin(), E1.b.end(), E2.b.begin(), E3.b.begin(), std::plus<T>());
		return E3;
	}
};

template<class T> class CellField{
public:
	vector<T> inner;
	vector<T> boundary;
	vector<vector<T>*> boundaryCondition;


	void setUniformValue(Mesh<T> mesh, T value){
		inner.resize(mesh.cellNum);
		boundary.resize(mesh.boundaryFaceNum);
		inner.assign(inner.size(), value);
		boundary.assign(boundary.size(), value);

	}
	void assignBoundary(Mesh<T> mesh) {
		for (int p = 0; p<mesh.bFbegin.size() - 1; p++){
			for (int f = mesh.bFbegin[p]; f<mesh.bFbegin[p + 1]; f++){
				int boundaryF = mesh.bF_1[f];
				double a = boundaryCondition[p]->at(0);
				double b = boundaryCondition[p]->at(1);
				double c = boundaryCondition[p]->at(2);
				double d = mesh.F_d[boundaryF];
				int ownner = mesh.FC[boundaryF * 2 + 0];
				double phi = inner[ownner];
				boundary[f] = (a / d*phi + c) / (a / d + b);
			}
		}
	}
};