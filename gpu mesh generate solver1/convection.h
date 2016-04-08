#include "Tool.h"
#include "faceReconstruct.h"
template<class T> class I_Convection{
public:
	virtual CSR<T> apply(vector<T> u_f, vector<T> v_f, CellField<T> phi, Mesh<T> mesh) = 0;
};

template<class T> class UpwindConvection:public I_Convection<T>{
public:
	virtual CSR<T> apply(vector<T> u_f, vector<T> v_f, CellField<T> phi, Mesh<T> mesh){
		CSR<T> eq(mesh);
		//copy(eq->A.begin(), eq->A.end(), ostream_iterator<T>(cout, " "));

		for (int f = 0; f < mesh.iF.size(); f++) {
			int inerF = mesh.iF[f];
			int owner = mesh.FC[inerF * 2 + 0];
			int neighbor = mesh.FC[inerF * 2 + 1];

			T Sf_x = mesh.F_area[inerF] * mesh.F_normal[inerF * 2 + 0];
			T Sf_y = mesh.F_area[inerF] * mesh.F_normal[inerF * 2 + 1];
			T mdot = u_f[inerF] * Sf_x + v_f[inerF] * Sf_y;
			//cout << mdot << endl;
			eq.A[mesh.ON[f * 2 + 0]] = VectorMath<T>::myMin(mdot, 0);
			eq.A[mesh.IA[owner]] += VectorMath<T>::myMax(mdot, 0);

			eq.A[mesh.ON[f * 2 + 1]] = VectorMath<T>::myMin(-mdot, 0);
			eq.A[mesh.IA[neighbor]] += VectorMath<T>::myMax(-mdot, 0);
			//copy(eq->A.begin(), eq->A.end(), ostream_iterator<T>(cout, " "));

		}
		for (int pat = 0; pat < mesh.bFbegin.size() - 1; pat++) {
			for (int f = mesh.bFbegin[pat]; f < mesh.bFbegin[pat + 1]; f++) {
				int boundaryF = mesh.bF_1[f];
				int owner = mesh.FC[boundaryF * 2 + 0];

				T Sf_x = mesh.F_area[boundaryF] * mesh.F_normal[boundaryF * 2 + 0];
				T Sf_y = mesh.F_area[boundaryF] * mesh.F_normal[boundaryF * 2 + 1];
				T mdot = u_f[boundaryF] * Sf_x + v_f[boundaryF] * Sf_y;

				T convection_oo = VectorMath<T>::myMax(mdot, 0);
				T convection_on = VectorMath<T>::myMin(mdot, 0);


				T a = phi.boundaryCondition[pat]->at(0);
				T b1 = phi.boundaryCondition[pat]->at(1);
				T c = phi.boundaryCondition[pat]->at(2);
				T d = mesh.F_d[boundaryF];

				eq.A[mesh.IA[owner]] += convection_oo + a / d / (a / d + b1)*convection_on;
				eq.b[owner] -= c / (a / d + b1)*convection_on;
			}
		}
		
		return eq;
	}
};