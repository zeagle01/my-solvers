
#include "dataStructure.h"
template<class T> class I_Laplace{
	public:
		virtual CSR<T> apply(vector<T> gama_f, CellField<T> phi,Mesh<T> mesh)=0;
		virtual CSR<T> apply( T gama_f, CellField<T> phi, Mesh<T> mesh) = 0;
};

template<class T> class Laplace:public I_Laplace<T>{
public:


	virtual CSR<T> apply(T gama_f,CellField<T> phi, Mesh<T> mesh ){
		vector<T> gama_field(mesh.faceNum, gama_f);
		return apply(gama_field, phi, mesh);
	}

	virtual CSR<T> apply(vector<T> gama_f, CellField<T> phi,Mesh<T> mesh){
		CSR<T> *eq = new CSR<T>();
		eq-> A .resize(mesh.JA.size());
		eq-> b.resize(mesh.cellNum);

		for (int f = 0; f < mesh.iF.size(); f++){
			int inerF = mesh.iF[f];
			int owner = mesh.FC[inerF * 2 + 0];
			int neighbor = mesh.FC[inerF * 2 + 1];
			T S_f = mesh.F_area[inerF];
			T d = mesh.F_d[inerF];
			T diffusion_on = gama_f[inerF] * S_f / d;
			T diffusion_oo = -diffusion_on;
			T diffusion_no = diffusion_on;
			T diffusion_nn = -diffusion_on;

			eq->A[mesh.IA[owner]] += diffusion_oo;
			eq->A[mesh.ON[f * 2 + 0]] += diffusion_on;
			eq->A[mesh.IA[neighbor]] += diffusion_nn;
			eq->A[mesh.ON[f * 2 + 1]] += diffusion_no;
		}
		
		for (int pat = 0; pat < mesh.bF.size(); pat++) {
			for (int f = 0; f < mesh.bF[pat]->size(); f++) {
				int boundaryF=mesh.bF[pat]->at(f);
				int owner = mesh.FC[boundaryF * 2 + 0];
				T S_f = mesh.F_area[boundaryF];
				T d = mesh.F_d[boundaryF];

				T diffusion_on = gama_f[boundaryF] * S_f / d;
				T diffusion_oo = -diffusion_on;
				phi.boundaryCondition[pat][1];
				T a = phi.boundaryCondition[pat ]->at(0);
				T b1 = phi.boundaryCondition[pat]->at(1);
				T c = phi.boundaryCondition[pat]->at(2);
				eq->A[mesh.IA[owner]] += diffusion_oo + a / d / (a / d + b1)*diffusion_on;
				eq->b[owner] -= c / (a / d + b1)*diffusion_on;
			}
		}
		return *eq;
	}


};