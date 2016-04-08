template<class T> class I_Temporal{
public :
	virtual CSR<T> apply(vector<CellField<T>> phi, T dt, Mesh<T> mesh) = 0;
};





template<class T> class ImplictEuler:public I_Temporal<T>{
public:

	//ImplictEuler<T>(){}
	virtual CSR<T> apply(vector<CellField<T>> phi, T dt,Mesh<T> mesh) {
		CSR<T> eq (mesh);
		CellField<T> phi0 = phi[0];
		for (int c = 0; c<mesh.cellNum; c++){
			eq.A[mesh.IA[c]] = mesh.C_volumn[c] / dt;
			eq.b[c] = phi0.inner[c] * mesh.C_volumn[c] / dt;
		}
		return eq;
	}
};