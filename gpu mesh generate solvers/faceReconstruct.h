
template <class T> class I_FaceReconstruct{
public :
	virtual vector<T> apply(CellField<T> phi,Mesh<T> mesh)=0;
};

template <class T> class FaceReconstruct_BEB :public I_FaceReconstruct<T>{
public:
	virtual vector<T> apply(CellField<T> gama, Mesh<T> mesh) {
			vector<T> faceValue(mesh.faceNum);
			for (int f = 0; f<mesh.iF.size(); f++){
				int innerF = mesh.iF[f];
				int ownner = mesh.FC[innerF * 2 + 0];
				int neighbor = mesh.FC[innerF * 2 + 1];
				T x0 = mesh.C_centroid[ownner * 2 + 0];
				T y0 = mesh.C_centroid[ownner * 2 + 1];
				T x1 = mesh.F_centroid[innerF * 2 + 0];
				T y1 = mesh.F_centroid[innerF * 2 + 1];
				T x2 = mesh.C_centroid[neighbor * 2 + 0];
				double y2 = mesh.C_centroid[neighbor * 2 + 1];
				double g1 = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
				double g2 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
				(faceValue)[innerF] = g1 / (g1 + g2)*gama.inner[neighbor] + g2 / (g1 + g2)*gama.inner[ownner];
			}
			boundaryHandle(gama, mesh,faceValue);
			return faceValue;		
	}

	virtual void boundaryHandle(CellField<T> gama, Mesh<T> mesh,vector<T> &faceValue) {
		for (int f = 0; f<mesh.boundaryFaceNum; f++){
			int boundaryF = mesh.bF_1[f];
			faceValue[boundaryF] = gama.boundary[f];
		}
	}
};

template<class T> class LinearFaceReconstruct_BEO :public FaceReconstruct_BEB<T>{

public:
	virtual  void boundaryHandle(CellField<T> gama, Mesh<T> mesh,vector<T> &faceValue) {
		for (int f = 0; f<mesh.boundaryFaceNum; f++){
			int boundaryF = mesh.bF_1[f];
			int owner = mesh.FC[boundaryF * 2 + 0];
			faceValue[boundaryF] = gama.inner[owner];
		}
	}
};



template<class T> class B_RhieChow{
	LinearFaceReconstruct_BEO<T> beo;
	FaceReconstruct_BEB<T> beb;
public:
	vector<vector<T>> calculateD(CSR<T> u_momentum, CSR<T> v_momentum,Mesh<T> mesh){
		vector<vector<T>> D_f(2, vector<T>());
		D_f[0].resize(mesh.cellNum);
		D_f[1].resize(mesh.cellNum);
		vector<CellField<T>> D(2,CellField<T>());
		D[0].setUniformValue(mesh,0);
		D[1].setUniformValue(mesh,0);
		for (int c = 0; c<mesh.cellNum; c++){
			D[0].inner[c] = mesh.C_volumn[c] / u_momentum.A[mesh.IA[c]];
			D[1].inner[c] = mesh.C_volumn[c] / v_momentum.A[mesh.IA[c]];
		}
		D_f[0] = beo.apply(D[0],mesh);
		D_f[1] = beo.apply(D[1], mesh);
		return D_f;
	}
	vector<vector<T>> reconstructVelocity(vector<vector<T>> D_f, vector< CellField<T>>& p_cellGrad, vector<vector<T>> p_faceGrad, vector<CellField<T>> velocity, Mesh<T> mesh){
		
		vector<vector<T>> v_RC(2, vector<T>());
		v_RC[0].resize(mesh.faceNum);
		v_RC[1].resize(mesh.faceNum);
		vector<T> p_avrgFaceGrad_x = beo.apply(p_cellGrad[0],mesh);
		vector<T> p_avrgFaceGrad_y = beo.apply(p_cellGrad[1],mesh);
		vector<T> u_avrg = beb.apply(velocity[0],mesh);
		vector<T> v_avrg = beb.apply(velocity[1],mesh);
		for (int f = 0; f<mesh.faceNum; f++){
			v_RC[0][f] = u_avrg[f] + D_f[0][f] * (p_avrgFaceGrad_x[f] - p_faceGrad[0][f]);
			v_RC[1][f] = v_avrg[f] + D_f[1][f] * (p_avrgFaceGrad_y[f] - p_faceGrad[1][f]);
		}
		return v_RC;
	}
};