
template<class T> class I_GaussCellGradient {
public:
	virtual vector<CellField<T>> apply(vector<T> faceValue, Mesh<T> mesh) = 0;
};


template<class T> class GaussCellGradient :public I_GaussCellGradient<T>{
	virtual vector<CellField<T>> apply(vector<T> faceValue, Mesh<T> mesh){

		vector<CellField<T>> grad(2, CellField<T>());
		grad[0].setUniformValue(mesh,0);
		grad[1].setUniformValue(mesh, 0);
		
		for (int f = 0; f<mesh.iF.size(); f++){
			int innerF = mesh.iF[f];
			int ownner = mesh.FC[innerF * 2 + 0];
			int neighbor = mesh.FC[innerF * 2 + 1];
			grad[0].inner[ownner] += mesh.F_normal[innerF * 2 + 0] * mesh.F_area[innerF] / mesh.C_volumn[ownner] * faceValue[innerF];
			grad[1].inner[ownner] += mesh.F_normal[innerF * 2 + 1] * mesh.F_area[innerF] / mesh.C_volumn[ownner] * faceValue[innerF];

			grad[0].inner[neighbor] -= mesh.F_normal[innerF * 2 + 0] * mesh.F_area[innerF] / mesh.C_volumn[neighbor] * faceValue[innerF];
			grad[1].inner[neighbor] -= mesh.F_normal[innerF * 2 + 1] * mesh.F_area[innerF] / mesh.C_volumn[neighbor] * faceValue[innerF];
		}

		for (int f = 0; f<mesh.boundaryFaceNum; f++){
			int boundaryF = mesh.bF_1[f];
			int ownner = mesh.FC[boundaryF * 2 + 0];
			grad[0].inner[ownner] += mesh.F_normal[boundaryF * 2 + 0] * mesh.F_area[boundaryF] / mesh.C_volumn[ownner] * faceValue[boundaryF];
			grad[1].inner[ownner] += mesh.F_normal[boundaryF * 2 + 1] * mesh.F_area[boundaryF] / mesh.C_volumn[ownner] * faceValue[boundaryF];
		}
		return grad;
	}
	
};



template<class T> class I_FaceGradient {
public:
	virtual vector<vector<T>> apply(CellField<T> phi , Mesh<T> mesh) = 0;
};

template<class T> class PlainFaceGradient :public  I_FaceGradient<T>{

public :
	vector<vector<T>> apply(CellField<T> phi,Mesh<T> mesh) {
		vector<vector<T>> grad (2,vector<T>());
		grad[0].resize(mesh.faceNum);
		grad[1].resize(mesh.faceNum);
		for (int f = 0; f<mesh.iF.size(); f++){
			int innerF = mesh.iF[f];
			int owner = mesh.FC[innerF * 2 + 0];
			int neighbor = mesh.FC[innerF * 2 + 1];
			T phi2 = phi.inner[neighbor];
			T phi1 = phi.inner[owner];
			T temp = (phi2 - phi1) / mesh.F_d[innerF];
			grad[0][innerF] = temp*mesh.F_normal[innerF * 2 + 0];
			grad[1][innerF] = temp*mesh.F_normal[innerF * 2 + 1];
		}
		for (int f = 0; f<mesh.boundaryFaceNum; f++){
			int boundaryF = mesh.bF_1[f];
			int owner = mesh.FC[boundaryF * 2 + 0];
			T phi2 = phi.boundary[f];
			T phi1 = phi.inner[owner];
			T temp = (phi2 - phi1) / mesh.F_d[boundaryF];
			grad[0][boundaryF] = temp*mesh.F_normal[boundaryF * 2 + 0];
			grad[1][boundaryF] = temp*mesh.F_normal[boundaryF * 2 + 1];
		}
		return grad;
	}

};
