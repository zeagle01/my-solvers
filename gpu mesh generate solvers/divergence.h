template<class T> class ExplicitDivergence {

public:
	vector<T> apply(vector<vector<T>> V_f,Mesh<T> mesh){
		vector<T> div (mesh.cellNum);
		for (int f = 0; f<mesh.iF.size(); f++){
			int innerF = mesh.iF[f];
			int ownner = mesh.FC[innerF * 2 + 0];
			int neighbor = mesh.FC[innerF * 2 + 1];
			double face_x = mesh.F_area[innerF] * mesh.F_normal[innerF * 2 + 0];
			double face_y = mesh.F_area[innerF] * mesh.F_normal[innerF * 2 + 1];

			double flux = face_x*V_f[0][innerF] + face_y*V_f[1][innerF];
			div[ownner] += flux / mesh.C_volumn[ownner];
			div[neighbor] -= flux / mesh.C_volumn[neighbor];
		}
		for (int f = 0; f<mesh.boundaryFaceNum; f++){
			int boundaryF = mesh.bF_1[f];
			int ownner = mesh.FC[boundaryF * 2 + 0];
			double face_x = mesh.F_area[boundaryF] * mesh.F_normal[boundaryF * 2 + 0];
			double face_y = mesh.F_area[boundaryF] * mesh.F_normal[boundaryF * 2 + 1];
			double flux = face_x*V_f[0][boundaryF] + face_y*V_f[1][boundaryF];
			div[ownner] += flux / mesh.C_volumn[ownner];
		}
		return div;
	}
};
