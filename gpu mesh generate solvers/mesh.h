#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include "Function.h"
using namespace std;
template<class T>
class Mesh
{
public:
	int Nx;
	int Ny;
	T Lx;
	T Ly;
	vector < vector<int> > boundaryEdge;
	ExpansionFunction<T>* x_expantionRatio;
	ExpansionFunction<T>* y_expantionRatio;



	vector<T> X;//mesh points coodinate
	vector<int> FP;//face point
	vector<int> FC;//face cell
	vector<int> CF;//cell face
	vector<int> CP;//cell point
	vector<int> iF;//inner face
	vector<vector<int>*> bF;//boundary face
	vector<int> bF_1;
	vector<int> bFbegin;



	vector<T> F_area;
	vector<T> F_d;
	vector<T> C_centroid;
	vector<T> F_centroid;
	vector<T> F_normal;
	vector<T> C_volumn;


	int faceNum;
	int innerFaceNum;
	int cellNum;
	int boundaryFaceNum;



	//compress sparse row
	vector<int> IA;
	vector<int> JA;

	vector<int> ON;


	Mesh(){

	}

	Mesh(int Nx, int Ny, T Lx, T Ly, ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio, vector < vector<int> > boundaryEdge){
		this->Nx = Nx;
		this->Ny = Ny;
		this->Lx = Lx;
		this->Ly = Ly;
		this->x_expantionRatio = x_expantionRatio;
		this->y_expantionRatio = y_expantionRatio;
		this->boundaryEdge = boundaryEdge;
	}
	virtual ~Mesh(){
	}
	void infer(){
		faceNum = Nx*(Ny + 1) + (Nx + 1)*Ny;
		cellNum = Nx*Ny;
		innerFaceNum = Nx*(Ny - 1) + (Nx - 1)*Ny;
		boundaryFaceNum = faceNum - innerFaceNum;
	}
	void generateTopology(){
		infer();
		generatePoints();
		generateFaces();
		generateCells();
		generateFC();
		generateCompressedSparseRow();
		generateGeo();
	}
	void generateCompressedSparseRow(){
		IA.resize(cellNum + 1);
		for (int f = 0; f<innerFaceNum; f++){
			int innerF = iF[f];
			int owner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			IA[owner + 1]++;
			IA[neighbor + 1]++;
		}
		for (int i = 1; i<cellNum + 1; i++){
			IA[i]++;
			IA[i] += IA[i - 1];
		}
		JA .resize(IA[cellNum]);
		vector<int> count(cellNum);
		ON.resize(innerFaceNum * 2);
		for (int f = 0; f<innerFaceNum; f++){
			int innerF = iF[f];
			int owner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			JA[IA[owner]] = owner;
			count[owner]++;
			ON[f * 2 + 0] = IA[owner] + count[owner];
			JA[ON[f * 2 + 0]] = neighbor;

			JA[IA[neighbor]] = neighbor;
			count[neighbor]++;
			ON[f * 2 + 1] = IA[neighbor] + count[neighbor];
			JA[ON[f * 2 + 1]] = owner;
		}
	}
	void generatePoints(){
		X.resize((Nx + 1)*(Ny + 1) * 2);
		for (int i = 0; i<Nx + 1; i++){
			for (int j = 0; j<Ny + 1; j++){
				T x = Lx*x_expantionRatio->operator()(1.0*i / Nx, 0);
				T y = Ly*y_expantionRatio->operator()(0, 1.0*j / Ny);
				int p = i*(Ny + 1) + j;
				X[p * 2 + 0] = x;
				X[p * 2 + 1] = y;
			}
		}
	}
	void generateFaces(){
		FP.resize((Nx*(Ny + 1) + (Nx + 1)*Ny) * 2);

		for (int i = 0; i<Nx + 1; i++){
			for (int j = 0; j<Ny; j++){
				int f = i*Ny + j;
				int p1 = i*(Ny + 1) + j;
				int p2 = p1 + 1;
				FP[f * 2 + 0] = p1;
				FP[f * 2 + 1] = p2;
			}
		}


		int xsweep = (Nx + 1)*Ny;
		for (int i = 0; i<Nx; i++){
			for (int j = 0; j<Ny + 1; j++){
				int f = i*(Ny + 1) + j;
				int p1 = i*(Ny + 1) + j;
				int p2 = p1 + Ny + 1;
				FP[(xsweep + f) * 2 + 0] = p1;
				FP[(xsweep + f) * 2 + 1] = p2;
			}
		}
	}
	int findBoudaryEdge(int e) {


		for (int i = 0; i<boundaryEdge.size(); i++){
			auto found = find(boundaryEdge[i].begin(), boundaryEdge[i].end(), e);
			if (found != boundaryEdge[i].end()){
				return i;
			}

		}
		return -1;
	}
	void addPatchFace(vector<vector<int>*> bF, int i, int j, int f, int edge, int owner) {
		int pat = findBoudaryEdge(edge);
		bF[pat]->push_back(f);
		FC[f * 2 + 0] = owner;
		FC[f * 2 + 1] = -1;
	}
	void generateFC(){
		FC.resize((Nx*(Ny + 1) + (Nx + 1)*Ny) * 2);

		for (int i = 0; i < boundaryEdge.size(); i++){
			bF.push_back(new vector<int>());
		}


		for (int i = 0; i < Nx + 1; i++){
			for (int j = 0; j < Ny; j++){
				int f = i*Ny + j;
				if (i == 0){
					addPatchFace(bF, i, j, f, 0, i*Ny + j);
				}
				else if (i == Nx){
					addPatchFace(bF, i, j, f, 2, (i - 1)*Ny + j);
				}
				else{
					iF.push_back(f);
					int cleft = (i - 1)*Ny + j;
					int cright = i*Ny + j;
					FC[f * 2 + 0] = cleft;
					FC[f * 2 + 1] = cright;
				}
			}
		}

		int xsweep = (Nx + 1)*Ny;
		for (int i = 0; i < Nx; i++){
			for (int j = 0; j < Ny + 1; j++){
				int f = i*(Ny + 1) + j;
				if (j == 0){
					addPatchFace(bF, i, j, xsweep + f, 3, i*Ny + j);
				}
				else if (j == Ny){
					addPatchFace(bF, i, j, xsweep + f, 1, i*Ny + j - 1);
				}
				else{
					iF.push_back(xsweep + f);
					int cup = i*Ny + j;
					int cdown = i*Ny + j - 1;
					FC[(xsweep + f) * 2 + 0] = cdown;
					FC[(xsweep + f) * 2 + 1] = cup;
				}
			}
		}

		bFbegin.resize(bF.size() + 1);
		bFbegin[0] = 0;
		int sum = 0;
		for (int i = 0; i<bF.size(); i++){
			sum += bF[i]->size();
			bFbegin[i + 1] = sum;
		}

		for (int pat = 0; pat < bF.size(); ++pat){
			for (int f = 0; f < bF[pat]->size(); f++){
				bF_1.push_back(bF[pat]->at(f));
			}
		}
	}
	void findBoudaryEdge();
	void generateCells(){
		CF.resize(Nx*Ny * 4);
		CP.resize(Nx*Ny * 4);
		for (int i = 0; i < Nx; i++){
			for (int j = 0; j < Ny; j++){
				int c = j + i*Ny;
				int fv1 = j + i*Ny;
				int fv2 = j + (i + 1)*Ny;
				int fh1 = Ny*(Nx + 1) + i*(Ny + 1) + j;
				int fh2 = Ny*(Nx + 1) + i*(Ny + 1) + j + 1;
				CF[c * 4 + 0] = fv1;
				CF[c * 4 + 1] = fv2;
				CF[c * 4 + 2] = fh1;
				CF[c * 4 + 3] = fh2;

				int p0 = i*(Ny + 1) + j;
				int p1 = i*(Ny + 1) + j + 1;
				int p2 = (i + 1)*(Ny + 1) + j + 1;
				int p3 = (i + 1)*(Ny + 1) + j;
				CP[c * 4 + 0] = p0;
				CP[c * 4 + 1] = p1;
				CP[c * 4 + 2] = p2;
				CP[c * 4 + 3] = p3;
			}
		}
	}
	void generateGeo(){
		C_centroid.resize(cellNum * 2);
		C_volumn.resize(cellNum);
		//special for rectangle shape
		for (int c = 0; c<cellNum; c++){
			int p0 = CP[c * 4 + 0];
			int p1 = CP[c * 4 + 1];
			int p2 = CP[c * 4 + 2];
			int p3 = CP[c * 4 + 3];
			T x0 = X[p0 * 2 + 0];
			T x1 = X[p1 * 2 + 0];
			T x2 = X[p2 * 2 + 0];
			T x3 = X[p3 * 2 + 0];
			T xc = 0.25*(x0 + x1 + x2 + x3);
			C_centroid[c * 2 + 0] = xc;
			T y0 = X[p0 * 2 + 1];
			T y1 = X[p1 * 2 + 1];
			T y2 = X[p2 * 2 + 1];
			T y3 = X[p3 * 2 + 1];
			T yc = 0.25*(y0 + y1 + y2 + y3);
			C_centroid[c * 2 + 1] = yc;

			T t_area0 = VectorMath<T>::triangleArea(xc - x0, yc - y0, xc - x1, yc - y1);
			T t_area1 = VectorMath<T>::triangleArea(xc - x1, yc - y1, xc - x2, yc - y2);
			T t_area2 = VectorMath<T>::triangleArea(xc - x2, yc - y2, xc - x3, yc - y3);
			T t_area3 = VectorMath<T>::triangleArea(xc - x3, yc - y3, xc - x0, yc - y0);

			C_volumn[c] = t_area0 + t_area1 + t_area2 + t_area3;

		}
		//Face
		F_centroid.resize(faceNum * 2);
		F_area.resize(faceNum);
		F_normal.resize(faceNum * 2);

		for (int f = 0; f<faceNum; f++){
			int p1 = FP[f * 2 + 0];
			int p2 = FP[f * 2 + 1];
			T x0 = X[p1 * 2 + 0];
			T x1 = X[p2 * 2 + 0];
			F_centroid[f * 2 + 0] = 0.5*(x0 + x1);
			T y0 = X[p1 * 2 + 1];
			T y1 = X[p2 * 2 + 1];
			F_centroid[f * 2 + 1] = 0.5*(y0 + y1);
			F_area[f] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));

			T* orth = VectorMath<T>::coorthogonal(x0 - x1, y0 - y1);
			int ownner = FC[f * 2 + 0];
			//int neighbor=FC[f*2+1];
			T cx0 = C_centroid[ownner * 2 + 0];
			T cy0 = C_centroid[ownner * 2 + 1];
			T cx1 = F_centroid[f * 2 + 0];
			T cy1 = F_centroid[f * 2 + 1];
			T* norm = VectorMath<T>::normalize(orth[0], orth[1]);
			if (VectorMath<T>::dot(cx1 - cx0, cy1 - cy0, orth[0], orth[1])<0){
				norm = VectorMath<T>::reverse(norm[0], norm[1]);
			}
			F_normal[f * 2 + 0] = norm[0];
			F_normal[f * 2 + 1] = norm[1];
			delete[] norm;
			delete[] orth;
		}

		//F_d
		F_d.resize(faceNum);
		for (int f = 0; f<iF.size(); f++){
			int innerF = iF[f];
			int ownner = FC[innerF * 2 + 0];
			int neighbor = FC[innerF * 2 + 1];
			T x0 = C_centroid[ownner * 2 + 0];
			T y0 = C_centroid[ownner * 2 + 1];
			T x1 = C_centroid[neighbor * 2 + 0];
			T y1 = C_centroid[neighbor * 2 + 1];
			F_d[innerF] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
		}
		for (auto pat : bF){
			for (auto boundaryF : (*pat)){
				int ownner = FC[boundaryF * 2 + 0];
				T x0 = C_centroid[ownner * 2 + 0];
				T y0 = C_centroid[ownner * 2 + 1];
				T x1 = F_centroid[boundaryF * 2 + 0];
				T y1 = F_centroid[boundaryF * 2 + 1];
				F_d[boundaryF] = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
			}
		}
	}

};