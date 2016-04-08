#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include "Function.h"
#include <fstream>
#include <tuple>
using namespace std;
template<class T>
class Mesh
{
public:
	//int Nx;
	//int Ny;
	//T Lx;
	//T Ly;
	//vector < vector<int> > boundaryEdge;
	//ExpansionFunction<T>* x_expantionRatio;
	//ExpansionFunction<T>* y_expantionRatio;



	vector<vector<T>> X_c;//mesh points coodinate ,colummn based;
	vector < vector<T>> X_r;//mesh points coodinates, row based;
	vector<vector<int>> FP;//face point
	//vector<int> FC;//face cell
	//vector<int> CF;//cell face
	vector<vector<int>> CP;//cell point
	//vector<int> iF;//inner face
	//vector<vector<int>*> bF;//boundary face
	//vector<int> bF_1;
	//vector<int> bFbegin;



	vector<T> F_area;
	vector<T> F_d;
	vector<T> C_centroid;
	vector<T> F_centroid;
	vector<T> F_normal;
	vector<T> C_volumn;
	vector<tuple<int, int, string>> physicalSet;

	int faceNum;
	int innerFaceNum;
	int cellNum;
	int boundaryFaceNum;



	//compress sparse row
	vector<int> IA;
	vector<int> JA;

	vector<int> ON;
	int nodeNum;

	Mesh(){

	}
	

	void ReadFromGshFile(string fileName){
		ifstream input(fileName);
		string line;
		int level = 0;
		while (true){
			getline(input, line);
			if (input.eof())
			{
				break;
			}
			if (line == "$PhysicalNames"){
				int physicalNum;
				input >> physicalNum;
				for (int p = 0; p < physicalNum; p++){
					int d;
					int numberId;
					string name;	
					input >> d >> numberId >> name;
					trim_quote(name);
					tuple<int, int, string> t(d,numberId,name);
					physicalSet.push_back(t);
					//cout << d << "  " << numberId << " " << name << endl;
					//cout << get<0>(t) <<" "<<get<1>(t)<<" "<<get<2>(t)<< endl;
				}				
			}
			if (line == "$Nodes"){
				input >> nodeNum;
				X_c.resize(nodeNum);
				//X_r.push_back(vector<T>());
				//X_r.push_back(vector<T>());
				X_r.resize(2);
				X_r[0].resize(nodeNum);
				X_r[1].resize(nodeNum);
				for (int n = 0; n < nodeNum; n++){
					int orderNum;
					double x;
					double y;
					double z;
					input >> orderNum >> x >> y >> z;
					X_c[n].resize(2);
					X_c[n][0] = x;
					X_c[n][1] = y;
					X_r[0][n] = x;
					X_r[1][n] = y;					
				}
			}
			if (line == "$Elements"){
				int elementNum;
				input >> elementNum;
				for (int e = 0; e < elementNum; e++){
					int orderNum, elementType, reserveWord, physicalNumId, baseUnitId;
					input >> orderNum >> elementType >> reserveWord >> physicalNumId>> baseUnitId;
					if (elementType == 1){
						FP.push_back(vector<int>());
						FP[FP.size() - 1].resize(2);
						int p0, p1;
						input >> p0 >> p1;
						FP[FP.size() - 1][0]=p0;
						FP[FP.size() - 1][1] = p1;
					}
					else if (elementType == 2){
						CP.push_back(vector<int>());
						CP[CP.size() - 1].resize(3);
						int p0, p1, p2;
						input >> p0 >> p1>>p2;
						CP[CP.size() - 1][0] = p0;
						CP[CP.size() - 1][1] = p1;
						CP[CP.size() - 1][2] = p2;
					}
					else if (elementType == 3){
						CP.push_back(vector<int>());
						CP[CP.size() - 1].resize(4);
						int p0, p1, p2,p3;
						input >> p0 >> p1 >> p2>>p3;
						CP[CP.size() - 1][0] = p0;
						CP[CP.size() - 1][1] = p1;
						CP[CP.size() - 1][2] = p2;
						CP[CP.size() - 1][3] = p3;
					}
				}
			}

		}
		
	}


	void trim_quote(string &s){		
		s.erase(0, 1);
		s.erase(s.size() - 1, 1);
	}

	virtual ~Mesh(){
	}


};