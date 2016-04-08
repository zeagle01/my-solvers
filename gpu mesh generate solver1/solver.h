




#pragma once
#include "mesh.h"
#include "diffusion.h"
#include "convection.h"
#include "LinearEquationSolver.h"
#include "gpu_solver.h"
#include <fstream>
#include "printer.h"
template<class T> class Solver{
public:
	Mesh<T> mesh;
	int Nx;
	int Ny;
	T Lx;
	T Ly;
	vector < vector<int> > boundaryEdge;
	ExpansionFunction<T>* x_expantionRatio;
	ExpansionFunction<T>* y_expantionRatio;
	int variableNum;
	T convergenceThrehold;
	vector<T> lesConvergenceThrehold;

	A_Printer<T>* printer;
	vector<CellField<T>*> phi;
	vector<I_LinearEquationSolver<T>*> les;

	vector<vector<vector<T>>>* boundaryConditions;
	void buildVariable(){
		for (int i = 0; i < variableNum; i++){
			phi.push_back(new CellField<T>());
			phi[i]->inner.resize(mesh.cellNum);
			phi[i]->boundary.resize(mesh.boundaryFaceNum);
			for (int j = 0; j < boundaryEdge.size(); j++){
				phi[i]->boundaryCondition.push_back(  &(*boundaryConditions)[j][i]);
			}
			
		}
	}
	void setLinearEquationSolverConvergenceThrehold(){
		for (int i = 0; i < les.size(); ++i){
			les[i]->set_convergedThrehold(lesConvergenceThrehold[i]);
		}
	}


	void config(){
		buildMesh();
		buildVariable();
		setLinearEquationSolverConvergenceThrehold();
	}
	Solver(string configFile){
		boost::property_tree::ptree pt;
		boost::property_tree::json_parser::read_json(configFile, pt);
		Nx=pt.get_child("Nx").get_value<int>();
		Ny=pt.get_child("Ny").get_value<int>();
		Lx = pt.get_child("Lx").get_value<T>();
		Ly = pt.get_child("Ly").get_value<T>();
		variableNum = pt.get_child("variableNum").get_value<int>();
		//x_expantionRatio = (ExpansionFunction<T>)(services.at(pt.get_child("x_expansion").get_value<string>()));
		
	}
	Solver(int Nx, int Ny, T Lx, T Ly, int variableNum,
		ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio, 
		vector < vector<int> > boundaryEdge,
		vector<vector<vector<T>>>* boundaryConditions,
		vector<I_LinearEquationSolver<T>*> les, vector<T> lesConvergenceThrehold,
		A_Printer<T>* printer,
		T convergenceThrehold
		){
		this->Nx = Nx;
		this->Ny = Ny;
		this->Lx = Lx;
		this->Ly = Ly;
		this->variableNum = variableNum;
		this->x_expantionRatio = x_expantionRatio;
		this->y_expantionRatio = y_expantionRatio;
		this->boundaryEdge = boundaryEdge;
		this->boundaryConditions = boundaryConditions;
		this->les = les;
		this->lesConvergenceThrehold = lesConvergenceThrehold;
		this->printer = printer;
		this->convergenceThrehold = convergenceThrehold;
	}

	void buildMesh(){
		mesh.x_expantionRatio = x_expantionRatio;
		mesh.y_expantionRatio = y_expantionRatio;
		mesh.Nx = Nx;
		mesh.Ny = Ny;
		mesh.Lx = Lx;
		mesh.Ly = Ly;
		mesh.boundaryEdge = boundaryEdge;
		mesh.generateTopology();
	}


	virtual void solve()=0;
};










template<class T> class DiffusionSolver :public Solver<T>{
public:

	I_Laplace<T>* diffusion;
	T gama;
	DiffusionSolver(string configFile) :Solver(configFile){

	}
	DiffusionSolver(int Nx, int Ny, T Lx, T Ly, T gama,
		ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio,
		vector < vector<int> > boundaryEdge, vector<vector<vector<T>>>* boundaryConditions,
		I_Laplace<T>* diffusion,
		vector<I_LinearEquationSolver<T>*> les, vector<T> lesConvergenceThrehold,
		A_Printer<T>* printer,
		T convergenceThrehold
		)

		:Solver(Nx, Ny, Lx, Ly, 1, x_expantionRatio, y_expantionRatio, boundaryEdge, boundaryConditions, les, lesConvergenceThrehold, printer ,convergenceThrehold),
		gama(gama), diffusion(diffusion){
	
	}

	virtual void solve(){
		config();
		CSR<T> eq = diffusion->apply(gama, *(phi[0]), mesh);
		les[0]->solve(eq, phi[0], mesh);
		phi[0]->assignBoundary(mesh);	
		printer->print("out.dat", mesh, phi);

		//std::ofstream fout("a.dat");
		//copy(phi[0]->inner.cbegin(), phi[0]->inner.cend(), ostream_iterator<T>(fout, "\n") );
		//for_each(phi[0]->inner.begin(), phi[0]->inner.end(), print<T>());
	}
};



template<class T> class ConvectionDiffusionSolver :public Solver<T>{
public:

	I_Laplace<T>* diffusion;
	I_Convection<T>* convection;
	vector<I_FaceReconstruct<T>*> faceReconstruct;
	vector<T> uniformVelocity;
	T gama;

	ConvectionDiffusionSolver(int Nx, int Ny, T Lx, T Ly, T gama,
		ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio,
		vector < vector<int> > boundaryEdge, vector<vector<vector<T>>>* boundaryConditions,
		vector<T> uniformVelocity,
		I_Laplace<T>* diffusion, I_Convection<T>* convection,
		vector<I_FaceReconstruct<T>*> faceReconstruct,
		vector<I_LinearEquationSolver<T>*> les, vector<T> lesConvergenceThrehold,
		A_Printer<T>* printer,
		T convergenceThrehold
		)

		:Solver(Nx, Ny, Lx, Ly, 1, x_expantionRatio, y_expantionRatio, boundaryEdge, boundaryConditions, les, lesConvergenceThrehold, printer, convergenceThrehold),
		gama(gama), diffusion(diffusion), convection(convection), faceReconstruct(faceReconstruct), uniformVelocity(uniformVelocity){

	}

	virtual void solve(){
		config();
		CSR<T> diffusionTerm = diffusion->apply(gama, *(phi[0]), mesh);
		CellField<T> u;
		u.setUniformValue(mesh, uniformVelocity[0]);
		CellField<T> v;
		v.setUniformValue(mesh, uniformVelocity[1]);

		vector<T> u_f=faceReconstruct[0]->apply(u,mesh);
		vector<T> v_f = faceReconstruct[0]->apply(v, mesh);
		CSR<T> convectionTerm = convection->apply(u_f,v_f, *phi[0], mesh);
		CSR<T> eq = CSR<T>::minus(convectionTerm, diffusionTerm);
		les[0]->solve(eq, phi[0], mesh);
		phi[0]->assignBoundary(mesh);
		printer->print("out.dat", mesh, phi);

	}
};







