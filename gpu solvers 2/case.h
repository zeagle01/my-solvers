

#ifndef CASE_H_H
#define CASE_H_H


#include "mesh.h"
#include "printer.h"
#include "serviceLocator.h"
#include "configReader.h"
#include "linearEquationSolver.h"
class Case{
public:
	Mesh mesh;
	int Nx;
	int Ny;
	double Lx;
	double Ly;
	vector < vector<int>> boundaryEdge;

	int variableNum;
	double convergenceThrehold;
	vector<double> lesConvergenceThrehold;


	vector<A_Printer*> printers;
	vector<ExpansionFunction*> expantionFunctions;

	A_ConfigReader* configReader;
	vector<CellField> phi;

	vector<I_LinearEquationSolver*> les;

	vector<vector<vector<double>>> boundaryConditions;

	void buildVariable(){
		for (int i = 0; i < variableNum; i++){
			phi.push_back( CellField(mesh));
			for (int j = 0; j < boundaryEdge.size(); j++){
				phi[i].boundaryCondition.push_back(boundaryConditions[j][i]);
			}

		}
	}

	void config(){
		readConfig();
		buildMesh();
		buildVariable();

	}
	

	void readConfig(){
		Nx = configReader->readInt("Nx");
		Ny = configReader->readInt("Ny");
		Lx = configReader->readDouble("Lx");
		Ly = configReader->readDouble("Ly");
		variableNum = configReader->readInt("variable_num");
		expantionFunctions = configReader->readExpantionFunction("expansion_functions");
		printers = configReader->readPrinter("printers");	
		boundaryEdge = configReader->read_2D_array<int>("boundarySet");
		boundaryConditions = configReader->read_3D_array<double>("boundary_conditions");
		les = configReader->readLinearEquationSolver("linear_equation_solvers");

		/*
		boost::property_tree::ptree les= pt.get_child("linear_equation_solvers");
		*/
	}

	

	Case(A_ConfigReader* configReader){
		this->configReader = configReader;
	}

	void buildMesh(){
		mesh.x_expantionRatio = expantionFunctions[0];
		mesh.y_expantionRatio = expantionFunctions[1];
		mesh.Nx = Nx;
		mesh.Ny = Ny;
		mesh.Lx = Lx;
		mesh.Ly = Ly;
		mesh.boundaryEdge = boundaryEdge;
		mesh.generateTopology();
	}

};


#endif