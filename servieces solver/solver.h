#ifndef SOLVER_H_H
#define SOLVER_H_H


#include "mesh.h"
#include "servicesLocator.h"
#include "dataStructure.h"
#include "printer.h"
#include "linearEquationSolver.h"
#include "configReader.h"
#include "configReader.h"
class Solver :public interface_t{


public:
	Mesh* mesh;
	ConfigReader configReader;


	int variableNum;
	double convergenceThrehold;
	vector<double> lesConvergenceThrehold;


	vector<A_Printer*> printers;

	vector<CellField> phi;

	vector<I_LinearEquationSolver*> les;

	vector<double> boundaryConditions;

	void buildVariable(){
		for (int i = 0; i < variableNum; i++){
			phi.push_back(CellField(mesh));
			for (int j = 0; j < mesh->boundaryPatchNum; j++){
				phi[i].boundaryCondition.push_back(boundaryConditions[(j*variableNum + i) * 3 + 0]);
				phi[i].boundaryCondition.push_back(boundaryConditions[(j*variableNum + i) * 3 + 1]);
				phi[i].boundaryCondition.push_back(boundaryConditions[(j*variableNum + i) * 3 + 2]);
			}
		}
	}

	virtual void config(){
		variableNum = configReader.readInt("variable_num");
		printers = configReader.readPrinter("printers");

		boundaryConditions = configReader.read_1D_array<double>("boundary_conditions");

		les = configReader.readLinearEquationSolver("linear_equation_solvers");
		mesh = configReader.readMeshClass("mesh");

		mesh->buildMeshData();
		buildVariable();
	}

	virtual void solve() = 0;
};
#endif