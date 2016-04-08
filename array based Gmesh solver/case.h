

#ifndef CASE_H_H
#define CASE_H_H


#include "mesh.h"
#include "printer.h"
#include "serviceLocator.h"
#include "configReader.h"
#include "linearEquationSolver.h"
#include "dataStructure.h"
class Case{
public:
	Mesh* mesh;
	A_ConfigReader* configReader;


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
	

	void config(){
		readConfig();
		buildMesh();
		buildVariable();

	}

	
	virtual void readConfig() = 0;
	

	

	Case(A_ConfigReader* configReader){
		this->configReader = configReader;
	}

	virtual void buildMesh() = 0; 

};




class GmeshLaplaceCase:public Case{
public:
	string meshFile;

	GmeshLaplaceCase(A_ConfigReader* configReader) :Case(configReader){
	}


	virtual void buildMesh(){
		mesh->buildMeshData();
	}

	virtual void readConfig(){
		variableNum = configReader->readInt("variable_num");
		printers = configReader->readPrinter("printers");

		boundaryConditions = configReader->read_1D_array<double>("boundary_conditions");

		les = configReader->readLinearEquationSolver("linear_equation_solvers");
		mesh = configReader->readMeshClass("mesh");

	}
};





#endif