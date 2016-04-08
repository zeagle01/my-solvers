#include "case.h"
#include "mesh.h"
#include "laplaceSolver.h"

#include <boost/timer.hpp>
int main(){


	boost::timer t;  //����һ����ʱ�࣬��ʼ��ʱ

	/*
	Mesh mesh;
	mesh.buildMeshData("unstructure_2d.msh");
	*/
	string s = "diffusionSolver.json";
	Case* mycase = new GmeshLaplaceCase(new A_ConfigReader(s));
	//Solver* solver = new UnorthogonalCorrectionDiffusionSolver(mycase);
	Solver* solver = new DiffusionSolver(mycase);
	solver->solve();

	std::cout << "����ʱ�䣺" << t.elapsed() << std::endl;
	system("pause");
}