#include "case.h"
#include "mesh.h"
#include "laplaceSolver.h"

#include <boost/timer.hpp>
int main(){


	boost::timer t;  //定义一个计时类，开始计时

	/*
	Mesh mesh;
	mesh.buildMeshData("unstructure_2d.msh");
	*/
	string s = "diffusionSolver.json";
	Case* mycase = new GmeshLaplaceCase(new A_ConfigReader(s));
	//Solver* solver = new UnorthogonalCorrectionDiffusionSolver(mycase);
	Solver* solver = new DiffusionSolver(mycase);
	solver->solve();

	std::cout << "运行时间：" << t.elapsed() << std::endl;
	system("pause");
}