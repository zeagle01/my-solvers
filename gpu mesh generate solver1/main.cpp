#include<iostream>
#include "pressureVelocityCoupleSolver.h"
#include "gpu_solver.h"
#include <boost\property_tree\json_parser.hpp>
#include "VelocityCorrectionIBM.h"

using namespace std;

//cases

template<class T> struct run_diffusion{
	run_diffusion(){
		DiffusionSolver<T> solver
			(200, 200, 1, 1, 1.0,
			new UniformExpansion<T>(), new UniformExpansion<T>(),
			vector<vector<int>>{{ 1, 2, 3 }, { 0 }, },
			new vector<vector<vector<T>>>{
				{ { 0, 1, 0 } },
				{ { 0, 1, 1 } }
		},
			new Laplace<T>(),
			vector < I_LinearEquationSolver < T >*>{new GS_Solver<T>()}, vector<T>{1e-6},
			new Printer<T>(),
			1e-6
			);

		solver.solve();
	}
};

template<class T> struct run_gpu_diffusion{
	run_gpu_diffusion(){
		GPU_DiffusionSolver<T> solver
			(200, 200, 1, 1, 1.0,
			new UniformExpansion<T>(), new UniformExpansion<T>(),
			vector<vector<int>>{{ 1, 2, 3 }, { 0 }, },
			new vector<vector<vector<T>>>{
				{ { 0, 1, 0 } },
				{ { 0, 1, 1 } }
		},
			new Laplace<T>(),
			vector < I_LinearEquationSolver < T >*>{new GS_Solver<T>()}, vector<T>{(T)1e-6},
			new Printer<T>(),
			1e-6
			);

		solver.solve();
	}
};


template<class T> struct run_convection_diffusion{
	run_convection_diffusion(){
		ConvectionDiffusionSolver<T> solver
			(100, 100, 1, 1, 0,											//Nx,Ny,lx,Ly,gama
			new UniformExpansion<T>(), new UniformExpansion<T>(),		//x_expansion,y_expansion											
			vector<vector<int>>{{ 1,2}, { 0 }, { 3 } },					//boundaryEdge			
			new vector<vector<vector<T>>>{								//boundaryConditions
				{ 
					{ 1, 0, 0 }

				},
				{
					{ 0, 1, 1 } 
	
				},
				{
					{ 0, 1, 0 }
				}
		},
			vector<T>{1, 1},											//velocity field
			new Laplace<T>(), new UpwindConvection<T>(),
			vector<I_FaceReconstruct<T>*>{new FaceResonctruct_BEB<T>()},
			vector < I_LinearEquationSolver < T >*>{new GS_Solver<T>()}, vector<T>{1e-6},
			new Printer<T>(),
			1e-6
			);

		solver.solve();
	}
};

template<class T> struct run_pressure_correction{
	run_pressure_correction(){
		PressureCorrectionSolver<T> solver
			(100, 100, 1, 1, 0.01,1e6,20000,100,									//Nx,Ny,lx,Ly,gama,dt,maxStep,checkStep
			new UniformExpansion<T>(), new UniformExpansion<T>(),		//x_expansion,y_expansion											
			vector<vector<int>>{{ 1 }, { 0, 2 }, { 3 } },					//boundaryEdge			
			new vector<vector<vector<T>>>{								//boundaryConditions
				{
					{ 0, 1, 1 },
					{ 0, 1, 0 },
					{ 1, 0, 0 }

				},
				{
					{ 0, 1, 0 },
					{ 0, 1, 0 },
					{ 1, 0, 0 }

				},
				{
					{ 0, 1, 0 },
					{ 0, 1, 0 },
					{ 1, 0, 0 }

				}
		},
			new Laplace<T>(), new UpwindConvection<T>(), new ImplictEuler<T>(),
			vector<I_FaceReconstruct<T>*>{new FaceReconstruct_BEB<T>()},
			new GaussCellGradient<T>(), new PlainFaceGradient<T>(), new B_RhieChow<T>(),
			vector < I_LinearEquationSolver < T >*>{new GS_Solver<T>()}, vector<T>{1e-6},
			new Printer<T>(),
			1e-6
			);

		solver.solve();
		//cout << endl;
		//copy(solver.phi[0]->inner.begin(), solver.phi[0]->inner.end(), ostream_iterator<T>(cout, "\n"));
	}
};


//end cases


#include <boost/foreach.hpp>

int main(){

	/*
	boost::property_tree::ptree pt;
	boost::property_tree::json_parser::read_json("jsonData.json", pt);

	boost::property_tree::ptree child_linktype = pt.get_child("rate.linktype");

	BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child_linktype) {
		cout << vt.second.get_value<double>() << "\t";
	}
	cout << endl;
	*/
	//DiffusionSolver<double> solver("solverConfig.json");
	boost::property_tree::ptree pt;
	boost::property_tree::json_parser::read_json("solverConfig.json", pt);
	cout << pt.get_child("Nx").get_value<int>() << endl;
	cout << pt.get_child("Ny").get_value<int>() << endl;
	//run_pressure_correction<double>();
	I_GPU_LinearEquationSolver<float> *gpu_les = new GPU_GS_Solver<float>();
	gpu_les->solve(NULL,NULL,NULL,NULL);
	//run_gpu_diffusion<float>();
	cout << "OK!" << endl;

	system("pause");

	return 0;
}


