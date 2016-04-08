#ifndef CONFIGREADER_H_H
#define CONFIGREADER_H_H


#include <boost/foreach.hpp>
#include <boost\property_tree\json_parser.hpp>
#include <string>
#include "servicesLocator.h"
#include "solver.h"
#include "mesh.h"
#include "printer.h"
#include "laplace.h"
using namespace std;

class ConfigReader{
	string configFile;
	boost::property_tree::ptree pt;
	servicelocator_t services;
public:
	ConfigReader(string configFile){
		this->configFile = configFile;
		boost::property_tree::json_parser::read_json(configFile, pt);
		buildService();
	}
	//register the class
	void buildService(){

	
		services.register_class<GmeshTriangleMesh>("Gmesh_triangle_mesh");
		services.register_class<VersatileGmeshMesh>("versatile_Gmesh_mesh");	
		
		services.register_class<Printer>("tecplot_without_boundary");
		services.register_class<UnstructurePrinter>("tecplot_unstructure");
		services.register_class<GS_Solver>("GS_iteration");
		//services.register_class<OpenMP_Jacobi_Solver>("OpenMP_Jacobi_iteration");
		//services.register_class<GPU_Jacobi_Solver>("GPU_Jacobi_iteration");
		//services.register_class<GPU_Jacobi_Solver_vector>("GPU_Jacobi_Solver_vector");
		services.register_class<Laplace>("Laplace");
		services.register_class<UnorthogonalCorrectionLaplace>("unstructure_Laplace");
		//services.register_class<GaussCellGradient>("GS_gradient");
		//services.register_class<FaceReconstruct_BEB>("linear_average_beb");
		//services.register_class<LinearFaceReconstruct_BEO>("linear_average_beo");
		

	}

	int readInt(string s){
		return pt.get_child(s).get_value<int>();
	}
	double readDouble(string s){
		return pt.get_child(s).get_value<double>();
	}

	string readString(string s){
		Solver* s;
		return pt.get_child(s).get_value<string>();
		
	}

	/*
	Solver* readSolverClass(string s){
		//Solver* r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_value<string>();
		//r = services.get_single_instance<Solver>(type);
		return 0;
	}
	*/

	
	
	Mesh* readMeshClass(string s){
		Mesh* r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_child("class").get_value<string>();
		string meshFile = child.get_child("mesh_file").get_value<string>();
		r = services.get_single_instance<Mesh>(type);
		r->meshFile = meshFile;
		return r;
	}

	/*
	I_GaussCellGradient* readCellGradient(string s){
		I_GaussCellGradient *r;
		boost::property_tree::ptree child = pt.get_child(s);
		string type = child.get_value<string>();
		r = services.get_single_instance<I_GaussCellGradient>(type);
		return r;
	}

	vector<I_FaceReconstruct*> readFaceReconstructor(string s){
		vector<I_FaceReconstruct*> r;
		boost::property_tree::ptree child = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child) {
			string type = vt.second.get_child("type").get_value<string>();
			r.push_back(services.get_single_instance<I_FaceReconstruct>(type));
		}
		return r;
	}
	*/

	vector<A_Printer*> readPrinter(string s){
		vector<A_Printer*> r;
		boost::property_tree::ptree child = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child) {
			string type = vt.second.get_child("type").get_value<string>();
			string zone_type = vt.second.get_child("zone_type").get_value<string>();

			r.push_back(services.get_single_instance<A_Printer>(type));
			r.front()->zoneType = zone_type;
		}
		return r;
	}
	


	vector<I_LinearEquationSolver*> readLinearEquationSolver(string s){
		vector<I_LinearEquationSolver*> r;
		boost::property_tree::ptree child = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child) {
			string type = vt.second.get_child("type").get_value<string>();
			r.push_back(services.get_single_instance<I_LinearEquationSolver>(type));
			r.back()->converge_threhold = vt.second.get_child("converge_threhold").get_value<double>();
			r.back()->max_step = vt.second.get_child("max_step").get_value<int>();
			r.back()->check_step = vt.second.get_child("check_step").get_value<int >();
		}
		return r;
	}




	I_Laplace* readLaplace(string s){
		I_Laplace* r;
		//boost::property_tree::ptree child = ;
		string type = pt.get_child(s).get_value<string>();
		r = services.get_single_instance<I_Laplace>(type);
		return r;
	}


	/*
	template<class T> vector<vector<vector<T>>> read_3D_array(string arrayName){
	boost::property_tree::ptree child = pt.get_child(arrayName);
	vector<vector<vector<T>>> r;
	BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child) {
	r.push_back(vector<vector<T>>());
	BOOST_FOREACH(boost::property_tree::ptree::value_type &vt1, vt.second) {
	r.back().push_back(vector<T>());
	BOOST_FOREACH(boost::property_tree::ptree::value_type &vt2, vt1.second) {
	r.back().back().push_back(vt2.second.get_value<T>());
	}
	}
	}
	return r;
	}

	template<class T> vector<vector<T>>read_2D_array(string arrayName){
		boost::property_tree::ptree boundarySet = pt.get_child(arrayName);
		vector<vector<T>> r;
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, boundarySet) {
			r.push_back(vector<T>());
			BOOST_FOREACH(boost::property_tree::ptree::value_type &vt1, vt.second) {
				r.back().push_back(vt1.second.get_value<T>());
			}
		}
		return r;
	}
	*/
	
	template<class T> vector<T> read_1D_array(string s){

		vector<T> r;
		boost::property_tree::ptree child_array = pt.get_child(s);
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child_array) {
			r.push_back(vt.second.get_value<T>());
		}
		return r;
	}
	
	/*
	template<class T> vector<T> read_1D_child_array(boost::property_tree::ptree child_array){
		//vector<T> r = child_array.get_value < vector<T>>();
		vector<T> r;
		BOOST_FOREACH(boost::property_tree::ptree::value_type &vt, child_array) {
			r.push_back(vt.second.get_value<T>());
		}
		return r;
	}


	*/
};







#endif