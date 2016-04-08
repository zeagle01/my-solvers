#include<iostream>
#include "mesh.h"
#include <boost\property_tree\json_parser.hpp>
using namespace std;







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
	boost::property_tree::ptree pt;
	boost::property_tree::json_parser::read_json("configData.json", pt);
	cout << pt.get_child("Nx").get_value<int>() << endl;
	cout << pt.get_child("Ny").get_value<int>() << endl;
	//run_pressure_correction<double>();

	Mesh<double> mesh;
	mesh.ReadFromGshFile("unstructure_2d.msh");



	//int index_file;
	//cg_open("grid_c.cgns", CG_MODE_WRITE, &index_file);

	cout << "OK!" << endl;


	system("pause");

	return 0;
}


