#ifndef RUNNER_H_H
#define RUNNER_H_H

#include <string>
#include "solver.h"
#include "configReader.h"
using namespace std;
class Runner{

public :

	Solver* solver;
	
	ConfigReader configReader;
	Runner(string s) :configReader(s){
		//solver = configReader.readSolverClass("solver");
		solver->configReader= configReader;
	}

	void run(){
		solver->config();
		solver->solve();
	}

};



#endif