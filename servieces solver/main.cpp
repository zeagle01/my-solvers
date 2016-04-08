

#include <iostream>
#include "runner.h"
using namespace std;


int main(int argc, char* argv[]){

	cout << argv[0] << endl;
	string s="diffusionSolver.json";
	Runner runner(s);
	runner.run();

	return 0;
}