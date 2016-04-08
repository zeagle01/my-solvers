#ifndef FUNCTION_H_H
#define FUNCTION_H_H


#include <iostream>
#include "serviceLocator.h"
#include <vector>
using namespace std;

class ExpansionFunction :public interface_t{
	vector<double> parameters;
public:
	virtual double operator()(double x, double y) = 0;
	void setParameters(vector<double> para){
		parameters = para;
	}
};


class UniformExpansion : public ExpansionFunction, public member_t<UniformExpansion>{
public:
	virtual double operator()(double x, double y){
		return x + y;
	}
};


class print{

public:
	void operator()(double e){
		cout << e << endl;
	}
};


#endif