#ifndef SOLVER_H_H
#define SOLVER_H_H

class Solver{
public:
	Case* cas;
	Solver(Case* cas){
		this->cas = cas;
	}
	virtual void solve() = 0;
};
#endif