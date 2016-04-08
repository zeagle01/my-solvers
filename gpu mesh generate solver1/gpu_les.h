//gpu_les

template<class T> class I_GPU_LinearEquationSolver {
public :
	virtual void solve(T* A, T* IA, T* JA, T* X) = 0;
	virtual void iterate(T* A, T* IA, T* JA, T* X, int iterateTimes) = 0;
	virtual void set_convergedThrehold(T d) = 0;
};




template<class T> class GPU_GS_Solver :public I_GPU_LinearEquationSolver<T> {
public :
	virtual void solve(T* A, T* IA, T* JA, T* X);
	virtual void iterate(T* A, T* IA, T* JA, T* X, int iterateTimes){}
	virtual void set_convergedThrehold(T d) {}
};