
template<class T> class I_LinearEquationSolver {
public:
	virtual void solve(CSR<T> eq, CellField<T>* phi, Mesh<T> mesh)=0;
	virtual void iterate(CSR<T> eq, CellField<T>* phi, Mesh<T> mesh,int iterateTimes) = 0;
	virtual void set_convergedThrehold(T d) = 0;
};


template<class T> class  GS_Solver: public I_LinearEquationSolver<T>{

private:
	T convergedThrehold;

public:
	virtual void set_convergedThrehold(T d){
		convergedThrehold = d;
	}

	virtual void iterate(CSR<T> eq, CellField<T>* phi, Mesh<T> mesh,int iterateTimes) {

		
		for (int it = 0; it < iterateTimes; ++it){
			for (int c = 0; c < mesh.cellNum; c++) {
				double temp = 0;
				for (int i = mesh.IA[c] + 1; i < mesh.IA[c + 1]; i++) {
					temp += eq.A[i] * phi->inner[mesh.JA[i]];
				}
				phi->inner[c] = (eq.b[c] - temp) / eq.A[mesh.IA[c]];
			}
		}
	}

	virtual void solve(CSR<T> eq, CellField<T>* phi, Mesh<T> mesh){
		T error = 0;
		do{
			CellField<T> pre=*(phi);
			iterate(eq, phi, mesh, 1);
			error=VectorMath<T>::rootOfSquareSum(phi->inner, pre.inner);
			/*
			transform(phi->inner.cbegin(), phi->inner.cend(), pre.inner.cbegin(),pre.inner.begin(), minus<T>());
			transform(pre.inner.begin(), pre.inner.end(), pre.inner.begin(), pre.inner.begin(), multiplies<T>());
			error = accumulate(pre.inner.begin(), pre.inner.end(), (T)0);//have to force convert to T type unless the result is always 0
			//error = inner_product(pre.inner.cbegin(), pre.inner.cend(), pre.inner.cbegin(),
			//	1, plus<T>(), plus<T>());
			error = sqrt(error);
			*/
		} while (error>convergedThrehold);
	}

};




