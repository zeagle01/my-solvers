
//

template <class T> class GPU_DiffusionSolver :public DiffusionSolver<T>{
public:
	int* IA;
	int * JA;
	T* X;
	T* A;


	int* d_IA;
	int * d_JA;
	T* d_X;
	T* d_A;
	GPU_DiffusionSolver(int Nx, int Ny, T Lx, T Ly, T gama,
		ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio,
		vector < vector<int> > boundaryEdge, vector<vector<vector<T>>>* boundaryConditions,
		I_Laplace<T>* diffusion,
		vector<I_LinearEquationSolver<T>*> les,vector<T> lesConvergenceThrehold,
		A_Printer<T>* printer,
		T convergenceThrehold
		) :DiffusionSolver(Nx, Ny, Lx, Ly, gama,
		x_expantionRatio, y_expantionRatio,
		boundaryEdge, boundaryConditions,
		diffusion,
		les, lesConvergenceThrehold,
		printer,
		convergenceThrehold
		){}


	void callGPU();
	void copyH2D();
	void copyD2H();
	virtual void solve(){
		config();
		CSR<T> eq = diffusion->apply(gama, *(phi[0]), mesh);
		int n = mesh.IA.size();
		IA = new int[n];
		for (int i = 0; i < n; i++){
			IA[i] = mesh.IA[i];
		}
		n = mesh.JA.size();
		JA = new int[n];
		A = new T[n];
		X = new T[mesh.cellNum];
		for (int i = 0; i < n; i++){
			JA[i] = mesh.JA[i];
			A[i] = eq.A[i];
		}
		copyH2D();
		callGPU();
		copyD2H();
		
	}


};