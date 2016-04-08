template<class T> class VelocityCorrectionIBM{
	PressureCorrectionSolver<T> solver;
	VelocityCorrectionIBM(
		int Nx, int Ny, T Lx, T Ly, T gama, T dt, int maxStep, int checkStep,
		ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio,
		vector < vector<int> > boundaryEdge, vector<vector<vector<T>>>* boundaryConditions,
		I_Laplace<T>* diffusion, I_Convection<T>* convection, I_Temporal<T>* temporal,
		vector<I_FaceReconstruct<T>*> faceReconstruct, I_GaussCellGradient<T>* cGrad, I_FaceGradient<T>* fGrad, B_RhieChow<T>* rc,
		vector<I_LinearEquationSolver<T>*> les, vector<T> lesConvergenceThrehold,
		A_Printer<T>* printer,
		T convergenceThrehold
		) :solver( Nx,  Ny,  Lx,  Ly,  gama,  dt,  maxStep,  checkStep,
		 x_expantionRatio,  y_expantionRatio,
		 boundaryEdge,  boundaryConditions,
		 diffusion,  convection,  temporal,
		 faceReconstruct, cGrad,  fGrad,  rc,
		 les, lesConvergenceThrehold,
		 printer,
		 convergenceThrehold){
	}

	vector<int> innerArea;
	

};