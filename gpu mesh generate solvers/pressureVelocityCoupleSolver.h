#include "solver.h"
#include "temporalTerm.h"
#include "gradient.h"
#include "divergence.h"
template<class T> class PressureVelocityCoupleSolver :public Solver<T>{
public:
	I_Laplace<T>* diffusion;
	I_Convection<T>* convection;
	vector<I_FaceReconstruct<T>*> faceReconstruct;
	T gama;
	T dt;
	int step;
	int maxStep;
	int checkStep;
	T error;
	PressureVelocityCoupleSolver(int Nx, int Ny, T Lx, T Ly, T gama, T dt, int maxStep,int checkStep,
		ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio,
		vector < vector<int> > boundaryEdge, vector<vector<vector<T>>>* boundaryConditions,
		I_Laplace<T>* diffusion, I_Convection<T>* convection,
		vector<I_FaceReconstruct<T>*> faceReconstruct,
		vector<I_LinearEquationSolver<T>*> les, vector<T> lesConvergenceThrehold,
		A_Printer<T>* printer,
		T convergenceThrehold
		)

		:Solver(Nx, Ny, Lx, Ly, 3,x_expantionRatio, y_expantionRatio, boundaryEdge, boundaryConditions, les, lesConvergenceThrehold, printer, convergenceThrehold),
		gama(gama), dt(dt), maxStep(maxStep), checkStep(checkStep),diffusion(diffusion), convection(convection), faceReconstruct(faceReconstruct){

	}	
};




template<class T> class PressureCorrectionSolver :public  PressureVelocityCoupleSolver<T>{
public:



	PressureCorrectionSolver(int Nx, int Ny, T Lx, T Ly, T gama, T dt, int maxStep,int checkStep,
		ExpansionFunction<T>* x_expantionRatio, ExpansionFunction<T>* y_expantionRatio,
		vector < vector<int> > boundaryEdge, vector<vector<vector<T>>>* boundaryConditions,
		I_Laplace<T>* diffusion, I_Convection<T>* convection, I_Temporal<T>* temporal,
		vector<I_FaceReconstruct<T>*> faceReconstruct, I_GaussCellGradient<T>* cGrad, I_FaceGradient<T>* fGrad, B_RhieChow<T>* rc,
		vector<I_LinearEquationSolver<T>*> les, vector<T> lesConvergenceThrehold,
		A_Printer<T>* printer,
		T convergenceThrehold
		)

		:PressureVelocityCoupleSolver(Nx, Ny, Lx, Ly, gama, dt,  maxStep, checkStep,
		x_expantionRatio,y_expantionRatio,
		boundaryEdge,  boundaryConditions,
		diffusion, convection,
		faceReconstruct,
		les, lesConvergenceThrehold,
		printer,
		convergenceThrehold
		), temporal(temporal), cGrad(cGrad), fGrad(fGrad),rc(rc){}

	I_Temporal<T>* temporal;
	I_GaussCellGradient<T>* cGrad;
	I_FaceGradient<T>* fGrad;
	B_RhieChow<T>* rc;

	ExplicitDivergence<T> div;
	virtual void solve(){
		config();


		CSR<T> u_diffusion = diffusion->apply(gama, *(phi[0]), mesh);
		CSR<T> v_diffusion = diffusion->apply(gama, *(phi[1]), mesh);
		CSR<T>  p_eq= diffusion->apply((T)1, *(phi[2]), mesh);
		CellField<T> p_prime = *(phi[2]);
		//vector<T> u_f()
		do{
			vector<CellField<T>> phi0(3, CellField<T>());
			phi0[0] = *(phi[0]);
			phi0[1] = *(phi[1]);
			phi0[2] = *(phi[2]);

			vector<T> u_f = faceReconstruct[0]->apply(*(phi[0]), mesh);
			vector<T> v_f = faceReconstruct[0]->apply(*(phi[1]), mesh);
			vector<T> p_f = faceReconstruct[0]->apply(*(phi[2]), mesh);

			CSR<T> u_convection = convection->apply(u_f, v_f, *(phi[0]), mesh);
			vector<CellField<T>> u(1, *(phi[0]));
			CSR<T> u_timeTerm = temporal->apply(u, dt, mesh);
			CSR<T> u_momentum = CSR<T>::minus(u_convection, u_diffusion);
			u_momentum = CSR<T>::plus(u_momentum, u_timeTerm);

			CSR<T> v_convection = convection->apply(u_f, v_f, *(phi[1]), mesh);
			vector<CellField<T>> v(1, *(phi[1]));
			CSR<T> v_timeTerm = temporal->apply(v, dt, mesh);
			CSR<T> v_momentum = CSR<T>::minus(v_convection, v_diffusion);
			v_momentum = CSR<T>::plus(v_momentum, v_timeTerm);


			vector<CellField<T>> p_cellGrad = cGrad->apply(p_f,mesh);
			for (int c = 0; c<mesh.cellNum; c++){
				u_momentum.b[c] -= mesh.C_volumn[c] * p_cellGrad[0].inner[c];
				v_momentum.b[c] -= mesh.C_volumn[c] * p_cellGrad[1].inner[c];
			}

			les[0]->iterate(u_momentum, phi[0], mesh,1);
			les[0]->iterate(v_momentum, phi[1], mesh,1);
			//les[0]->solve(u_momentum, phi[0],mesh);
			//les[0]->solve(v_momentum, phi[1], mesh);
			phi[0]->assignBoundary(mesh);
			phi[1]->assignBoundary(mesh);



			vector<vector<T>> fgrad = fGrad->apply(*(phi[2]), mesh);

			vector<vector<T>> D_f = rc->calculateD(u_momentum, v_momentum, mesh);
			vector<CellField<T>> V(2, CellField<T>());
			V[0] = *(phi[0]); V[1] = *(phi[1]);
			vector<vector<T>> v_Rc = rc->reconstructVelocity(D_f, p_cellGrad, fgrad, V, mesh);
			vector<T> mdot = div.apply(v_Rc, mesh);

			for (int c = 0; c<mesh.cellNum; c++){
				p_eq.b[c] = mdot[c] * mesh.C_volumn[c] / dt;
			}

			
			 les[0]->iterate(p_eq, &p_prime, mesh, 1);
			//les[0]->solve(p_eq, &p_prime, mesh);
			setRelativePressure(p_prime);
			p_prime.assignBoundary(mesh);

			vector<T> p_prime_f = faceReconstruct[0]->apply(p_prime, mesh);
			vector<CellField<T>> p_prime_cgrad = cGrad->apply(p_prime_f, mesh);
			for (int c = 0; c < mesh.cellNum; c++) {
				phi[0]->inner[c] = phi[0]->inner[c] - p_prime_cgrad[0].inner[c] * dt;
				phi[1]->inner[c] = phi[1]->inner[c] - p_prime_cgrad[1].inner[c] * dt;
				phi[2]->inner[c] += p_prime.inner[c];
			}
			phi[0]->assignBoundary(mesh);
			phi[1]->assignBoundary(mesh);
			phi[2]->assignBoundary(mesh);

			if (step%checkStep == 0){
				error = VectorMath<T>::rootOfSquareSum(phi[0]->inner, phi0[0].inner);
				error += VectorMath<T>::rootOfSquareSum(phi[1]->inner, phi0[1].inner);
				error += VectorMath<T>::rootOfSquareSum(phi[2]->inner, phi0[2].inner);
				printer->print("out.dat", mesh, phi);
				cout << step <<"	"<< error << endl;
			}
			step++;
		} while (step<maxStep&&error>convergenceThrehold);

		printer->print("out.dat", mesh, phi);
	}

	void setRelativePressure(CellField<T> & phi){
		for (int c = 0; c < mesh.cellNum; c++){
			phi.inner[c] = phi.inner[c] - phi.inner[0];
		}
	}
};