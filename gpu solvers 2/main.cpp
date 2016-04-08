

#include "case.h"
#include "Solver.h"
#include "configReader.h"
#include "laplaceSolver.h"
//#include <cuda_runtime.h>
//#include "device_launch_parameters.h"

#include <boost/timer.hpp>
//#define BOOST_DATE_TIME_SOURCE  



int main(){

	boost::timer t;  //����һ����ʱ�࣬��ʼ��ʱ



	//boost::posix_time::ptime time_now, time_now1;
	//boost::posix_time::millisec_posix_time_system_config::time_duration_type time_elapse;

	// ����Ϊ΢��Ϊ��λ;������Խ�microsec_clock�滻��second_clock����Ϊ��λ;  
	//time_now = boost::posix_time::microsec_clock::universal_time();


	string s = "solverConfig.json";
	Case* mycase = new Case(new A_ConfigReader(s));
	Solver* solver = new DiffusionSolver(mycase);
	solver->solve();

	//std::cout << "�ɶ��������ʱ��:" << t.elapsed_max() / 3600 << "h" << std::endl;

	//std::cout << "�ɶ��������ʱ��:" << t.elapsed_min() << "s" << std::endl;

	std::cout << "ʹ��ʱ��Ϊ��" << t.elapsed() << std::endl;

	/*
	// sleep 100����;  
	boost::this_thread::sleep(boost::posix_time::millisec(100));
	*/

	//time_now1 = boost::posix_time::microsec_clock::universal_time();

	//time_elapse = time_now1 - time_now;

	// ����GetTickCount��ֻ����ߵõ�����2��ʱ���ticketֵ�Ĳ��΢��Ϊ��λ;  
	//int ticks = time_elapse.ticks();
	
	// �õ�����ʱ����������;  
	//int sec = time_elapse.total_seconds();
	//cout << sec <<" "<<ticks<< endl;

	/* //cudaTimer
	cudaEvent_t begin, stop;
	cudaEventCreate(&begin);
	cudaEventCreate(&stop);
	cudaEventRecord(begin, 0);

	string s = "solverConfig.json";
	Case* mycase = new Case(new A_ConfigReader(s));
	Solver* solver = new DiffusionSolver(mycase);
	solver->solve();

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, begin, stop);
	std::cout <<"����ʱ��Ϊ" << elapsedTime << std::endl;
	cudaEventDestroy(begin);
	cudaEventDestroy(stop);
	*/
	system("pause");
	return 0;
}