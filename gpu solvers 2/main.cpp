

#include "case.h"
#include "Solver.h"
#include "configReader.h"
#include "laplaceSolver.h"
//#include <cuda_runtime.h>
//#include "device_launch_parameters.h"

#include <boost/timer.hpp>
//#define BOOST_DATE_TIME_SOURCE  



int main(){

	boost::timer t;  //定义一个计时类，开始计时



	//boost::posix_time::ptime time_now, time_now1;
	//boost::posix_time::millisec_posix_time_system_config::time_duration_type time_elapse;

	// 这里为微秒为单位;这里可以将microsec_clock替换成second_clock以秒为单位;  
	//time_now = boost::posix_time::microsec_clock::universal_time();


	string s = "solverConfig.json";
	Case* mycase = new Case(new A_ConfigReader(s));
	Solver* solver = new DiffusionSolver(mycase);
	solver->solve();

	//std::cout << "可度量的最大时间:" << t.elapsed_max() / 3600 << "h" << std::endl;

	//std::cout << "可度量的最大时间:" << t.elapsed_min() << "s" << std::endl;

	std::cout << "使用时间为：" << t.elapsed() << std::endl;

	/*
	// sleep 100毫秒;  
	boost::this_thread::sleep(boost::posix_time::millisec(100));
	*/

	//time_now1 = boost::posix_time::microsec_clock::universal_time();

	//time_elapse = time_now1 - time_now;

	// 类似GetTickCount，只是这边得到的是2个时间的ticket值的差，以微秒为单位;  
	//int ticks = time_elapse.ticks();
	
	// 得到两个时间间隔的秒数;  
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
	std::cout <<"运行时间为" << elapsedTime << std::endl;
	cudaEventDestroy(begin);
	cudaEventDestroy(stop);
	*/
	system("pause");
	return 0;
}