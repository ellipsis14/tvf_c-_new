#ifndef __INPUTOUTPUT_H
#define __INPUTOUTPUT_H

#include <dolfin.h>

#include "TvfGrad.h"
#include "Tvf4.h"
#include "converge.h"

#include <random>
#include <cmath>
#include <fstream>
#include <csignal>
#include <iostream>

// //NEWTON SOLVER FILES RE-WRITTEN:
// #include "./jdolfin/jNonlinearVariationalSolver.cpp"
// #include "./jdolfin/jNewtonSolver.cpp"

using namespace dolfin;

struct initParams{
	double TimeStepSize;
	std::array<double,2> gridSizeMicrons;
	size_t numNodes;
	
	double Alpha;
	double Beta;
	double Delta;
	double MaximumTimeSteps;
	double ScaledTimeStepSize;

};

// class dolfinVersionContents;//dummy class
// // {
// // 	friend std::ostream &operator<<(std::ostream &os, const dolfinVersionContents &io);
// // };

class inputoutput
{
	public:

	int parseInputLine(int argc, char* argv[]);

    std::string 
    // initOutputFiles(struct initParams &p, std::shared_ptr<jNonlinearVariationalSolver> solver);
    initOutputFiles(struct initParams &p);
	void 
	finalizeLogFile(double timingData[][5]);
	void 
	logTimingData(double t, double steps, double tstep, double dx0, double contraction1,
  			double timingData[][5]);

	public:
	bool isLocalComputer = false;
	bool isAugsbergCluster = false;
	bool isOpuntiaCluster = false;
	int  slurmArrayIndex = -1;

	private:
	// dolfinVersionContents dolfinVC;
	std::stringstream sstream;
	std::ofstream logFile;
	std::string fpath;

	int timingCounter = 0;
};
/*
std::ostream &operator<<(std::ostream &os, const dolfinVersionContents &io)//dummy comm, not needed
{
	//c++ lambda function (allows function def w/in function)
    auto outputMaps = [&](std::map <std::string, std::string> data)
	{
		std::map <std::string, std::string> :: iterator itr;
	    for (itr = data.begin(); itr != data.end(); ++itr)
	    {
	        os  <<  '\t' << itr->first 
	            <<  '\t' << itr->second << '\n';
	    }
	};
	os << "dolfin_version(): " << dolfin_version() << std::endl;
	os << "ufc_signature(): " << ufc_signature() << std::endl;
	os << "git_commit_hash(): " << git_commit_hash() << std::endl;

	os << "sizeof_la_index(): " << sizeof_la_index() << std::endl;
	os << "has_openmp(): " << has_openmp() << std::endl;
	os << "has_mpi(): " << has_mpi() << std::endl;
	os << "has_petsc(): " << has_petsc() << std::endl;
	os << "has_slepc(): " << has_slepc() << std::endl;
	os << "has_scotch(): " << has_scotch() << std::endl;
	os << "has_umfpack(): " << has_umfpack() << std::endl;
	os << "has_cholmod(): " << has_cholmod() << std::endl;
	os << "has_parmetis(): " << has_parmetis() << std::endl;
	os << "has_zlib(): " << has_zlib() << std::endl;
	os << "has_hdf5(): " << has_hdf5() << std::endl;
	os << "has_vtk(): " << has_vtk() << std::endl;

    // printing map gquiz1
    std::map <std::string, std::string> :: iterator itr;

	os << "linear_algebra_backends(): " << std::endl;    		outputMaps(linear_algebra_backends());
	os << "linear_solver_methods(): " << std::endl;    			outputMaps(linear_solver_methods());
	os << "lu_solver_methods(): " << std::endl;    				outputMaps(lu_solver_methods());
	os << "krylov_solver_methods(): " << std::endl;    			outputMaps(krylov_solver_methods());
	os << "krylov_solver_preconditioners(): " << std::endl;    	outputMaps(krylov_solver_preconditioners());

    return os;
}
*/
#endif
