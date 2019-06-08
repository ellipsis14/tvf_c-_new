#ifndef __TVF_H
#define __TVF_H

#include <dolfin.h>

#include "Tvf4.h"
#include "converge.h"
#include "TvfGrad.h"
#include "ErrorIPDG.h"

#include <random>
#include <cmath>
#include <fstream>
#include <csignal>
#include <iostream>

//NEWTON SOLVER FILES RE-WRITTEN:
#include "../rdolfin/jNonlinearVariationalSolver.h"
#include "../rdolfin/jNewtonSolver.h"

class PeriodicBC;

class TVF
{
public:
	TVF(struct initParams &);
			 
	~TVF();

	std::string printMeshData(void);

	void stepDiffusion();

	std::tuple<std::size_t, bool, double, double>
	step_u();
	void resetSolutions();
	void copySolutions();
	void changeTimeStep(double scale);
private:
	std::stringstream sstream;
	std::size_t numMeshNodes = 0;

	unsigned int **dof_lookup;
	double gridWidthMicrons;
	double gridHeightMicrons;


public:
	std::shared_ptr<dolfin::RectangleMesh> mesh;
	
	std::shared_ptr<PeriodicBC> PBC; 	

	std::shared_ptr<Tvf4::FunctionSpace> 	V_tvf;

	std::shared_ptr<dolfin::Function> 		u,u0;
	std::shared_ptr<Tvf4::ResidualForm> 	F;
	std::shared_ptr<Tvf4::JacobianForm> 	J;

	//GRADIENTS:
	std::shared_ptr<TvfGrad::FunctionSpace> 	G_tvf;
    std::shared_ptr<dolfin::Function>		    gradU;

	std::shared_ptr<TvfGrad::LinearForm>        L_gradU;
	std::shared_ptr<TvfGrad::BilinearForm>      a_gradU;

	std::shared_ptr<TvfGrad::FunctionSpace> 	G_tvf2;
	std::shared_ptr<dolfin::Function>		    gradU2;

    std::shared_ptr<TvfGrad::LinearForm>       L_gradU2;
	std::shared_ptr<TvfGrad::BilinearForm>     a_gradU2;



	std::shared_ptr<dolfin::NonlinearVariationalProblem> Problem_Tvf;
    std::shared_ptr<dolfin::jNonlinearVariationalSolver> Solver_Tvf;


  //extra buffers of prev. soln. and differences:
    std::shared_ptr<dolfin::Function> 		u_prev;
	std::shared_ptr<dolfin::Function> 		uDiff;


    std::shared_ptr<dolfin::LinearVariationalProblem> guLVP;
	std::shared_ptr<dolfin::LinearVariationalSolver> guLVS;


    std::shared_ptr<dolfin::Constant> newScaledTimeStepSize;
  // double currentTime, currentTimeStep;
  	double currentTimeStep;


	std::vector<std::array<double,2>> coords;
	std::vector<double> mesh_coords;
	std::vector<std::size_t> vertex_from_dof;	//dof_to_vertex_map(const FunctionSpace& space);
  	std::vector<dolfin::la_index> dof_from_vertex;//vertex_to_dof_map(const FunctionSpace& space);

};


#endif
