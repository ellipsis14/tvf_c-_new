// Copyright (C) 2005-2008 Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Anders Logg 2005-2009
// Modified by Martin Alnes 2008
// Modified by Johan Hake 2010
//
// First added:  2005-10-23
// Last changed: 2014-05-27

#include <cmath>
#include <string>

#include <dolfin/common/constants.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/la/DefaultFactory.h>
#include <dolfin/la/LinearSolver.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/LUSolver.h>
#include <dolfin/la/KrylovSolver.h>
#include <dolfin/log/log.h>
#include <dolfin/common/MPI.h>
#include "NonlinearProblem.h"
#include "jNewtonSolver.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Parameters jNewtonSolver::default_parameters()
{
  Parameters p("newton_solver");

  p.add("linear_solver",           "default");
  p.add("preconditioner",          "default");
  p.add("maximum_iterations",      50);
  p.add("relative_tolerance",      1e-9);
  p.add("absolute_tolerance",      1e-10);
  p.add("convergence_criterion",   "residual");
  p.add("report",                  true);
  p.add("error_on_nonconvergence", true);
  p.add<double>("relaxation_parameter");

  //p.add("reuse_preconditioner", false);

  p.add(LUSolver::default_parameters());
  p.add(KrylovSolver::default_parameters());

  return p;
}
//-----------------------------------------------------------------------------
jNewtonSolver::jNewtonSolver(MPI_Comm comm,
                           std::shared_ptr<GenericLinearSolver> solver,
                           GenericLinearAlgebraFactory& factory)
  : Variable("Newton solver", "unnamed"), _newton_iteration(0),
    _krylov_iterations(0), _relaxation_parameter(1.0), _residual(0.0),
    _residual0(0.0), _solver(solver), _matA(factory.create_matrix(comm)),
    _matP(factory.create_matrix(comm)), _dx(factory.create_vector(comm)),
    _b(factory.create_vector(comm)), _mpi_comm(comm)
{
  // Set default parameters
  parameters = default_parameters();

  // Override linear solver type if solver passed in
  if (solver)
  {
    parameters["linear_solver"] = "user_defined";
    parameters["preconditioner"] = "user_defined";
  }
}
//-----------------------------------------------------------------------------
jNewtonSolver::jNewtonSolver(MPI_Comm comm)
  : jNewtonSolver(comm, nullptr, DefaultFactory::factory()) { }
//-----------------------------------------------------------------------------
jNewtonSolver::~jNewtonSolver()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
//jMODS
std::tuple<std::size_t, bool, double, double>
// std::pair<std::size_t, bool>
jNewtonSolver::solve(NonlinearProblem& nonlinear_problem,
                    GenericVector& x)
{

  dolfin_assert(_matA);
  dolfin_assert(_b);
  dolfin_assert(_dx);

  // Extract parameters
  const std::string convergence_criterion = parameters["convergence_criterion"];
  const std::size_t maxiter = parameters["maximum_iterations"];
  if (parameters["relaxation_parameter"].is_set())
    set_relaxation_parameter(parameters["relaxation_parameter"]);

  // Create linear solver if not already created
  const std::string solver_type = parameters["linear_solver"];
  const std::string pc_type = parameters["preconditioner"];
  if (!_solver)
  {
    _solver = std::make_shared<LinearSolver>(x.mpi_comm(), solver_type, pc_type);
  }
  dolfin_assert(_solver);

  // Set parameters for linear solver
  _solver->update_parameters(parameters(_solver->parameter_type()));

  // Reset iteration counts
  _newton_iteration = 0;
  _krylov_iterations = 0;

  // Compute F(u)
  nonlinear_problem.form(*_matA, *_matP, *_b, x);
  nonlinear_problem.F(*_b, x);

  // Check convergence
  bool newton_converged = false;
  if (convergence_criterion == "residual")
    newton_converged = converged(*_b, nonlinear_problem, 0);
  else if (convergence_criterion == "incremental")
  {
    // We need to do at least one Newton step with the ||dx||-stopping
    // criterion.
    newton_converged = false;
  }
  else
  {
    dolfin_error("NewtonSolver.cpp",
                 "check for convergence",
                 "The convergence criterion %s is unknown, known criteria are 'residual' or 'incremental'",
                 convergence_criterion.c_str());
  }

  // Start iterations
  while (!newton_converged && _newton_iteration < maxiter)
  {
    // Compute Jacobian
    nonlinear_problem.J(*_matA, x);
    nonlinear_problem.J_pc(*_matP, x);

    // Setup (linear) solver (including set operators)
    solver_setup(_matA, _matP, nonlinear_problem, _newton_iteration);

    // Perform linear solve and update total number of Krylov
    // iterations
    if (!_dx->empty())
      _dx->zero();
    //jMODS:  CALL TO SOLVER:
////////////////////////////////////////////////////////////////////////    
    _krylov_iterations += _solver->solve(*_dx, *_b);
////////////////////////////////////////////////////////////////////////    

    // Update solution
    update_solution(x, *_dx, _relaxation_parameter,
                    nonlinear_problem, _newton_iteration);

    // Increment iteration count
    _newton_iteration++;

    // FIXME: This step is not needed if residual is based on dx and
    //        this has converged.
    // FIXME: But, this function call may update internal variable, etc.
    // Compute F
    nonlinear_problem.form(*_matA, *_matP, *_b, x);
    nonlinear_problem.F(*_b, x);


    //jMODS:
     double  r_dx = (*_dx).norm("l2");
     double  r_b = (*_b).norm("l2");

    if (_newton_iteration == 1)//j: already incremented above
    {//save the first residuals:
      _r_dx_nprev = r_dx;
      _r_b_nprev = r_b;
      _dx0 = r_dx;
      _b0 = r_b;
    }
    // jMOD: iterative Relative residual
     double r_dx_contraction_factor = r_dx/_r_dx_nprev;
     double r_b_contraction_factor = r_b/_r_b_nprev;

     if(_newton_iteration <= 2)
     {//save the first contraction factor:
        _r1_dx_contraction_factor = r_dx_contraction_factor;
        _r1_b_contraction_factor = r_b_contraction_factor;
     }
    //jMOD: update prev. residual buffer
    _r_dx_nprev = r_dx;
    _r_b_nprev = r_b;
    //jMOD:
    info("jTheta_dx = %.3e jTheta_b = %.3e [norm(l2)](dx,b) = (%.3e, %.3e)",
         r_dx_contraction_factor, r_b_contraction_factor, r_dx, r_b);

    //JW: RETURN IF CONTRACTION FACTOR > 1.0 HERE:
    if( (r_dx_contraction_factor > 1.0) )
    // if( (r_dx_contraction_factor > 1.0) || (r_b_contraction_factor > 1.0))
    // if(0)
    {
      double correction_dx = (sqrt(2.0) - 1.0) / (sqrt(4.0*r_dx_contraction_factor + 1.0) - 1.0);
      double correction_b = (sqrt(2.0) - 1.0) / (sqrt(4.0*r_b_contraction_factor + 1.0) - 1.0);
      info("CONTRACTION > 1.0: COMPUTED CORRECTIONS: %.3e, %.3e");
      // info("CONTRACTION > 0.85: COMPUTED CORRECTIONS: %.3e, %.3e", correction_dx, correction_b);
      return std::make_tuple(_newton_iteration, false, _dx0, correction_dx);
    }

    // Test for convergence
    if (convergence_criterion == "residual")
    {
      newton_converged = converged(*_b, nonlinear_problem, 
                                    _newton_iteration);
    }
    else if (convergence_criterion == "incremental")
    {
      // Subtract 1 to make sure that the initial residual0 is
      // properly set.
      newton_converged = converged(*_dx, nonlinear_problem,
                                    _newton_iteration - 1);
    }
    else
    {
      dolfin_error("NewtonSolver.cpp",
                   "check for convergence",
                   "The convergence criterion %s is unknown, known criteria are 'residual' or 'incremental'",
                   convergence_criterion.c_str());
    }
  }

  if (newton_converged)
  {
    if (dolfin::MPI::rank(_mpi_comm) == 0)
    {
      info("Newton solver finished in %d iterations and %d linear solver iterations.",
           _newton_iteration, _krylov_iterations);
    }
  }
  else
  {
    const bool error_on_nonconvergence = parameters["error_on_nonconvergence"];
    if (error_on_nonconvergence)
    {
      if (_newton_iteration == maxiter)
      {
        dolfin_error("NewtonSolver.cpp",
                     "solve nonlinear system with NewtonSolver",
                     "Newton solver did not converge because maximum number of iterations reached");
      }
      else
      {
        dolfin_error("NewtonSolver.cpp",
                     "solve nonlinear system with NewtonSolver",
                     "Newton solver did not converge");
      }
    }
    else
      warning("Newton solver did not converge.");
  }

//jMODS
  return std::make_tuple(_newton_iteration, newton_converged, _dx0, _r1_dx_contraction_factor);
  // return std::make_pair(_newton_iteration, newton_converged);
}
//-----------------------------------------------------------------------------
std::size_t jNewtonSolver::iteration() const
{
  return _newton_iteration;
}
//-----------------------------------------------------------------------------
std::size_t jNewtonSolver::krylov_iterations() const
{
  return _krylov_iterations;
}
//-----------------------------------------------------------------------------
double jNewtonSolver::residual() const
{
  return _residual;
}
//-----------------------------------------------------------------------------
double jNewtonSolver::residual0() const
{
  return _residual0;
}
//-----------------------------------------------------------------------------
double jNewtonSolver::relative_residual() const
{
  return _residual/_residual0;
}
//-----------------------------------------------------------------------------
GenericLinearSolver& jNewtonSolver::linear_solver() const
{
  if (!_solver)
  {
    dolfin_error("NewtonSolver.cpp",
                 "access linear solver for Newton solver",
                 "The linear solver will not be initialized until solve() "
                 "has been called. For control of the linear solver, pass "
                 "a linear solver to the constructor of NewtonSolver");
  }

  return *_solver;
}
//-----------------------------------------------------------------------------
bool jNewtonSolver::converged(const GenericVector& r,
                             const NonlinearProblem& nonlinear_problem,
                             std::size_t newton_iteration)
{
  const double rtol = parameters["relative_tolerance"];
  const double atol = parameters["absolute_tolerance"];
  const bool report = parameters["report"];

  _residual = r.norm("l2");

  // If this is the first iteration step, set initial residual
  if (newton_iteration == 0)
    _residual0 = _residual;

  // Relative residual
  const double relative_residual = _residual/_residual0;

  // Output iteration number and residual
  if (report && dolfin::MPI::rank(_mpi_comm) == 0)
  {
    info("jNewton iteration %d: r (abs) = %.3e (tol = %.3e) r (rel) = %.3e (tol = %.3e)",
         newton_iteration, _residual, atol, relative_residual, rtol);
  }

  // Return true if convergence criterion is met
  if (relative_residual < rtol || _residual < atol)
    return true;
  else
    return false;
}
//-----------------------------------------------------------------------------
void jNewtonSolver::solver_setup(std::shared_ptr<const GenericMatrix> A,
                                std::shared_ptr<const GenericMatrix> P,
                                const NonlinearProblem& nonlinear_problem,
                                std::size_t interation)
{
  // Update Jacobian in linear solver (and preconditioner if given)
  if (_matP->empty())
  {
    _solver->set_operator(A);
    log(TRACE, "NewtonSolver: using Jacobian as preconditioner matrix");
  }
  else
  {
    _solver->set_operators(A, P);
  }
}
//-----------------------------------------------------------------------------
void jNewtonSolver::update_solution(GenericVector& x, const GenericVector& dx,
                                   double relaxation_parameter,
                                   const NonlinearProblem& nonlinear_problem,
                                   std::size_t interation)
{
  if (relaxation_parameter == 1.0)
    x -= dx;
  else
    x.axpy(-relaxation_parameter, dx);
}
//-----------------------------------------------------------------------------
