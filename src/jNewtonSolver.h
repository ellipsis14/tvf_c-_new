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
// Modified by Anders Logg 2006-2011
// Modified by Anders E. Johansen 2011
//
// First added:  2005-10-23
// Last changed: 2013-11-20

#ifndef __jNEWTON_SOLVER_H
#define __jNEWTON_SOLVER_H

#include <utility>
#include <memory>
#include <dolfin/common/MPI.h>
#include <dolfin/common/Variable.h>

namespace dolfin
{

  // Forward declarations
  class GenericLinearSolver;
  class GenericLinearAlgebraFactory;
  class GenericMatrix;
  class GenericVector;
  class NonlinearProblem;

  /// This class defines a Newton solver for nonlinear systems of
  /// equations of the form :math:`F(x) = 0`.

  class jNewtonSolver : public Variable
  {
  public:

    /// Create nonlinear solver
    explicit jNewtonSolver(MPI_Comm comm=MPI_COMM_WORLD);

    /// Create nonlinear solver using provided linear solver
    ///
    /// *Arguments*
    ///     comm (_MPI_Ccmm_)
    ///         The MPI communicator.
    ///     solver (_GenericLinearSolver_)
    ///         The linear solver.
    ///     factory (_GenericLinearAlgebraFactory_)
    ///         The factory.
    jNewtonSolver(MPI_Comm comm, std::shared_ptr<GenericLinearSolver> solver,
                 GenericLinearAlgebraFactory& factory);

    /// Destructor
    virtual ~jNewtonSolver();

    /// Solve abstract nonlinear problem :math:`F(x) = 0` for given
    /// :math:`F` and Jacobian :math:`\dfrac{\partial F}{\partial x}`.
    ///
    /// *Arguments*
    ///     nonlinear_function (_NonlinearProblem_)
    ///         The nonlinear problem.
    ///     x (_GenericVector_)
    ///         The vector.
    ///
    /// *Returns*
    ///     std::pair<std::size_t, bool>
    ///         Pair of number of Newton iterations, and whether
    ///         iteration converged)
    
    //jMODS
    std::tuple<std::size_t, bool, double, double>
    solve(NonlinearProblem& nonlinear_function,
                                       GenericVector& x);
    // std::pair<std::size_t, bool> solve(NonlinearProblem& nonlinear_function,
    //                                    GenericVector& x);

    /// Return current Newton iteration number
    ///
    /// *Returns*
    ///     std::size_t
    ///         The iteration number.
    std::size_t iteration() const;

    /// Return number of Krylov iterations elapsed since
    /// solve started
    ///
    /// *Returns*
    ///     std::size_t
    ///         The number of iterations.
    std::size_t krylov_iterations() const;

    /// Return current residual
    ///
    /// *Returns*
    ///     double
    ///         Current residual.
    double residual() const;

    /// Return initial residual
    ///
    /// *Returns*
    ///     double
    ///         Initial residual.
    double residual0() const;

    /// Return current relative residual
    ///
    /// *Returns*
    ///     double
    ///       Current relative residual.
    double relative_residual() const;

    /// Return the linear solver
    ///
    /// *Returns*
    ///     _GenericLinearSolver_
    ///         The linear solver.
    GenericLinearSolver& linear_solver() const;

    /// Default parameter values
    ///
    /// *Returns*
    ///     _Parameters_
    ///         Parameter values.
    static Parameters default_parameters();

    /// Set relaxation parameter. Default value 1.0 means full
    /// Newton method, value smaller than 1.0 relaxes the method
    /// by shrinking effective Newton step size by the given factor.
    ///
    /// *Arguments*
    ///     relaxation_parameter(double)
    ///         Relaxation parameter value.
    void set_relaxation_parameter(double relaxation_parameter)
    { _relaxation_parameter = relaxation_parameter; }

    /// Get relaxation parameter
    ///
    /// *Returns*
    ///     double
    ///         Relaxation parameter value.
    double get_relaxation_parameter()
    { return _relaxation_parameter; }

  protected:

    /// Convergence test. It may be overloaded using virtual inheritance and
    /// this base criterion may be called from derived, both in C++ and Python.
    ///
    /// *Arguments*
    ///     r (_GenericVector_)
    ///         Residual for criterion evaluation.
    ///     nonlinear_problem (_NonlinearProblem_)
    ///         The nonlinear problem.
    ///     iteration (std::size_t)
    ///         Newton iteration number.
    ///
    /// *Returns*
    ///     bool
    ///         Whether convergence occurred.
    virtual bool converged(const GenericVector& r,
                           const NonlinearProblem& nonlinear_problem,
                           std::size_t iteration);

    /// Setup solver to be used with system matrix A and preconditioner
    /// matrix P. It may be overloaded to get finer control over linear
    /// solver setup, various linesearch tricks, etc. Note that minimal
    /// implementation should call *set_operators* method of the linear
    /// solver.
    ///
    /// *Arguments*
    ///     A (_std::shared_ptr<const GenericMatrix>_)
    ///         System Jacobian matrix.
    ///     J (_std::shared_ptr<const GenericMatrix>_)
    ///         System preconditioner matrix.
    ///     nonlinear_problem (_NonlinearProblem_)
    ///         The nonlinear problem.
    ///     iteration (std::size_t)
    ///         Newton iteration number.
    virtual void solver_setup(std::shared_ptr<const GenericMatrix> A,
                              std::shared_ptr<const GenericMatrix> P,
                              const NonlinearProblem& nonlinear_problem,
                              std::size_t iteration);

    /// Update solution vector by computed Newton step. Default
    /// update is given by formula::
    ///
    ///   x -= relaxation_parameter*dx
    ///
    /// *Arguments*
    ///     x (_GenericVector>_)
    ///         The solution vector to be updated.
    ///     dx (_GenericVector>_)
    ///         The update vector computed by Newton step.
    ///     relaxation_parameter (double)
    ///         Newton relaxation parameter.
    ///     nonlinear_problem (_NonlinearProblem_)
    ///         The nonlinear problem.
    ///     iteration (std::size_t)
    ///         Newton iteration number.
    virtual void update_solution(GenericVector& x,
                                 const GenericVector& dx,
                                 double relaxation_parameter,
                                 const NonlinearProblem& nonlinear_problem,
                                 std::size_t iteration);

  private:

    // Current number of Newton iterations
    std::size_t _newton_iteration;

    // Accumulated number of Krylov iterations since solve began
    std::size_t _krylov_iterations;

    // Relaxation parameter
    double _relaxation_parameter;

    // Most recent residual and initial residual
    double _residual, _residual0;
    //jMOD:
    double _r_dx_nprev, _r_b_nprev, _dx0, _b0;
    double _r1_dx_contraction_factor, _r1_b_contraction_factor;
    // Solver
    std::shared_ptr<GenericLinearSolver> _solver;

    // Jacobian matrix
    std::shared_ptr<GenericMatrix> _matA;

    // Preconditioner matrix
    std::shared_ptr<GenericMatrix> _matP;

    // Solution vector
    std::shared_ptr<GenericVector> _dx;

    // Residual vector
    std::shared_ptr<GenericVector> _b;

    // MPI communicator
    MPI_Comm _mpi_comm;

  };

}

#endif
