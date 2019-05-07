// Copyright (C) 2007-2012 Anders Logg and Garth N. Wells
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
// Modified by Kristian Oelgaard, 2007
// Modified by Johan Hake, 2009
// Modified by Joachim B Haga, 2012
// Modified by Mikael Mortensen, 2014
//
// First added:  2007-04-10
// Last changed: 2014-01-23
//
// FIXME: This class needs some cleanup, in particular collecting
//        all data from different representations into a common
//        data structure (perhaps an std::vector<std::size_t> with
//        facet indices).

#ifndef __DIRICHLET_BC_H
#define __DIRICHLET_BC_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>
#include <memory>
#include <unordered_map>

#include <dolfin/common/types.h>
#include <dolfin/common/Hierarchical.h>
#include <dolfin/common/MPI.h>
#include <dolfin/common/Variable.h>

namespace dolfin
{

  class GenericFunction;
  class FunctionSpace;
  class Facet;
  class GenericMatrix;
  class GenericVector;
  class SubDomain;
  template<typename T> class MeshFunction;

  /// Interface for setting (strong) Dirichlet boundary conditions.

  ///
  ///     u = g on G,
  ///
  /// where u is the solution to be computed, g is a function
  /// and G is a sub domain of the mesh.
  ///
  /// A DirichletBC is specified by the function g, the function space
  /// (trial space) and boundary indicators on (a subset of) the mesh
  /// boundary.
  ///
  /// The boundary indicators may be specified in a number of
  /// different ways.
  ///
  /// The simplest approach is to specify the boundary by a _SubDomain_
  /// object, using the inside() function to specify on which facets
  /// the boundary conditions should be applied. The boundary facets
  /// will then be searched for and marked *only* on the first call to
  /// apply. This means that the mesh could be moved after the first
  /// apply and the boundary markers would still remain intact.
  ///
  /// Alternatively, the boundary may be specified by a _MeshFunction_
  /// over facets labeling all mesh facets together with a number that
  /// specifies which facets should be included in the boundary.
  ///
  /// The third option is to attach the boundary information to the
  /// mesh. This is handled automatically when exporting a mesh from
  /// for example VMTK.
  ///
  /// The 'method' variable may be used to specify the type of method
  /// used to identify degrees of freedom on the boundary. Available
  /// methods are: topological approach (default), geometric approach,
  /// and pointwise approach. The topological approach is faster, but
  /// will only identify degrees of freedom that are located on a
  /// facet that is entirely on the boundary. In particular, the
  /// topological approach will not identify degrees of freedom for
  /// discontinuous elements (which are all internal to the cell). A
  /// remedy for this is to use the geometric approach. In the
  /// geometric approach, each dof on each facet that matches the
  /// boundary condition will be checked. To apply pointwise boundary
  /// conditions e.g. pointloads, one will have to use the pointwise
  /// approach. The three possibilities are "topological", "geometric"
  /// and "pointwise".
  ///
  /// Note: when using "pointwise", the boolean argument `on_boundary`
  /// in SubDomain::inside will always be false.
  ///
  /// The 'check_midpoint' variable can be used to decide whether or
  /// not the midpoint of each facet should be checked when a
  /// user-defined _SubDomain_ is used to define the domain of the
  /// boundary condition. By default, midpoints are always checked.
  /// Note that this variable may be of importance close to corners,
  /// in which case it is sometimes important to check the midpoint to
  /// avoid including facets "on the diagonal close" to a corner. This
  /// variable is also of importance for curved boundaries (like on a
  /// sphere or cylinder), in which case it is important *not* to
  /// check the midpoint which will be located in the interior of a
  /// domain defined relative to a radius.
  ///
  /// Note that there may be caching employed in BC computation for
  /// performance reasons. In particular, applicable DOFs are cached
  /// by some methods on a first apply(). This means that changing a
  /// supplied object (defining boundary subdomain) after first use may
  /// have no effect. But this is implementation and method specific.

  class DirichletBC : public Hierarchical<DirichletBC>, public Variable
  {

  public:

    /// map type used by DirichletBC
    typedef std::unordered_map<std::size_t, double> Map;

    /// Create boundary condition for subdomain
    ///
    /// @param[in] V (FunctionSpace)
    ///         The function space
    /// @param[in] g (GenericFunction)
    ///         The value
    /// @param[in] sub_domain (SubDomain)
    ///         The subdomain
    /// @param[in] method (std::string)
    ///         Optional argument: A string specifying
    ///         the method to identify dofs
    /// @param[in] check_midpoint (bool)
    DirichletBC(std::shared_ptr<const FunctionSpace> V,
                std::shared_ptr<const GenericFunction> g,
                std::shared_ptr<const SubDomain> sub_domain,
                std::string method="topological",
                bool check_midpoint=true);

    /// Create boundary condition for subdomain specified by index
    ///
    /// @param[in] V (FunctionSpace)
    ///         The function space.
    /// @param[in] g (GenericFunction)
    ///         The value.
    /// @param[in] sub_domains (MeshFnunction<std::size_t>)
    ///         Subdomain markers
    /// @param[in] sub_domain (std::size_t)
    ///         The subdomain index (number)
    /// @param[in] method (std::string)
    ///         Optional argument: A string specifying the
    ///         method to identify dofs.
    DirichletBC(std::shared_ptr<const FunctionSpace> V,
                std::shared_ptr<const GenericFunction> g,
                std::shared_ptr<const MeshFunction<std::size_t>> sub_domains,
                std::size_t sub_domain,
                std::string method="topological");

    // TODO: Remove/deprecate this function
    /// Create boundary condition for boundary data included in the mesh
    ///
    /// @param[in] V (FunctionSpace)
    ///         The function space.
    /// @param[in]  g (GenericFunction)
    ///         The value.
    /// @param[in] sub_domain (std::size_t)
    ///         The subdomain index (number)
    /// @param[in] method (std::string)
    ///         Optional argument: A string specifying the
    ///         method to identify dofs.
    DirichletBC(std::shared_ptr<const FunctionSpace> V,
                std::shared_ptr<const GenericFunction> g,
                std::size_t sub_domain,
                std::string method="topological");

    /// Create boundary condition for subdomain by boundary markers
    /// (cells, local facet numbers)
    ///
    /// @param[in] V (FunctionSpace)
    ///         The function space.
    /// @param[in] g (GenericFunction)
    ///         The value.
    /// @param[in] markers (std::vector<std:size_t>&)
    ///         Subdomain markers (facet index local to process)
    /// @param[in] method (std::string)
    ///         Optional argument: A string specifying the
    ///         method to identify dofs.
    DirichletBC(std::shared_ptr<const FunctionSpace> V,
                std::shared_ptr<const GenericFunction> g,
                const std::vector<std::size_t>& markers,
                std::string method="topological");

    /// Copy constructor. Either cached DOF data are copied.
    ///
    /// @param[in] bc (DirichletBC&)
    ///         The object to be copied.
    DirichletBC(const DirichletBC& bc);

    /// Destructor
    ~DirichletBC();

    /// Assignment operator. Either cached DOF data are assigned.
    ///
    /// @param[in] bc (DirichletBC)
    ///         Another DirichletBC object.
    const DirichletBC& operator= (const DirichletBC& bc);

    /// Apply boundary condition to a matrix
    ///
    /// @param[in,out] A (GenericMatrix)
    ///         The matrix to apply boundary condition to.
    void apply(GenericMatrix& A) const;

    /// Apply boundary condition to a vector
    ///
    /// @param[in,out] b (GenericVector)
    ///         The vector to apply boundary condition to.
    void apply(GenericVector& b) const;

    /// Apply boundary condition to a linear system
    ///
    /// @param[in,out] A (GenericMatrix)
    ///         The matrix to apply boundary condition to.
    /// @param[in,out] b (GenericVector)
    ///         The vector to apply boundary condition to.
    void apply(GenericMatrix& A, GenericVector& b) const;

    /// Apply boundary condition to vectors for a nonlinear problem
    ///
    /// @param[in,out] b (GenericVector)
    ///         The vector to apply boundary conditions to.
    /// @param[in] x (GenericVector)
    ///         Another vector (nonlinear problem).
    void apply(GenericVector& b, const GenericVector& x) const;

    /// Apply boundary condition to a linear system for a nonlinear problem
    ///
    /// @param[in,out] A (GenericMatrix)
    ///         The matrix to apply boundary conditions to.
    /// @param[in,out] b (GenericVector)
    ///         The vector to apply boundary conditions to.
    /// @param[in] x (GenericVector)
    ///         Another vector (nonlinear problem).
    void apply(GenericMatrix& A, GenericVector& b,
               const GenericVector& x) const;

    /// Get Dirichlet dofs and values. If a method other than 'pointwise' is
    /// used in parallel, the map may not be complete for local vertices since
    /// a vertex can have a bc applied, but the partition might not have a
    /// facet on the boundary. To ensure all local boundary dofs are marked,
    /// it is necessary to call gather() on the returned boundary values.
    ///
    /// @param[in,out] boundary_values (Map&)
    ///         Map from dof to boundary value.
    void get_boundary_values(Map& boundary_values) const;

    /// Get boundary values from neighbour processes. If a method other than
    /// "pointwise" is used, this is necessary to ensure all boundary dofs are
    /// marked on all processes.
    ///
    /// @param[in,out] boundary_values (Map&)
    ///         Map from dof to boundary value.
    void gather(Map& boundary_values) const;

    /// Make rows of matrix associated with boundary condition zero,
    /// useful for non-diagonal matrices in a block matrix.
    ///
    /// @param[in,out] A (GenericMatrix&)
    ///         The matrix
    void zero(GenericMatrix& A) const;

    /// Make columns of matrix associated with boundary condition
    /// zero, and update a (right-hand side) vector to reflect the
    /// changes. Useful for non-diagonals.
    ///
    /// @param[in,out] A (GenericMatrix&)
    ///         The matrix
    /// @param[in,out] b (GenericVector&)
    ///         The vector
    /// @param[in] diag_val (double)
    ///         This parameter would normally be -1, 0 or 1.
    void zero_columns(GenericMatrix& A, GenericVector& b,
                      double diag_val=0) const;

    /// Return boundary markers
    ///
    /// @return std::vector<std::size_t>&
    ///         Boundary markers (facets stored as pairs of cells and
    ///         local facet numbers).
    const std::vector<std::size_t>& markers() const;

    /// Return function space V
    ///
    /// @return FunctionSpace
    ///         The function space to which boundary conditions are applied.
    std::shared_ptr<const FunctionSpace> function_space() const
    { return _function_space; }

    /// Return boundary value g
    ///
    /// @return GenericFunction
    ///         The boundary values.
    std::shared_ptr<const GenericFunction> value() const;

    /// Return shared pointer to subdomain
    ///
    /// @return SubDomain
    ///         Shared pointer to subdomain.
    std::shared_ptr<const SubDomain> user_sub_domain() const;

    /// Set value g for boundary condition, domain remains unchanged
    ///
    /// @param[in] g (GenericFucntion)
    ///         The value.
    void set_value(std::shared_ptr<const GenericFunction> g);

    /// Set value to 0.0
    void homogenize();

    /// Return method used for computing Dirichlet dofs
    ///
    /// @return std::string
    ///         Method used for computing Dirichlet dofs ("topological",
    ///         "geometric" or "pointwise").
    std::string method() const;

    /// Default parameter values
    /// @return Parameters
    static Parameters default_parameters()
    {
      Parameters p("dirichlet_bc");
      p.add("use_ident", true);
      p.add("check_dofmap_range", true);
      return p;
    }

  private:

    class LocalData;

    // Apply boundary conditions, common method
    void apply(GenericMatrix* A, GenericVector* b,
               const GenericVector* x) const;

    // Check input data to constructor
    void check() const;

    // Initialize facets (from sub domain, mesh, etc)
    void init_facets(const MPI_Comm mpi_comm) const;

    // Initialize sub domain markers from sub domain
    void
      init_from_sub_domain(std::shared_ptr<const SubDomain> sub_domain) const;

    // Initialize sub domain markers from MeshFunction
    void init_from_mesh_function(const MeshFunction<std::size_t>& sub_domains,
                                 std::size_t sub_domain) const;

    // Initialize sub domain markers from mesh
    void init_from_mesh(std::size_t sub_domain) const;

    // Compute dofs and values for application of boundary conditions
    // using given method
    void compute_bc(Map& boundary_values, LocalData& data,
                    std::string method) const;

    // Compute boundary values for facet (topological approach)
    void compute_bc_topological(Map& boundary_values,
                                LocalData& data) const;

    // Compute boundary values for facet (geometrical approach)
    void compute_bc_geometric(Map& boundary_values,
                              LocalData& data) const;

    // Compute boundary values for facet (pointwise approach)
    void compute_bc_pointwise(Map& boundary_values,
                              LocalData& data) const;

    // Check if the point is in the same plane as the given facet
    bool on_facet(const double* coordinates, const Facet& facet) const;

    // Check arguments for compatibility of tensors and dofmap,
    // dim is means an axis to which bc applies
    void check_arguments(GenericMatrix* A, GenericVector* b,
                         const GenericVector* x,
                         std::size_t dim) const;

    // The function space (possibly a sub function space)
    std::shared_ptr<const FunctionSpace> _function_space;

    // The function
    std::shared_ptr<const GenericFunction> _g;

    // Search method
    std::string _method;

    // Possible search methods
    static const std::set<std::string> methods;

  public:

    // User defined sub domain
    std::shared_ptr<const SubDomain> _user_sub_domain;

  private:

    // Cached number of bc dofs, used for memory allocation on second use
    mutable std::size_t _num_dofs;

    // Boundary facets, stored by facet index (local to process)
    mutable std::vector<std::size_t> _facets;

    // Cells attached to boundary, stored by cell index with map to
    // local dof number
    mutable std::map<std::size_t, std::vector<std::size_t>>
      _cells_to_localdofs;

    // User defined mesh function
    std::shared_ptr<const MeshFunction<std::size_t>> _user_mesh_function;

    // User defined sub domain marker for mesh or mesh function
    std::size_t _user_sub_domain_marker;

    // Flag for whether midpoints should be checked
    bool _check_midpoint;

    // Local data for application of boundary conditions
    class LocalData
    {
    public:

      // Constructor
      LocalData(const FunctionSpace& V);

      // Coefficients
      std::vector<double> w;

      // Facet dofs
      std::vector<std::size_t> facet_dofs;

      // Coordinates for dofs
      boost::multi_array<double, 2> coordinates;

    };


  };

}

#endif
