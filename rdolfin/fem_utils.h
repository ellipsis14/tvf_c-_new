// Copyright (C) 2013, 2015, 2016 Johan Hake, Jan Blechta
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

#ifndef __FEM_UTILS_H
#define __FEM_UTILS_H

#include <vector>

#include <dolfin/common/types.h>

namespace dolfin
{
  // Forward declarations
  class Mesh;
  class FunctionSpace;
  class Function;
  class MeshGeometry;

  /// Return a map between dof indices and vertex indices
  ///
  /// Only works for FunctionSpace with dofs exclusively on vertices.
  /// For mixed FunctionSpaces vertex index is offset with the number
  /// of dofs per vertex.
  ///
  /// In parallel the returned map maps both owned and unowned dofs
  /// (using local indices) thus covering all the vertices. Hence the
  /// returned map is an inversion of _vertex_to_dof_map_.
  ///
  /// @param    space (_FunctionSpace_)
  ///         The FunctionSpace for what the dof to vertex map should
  ///         be computed for
  ///
  /// @return   std::vector<std::size_t>
  ///         The dof to vertex map
  std::vector<std::size_t> dof_to_vertex_map(const FunctionSpace& space);

  /// Return a map between vertex indices and dof indices
  ///
  /// Only works for FunctionSpace with dofs exclusively on vertices.
  /// For mixed FunctionSpaces dof index is offset with the number of
  /// dofs per vertex.
  ///
  /// @param    space (_FunctionSpace_)
  ///         The FunctionSpace for what the vertex to dof map should
  ///         be computed for
  ///
  /// @return    std::vector<dolfin::la_index>
  ///         The vertex to dof map
  std::vector<dolfin::la_index> vertex_to_dof_map(const FunctionSpace& space);

  /// Sets mesh coordinates from function
  ///
  /// Mesh connectivities d-0, d-1, ..., d-r are built on function mesh
  /// (where d is topological dimension of the mesh and r is maximal
  /// dimension of entity associated with any coordinate node). Consider
  /// clearing unneeded connectivities when finished.
  ///
  /// @param   geometry (_MeshGeometry_)
  ///         Mesh geometry to be set
  /// @param    position (_Function_)
  ///         Vectorial Lagrange function with matching degree and mesh
  void set_coordinates(MeshGeometry& geometry, const Function& position);

  /// Stores mesh coordinates into function
  ///
  /// Mesh connectivities d-0, d-1, ..., d-r are built on function mesh
  /// (where d is topological dimension of the mesh and r is maximal
  /// dimension of entity associated with any coordinate node). Consider
  /// clearing unneeded connectivities when finished.
  ///
  /// @param   position (_Function_)
  ///         Vectorial Lagrange function with matching degree and mesh
  /// @param    geometry (_MeshGeometry_)
  ///         Mesh geometry to be stored
  void get_coordinates(Function& position, const MeshGeometry& geometry);

  /// Creates mesh from coordinate function
  ///
  /// Topology is given by underlying mesh of the function space and
  /// geometry is given by function values. Hence resulting mesh
  /// geometry has a degree of the function space degree. Geometry of
  /// function mesh is ignored.
  ///
  /// Mesh connectivities d-0, d-1, ..., d-r are built on function mesh
  /// (where d is topological dimension of the mesh and r is maximal
  /// dimension of entity associated with any coordinate node). Consider
  /// clearing unneeded connectivities when finished.
  ///
  /// @param coordinates
  /// (_Function_)
  ///         Vector Lagrange function of any degree
  ///
  /// @return Mesh
  ///         The mesh
  Mesh create_mesh(Function& coordinates);
}

#endif
