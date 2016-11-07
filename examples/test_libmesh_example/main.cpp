// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1> Systems Example 4 - Linear Elastic Cantilever </h1>
// \author David Knezevic
// \date 2012
//
// In this example we model a homogeneous isotropic cantilever
// using the equations of linear elasticity. We set the Poisson ratio to
// \nu = 0.3 and clamp the left boundary and apply a vertical load at the
// right boundary.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include "libmesh/linear_solver.h"
#include "libmesh/enum_preconditioner_type.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
// Bring in everything from the libMesh namespace
using namespace libMesh;



// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries
  LibMeshInit init (argc, argv);

  // Initialize the cantilever mesh
  const unsigned int dim = 2;

  // Create a 2D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);

  MeshTools::Generation::build_square (mesh,
                                       2, 2,
                                       0., 1.,
                                       0., 0.2,
                                       TRI3);
  // Print information about the mesh to the screen.
//  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Add two displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
//  unsigned int v_var = system.add_variable("uy", FIRST, LAGRANGE);
  system.add_vector("vec");

  equation_systems.init();
   equation_systems.parameters.set<double>("nu") = 0.5;
   equation_systems.parameters.set<double>("E") = 250.0;
  equation_systems.print_info();

  equation_systems.write("test.dat",libMesh::WRITE);

    ExodusII_IO exo(mesh) ;
    exo.write_equation_systems ("test.exo", equation_systems);
    exo.write_timestep ("test.exo", equation_systems, 1, 0.1);
    system.solution->set(2,2);
    exo.write_timestep ("test.exo", equation_systems, 2, 0.2);

  Mesh mesh2(init.comm(), dim);
    MeshTools::Generation::build_square (mesh2,
                                       2, 3,
                                       0., 1.,
                                       0., 0.2,
                                       TRI3);
  // We may need XDR support compiled in to read binary .xdr files
  std::string meshfile = "test.exo";
      // Read the input mesh.
  std::cout << "Reading Mesh " << std::endl;
    ExodusII_IO exo2(mesh2) ;
    exo2.read (&meshfile[0]);
  std::cout << "Done " << std::endl;


//  mesh2.print_info();
  // Create an equation systems object.
  EquationSystems equation_systems2 (mesh2);
  LinearImplicitSystem & system2 =
    equation_systems2.add_system<LinearImplicitSystem> ("Elasticity");
    // Add two displacement variables, u and v, to the system
  unsigned int u2_var = system2.add_variable("u", FIRST, LAGRANGE);
//  unsigned int v2_var = system2.add_variable("uy", FIRST, LAGRANGE);
    equation_systems2.init();
  equation_systems2.print_info();

  exo2.copy_nodal_solution(system2, "u", "u", 2);

    ExodusII_IO exo_2(mesh) ;
    exo_2.write_equation_systems ("test2.exo", equation_systems);
    exo_2.write_timestep ("test2.exo", equation_systems, 1, 0.0);


  // All done.
  return 0;
}



