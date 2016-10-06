/*
============================================================================

	.______    _______     ___   .___________.    __  .___________.
    |   _  \  |   ____|   /   \  |           |   |  | |           |
    |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
    |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |     
    |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |     
    |______/  |_______/__/     \__\  |__|        |__|     |__|     
                                                               
    BeatIt - code for cardiovascular simulations
    Copyright (C) 2016 Simone Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
*/

/**
 * \file main.cpp
 *
 * \class main
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 11, 2016
 *
 */


// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/exodusII_io.h"


int main (int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace

    using namespace libMesh;
    MPI_Init (&argc, &argv);
    std::cout << "Starting ... "  << std::endl;

      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);
      std::cout << "Ending ... "  << std::endl;

      // Create a mesh, with dimension to be overridden later, distributed
      // across the default MPI communicator.
      Mesh mesh(init.comm());

      // Use the MeshTools::Generation mesh generator to create a uniform
      // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
      // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
      // elements instead of the default QUAD4's we used in example 2
      // allow us to use higher-order approximation.
      MeshTools::Generation::build_cube (mesh,
                                         2, 2, 2,
                                         0., 1.,
                                         0., 1.,
                                         0., 1.,
                                         TET4);

      // Create an equation systems object.
      EquationSystems equation_systems (mesh);
      // Declare the Poisson system and its variables.
      // The Poisson system is another example of a steady system.
      equation_systems.add_system<LinearImplicitSystem> ("Poisson");

      // Adds the variable "u" to "Poisson".  "u"
      // will be approximated using second-order approximation.
      equation_systems.get_system("Poisson").add_variable("u", FIRST);
      equation_systems.get_system("Poisson").add_variable("p", FIRST);
      equation_systems.get_system("Poisson").add_vector("prova");
      equation_systems.init();

      auto my_length = equation_systems.get_system("Poisson").get_vector(0).local_size();

      std::cout << "Comm size = " << init.comm().size() << std::endl;


      if( 0 == init.comm().rank() ) std::cout << "My Length on 0 is " << my_length << std::endl;
      if( 1 == init.comm().rank() ) std::cout << "My Length on 1 is " << my_length << std::endl;


      auto first_local_index = equation_systems.get_system("Poisson").get_vector(0).first_local_index();
      auto last_local_index = equation_systems.get_system("Poisson").get_vector(0).last_local_index();

      equation_systems.get_system("Poisson").get_vector(0).set(first_local_index+1, 3.0);
      equation_systems.get_system("Poisson").get_vector(0).print(std::cout);

      auto rank = equation_systems.comm().rank();
      for(int i = first_local_index; i < last_local_index; ++i)
      {
             if(rank == 0 )  equation_systems.get_system("Poisson").solution->add(i, 1.0);
             if(rank == 1 )  equation_systems.get_system("Poisson").solution->add(i, -1.0);
      }
      if(rank == 1 )  equation_systems.get_system("Poisson").solution->add(first_local_index-1, -1.0);
      equation_systems.get_system("Poisson").solution->close();
      equation_systems.get_system("Poisson").solution->print(std::cout);

      if(rank == 1 )  equation_systems.get_system("Poisson").solution->set(last_local_index-1, 0.0);
      equation_systems.comm().barrier();
      equation_systems.get_system("Poisson").solution->close();
      equation_systems.get_system("Poisson").solution->print(std::cout);
      if(rank == 0 )  equation_systems.get_system("Poisson").solution->add(last_local_index, 2.0);
      equation_systems.get_system("Poisson").solution->close();
      equation_systems.get_system("Poisson").solution->print(std::cout);

      libMesh::ExodusII_IO exo(mesh);
      exo.write_equation_systems ("meshfile.e", equation_systems);

	  auto init_val = equation_systems.get_system("Poisson").solution->operator()(first_local_index);
	  auto last_val = equation_systems.get_system("Poisson").solution->operator()(last_local_index-1);

	  std::cout << "init_val: " << init_val << ", last_val " << last_val << std::endl;
	  if(init_val == 1.0 && last_val == 0.0) return EXIT_SUCCESS;
	  else return EXIT_FAILURE;
}


