// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// <h1>Introduction Example 3 - Solving a Poisson Problem</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This is the third example program.  It builds on
// the second example program by showing how to solve a simple
// Poisson system.  This example also introduces the notion
// of customized matrix assembly functions, working with an
// exact solution, and using element iterators.
// We will not comment on things that
// were already explained in the second example.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>


// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/vtk_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/enum_solver_package.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include <libmesh/point_locator_tree.h>

#include "Util/Noise.hpp"
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void setup_mesh_subdomains(Mesh& mesh, GetPot& data)
{
    std::cout << "setup_mesh_subdomains: " << std::endl;
    double x_san_interface = data("x_san_interface", 0.45);
    double x_atria_interface = data("x_atria_interface", 0.55);
    double y_interface = data("y_interface", 0.5);

    for(auto & elem : mesh.active_element_ptr_range() )
    {
        Point c = elem->centroid();
        double x = c(0);
        double y = c(1);
        if( y > y_interface ) elem->subdomain_id() = 4;
        else
        {
            if(x < x_san_interface) elem->subdomain_id() = 1;
            else if(x >= x_atria_interface) elem->subdomain_id() = 2;
            else elem->subdomain_id() = 3;
        }

    }
}
void setup_fibrosis_subdomains(EquationSystems& es, GetPot& data)
{
    std::cout << "setup_fibrosis_subdomains: " << std::endl;

    auto& noise_system = es.get_system<libMesh::ExplicitSystem>("Noise");
    auto& dof_map = noise_system.get_dof_map();
    std::vector<libMesh::dof_id_type> dof_indices;
    double threshold = data("threshold", 0.5);

    for(auto & elem : es.get_mesh().active_element_ptr_range() )
    {
        dof_map.dof_indices(elem, dof_indices);
        double noise = (*noise_system.solution)(dof_indices[0]);
        if(noise > threshold && elem->subdomain_id() == 1)
        {
            elem->subdomain_id() = 3;
        }

    }
}


int main (int argc, char ** argv)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init (argc, argv);


  GetPot commandLine ( argc, argv );
  std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
  GetPot data(datafile_name);

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
  // elements instead of the default QUAD4's we used in example 2
  // allow us to use higher-order approximation.
  std::string mesh_file = data("mesh_file","NONE");
  std::cout << std::endl << mesh_file << std::endl;

  Mesh mesh(init.comm());

  libMesh::ExodusII_IO importer(mesh);

  if(mesh_file != "NONE")
  {
      std::cout << "Importing mesh" << std::endl;
      importer.read(mesh_file);
      mesh.prepare_for_use();
  }
  else
  {
      std::cout << "Creating square" << std::endl;
      double mesh_min_x = data("mesh_min_x",-1.0);
      double mesh_min_y = data("mesh_min_y",-1.0);
      double mesh_min_z = data("mesh_min_z",-1.0);
      double mesh_max_x = data("mesh_max_x",1.0);
      double mesh_max_y = data("mesh_max_y",1.0);
      double mesh_max_z = data("mesh_max_z",1.0);
      int mesh_el_x = data("mesh_el_x",64);
      int mesh_el_y = data("mesh_el_y",64);
      int mesh_el_z = data("mesh_el_z",0);

      if(mesh_el_z <= 0)
      {
          MeshTools::Generation::build_square (mesh,
                  mesh_el_x, mesh_el_y,
                  mesh_min_x, mesh_max_x,
                  mesh_min_y, mesh_max_y,
                   TRI3);
      }
      else
      {
          MeshTools::Generation::build_cube (mesh,
                  mesh_el_x, mesh_el_y, mesh_el_z,
                  mesh_min_x, mesh_max_x,
                  mesh_min_y, mesh_max_y,
                  mesh_min_z, mesh_max_z,
                   TET4);
      }
  }


  // Print information about the mesh to the screen.
  // Note that 5x5 QUAD9 elements actually has 11x11 nodes,
  // so this mesh is significantly larger than the one in example 2.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  Noise noise;
  noise.setup(data, "");
  noise.generate_noise_field(equation_systems);

  setup_mesh_subdomains(mesh, data);
  setup_fibrosis_subdomains(equation_systems, data);

  std::cout << "Create EP ..." << std::endl;

  std::string ep_data_filename = data("ep_data","ep.beat");
  GetPot data_ep(ep_data_filename);
  // Create the TimeData object
  BeatIt::TimeData datatime;
  // Set it up using the input file
  datatime.setup(data, "");
  // Output on screen the stored variables
  datatime.print();
  // Read the model from the input file
  std::string model = data_ep("model", "monowave");
  std::cout << "Create Electrophysiology model ..." << std::endl;
  // Create the model
  BeatIt::ElectroSolver* solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, equation_systems);
  // Setup the EP model using the input file
  std::cout << "Calling setup..." << std::endl;
  solver->setup(data_ep, model);
  // Initialize systems
  std::cout << "Calling init ..." << std::endl;
  equation_systems.print_info();
  // Set up initial conditions at time
  solver->init(datatime.M_startTime);
  std::cout << "Assembling matrices" << std::endl;
  solver->assemble_matrices(datatime.M_dt);
  // output file counter
  int save_iter = 0;
  // Export initial condition at time
//    solver->save_exo_timestep(save_iter, datatime.M_time);
  solver->save_potential(save_iter, datatime.M_startTime);
  // Parameters to save the activation times
  // A node is activated if the transmembrane potential > threshold
  double threshold = data("threshold", -10.0);
  // We export the activation times every at_save_iter iterations
  int at_save_iter = data("at_save_iter", 25);

  // Old parameters that define the method, not to be changed
  // TODO: clean this part
  std::string system_mass = data_ep(model + "/diffusion_mass", "mass");
  std::string iion_mass = data_ep(model + "/reaction_mass", "lumped_mass");
  bool useMidpointMethod = false;
  int step0 = 0;
  int step1 = 1;

  // Start loop in time
  std::cout << "Time loop starts:" << std::endl;
  // Control the time loop using the TimeData object
  for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
  {
      // We are doing a new iteration
      // let's update first the information in the TimeData object
      datatime.advance();
      // Ouput to screen the current time and the current iteration
      std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter << std::endl;
      // Advance the solution in time: u_n <- u_n+1
      solver->advance();
      // Solve ionic model and evaluate ionic currents
      solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
      // Solve monodomain model
      solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
      // Update the activation times
      solver->update_activation_time(datatime.M_time, threshold);
      //Export the solution if at the right timestep
      if (0 == datatime.M_iter % datatime.M_saveIter)
      {
          // update output file counter
          save_iter++;
          // export current solution
          solver->save_potential(save_iter, datatime.M_time);
//            solver->save_exo_timestep(save_iter, datatime.M_time);
      }
      // export the activation times if at the corresponding timestep
      if (0 == datatime.M_iter % (at_save_iter * datatime.M_saveIter))
      {
          // export activation times
          solver->save_activation_times(save_iter);
      }

  }

  // // export activation times again
  solver->save_activation_times(save_iter);
  // delete solver before ending the simulation
  // avoiding memory leaks
  delete solver;

//  std::string output_file = data("output","out.e");
//  ExodusII_IO exporter_lattice(mesh);
//  exporter_lattice.write_equation_systems(output_file,
//          equation_systems);
//  exporter_lattice.write_element_data(equation_systems);

  // All done.
  return 0;
}
