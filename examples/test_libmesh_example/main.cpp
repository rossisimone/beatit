// C++ include files that we need
#include <iostream>

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


int main (int argc, char ** argv)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init (argc, argv);
  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());
  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.
  MeshTools::Generation::build_square (mesh,
                                        1,   1,
                                       -1., 1.,
                                       -1., 1.,
                                       TRI3);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  // Initialize the data structures for the equation system.
  equation_systems.add_system<LinearImplicitSystem> ("sys");

  // Adds the variable "u" to "DG".  "u"
  // will be approximated using first-order approximation.
  equation_systems.get_system("sys").add_variable("u", FIRST);
  equation_systems.get_system("sys").add_variable("v", FIRST);
  equation_systems.get_system("sys").add_variable("p", FIRST);

  equation_systems.add_system<LinearImplicitSystem> ("disp");
  equation_systems.get_system("disp").add_variable("uu", FIRST);
  equation_systems.get_system("disp").add_variable("vv", FIRST);

  // Prints information about the system to the screen.
  equation_systems.init();
  equation_systems.print_info();

  equation_systems.get_system("sys").solution->print(std::cout);
  auto it = mesh.active_nodes_begin();
  auto it_end = mesh.active_nodes_end();
  auto sys_num = equation_systems.get_system("sys").number();
  for (; it != it_end; it++)
  {
      auto * node = *it;
      auto n_vars = node->n_vars(sys_num);
      unsigned int c = 0;
      for(auto v = 0; v != n_vars; ++v)
      {
          auto dof_id = node -> dof_number (sys_num, v, c);
          double value = v+1;
          equation_systems.get_system("sys").solution->set(dof_id, value);
      }
  }
  equation_systems.get_system("sys").solution->close();
  equation_systems.get_system("sys").solution->print(std::cout);

  equation_systems.get_system("disp").solution->print(std::cout);


  it = mesh.active_nodes_begin();
  auto disp_sys_num = equation_systems.get_system("disp").number();
  for (; it != it_end; it++)
  {
      auto * node = *it;
      auto n_vars = node->n_vars(sys_num);
      auto disp_n_vars = node->n_vars(disp_sys_num);
      //std::cout << "n_vars: " << n_vars << ", disp_n_vars: " << disp_n_vars << std::endl;
      for(auto v = 0; v != disp_n_vars; ++v)
      {
          auto nc_disp = node->n_comp(disp_sys_num, v);
          auto nc = node->n_comp(sys_num, v);
          //std::cout << "nc disp: " << nc_disp << ", nc: " << nc << std::endl;
          for(auto c = 0; c != nc_disp; ++c)
          {
              auto dof_id_disp = node -> dof_number (disp_sys_num, v, c);
              auto dof_id = node -> dof_number (sys_num, v, c);
              //std::cout << "dof_id_disp: " << dof_id_disp << ", dof_id: " << dof_id << std::endl;
              double value = (*equation_systems.get_system("sys").solution)(dof_id);
              //std::cout << "value: " << value << std::endl;
              equation_systems.get_system("disp").solution->set(dof_id_disp, value);
          }
      }
  }
  equation_systems.get_system("disp").solution->print(std::cout);

  return 0;
}



