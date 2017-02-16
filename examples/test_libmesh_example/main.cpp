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
  equation_systems.add_system<LinearImplicitSystem> ("DG");

  // Adds the variable "u" to "DG".  "u"
  // will be approximated using first-order approximation.
  equation_systems.get_system("DG").add_variable("u", FIRST, L2_LAGRANGE);
  equation_systems.init();
  // Prints information about the system to the screen.
  equation_systems.print_info();

  //INITIALIZE u VECTOR:
  libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  std::vector<libMesh::dof_id_type> dof_indices;
  libMesh::DenseVector<libMesh::Number> Fe;
  const libMesh::DofMap & dof_map = equation_systems.get_system("DG").get_dof_map();
  libMesh::FEType fe_disp = dof_map.variable_type(0);
  const unsigned int dim = mesh.mesh_dimension();

  for (; el != end_el; ++el)
  {
      const libMesh::Elem * elem = *el;
      auto elID = elem->id();

      dof_map.dof_indices(elem, dof_indices);
      Fe.resize( dof_indices.size() );
      Fe(0) = dof_indices[0];
      Fe(1) = dof_indices[1];
      Fe(2) = dof_indices[2];

      equation_systems.get_system("DG").solution->insert(Fe, dof_indices);
  }
  equation_systems.get_system("DG").solution->close();
  equation_systems.get_system("DG").solution->print(std::cout);

  std::cout << "Exporting timestep 0 ... " << std::flush;
  ExodusII_IO exp(mesh);
  exp._timestep = 1;
  exp.write_discontinuous_exodusII("dg.e-s.0000", equation_systems);
  VTKIO vtk0(mesh);
  vtk0.write_equation_systems("dg.0.vtk", equation_systems);
  std::cout << " done. " << std::endl;

  el = mesh.active_local_elements_begin();
  for (; el != end_el; ++el)
  {
      const libMesh::Elem * elem = *el;
      auto elID = elem->id();

      dof_map.dof_indices(elem, dof_indices);
      Fe.resize( dof_indices.size() );
      Fe(0) = dof_indices[0]+1;
      Fe(1) = dof_indices[1]+1;
      Fe(2) = dof_indices[2]+1;

      equation_systems.get_system("DG").solution->insert(Fe, dof_indices);
  }
  equation_systems.get_system("DG").solution->close();
  std::cout << "Exporting timestep 1  ... " << std::flush;
  ExodusII_IO exp2(mesh);
  exp2._timestep = 2;
  exp2.write_discontinuous_exodusII("dg.e-s.0001", equation_systems);
  vtk0.write_equation_systems("dg.1.vtk", equation_systems);

  std::cout << " done. " << std::endl;

  el = mesh.active_local_elements_begin();
  for (; el != end_el; ++el)
  {
      const libMesh::Elem * elem = *el;
      auto elID = elem->id();

      dof_map.dof_indices(elem, dof_indices);
      Fe.resize( dof_indices.size() );
      Fe(0) = dof_indices[0]+2;
      Fe(1) = dof_indices[1]+2;
      Fe(2) = dof_indices[2]+2;

      equation_systems.get_system("DG").solution->insert(Fe, dof_indices);
  }
  equation_systems.get_system("DG").solution->close();
  std::cout << "Exporting timestep 2 ... " << std::flush;
  ExodusII_IO exp3(mesh);
  exp3._timestep = 3;
  exp3.write_discontinuous_exodusII("dg.e-s.0002", equation_systems);
  vtk0.write_equation_systems("dg.2.vtk", equation_systems);

  std::cout << " done. " << std::endl;
  return 0;
}



