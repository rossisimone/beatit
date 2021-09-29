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

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void assemble_poisson(EquationSystems & es,
                      const std::string & system_name);

// Function prototype for the exact solution.
Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.)
{
    return 0;
}

namespace Parameters
{
static int Nf = 1;
static double gamma_noise = 0.5;
static double s = 15.0;
static double L_fibers = 10.0;
static double PHI_fibers = 10.0*M_PI;

static double sigma = 1.0;
static double gamma = 1.0;
static double beta = 0.0;
static double kappa = 1.0;
static double L = 4.0;
static double Ly = 4.0;
static libMesh::ExplicitSystem * g_system = nullptr;
static libMesh::ExplicitSystem * f_system = nullptr;
static libMesh::PointLocatorTree *  locator = nullptr;

static double grid_min_x = 0;
static double grid_min_y = 0;
static double grid_min_z = 0;
static double grid_max_x = 1;
static double grid_max_y = 1;
static double grid_max_z = 1;
static int grid_el_x = 4;
static int grid_el_y = 4;
static int grid_el_z = 0;

bool use_white_noise = false;
bool use_fiber_noise = false;

}

static std::string mesh_file = "NONE";

// NOISE
static std::random_device rd;
static std::mt19937 mt(rd());
static std::uniform_real_distribution<double> uniform_dist(0.0, 2*M_PI);
static std::normal_distribution<double> gauss_dist(0.0, 1.0);


int main (int argc, char ** argv)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  GetPot commandLine ( argc, argv );
  std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
  GetPot data(datafile_name);


  Parameters::grid_min_x = data("grid_min_x",0.0);
  Parameters::grid_min_y = data("grid_min_y",0.0);
  Parameters::grid_min_z = data("grid_min_z",0.0);
  Parameters::grid_max_x = data("grid_max_x",1.0);
  Parameters::grid_max_y = data("grid_max_y",1.0);
  Parameters::grid_max_z = data("grid_max_z",1.0);
  Parameters::grid_el_x = data("grid_el_x",4);
  Parameters::grid_el_y = data("grid_el_y",4);
  Parameters::grid_el_z = data("grid_el_z",4);
  Parameters::use_white_noise = data("use_white_noise",false);
  Parameters::use_fiber_noise = data("use_fiber_noise",false);

  Parameters::gamma_noise = data("gamma_noise",0.5);
  Parameters::Nf = data("Nf",1);
  Parameters::s = data("s",15);
  Parameters::L_fibers = data("L_fibers",10.0);
  Parameters::PHI_fibers = data("PHI_fibers",10.0*M_PI);

  Parameters::sigma = data("sigma",1.0);
  Parameters::gamma = data("gamma",1.0);
  Parameters::beta = data("beta",0.0);
  Parameters::kappa = data("kappa",1.0);
  int N  = data("nel",15);
  Parameters::L  = data("L", 4.0);

  Parameters::Ly  = data("Ly", -1.0);
  if(Parameters::Ly < 0) Parameters::Ly = Parameters::L;
  // Brief message to the user regarding the program name
  // and command line arguments.
  libMesh::out << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << "L: " << Parameters::L << ", Ly: " << Parameters::Ly << std::endl << std::endl;

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh lattice(init.comm());
  double L_lattice =ceil(Parameters::L);

  if(Parameters::grid_el_z <= 0)
  {
      MeshTools::Generation::build_square (lattice,
              Parameters::grid_el_x, Parameters::grid_el_y,
              Parameters::grid_min_x, Parameters::grid_max_x,
              Parameters::grid_min_y, Parameters::grid_max_y,
                                           QUAD4);
  }
  else
  {
      MeshTools::Generation::build_cube (lattice,
              Parameters::grid_el_x, Parameters::grid_el_y, Parameters::grid_el_z ,
              Parameters::grid_min_x, Parameters::grid_max_x,
              Parameters::grid_min_y, Parameters::grid_max_y,
              Parameters::grid_min_z, Parameters::grid_max_z,
                                           HEX8);
  }

  Parameters::locator = new libMesh::PointLocatorTree(lattice);

  EquationSystems lattice_equation_systems(lattice);
  ExplicitSystem & g_system = lattice_equation_systems.add_system<ExplicitSystem>("Perlin");
  Parameters::g_system = &g_system;
  // Add the variables "u" & "v" to "Stokes".  They
  // will be approximated using FIRST-order approximation.
  g_system.add_variable("ux", FIRST);
  g_system.add_variable("uy", FIRST);
  g_system.add_variable("uz", FIRST);
  lattice_equation_systems.init();
  std::vector<dof_id_type> dof_indices;
  // define lattice gradients;
  for (auto & node : lattice.local_node_ptr_range())
  {
      g_system.get_dof_map().dof_indices(node, dof_indices);
      double angle1 = uniform_dist(mt);
      double angle2 = 0.0;
      if(Parameters::grid_el_z > 0 ) angle2 = uniform_dist(mt);


      double ux = std::cos(angle1) * std::cos(angle2);
      double uy = std::sin(angle1) * std::cos(angle2);
      double uz = std::sin(angle2);

      libMesh::Point u(ux,uy, uz);
      double unorm = u.norm();
      u = u / u.norm();

      g_system.solution->set(dof_indices[0], u(0));
      g_system.solution->set(dof_indices[1], u(1));
      g_system.solution->set(dof_indices[2], u(2));
  }

  libMesh::ExodusII_IO exporter(lattice);
  exporter.write_equation_systems("lattice.e",
          lattice_equation_systems);
  exporter.write_element_data(lattice_equation_systems);

  // Use the MeshTools::Generation mesh generator to create a uniform
  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
  // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
  // elements instead of the default QUAD4's we used in example 2
  // allow us to use higher-order approximation.
  mesh_file = data("mesh_file","NONE");
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
  double mesh_max_x = data("grid_max_x",1.0);
  double mesh_max_y = data("grid_max_y",1.0);
  double mesh_max_z = data("grid_max_z",1.0);
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

  // Declare the Poisson system and its variables.
  // The Poisson system is another example of a steady system.
  equation_systems.add_system<LinearImplicitSystem> ("Poisson");

  // Adds the variable "u" to "Poisson".  "u"
  // will be approximated using second-order approximation.
  equation_systems.get_system("Poisson").add_variable("u", FIRST);

  // Give the system a pointer to the matrix assembly
  // function.  This will be called when needed by the
  // library.
  equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);


  // Declare the Poisson system and its variables.
  // The Poisson system is another example of a steady system.
  equation_systems.add_system<LinearImplicitSystem> ("Noise");
  // Adds the variable "u" to "Poisson".  "u"
  // will be approximated using second-order approximation.
  equation_systems.get_system("Noise").add_variable("w", CONSTANT, MONOMIAL);
  equation_systems.get_system("Noise").add_variable("p", CONSTANT, MONOMIAL);
  equation_systems.get_system("Noise").add_variable("F", CONSTANT, MONOMIAL);


  ExplicitSystem & f_system = equation_systems.add_system<ExplicitSystem>("fibers");
  Parameters::f_system = &f_system;
  // Add the variables "u" & "v" to "Stokes".  They
  // will be approximated using FIRST-order approximation.
  f_system.add_variable("fx", CONSTANT, MONOMIAL);
  f_system.add_variable("fy", CONSTANT, MONOMIAL);
  f_system.add_variable("fz", CONSTANT, MONOMIAL);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  if("NONE" != mesh_file)
  {
      importer.copy_elemental_solution (f_system, "fx", "fibersx");
      importer.copy_elemental_solution (f_system, "fy", "fibersy");
      importer.copy_elemental_solution (f_system, "fz", "fibersz");
  }

  int nrefs = data("nr" , 0);
  if ( nrefs > 0 )
  {
      std::cout << "refining: " << nrefs << " times" << std::endl;
      for(int nr = 0; nr < nrefs; nr++)
      {
      libMesh::MeshRefinement refine(mesh);
      std::cout << "refining uniformly #" << nr + 1 << std::endl;
      refine.uniformly_refine();
      std::cout << "mesh prepare for use"<< std::endl;
      mesh.prepare_for_use(false);
      std::cout << "reinit systems" << std::endl;
      equation_systems.reinit();
      }
  }
  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Solve the system "Poisson".  Note that calling this
  // member will assemble the linear system and invoke
  // the default numerical solver.  With PETSc the solver can be
  // controlled from the command line.  For example,
  // you can invoke conjugate gradient with:
  //
  // ./introduction_ex3 -ksp_type cg
  //
  // You can also get a nice X-window that monitors the solver
  // convergence with:
  //
  // ./introduction-ex3 -ksp_xmonitor
  //
  // if you linked against the appropriate X libraries when you
  // built PETSc.
  equation_systems.get_system("Poisson").solve();

  double u_max =   equation_systems.get_system("Poisson").solution->max();
  double u_min =   equation_systems.get_system("Poisson").solution->min();
  double du = u_max - u_min;
  equation_systems.get_system("Poisson").solution->add(-u_min);
  (*equation_systems.get_system("Poisson").solution) /= du;
	
#if defined(LIBMESH_HAVE_VTK) && !defined(LIBMESH_ENABLE_PARMESH)

  // After solving the system write the solution
  // to a VTK-formatted plot file.
//  VTKIO (mesh).write_equation_systems ("out.pvtu", equation_systems);
  std::string output_file = data("output","out.e");
  ExodusII_IO exporter_lattice(mesh);
  exporter_lattice.write_equation_systems(output_file,
          equation_systems);
  exporter_lattice.write_element_data(equation_systems);

#endif // #ifdef LIBMESH_HAVE_VTK

  delete Parameters::locator;
  // All done.
  return 0;
}



// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_poisson(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Poisson");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Poisson");
  LinearImplicitSystem & w_system = es.get_system<LinearImplicitSystem> ("Noise");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();
  const DofMap & w_dof_map = w_system.get_dof_map();
  const DofMap & f_dof_map = Parameters::f_system->get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.  Introduction Example 4
  // describes some advantages of  std::unique_ptr's in the context of
  // quadrature rules.
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, SECOND);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim-1, FIFTH);

  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point> & q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> w_dof_indices;
  std::vector<dof_id_type> g_dof_indices;
  std::vector<dof_id_type> f_dof_indices;

  // find size
  double hx = (Parameters::grid_max_x - Parameters::grid_min_x) / Parameters::grid_el_x;
  double hy = (Parameters::grid_max_y - Parameters::grid_min_y) / Parameters::grid_el_y;
  double hz = (Parameters::grid_max_z - Parameters::grid_min_z) / Parameters::grid_el_z;
  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element iterators are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // iterator el will iterate from the first to the last element on
  // the local processor.  The iterator end_el tells us when to stop.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      w_dof_map.dof_indices (elem, w_dof_indices);
      f_dof_map.dof_indices (elem, f_dof_indices);

      // white noise
      double W =  gauss_dist(mt);
      w_system.solution->set(w_dof_indices[0], W);
      // Perlin noise
      double P = 0.;
      double FP = 0.0;
      libMesh::VectorValue<double> f0;

      {
          libMesh::Point x = elem->centroid();


          const Elem* mapped_element = (*Parameters::locator)(x);


          P = 0.5;
          if(mapped_element)
          {
              for( auto & node : mapped_element->node_ref_range())
              {
                  Parameters::g_system->get_dof_map().dof_indices(&node, g_dof_indices);
                  libMesh::Point xi(node);
                  libMesh::Point gi( (*Parameters::g_system->solution)(g_dof_indices[0]),
                                    (*Parameters::g_system->solution)(g_dof_indices[1]),
                                    (*Parameters::g_system->solution)(g_dof_indices[2]) );
                  double gc = 1.0;
                  double PNf = 0.0;
                  for(int cf = 0; cf < Parameters::Nf; ++cf)
                  {
                  libMesh::Point ri(x-xi);
                  double dX = 1-std::abs(ri(0)) / hx;
                  double dY = 1-std::abs(ri(1)) / hy;
                  double dZ = 1-std::abs(ri(2)) / hz;

                  double aX = ((6*dX - 15)*dX + 10)*dX*dX*dX;
                  double aY = ((6*dY - 15)*dY + 10)*dY*dY*dY;
                  double aZ = ((6*dZ - 15)*dZ + 10)*dZ*dZ*dZ;


                  //if(Parameters::grid_el_z <= 0) P += std::sqrt(2.0) * 0.5 * aX * aY * gi.contract(ri);
                  //else P += 1.0 / std::sqrt(3.0) * aX * aY * aZ* gi.contract(ri);
                  if(Parameters::grid_el_z <= 0) PNf += gc * aX * aY * gi.contract(ri);
                  else PNf += gc * aX * aY * aZ* gi.contract(ri);
                  gc *= Parameters::gamma_noise;
                  }
                  P += (1.0 - Parameters::gamma_noise) / (1.0 - std::pow(Parameters::gamma_noise, Parameters::Nf) ) /std::sqrt(dim) * PNf;
              }
          }
          w_system.solution->set(w_dof_indices[1], P);

          double Fcos = 0.5 + 0.5 * std::cos(2.0*M_PI/Parameters::L_fibers * x(1) + Parameters::PHI_fibers * P);
          FP = std::pow(Fcos, Parameters::s);
          w_system.solution->set(w_dof_indices[2], FP);


          if(mesh_file == "NONE")
          {
              // fibers
              x = elem->centroid();
              f0(0) = 1.0;
              f0(1) = 0.0;
              f0(2) = 0.0;
              //f0(0) = -x(1);
              //f0(1) = x(0);
              //f0(2) = 0.0;
              if(f0.norm() > 0) f0 /= f0.norm();
              else
              {
                  f0(0) = 1.0;
                  f0(1) = 0.0;
                  f0(2) = 0.0;
              }
              Parameters::f_system->solution->set(f_dof_indices[0], f0(0));
              Parameters::f_system->solution->set(f_dof_indices[1], f0(1));
              Parameters::f_system->solution->set(f_dof_indices[2], f0(2));
          }
      }

      // Cache the number of degrees of freedom on this element, for
      // use as a loop bound later.  We use cast_int to explicitly
      // convert from size() (which may be 64-bit) to unsigned int
      // (which may be 32-bit but which is definitely enough to count
      // *local* degrees of freedom.
      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // With one variable, we should have the same number of degrees
      // of freedom as shape functions.
      libmesh_assert_equal_to (n_dofs, phi.size());

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).

      // The  DenseMatrix::resize() and the  DenseVector::resize()
      // members will automatically zero out the matrix  and vector.
      Ke.resize (n_dofs, n_dofs);

      Fe.resize (n_dofs);

      f0(0) = (*Parameters::f_system->solution)(f_dof_indices[0]);
      f0(1) = (*Parameters::f_system->solution)(f_dof_indices[1]);
      f0(2) = (*Parameters::f_system->solution)(f_dof_indices[2]);


      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          const Real x = q_point[qp](0);
          const Real y = q_point[qp](1);

          libMesh::TensorValue<double> H;
          libMesh::TensorValue<double> I;

//          f0(0) = -y;
//          f0(1) = x;
//          f0(2) = 0.0;
//          f0 /= f0.norm();
          for(int i = 0; i < 3; ++i)
          {
              I(i,i) = 1.0;
              for(int j = 0; j < 3; ++j)
              {
                  H(i,j) = Parameters::gamma * I(i,j) + Parameters::beta * f0(i) * f0(j);
              }
          }

          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test functions (i) against
          // the trial functions (j).
          for (unsigned int i=0; i != n_dofs; i++)
            for (unsigned int j=0; j != n_dofs; j++)
              {
                Ke(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
                Ke(i,j) += JxW[qp]*( dphi[i][qp] * ( H * dphi[j][qp] ) );
              }

          // This is the end of the matrix summation loop
          // Now we build the element right-hand-side contribution.
          // This involves a single loop in which we integrate the
          // "forcing function" in the PDE against the test functions.
          {


            Real fxy = Parameters::sigma * P;
            if(Parameters::use_white_noise) fxy = Parameters::sigma * W;
            else if(Parameters::use_fiber_noise) fxy = Parameters::sigma * FP;
            for (unsigned int i=0; i != n_dofs; i++)
              Fe(i) += JxW[qp]*fxy*phi[i][qp];
          }
        }

      // We have now reached the end of the RHS summation,
      // and the end of quadrature point loop, so
      // the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions.
      //
      // There are several ways Dirichlet boundary conditions
      // can be imposed.  A simple approach, which works for
      // interpolary bases like the standard Lagrange polynomials,
      // is to assign function values to the
      // degrees of freedom living on the domain boundary. This
      // works well for interpolary bases, but is more difficult
      // when non-interpolary (e.g Legendre or Hierarchic) bases
      // are used.
      //
      // Dirichlet boundary conditions can also be imposed with a
      // "penalty" method.  In this case essentially the L2 projection
      // of the boundary values are added to the matrix. The
      // projection is multiplied by some large factor so that, in
      // floating point arithmetic, the existing (smaller) entries
      // in the matrix and right-hand-side are effectively ignored.
      //
      // This amounts to adding a term of the form (in latex notation)
      //
      // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
      //
      // where
      //
      // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
      {

        // The following loop is over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (auto side : elem->side_index_range())
          if (elem->neighbor_ptr(side) == nullptr)
            {
              // The value of the shape functions at the quadrature
              // points.
              const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              // The XYZ locations (in physical space) of the
              // quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector<Point> & qface_point = fe_face->get_xyz();

              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, side);

              // Some shape functions will be 0 on the face, but for
              // ease of indexing and generality of code we loop over
              // them anyway
              libmesh_assert_equal_to (n_dofs, phi_face.size());

              // Loop over the face quadrature points for integration.
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // The location on the boundary of the current
                  // face quadrature point.
                  const Real xf = qface_point[qp](0);
                  const Real yf = qface_point[qp](1);

                  // The penalty value.  \frac{1}{\epsilon}
                  // in the discussion above.
                  const Real penalty = 0.e10;

                  // The boundary value.
                  const Real value = exact_solution(xf, yf);

                  // Matrix contribution of the L2 projection.
                  for (unsigned int i=0; i != n_dofs; i++)
                    for (unsigned int j=0; j != n_dofs; j++)
                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

                  // Right-hand-side contribution of the L2
                  // projection.
                  for (unsigned int i=0; i != n_dofs; i++)
                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                }
            }
      }

      // We have now finished the quadrature point loop,
      // and have therefore applied all the boundary conditions.

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The  SparseMatrix::add_matrix()
      // and  NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

  // All done!
}
