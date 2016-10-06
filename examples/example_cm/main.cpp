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



// <h1>Vector Finite Element Example 1 - Solving an uncoupled Poisson Problem</h1>
// \author Paul Bauman
// \date 2012
//
// This is the first vector FE example program.  It builds on
// the introduction_ex3 example program by showing how to solve a simple
// uncoupled Poisson system using vector Lagrange elements.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io_helper.h"
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




// Matrix and right-hand side assemble
void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name);


void assemble_pressure(EquationSystems & es,
                         const std::string & system_name);

// Define the elasticity tensor, which is a fourth-order tensor
// i.e. it has four indices i, j, k, l
Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l);

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries
  LibMeshInit init (argc, argv);
    GetPot data("data.pot");

  // Initialize the cantilever mesh
  const unsigned int dim = 2;

  // Skip this 2D example if libMesh was compiled as 1D-only.
//  libmesh_example_requires(dim <= LIBMESH_DIM, "2D support");

  // Create a 2D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
//  MeshTools::Generation::build_square (mesh,
//                                        80,  16,
//                                        0.,  1.,
//                                        0.,  0.2,
//                                       TRI6);
////
//

      // We may need XDR support compiled in to read binary .xdr files
    std::string meshfile = data("input_mesh_name", "cm_t6_split2_m1.e");

    mesh.read(&meshfile[0]);
    // Print information about the mesh to the screen.
    mesh.print_info();
    std::cout << "Reading mesh Done!" << std::endl;
  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Add two displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("ux", SECOND, LAGRANGE);
  unsigned int v_var = system.add_variable("uy",  SECOND, LAGRANGE);

  system.attach_assemble_function (assemble_elasticity);

    LinearImplicitSystem & system_p =
    equation_systems.add_system<LinearImplicitSystem> ("pressure");
  unsigned int p_var = system_p.add_variable("p",  SECOND, LAGRANGE);
  system_p.attach_assemble_function (assemble_pressure);
  // Construct a Dirichlet boundary condition object
  // We impose a "clamped" boundary condition on the
  // "left" boundary, i.e. bc_id = 3
  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(4);

  // Create a vector storing the variable numbers which the BC applies to
  std::vector<unsigned int> variables(2);
  variables[0] = u_var; variables[1] = v_var;

  // Create a ZeroFunction to initialize dirichlet_bc
  ZeroFunction<> zf;

  DirichletBoundary dirichlet_bc(boundary_ids,
                                 variables,
                                 &zf);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  // Initialize the data structures for the equation system.
    std::cout << "Init systems ... ";
  equation_systems.init();
  std::cout << " Done.  " << std::endl;

  // Print information about the system to the screen.
  //equation_systems.print_info();

  // Solve the system
  std::cout << "Solving system ... ";
  system.solve();
  system_p.solve();
  std::cout << " Done.  " << std::endl;

  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems("displacement.e", equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}


void assemble_pressure(EquationSystems & es,
                         const std::string & system_name)
{

	  // Define the Poisson ratio
   double nu = 0.3;

  // Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
//  const Real lambda_1 = nu / ((1. + nu) * (1. - 2.*nu));
//  const Real lambda_2 = 0.5 / (1 + nu);

   double E = 250.0;
   double mu = E / 3.0 / (1.-nu);
   double kappa = E / 2.0 / (1.- 2*nu);
   kappa = 1.0;
  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");
  const unsigned int ux_var = system.variable_number ("ux");
  const unsigned int uy_var = system.variable_number ("uy");
  const DofMap & dof_map = system.get_dof_map();


  LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem>("pressure");

  const unsigned int p_var = system_p.variable_number ("p");
    const DofMap & dof_map_p = system_p.get_dof_map();
  FEType fe_type_p = dof_map_p.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type_p));
  QGauss qrule (dim, fe_type_p.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<libMesh::Real> & JxW = fe->get_JxW();

    // The physical XY locations of the quadrature points on the element.
    // These might be useful for evaluating spatially varying material
    // properties at the quadrature points.
    const std::vector<libMesh::Point> & q_point = fe->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<libMesh::Real> > & phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<libMesh::RealGradient> > & dphi =
            fe->get_dphi();
      // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseVector<libMesh::Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_p;
    std::vector<libMesh::dof_id_type> dof_indices_ux;
    std::vector<libMesh::dof_id_type> dof_indices_uy;
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
    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();
    // Loop over the elements.  Note that  ++el is preferred to
    // el++ since the latter requires an unnecessary temporary
    // object.
    for (; el != end_el; ++el)
    {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const libMesh::Elem * elem = *el;

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_ux, ux_var);
        dof_map.dof_indices(elem, dof_indices_uy, uy_var);
        dof_map_p.dof_indices(elem, dof_indices_p);
                // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        Me.resize(dof_indices_p.size(), dof_indices_p.size());

        Fe.resize(dof_indices_p.size());


        libMesh::RealGradient dxu1(0.0,0.0);
        libMesh::RealGradient dxu2(0.0,0.0);

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for(unsigned int i = 0; i < phi.size(); ++i )
            {
                dxu1 += dphi[i][qp] *  (*system.solution)(dof_indices_ux[i]);
                dxu2 += dphi[i][qp] *  (*system.solution)(dof_indices_uy[i]);
            }
            double div_u = dxu1(0) + dxu2(1);


            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test funcions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i < phi.size(); i++)
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    // Mass term
                    Me(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                }

				for (unsigned int i = 0; i < phi.size(); i++)
                {
                    // Mass term
                    Fe(i) += JxW[qp] * kappa * div_u * phi[i][qp];
                }

        }

        system_p.matrix->add_matrix(Me, dof_indices_p);
        system_p.rhs->add_vector(Fe, dof_indices_p);
    }

}

void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name)
{
  libmesh_assert_equal_to (system_name, "Elasticity");

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

  const unsigned int u_var = system.variable_number ("ux");
  const unsigned int v_var = system.variable_number ("uy");

  const DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke),
    Kvu(Ke), Kvv(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe);

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();

      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);

      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=0, C_k=0;

                C_j=0, C_l=0;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=0, C_k=1;

                C_j=0, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=1, C_k=0;

                C_j=0, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=1, C_k=1;

                C_j=0, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
        }

      {
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor(side) == libmesh_nullptr)
            {
              const std::vector<std::vector<Real> > & phi_face = fe_face->get_phi();
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              fe_face->reinit(elem, side);

              if (mesh.get_boundary_info().has_boundary_id (elem, side, 2)) // Apply a traction on the right side
                {
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    for (unsigned int i=0; i<n_v_dofs; i++)
                      Fv(i) += JxW_face[qp] * (6.25) * phi_face[i][qp];
                }
            }
      }

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
}

Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l)
{
  // Define the Poisson ratio
   double nu = 0.49995;

  // Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
//  const Real lambda_1 = nu / ((1. + nu) * (1. - 2.*nu));
//  const Real lambda_2 = 0.5 / (1 + nu);

   double E = 250.0;
   double mu = E / 3.0 / (1.-nu);
   double kappa = E / 2.0 / (1.- 2*nu);
  // Define the Kronecker delta functions that we need here
  Real delta_ij = (i == j) ? 1. : 0.;
  Real delta_il = (i == l) ? 1. : 0.;
  Real delta_ik = (i == k) ? 1. : 0.;
  Real delta_jl = (j == l) ? 1. : 0.;
  Real delta_jk = (j == k) ? 1. : 0.;
  Real delta_kl = (k == l) ? 1. : 0.;

  return  kappa * delta_ij * delta_kl + mu * (delta_ik * delta_jl + delta_il * delta_jk - 2.0 / 3.0  * delta_ij * delta_kl);
}


//int main (int argc, char ** argv)
//{
//  // Initialize libraries.
//  LibMeshInit init (argc, argv);
//
//  // Brief message to the user regarding the program name
//  // and command line arguments.
//  libMesh::out << "Running " << argv[0];
//
//  for (int i=1; i<argc; i++)
//    libMesh::out << " " << argv[i];
//
//  libMesh::out << std::endl << std::endl;
//
//  // Skip this 2D example if libMesh was compiled as 1D-only.
//  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");
//
//  // Create a mesh, with dimension to be overridden later, on the
//  // default MPI communicator.
//  Mesh mesh(init.comm());
//
//  // Use the MeshTools::Generation mesh generator to create a uniform
//  // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
//  // to build a mesh of 15x15 QUAD9 elements.
//  MeshTools::Generation::build_square (mesh,
//                                       15, 15,
//                                       -1., 1.,
//                                       -1., 1.,
//                                       QUAD9);
//
//  // Print information about the mesh to the screen.
//  mesh.print_info();
//
//  // Create an equation systems object.
//  EquationSystems equation_systems (mesh);
//
//  // Declare the Poisson system and its variables.
//  // The Poisson system is another example of a steady system.
//  equation_systems.add_system<LinearImplicitSystem> ("Poisson");
//
//  // Adds the variable "u" to "Poisson".  "u"
//  // will be approximated using second-order approximation
//  // using vector Lagrange elements. Since the mesh is 2-D, "u" will
//  // have two components.
//  equation_systems.get_system("Poisson").add_variable("u", SECOND, libMesh::LAGRANGE_VEC);
//
//  // Give the system a pointer to the matrix assembly
//  // function.  This will be called when needed by the
//  // library.
//  equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
//
//  // Initialize the data structures for the equation system.
//  equation_systems.init();
//
//  // Prints information about the system to the screen.
//  equation_systems.print_info();
//
//  // Solve the system "Poisson".  Note that calling this
//  // member will assemble the linear system and invoke
//  // the default numerical solver.  With PETSc the solver can be
//  // controlled from the command line.  For example,
//  // you can invoke conjugate gradient with:
//  //
//  // ./vector_fe_ex1 -ksp_type cg
//  //
//  // You can also get a nice X-window that monitors the solver
//  // convergence with:
//  //
//  // ./vector_fe_ex1 -ksp_xmonitor
//  //
//  // if you linked against the appropriate X libraries when you
//  // built PETSc.
//  equation_systems.get_system("Poisson").solve();
//
//#ifdef LIBMESH_HAVE_EXODUS_API
//  ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);
//#endif
//
//#ifdef LIBMESH_HAVE_GMV
//  GMVIO(mesh).write_equation_systems("out.gmv", equation_systems);
//#endif
//
//  // All done.
//  return 0;
//}
//
//
//
//// We now define the matrix assembly function for the
//// Poisson system.  We need to first compute element
//// matrices and right-hand sides, and then take into
//// account the boundary conditions, which will be handled
//// via a penalty method.
//void assemble_poisson(EquationSystems & es,
//                      const std::string & system_name)
//{
//
//  // It is a good idea to make sure we are assembling
//  // the proper system.
//  libmesh_assert_equal_to (system_name, "Poisson");
//
//
//  // Get a constant reference to the mesh object.
//  const MeshBase & mesh = es.get_mesh();
//
//  // The dimension that we are running
//  const unsigned int dim = mesh.mesh_dimension();
//
//  // Get a reference to the LinearImplicitSystem we are solving
//  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Poisson");
//
//  // A reference to the  DofMap object for this system.  The  DofMap
//  // object handles the index translation from node and element numbers
//  // to degree of freedom numbers.  We will talk more about the  DofMap
//  // in future examples.
//  const DofMap & dof_map = system.get_dof_map();
//
//  // Get a constant reference to the Finite Element type
//  // for the first (and only) variable in the system.
//  FEType fe_type = dof_map.variable_type(0);
//
//  // Build a Finite Element object of the specified type.
//  // Note that FEVectorBase is a typedef for the templated FE
//  // class.
//  UniquePtr<FEVectorBase> fe (FEVectorBase::build(dim, fe_type));
//
//  // A 5th order Gauss quadrature rule for numerical integration.
//  QGauss qrule (dim, FIFTH);
//
//  // Tell the finite element object to use our quadrature rule.
//  fe->attach_quadrature_rule (&qrule);
//
//  // Declare a special finite element object for
//  // boundary integration.
//  UniquePtr<FEVectorBase> fe_face (FEVectorBase::build(dim, fe_type));
//
//  // Boundary integration requires one quadraure rule,
//  // with dimensionality one less than the dimensionality
//  // of the element.
//  QGauss qface(dim-1, FIFTH);
//
//  // Tell the finite element object to use our
//  // quadrature rule.
//  fe_face->attach_quadrature_rule (&qface);
//
//  // Here we define some references to cell-specific data that
//  // will be used to assemble the linear system.
//  //
//  // The element Jacobian * quadrature weight at each integration point.
//  const std::vector<Real> & JxW = fe->get_JxW();
//
//  // The physical XY locations of the quadrature points on the element.
//  // These might be useful for evaluating spatially varying material
//  // properties at the quadrature points.
//  const std::vector<Point> & q_point = fe->get_xyz();
//
//  // The element shape functions evaluated at the quadrature points.
//  // Notice the shape functions are a vector rather than a scalar.
//  const std::vector<std::vector<RealGradient> > & phi = fe->get_phi();
//
//  // The element shape function gradients evaluated at the quadrature
//  // points. Notice that the shape function gradients are a tensor.
//  const std::vector<std::vector<RealTensor> > & dphi = fe->get_dphi();
//
//  // Define data structures to contain the element matrix
//  // and right-hand-side vector contribution.  Following
//  // basic finite element terminology we will denote these
//  // "Ke" and "Fe".  These datatypes are templated on
//  //  Number, which allows the same code to work for real
//  // or complex numbers.
//  DenseMatrix<Number> Ke;
//  DenseVector<Number> Fe;
//
//  // This vector will hold the degree of freedom indices for
//  // the element.  These define where in the global system
//  // the element degrees of freedom get mapped.
//  std::vector<dof_id_type> dof_indices;
//
//  // Now we will loop over all the elements in the mesh.
//  // We will compute the element matrix and right-hand-side
//  // contribution.
//  //
//  // Element iterators are a nice way to iterate through all the
//  // elements, or all the elements that have some property.  The
//  // iterator el will iterate from the first to the last element on
//  // the local processor.  The iterator end_el tells us when to stop.
//  // It is smart to make this one const so that we don't accidentally
//  // mess it up!  In case users later modify this program to include
//  // refinement, we will be safe and will only consider the active
//  // elements; hence we use a variant of the active_elem_iterator.
//  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
//  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//
//  // Loop over the elements.  Note that  ++el is preferred to
//  // el++ since the latter requires an unnecessary temporary
//  // object.
//  for ( ; el != end_el ; ++el)
//    {
//      // Store a pointer to the element we are currently
//      // working on.  This allows for nicer syntax later.
//      const Elem * elem = *el;
//
//      // Get the degree of freedom indices for the
//      // current element.  These define where in the global
//      // matrix and right-hand-side this element will
//      // contribute to.
//      dof_map.dof_indices (elem, dof_indices);
//
//      // Compute the element-specific data for the current
//      // element.  This involves computing the location of the
//      // quadrature points (q_point) and the shape functions
//      // (phi, dphi) for the current element.
//      fe->reinit (elem);
//
//      // Zero the element matrix and right-hand side before
//      // summing them.  We use the resize member here because
//      // the number of degrees of freedom might have changed from
//      // the last element.  Note that this will be the case if the
//      // element type is different (i.e. the last element was a
//      // triangle, now we are on a quadrilateral).
//
//      // The  DenseMatrix::resize() and the  DenseVector::resize()
//      // members will automatically zero out the matrix  and vector.
//      Ke.resize (dof_indices.size(),
//                 dof_indices.size());
//
//      Fe.resize (dof_indices.size());
//
//      // Now loop over the quadrature points.  This handles
//      // the numeric integration.
//      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//        {
//          // Now we will build the element matrix.  This involves
//          // a double loop to integrate the test funcions (i) against
//          // the trial functions (j).
//          for (unsigned int i=0; i<phi.size(); i++)
//            for (unsigned int j=0; j<phi.size(); j++)
//              Ke(i,j) += JxW[qp] * dphi[i][qp].contract(dphi[j][qp]);
//
//          // This is the end of the matrix summation loop
//          // Now we build the element right-hand-side contribution.
//          // This involves a single loop in which we integrate the
//          // "forcing function" in the PDE against the test functions.
//          {
//            const Real x = q_point[qp](0);
//            const Real y = q_point[qp](1);
//            const Real eps = 1.e-3;
//
//            // "f" is the forcing function for the Poisson equation.
//            // In this case we set f to be a finite difference
//            // Laplacian approximation to the (known) exact solution.
//            //
//            // We will use the second-order accurate FD Laplacian
//            // approximation, which in 2D is
//            //
//            // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
//            //                u(i-1,j) + u(i+1,j) +
//            //                -4*u(i,j))/h^2
//            //
//            // Since the value of the forcing function depends only
//            // on the location of the quadrature point (q_point[qp])
//            // we will compute it here, outside of the i-loop
//            const Real fx = -(exact_solution(0, x, y-eps) +
//                              exact_solution(0, x, y+eps) +
//                              exact_solution(0, x-eps, y) +
//                              exact_solution(0, x+eps, y) -
//                              4.*exact_solution(0, x, y))/eps/eps;
//
//            const Real fy = -(exact_solution(1, x, y-eps) +
//                              exact_solution(1, x, y+eps) +
//                              exact_solution(1, x-eps, y) +
//                              exact_solution(1, x+eps, y) -
//                              4.*exact_solution(1, x, y))/eps/eps;
//
//            const RealGradient f(fx, fy);
//
//            for (unsigned int i=0; i<phi.size(); i++)
//              Fe(i) += JxW[qp]*f*phi[i][qp];
//          }
//        }
//
//      // We have now reached the end of the RHS summation,
//      // and the end of quadrature point loop, so
//      // the interior element integration has
//      // been completed.  However, we have not yet addressed
//      // boundary conditions.  For this example we will only
//      // consider simple Dirichlet boundary conditions.
//      //
//      // There are several ways Dirichlet boundary conditions
//      // can be imposed.  A simple approach, which works for
//      // interpolary bases like the standard Lagrange polynomials,
//      // is to assign function values to the
//      // degrees of freedom living on the domain boundary. This
//      // works well for interpolary bases, but is more difficult
//      // when non-interpolary (e.g Legendre or Hierarchic) bases
//      // are used.
//      //
//      // Dirichlet boundary conditions can also be imposed with a
//      // "penalty" method.  In this case essentially the L2 projection
//      // of the boundary values are added to the matrix. The
//      // projection is multiplied by some large factor so that, in
//      // floating point arithmetic, the existing (smaller) entries
//      // in the matrix and right-hand-side are effectively ignored.
//      //
//      // This amounts to adding a term of the form (in latex notation)
//      //
//      // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
//      //
//      // where
//      //
//      // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
//      {
//        // The following loop is over the sides of the element.
//        // If the element has no neighbor on a side then that
//        // side MUST live on a boundary of the domain.
//        for (unsigned int side=0; side<elem->n_sides(); side++)
//          if (elem->neighbor(side) == libmesh_nullptr)
//            {
//              // The value of the shape functions at the quadrature
//              // points.
//              const std::vector<std::vector<RealGradient> > & phi_face = fe_face->get_phi();
//
//              // The Jacobian * Quadrature Weight at the quadrature
//              // points on the face.
//              const std::vector<Real> & JxW_face = fe_face->get_JxW();
//
//              // The XYZ locations (in physical space) of the
//              // quadrature points on the face.  This is where
//              // we will interpolate the boundary value function.
//              const std::vector<Point> & qface_point = fe_face->get_xyz();
//
//              // Compute the shape function values on the element
//              // face.
//              fe_face->reinit(elem, side);
//
//              // Loop over the face quadrature points for integration.
//              for (unsigned int qp=0; qp<qface.n_points(); qp++)
//                {
//                  // The location on the boundary of the current
//                  // face quadrature point.
//                  const Real xf = qface_point[qp](0);
//                  const Real yf = qface_point[qp](1);
//
//                  // The penalty value.  \frac{1}{\epsilon}
//                  // in the discussion above.
//                  const Real penalty = 1.e10;
//
//                  // The boundary values.
//                  const RealGradient f(exact_solution(0, xf, yf),
//                                       exact_solution(1, xf, yf));
//
//                  // Matrix contribution of the L2 projection.
//                  for (unsigned int i=0; i<phi_face.size(); i++)
//                    for (unsigned int j=0; j<phi_face.size(); j++)
//                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
//
//                  // Right-hand-side contribution of the L2
//                  // projection.
//                  for (unsigned int i=0; i<phi_face.size(); i++)
//                    Fe(i) += JxW_face[qp]*penalty*f*phi_face[i][qp];
//                }
//            }
//      }
//
//      // We have now finished the quadrature point loop,
//      // and have therefore applied all the boundary conditions.
//
//      // If this assembly program were to be used on an adaptive mesh,
//      // we would have to apply any hanging node constraint equations
//      //dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
//
//      // The element matrix and right-hand-side are now built
//      // for this element.  Add them to the global matrix and
//      // right-hand-side vector.  The  SparseMatrix::add_matrix()
//      // and  NumericVector::add_vector() members do this for us.
//      system.matrix->add_matrix (Ke, dof_indices);
//      system.rhs->add_vector    (Fe, dof_indices);
//    }
//
//  // All done!
//}
//
//
///**
// * This is the exact solution that
// * we are trying to obtain.  We will solve
// *
// * - (u_xx + u_yy) = f
// *
// * and take a finite difference approximation using this
// * function to get f.  This is the well-known "method of
// * manufactured solutions".
// */
//Real exact_solution (const int component,
//                     const Real x,
//                     const Real y,
//                     const Real z = 0.)
//{
//  static const Real pi = acos(-1.);
//
//  switch (component)
//    {
//    case 0:
//      return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
//    case 1:
//      return sin(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z);
//    case 2:
//      return sin(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z)*cos(.5*pi*x*y*z);
//    default:
//      libmesh_error_msg("Invalid component = " << component);
//    }
//
//  // dummy
//  return 0.0;
//}
