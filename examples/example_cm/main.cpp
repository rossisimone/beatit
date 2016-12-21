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

// Matrix and right-hand side assemble
void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name);
void assemble_mass(EquationSystems & es,
                         const std::string & system_name);


// Define the elasticity tensor, which is a fourth-order tensor
// i.e. it has four indices i, j, k, l
Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                           unsigned int l, double nu , double E);

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries
  LibMeshInit init (argc, argv);

  // Initialize the cantilever mesh
  const unsigned int dim = 2;

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D support");

  GetPot commandLine ( argc, argv );
  std::string datafile_name = commandLine.follow ( "data.pot", 2, "-i", "--input" );
  GetPot data(datafile_name);
  // Create a 2D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  Mesh mesh2(init.comm(), dim);

  double nu = data("nu", 0.3);
  double E = data("E", 5e5);
  int elX = data("elX", 4);
  std::string elTypeName = data("elType", "TRI6");
  std::map<std::string, ElemType> orderMap;
  orderMap["TRI3"] = TRI3;
  orderMap["QUAD4"] = QUAD4;
  orderMap["TRI6"] = TRI6;
  orderMap["QUAD9"] = QUAD9;
  auto elType = orderMap.find(elTypeName)->second;
  auto order = FIRST;
  int elY = elX;
  if(elType == TRI6 || elType == QUAD9 ) order = SECOND;
  auto elType2 = TRI3;
  if( elType == QUAD9 ) elType2 = QUAD4;

  MeshTools::Generation::build_square (mesh,
                                       elX, elY,
                                       0., 1.,
                                       0., 0.2,
                                       elType);

  MeshTools::Generation::build_square (mesh2,
                                       elX, elY,
                                       0., 1.,
                                       0., 0.2,
                                       elType2);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  EquationSystems equation_systems2 (mesh2);
   auto & export_system_p =
    equation_systems2.add_system<ExplicitSystem> ("Pressure_Projection");
    export_system_p.add_variable("p",FIRST, LAGRANGE);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Add two displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("ux", order, LAGRANGE);
  unsigned int v_var = system.add_variable("uy", order, LAGRANGE);

  system.attach_assemble_function (assemble_elasticity);
  std::cout << "Add pressure system ... " << std::endl;

    LinearImplicitSystem & system_p =
    equation_systems.add_system<LinearImplicitSystem> ("Pressure_Projection");
    system_p.add_variable("p",FIRST, LAGRANGE);
  std::cout << "Add mass assembly ... " << std::endl;
  system_p.attach_assemble_function (assemble_mass);
  // Construct a Dirichlet boundary condition object
  // We impose a "clamped" boundary condition on the
  // "left" boundary, i.e. bc_id = 3
  std::cout << "Add BC... " << std::endl;
  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(3);

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

  std::cout << "Ready for init ... " << std::endl;
  // Initialize the data structures for the equation system.
  system.get_linear_solver()->set_preconditioner_type( LU_PRECOND );
  system_p.get_linear_solver()->set_preconditioner_type( SOR_PRECOND );
  equation_systems.init();
  equation_systems2.init();
  std::cout << "Ready print info  ... " << std::endl;
   equation_systems.parameters.set<double>("nu") = nu;
   equation_systems.parameters.set<double>("E") = E;

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  std::cout << "Solve system  ... " << std::endl;
  system.solve();
   system.solution->print(std::cout);
 system.rhs->print(std::cout);

  std::cout << "Export pressure  ... " << std::endl;
  system_p.solve();

  auto& p_proj = equation_systems.get_system<LinearImplicitSystem>("Pressure_Projection").solution;
  auto& p_proj_exp = equation_systems2.get_system<ExplicitSystem>("Pressure_Projection").solution;

  auto first = p_proj->first_local_index();
  auto last = p_proj->last_local_index();
  for(int i = first; i < last; ++i)
  {
	  p_proj_exp->set(i, (*p_proj)(i));
  }
  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems("displacement.e", equation_systems);
  ExodusII_IO (mesh2).write_equation_systems("pressure.e", equation_systems2);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}


void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name)
{
  libmesh_assert_equal_to (system_name, "Elasticity");

  const MeshBase & mesh = es.get_mesh();

  double nu = es.parameters.get<double>("nu");
  double E = es.parameters.get<double>("E");
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
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kuu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=0, C_k=1;

                C_j=0, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kuv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=1, C_k=0;

                C_j=0, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kvu(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
                // Tensor indices
                unsigned int C_i, C_j, C_k, C_l;
                C_i=1, C_k=1;

                C_j=0, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=0;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=0, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));

                C_j=1, C_l=1;
                Kvv(i,j) += JxW[qp]*(eval_elasticity_tensor(C_i, C_j, C_k, C_l, nu, E) * dphi[i][qp](C_j)*dphi[j][qp](C_l));
              }
        }
//      std::cout << "\n Ke: " << std::endl;
//      Ke.print(std::cout);
      {
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor(side) == libmesh_nullptr)
            {
              const std::vector<std::vector<Real> > & phi_face = fe_face->get_phi();
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              fe_face->reinit(elem, side);

              if (mesh.get_boundary_info().has_boundary_id (elem, side, 1)) // Apply a traction on the right side
                {
                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                    for (unsigned int i=0; i<n_v_dofs; i++)
                      Fv(i) += JxW_face[qp] * (-1.) * phi_face[i][qp];
                }
            }
      }

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
  system.matrix->close();
  	 system.matrix->print(std::cout);

}


void assemble_mass(EquationSystems & es,
                         const std::string & system_name)
{
 libmesh_assert_equal_to (system_name, "Elasticity");

  const MeshBase & mesh = es.get_mesh();
  double nu = es.parameters.get<double>("nu");
  double E = es.parameters.get<double>("E");

  double kappa = E / (3*(1-2*nu));
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");
  LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem>("Pressure_Projection");

  const unsigned int u_var = system.variable_number ("ux");
  const unsigned int v_var = system.variable_number ("uy");
 const unsigned int p_var = system_p.variable_number ("p");

  const DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, SECOND);
  fe->attach_quadrature_rule (&qrule);

  const DofMap & dof_map_p = system_p.get_dof_map();
  FEType fe_type_p = dof_map_p.variable_type(0);
  UniquePtr<FEBase> fe_p (FEBase::build(dim, fe_type_p));
  QGauss qrule_p (dim,SECOND);
  fe_p->attach_quadrature_rule (&qrule_p);


  const std::vector<Real> & JxW = fe_p->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe_p->get_dphi();
  const std::vector<std::vector<Real> > & phi = fe_p->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi_u = fe->get_dphi();
  const std::vector<std::vector<Real> > & phi_u = fe->get_phi();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;


  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_p;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map_p.dof_indices (elem, dof_indices_p, p_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      fe->reinit (elem);
      fe_p->reinit (elem);

         Ke.resize(dof_indices_p.size(), dof_indices_p.size());
          Fe.resize(dof_indices_p.size());
          Real du = 0.0;
          Real dv = 0.0;

//          std::cout << "Assembling Poisson ... " << std::endl;
         for (unsigned int qp = 0; qp < qrule_p.n_points(); qp++)
         {
//        	 std::cout << "\ngetting div " << std::endl;
             for (unsigned int i = 0; i < phi_u.size(); i++)
             {
            	 du += dphi_u[i][qp](0) * (*system.solution)(dof_indices_u[i]);
            	 dv += dphi_u[i][qp](1) * (*system.solution)(dof_indices_v[i]);
             }
//        	 std::cout << "Done " << std::endl;

             Real kdivu =kappa *  (du + dv);

             for (unsigned int i = 0; i < phi.size(); i++)
             {
                 Fe(i) +=  JxW[qp] * kdivu* phi[i][qp];

                 for (unsigned int j = 0; j < phi.size(); j++)
                 {
                     // mass term
                     Ke(i, i) += JxW[qp] * phi[i][qp] * phi[j][qp];
                 }
             }
         }
//          std::cout << "Add elemental contributions ... " << std::endl;
               system_p.matrix->add_matrix (Ke, dof_indices_p);
               system_p.rhs->add_vector    (Fe, dof_indices_p);

    }
}
Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l, double nu , double E)
{
  // Define the Poisson ratio

  // Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
  const Real lambda_1 = E*nu / ((1. + nu) * (1. - 2.*nu));
  const Real lambda_2 = E*0.5 / (1 + nu);

  // Define the Kronecker delta functions that we need here
  Real delta_ij = (i == j) ? 1. : 0.;
  Real delta_il = (i == l) ? 1. : 0.;
  Real delta_ik = (i == k) ? 1. : 0.;
  Real delta_jl = (j == l) ? 1. : 0.;
  Real delta_jk = (j == k) ? 1. : 0.;
  Real delta_kl = (k == l) ? 1. : 0.;
  double c = 1.0;
  return c *lambda_1 * delta_ij * delta_kl + lambda_2 * (delta_ik * delta_jl + delta_il * delta_jk);
}
