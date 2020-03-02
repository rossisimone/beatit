// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// <h1>Systems Example 1 - Stokes Equations</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example shows how a simple, linear system of equations
// can be solved in parallel.  The system of equations are the familiar
// Stokes equations for low-speed incompressible fluid flow.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include "libmesh/getpot.h"

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/petsc_linear_solver.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "Util/IO/io.hpp"
// The definition of a geometric element
#include "libmesh/elem.h"
#include "libmesh/zero_function.h"
#include "libmesh/mesh_refinement.h"
// Bring in everything from the libMesh namespace
using namespace libMesh;

std::set<boundary_id_type> inflows;
std::set<boundary_id_type> outflows;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_poisson(EquationSystems &es, const std::string &system_name);

void zf(DenseVector<Number> &output, const Point &p, const Real)
{
    output(0) = 0.0;
}
void zfv(DenseVector<Number> &output, const Point &p, const Real)
{
    output(0) = 0.0;
    output(1) = 0.0;
    output(2) = 0.0;
}
// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_stokes(EquationSystems &es, const std::string &system_name);

// The main program.
int main(int argc, char **argv)
{
    // Initialize libMesh.
    LibMeshInit init(argc, argv);
    GetPot data = BeatIt::readInputFile(argc, argv);
    // This example requires a linear solver package.
//    libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE, "--enable-petsc, --enable-trilinos, or --enable-eigen");

// This example NaNs with the Eigen sparse linear solvers and
// Trilinos solvers, but should work OK with either PETSc or
// Laspack.
//    libmesh_example_requires(libMesh::default_solver_package() != EIGEN_SOLVERS, "--enable-petsc or --enable-laspack");
//    libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc or --enable-laspack");

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    ParallelMesh mesh(init.comm());

    // Use the MeshTools::Generation mesh generator to create a uniform
    // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
    // to build a mesh of 8x8 Quad9 elements.  Building these
    // higher-order elements allows us to use higher-order
    // approximation, as in example 3.
//  MeshTools::Generation::build_square (mesh,
//                                       15, 15,
//                                       0., 1.,
//                                       0., 1.,
//                                       TRI3);

    std::string input_mesh = data("input", "NONE");
    std::cout << "Reading mesh" << std::endl;
    int n_refinements = data("nrefs", 0);
    BeatIt::serial_mesh_partition(init.comm(), input_mesh, &mesh, n_refinements);

//	mesh.read(input_mesh);
    std::cout << "Scaling mesh" << std::endl;
    double scale = data("scale", 1.0);
    MeshTools::Modification::scale(mesh, scale, scale, scale);
    std::cout << "Refining mesh" << std::endl;
//	MeshRefinement refinement(mesh);
//	refinement.uniformly_refine();
    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    std::cout << "Creating equation systems" << std::endl;
    EquationSystems equation_systems(mesh);

    std::cout << "Creating Poisson system" << std::endl;
    LinearImplicitSystem &system_profile = equation_systems.add_system < LinearImplicitSystem > ("Poisson");
    system_profile.add_variable("phi", FIRST);
    system_profile.attach_assemble_function(assemble_poisson);
    std::set < libMesh::boundary_id_type > dirichlet_poisson;
    // Set Dirichlet condition on lumen
    // HH16: sideset 14
    std::cout << "Adding Dirichlet BC for parabolic velocity profile" << std::endl;
    dirichlet_poisson.insert(14);
    std::vector<unsigned int> vars_p(1);
    vars_p[0] = 0;
    AnalyticFunction<> bcz(zf);
    libMesh::DirichletBoundary dirichlet_bc_poisson(dirichlet_poisson, vars_p, bcz);
    system_profile.get_dof_map().add_dirichlet_boundary(dirichlet_bc_poisson);

    // Declare the system and its variables.
    // Create a transient system named "Stokes"
    std::cout << "Creating Stokes system" << std::endl;
    LinearImplicitSystem &system = equation_systems.add_system < LinearImplicitSystem > ("Stokes");

    // Add the variables "u" & "v" to "Stokes".  They
    // will be approximated using second-order approximation.
    unsigned int ux = system.add_variable("ux", FIRST);
    unsigned int uy = system.add_variable("uy", FIRST);
    unsigned int uz = system.add_variable("uz", FIRST);

    // Add the variable "p" to "Stokes". This will
    // be approximated with a first-order basis,
    // providing an LBB-stable pressure-velocity pair.
    unsigned int p = system.add_variable("p", FIRST);

    // Give the system a pointer to the matrix assembly
    // function.
    system.attach_assemble_function(assemble_stokes);

    std::cout << "Adding Dirichlet BC for parabolic velocity profile" << std::endl;
    std::string in_bc = data("in_bc", "NONE");
    BeatIt::readList(in_bc, inflows);
    std::string out_bc = data("out_bc", "NONE");
    BeatIt::readList(out_bc, outflows);

    std::set < libMesh::boundary_id_type > dirichlet_stokes;

    std::string bcs = data("Dbc", "NONE");
    std::set < libMesh::boundary_id_type > bcs_vec;
    BeatIt::readList(bcs, dirichlet_stokes);
    std::vector<unsigned int> vars_vel_stokes(3);
    vars_vel_stokes[0] = ux;
    vars_vel_stokes[1] = uy;
    vars_vel_stokes[2] = uz;
    AnalyticFunction<> bczv(zfv);
    libMesh::ZeroFunction < Number > zero;
    libMesh::DirichletBoundary dirichlet_bc_stokes(dirichlet_stokes, vars_vel_stokes, &zero);
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc_stokes);

    // Initialize the data structures for the equation system.
    std::cout << "Initialize equation systems" << std::endl;
    equation_systems.init();

//	equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 250;
//	equation_systems.parameters.set < Real > ("linear solver tolerance") = TOLERANCE;

    // Prints information about the system to the screen.
    equation_systems.print_info();

    // Assemble & solve the linear system,
    // then write the solution.
    std::cout << "Solve Poisson system" << std::endl;
    typedef libMesh::PetscLinearSolver<libMesh::Number> PetscSolver;
    PetscSolver * linear_solver = dynamic_cast<PetscSolver*>(system_profile.get_linear_solver());
    KSPSetOptionsPrefix(linear_solver->ksp(),"poisson_");
    PCSetOptionsPrefix(linear_solver->pc(),"poisson_");
    KSPSetFromOptions(linear_solver->ksp());
    equation_systems.get_system("Poisson").solve();
    std::cout << "Solve Stokes system" << std::endl;
    equation_systems.get_system("Stokes").solve();

    std::string output = data("output", "out.e");
    std::cout << "Output" << std::endl;
#ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write_equation_systems (output,
            equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API

    // All done.
    std::cout << "Good luck with your simulation!" << std::endl;
    return 0;
}

void assemble_stokes(EquationSystems &es, const std::string& libmesh_dbg_var(system_name))
{
    std::cout << "Start Stokes Assembly" << std::endl;

    double scale = 1.0e-4;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Stokes");

    // Get a constant reference to the mesh object.
    const MeshBase &mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    LinearImplicitSystem &system = es.get_system < LinearImplicitSystem > ("Stokes");
    LinearImplicitSystem &system_phi = es.get_system < LinearImplicitSystem > ("Poisson");

    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = system.variable_number("ux");
    const unsigned int v_var = system.variable_number("uy");
    const unsigned int w_var = system.variable_number("uz");
    const unsigned int p_var = system.variable_number("p");

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_vel_type = system.variable_type(u_var);

    // Get the Finite Element type for "p".
    FEType fe_pres_type = system.variable_type(p_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    std::unique_ptr < FEBase > fe_vel(FEBase::build(dim, fe_vel_type));

    // Build a Finite Element object of the specified type for
    // the pressure variables.
    std::unique_ptr < FEBase > fe_pres(FEBase::build(dim, fe_pres_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    QGauss qrule(dim, fe_vel_type.default_quadrature_order());

    // Tell the finite element objects to use our quadrature rule.
    fe_vel->attach_quadrature_rule(&qrule);
    fe_pres->attach_quadrature_rule(&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe_vel->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

    // The element shape functions for the pressure variable
    // evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &psi = fe_pres->get_phi();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap &dof_map = system.get_dof_map();
    const DofMap &dof_map_phi = system_phi.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    DenseSubMatrix<Number> Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke), Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke), Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke), Kpu(Ke), Kpv(
            Ke), Kpw(Ke), Kpp(Ke);

    DenseSubVector<Number> Fu(Fe), Fv(Fe), Fw(Fe), Fp(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_u;
    std::vector < dof_id_type > dof_indices_v;
    std::vector < dof_id_type > dof_indices_w;
    std::vector < dof_id_type > dof_indices_p;
    std::vector < dof_id_type > dof_indices_phi;
    std::vector<double> phi_k;

    // BC stuff
    libMesh::UniquePtr < libMesh::FEBase > fe_face(libMesh::FEBase::build(dim, fe_vel_type));
    libMesh::QGauss qface(dim - 1, libMesh::FIRST);
    fe_face->attach_quadrature_rule(&qface);
    const std::vector<libMesh::Real> &JxW_face = fe_face->get_JxW();
    const std::vector<std::vector<libMesh::Real> > &phi_face = fe_face->get_phi();
    int n_phi = phi_face.size();
    const std::vector<libMesh::Point> &qface_point = fe_face->get_xyz();
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_face = fe_face->get_dphi();
    const std::vector<libMesh::Point> &normals = fe_face->get_normals();

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the active_elem_iterator.
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_u, u_var);
        dof_map.dof_indices(elem, dof_indices_v, v_var);
        dof_map.dof_indices(elem, dof_indices_w, w_var);
        dof_map.dof_indices(elem, dof_indices_p, p_var);
        dof_map_phi.dof_indices(elem, dof_indices_phi);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        const unsigned int n_w_dofs = dof_indices_w.size();
        const unsigned int n_p_dofs = dof_indices_p.size();

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_vel->reinit(elem);
        fe_pres->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);
        // Reposition the submatrices...  The idea is this:
        //
        //         -           -          -  -
        //        | Kuu Kuv Kup |        | Fu |
        //   Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
        //        | Kpu Kpv Kpp |        | Fp |
        //         -           -          -  -
        //
        // The DenseSubMatrix.reposition () member takes the
        // (row_offset, column_offset, row_size, column_size).
        //
        // Similarly, the DenseSubVector.reposition () member
        // takes the (row_offset, row_size)
        Kuu.reposition(u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition(u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
        Kup.reposition(u_var * n_u_dofs, p_var * n_u_dofs, n_u_dofs, n_p_dofs);

        Kvu.reposition(v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition(v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition(v_var * n_v_dofs, p_var * n_v_dofs, n_v_dofs, n_p_dofs);

        Kww.reposition(w_var * n_w_dofs, w_var * n_w_dofs, n_w_dofs, n_w_dofs);
        Kwp.reposition(w_var * n_w_dofs, p_var * n_w_dofs, n_w_dofs, n_p_dofs);

        Kpu.reposition(p_var * n_u_dofs, u_var * n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition(p_var * n_u_dofs, v_var * n_u_dofs, n_p_dofs, n_v_dofs);
        Kpw.reposition(p_var * n_u_dofs, w_var * n_u_dofs, n_p_dofs, n_w_dofs);
        Kpp.reposition(p_var * n_u_dofs, p_var * n_u_dofs, n_p_dofs, n_p_dofs);

        Fu.reposition(u_var * n_u_dofs, n_u_dofs);
        Fv.reposition(v_var * n_u_dofs, n_v_dofs);
        Fw.reposition(w_var * n_u_dofs, n_w_dofs);
        Fp.reposition(p_var * n_u_dofs, n_p_dofs);

        // Now we will build the element matrix.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            // Assemble the u-velocity row
            // uu coupling
            for (unsigned int i = 0; i < n_u_dofs; i++)
            {
                for (unsigned int j = 0; j < n_u_dofs; j++)
                {
                    Kuu(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    Kvv(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    Kww(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                }
            }

            // up coupling
            for (unsigned int i = 0; i < n_u_dofs; i++)
            {
                for (unsigned int j = 0; j < n_p_dofs; j++)
                {
                    Kup(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](0);
                    Kvp(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](1);
                    Kwp(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](2);
                }
            }

            // Assemble the pressure row
            // pu coupling
            for (unsigned int i = 0; i < n_p_dofs; i++)
            {
                for (unsigned int j = 0; j < n_u_dofs; j++)
                {
                    Kpu(i, j) += JxW[qp] * psi[i][qp] * dphi[j][qp](0);
                    Kpv(i, j) += JxW[qp] * psi[i][qp] * dphi[j][qp](1);
                    Kpw(i, j) += JxW[qp] * psi[i][qp] * dphi[j][qp](2);
                }
            }

            // pp coupling
            double h = elem->hmax();

            for (unsigned int i = 0; i < n_p_dofs; i++)
            {
                for (unsigned int j = 0; j < n_p_dofs; j++)
                {
                    Kpp(i, j) += 0.5 * h * h * JxW[qp] * dphi[i][qp] * dphi[j][qp];
                }
            }

        } // end of the quadrature point qp-loop

        // At this point the interior element integration has
        // been completed.  However, we have not yet addressed
        // boundary conditions.  For this example we will only
        // consider simple Dirichlet boundary conditions imposed
        // via the penalty method. The penalty method used here
        // is equivalent (for Lagrange basis functions) to lumping
        // the matrix resulting from the L2 projection penalty
        // approach introduced in example 3.

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.

        // loop over sides
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            // If on mesh boundary
            if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
//                    std::unique_ptr<const Elem> side_elem_ptr(elem->build_side_ptr(s));
                const unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, side);

                // if on inflow or lumen
                //if (boundary_id == 1 || boundary_id == 2 || boundary_id == 3 || boundary_id == 4 || boundary_id == 14)
                int sign = -1;

                if (outflows.find(boundary_id) != outflows.end())
                    sign *= -1;

                if (inflows.find(boundary_id) != inflows.end() || outflows.find(boundary_id) != outflows.end())
                //if ( boundary_id == 5 )
//				if ( boundary_id < 14 )
                {
                    fe_face->reinit(elem, side);
                    // Loop over face qp
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        const double xq = qface_point[qp](0);
                        const double yq = qface_point[qp](1);
                        const double zq = qface_point[qp](2);

                        double phi_value = 0.0;
                        // Evaluate parabolic profile
                        for (unsigned int l = 0; l < phi_face.size(); l++)
                        {
                            phi_value += phi_face[l][qp] * (*system_phi.current_local_solution)(dof_indices_phi[l]);
                        }

//						double inflow_magnitude = 1e3*phi_value;
                        double inflow_magnitude = 1e-2;
                        double ux = sign * inflow_magnitude * normals[qp](0);
                        double uy = sign * inflow_magnitude * normals[qp](1);
                        double uz = sign * inflow_magnitude * normals[qp](2);
                        double penalty = 1e10 * inflow_magnitude;
                        // Assemble RHS BC
                        for (unsigned int l = 0; l < phi_face.size(); l++)
                        {
                            Fu(l) += JxW_face[qp] * penalty * ux * phi_face[l][qp];
                            Fv(l) += JxW_face[qp] * penalty * uy * phi_face[l][qp];
                            Fw(l) += JxW_face[qp] * penalty * uz * phi_face[l][qp];
                        } // Assemble RHS BC

                        // Assemble Matrix BC
                        for (unsigned int l = 0; l < phi_face.size(); l++)
                        {
                            for (unsigned int k = 0; k < phi_face.size(); k++)
                            {
                                Kuu(l, k) += JxW_face[qp] * penalty * phi_face[l][qp] * phi_face[k][qp];
                                Kvv(l, k) += JxW_face[qp] * penalty * phi_face[l][qp] * phi_face[k][qp];
                                Kww(l, k) += JxW_face[qp] * penalty * phi_face[l][qp] * phi_face[k][qp];
                            }
                        } // Assemble Matrix BC

                    } // Loop over face qp
                } // if on inflow or lumen
                  //else if ( boundary_id == 1 )
                  //{
//
                //              }
//				else if  ( boundary_id == 14 )
//				{
//	                    fe_face->reinit(elem, side);
//	                    // Loop over face qp
//	                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
//	                    {
//	                        const double xq = qface_point[qp](0);
//	                        const double yq = qface_point[qp](1);
//	                        const double zq = qface_point[qp](2);
//	                        double penalty = 1e5;
//
//	                        // Assemble Matrix BC
//	                        for (unsigned int l = 0; l < phi_face.size(); l++)
//	                        {
//	                            for (unsigned int k = 0; k < phi_face.size(); k++)
//	                            {
//	                                Kuu(l, k) += JxW_face[qp] * penalty * phi_face[l][qp] * phi_face[k][qp];
//	                                Kvv(l, k) += JxW_face[qp] * penalty * phi_face[l][qp] * phi_face[k][qp];
//	                                Kww(l, k) += JxW_face[qp] * penalty * phi_face[l][qp] * phi_face[k][qp];
//	                            }
//	                        }// Assemble Matrix BC
//	                    } // Loop over face qp
//				}
//

            } // If on mesh boundary
        } // loop over sides

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations.
        //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The NumericMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    } // end of element loop
    std::cout << "End Stokes Assembly" << std::endl;

}

void assemble_poisson(EquationSystems &es, const std::string& libmesh_dbg_var(system_name))
{
    std::cout << "Start Poisson Assembly" << std::endl;
// It is a good idea to make sure we are assembling
// the proper system.

// Get a constant reference to the mesh object.
    const MeshBase &mesh = es.get_mesh();

// The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

// Get a reference to the Convection-Diffusion system object.
    LinearImplicitSystem &system = es.get_system < LinearImplicitSystem > ("Poisson");

// Numeric ids corresponding to each variable in the system
    const unsigned int u_var = system.variable_number("phi");

// Get the Finite Element type for "u".  Note this will be
// the same as the type for "v".
    FEType fe_vel_type = system.variable_type(u_var);

// Build a Finite Element object of the specified type for
// the velocity variables.
    std::unique_ptr < FEBase > fe_vel(FEBase::build(dim, fe_vel_type));

// A Gauss quadrature rule for numerical integration.
// Let the FEType object decide what order rule is appropriate.
//QGauss qrule(dim, fe_vel_type.default_quadrature_order());
    QGauss qrule(dim, SECOND);

// Tell the finite element objects to use our quadrature rule.
    fe_vel->attach_quadrature_rule(&qrule);

// Here we define some references to cell-specific data that
// will be used to assemble the linear system.
//
// The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe_vel->get_JxW();

// The element shape function gradients for the velocity
// variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();
    const std::vector<std::vector<Number>> &phi = fe_vel->get_phi();

// A reference to the DofMap object for this system.  The DofMap
// object handles the index translation from node and element numbers
// to degree of freedom numbers.  We will talk more about the DofMap
// in future examples.
    const DofMap &dof_map = system.get_dof_map();

// Define data structures to contain the element matrix
// and right-hand-side vector contribution.  Following
// basic finite element terminology we will denote these
// "Ke" and "Fe".
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

// This vector will hold the degree of freedom indices for
// the element.  These define where in the global system
// the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;

// Now we will loop over all the elements in the mesh that
// live on the local processor. We will compute the element
// matrix and right-hand-side contribution.  In case users later
// modify this program to include refinement, we will be safe and
// will only consider the active elements; hence we use a variant of
// the active_elem_iterator.
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        const unsigned int n_dofs = dof_indices.size();

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_vel->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);

        // Now we will build the element matrix.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            // Assemble the v-velocity row
            // vv coupling
            for (unsigned int i = 0; i < n_dofs; i++)
            {
                for (unsigned int j = 0; j < n_dofs; j++)
                {
                    Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);

                }
                Fe(i) += JxW[qp] * phi[i][qp];
            }
        } // end of the quadrature point qp-loop

        // At this point the interior element integration has
        // been completed.  However, we have not yet addressed
        // boundary conditions.  For this example we will only
        // consider simple Dirichlet boundary conditions imposed
        // via the penalty method. The penalty method used here
        // is equivalent (for Lagrange basis functions) to lumping
        // the matrix resulting from the L2 projection penalty
        // approach introduced in example 3.

//        {
//            // The following loops over the sides of the element.
//            // If the element has no neighbor on a side then that
//            // side MUST live on a boundary of the domain.
//            for (auto s : elem->side_index_range())
//            {
//                if (elem->neighbor_ptr(s) == libmesh_nullptr)
//                {
//                    const unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, s);
//                    if (boundary_id == 7)
//                    {
//                        std::unique_ptr<const Elem> side(elem->build_side_ptr(s));
//
//                        // Loop over the nodes on the side.
//                        for (auto ns : side->node_index_range())
//                        {
//                            // The penalty value.  \f$ \frac{1}{\epsilon \f$
//                            const Real penalty = 1.e5;
//
//                            // The boundary values.
//
//                            // Set u = 1 on the top boundary, 0 everywhere else
//                            // On sideset 1 & 2: u = 1, v = 0, w = 0;
//                            // On sidesets 3 & 4 & 5 & 6: sigma n = 0
//                            // else u = 0, v = 0, w = 0
//                            // Set v = 0 everywhere
//                            Real phi_value = 0.;
//
//                            // Find the node on the element matching this node on
//                            // the side.  That defined where in the element matrix
//                            // the boundary condition will be applied.
//                            for (auto n : elem->node_index_range())
//                                if (elem->node_id(n) == side->node_id(ns))
//                                {
//                                    // Matrix contribution.
//                                    Ke(n, n) += penalty;
//
//                                    // Right-hand-side contribution.
//                                    Fe(n) += penalty * phi_value;
//                                }
//                        } // end face node loop
//                    } // if on penalty boundaries
//                } // end if (elem->neighbor(side) == libmesh_nullptr)
//            }
//        } // end boundary condition section

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations.
        //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The NumericMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    } // end of element loop
    std::cout << "End Poisson Assembly" << std::endl;

}

