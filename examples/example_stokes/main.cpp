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

std::map<boundary_id_type, double> pressure_bc;

double tau;
double mu;
Order elOrder;

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

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    ParallelMesh mesh(init.comm());

    // Use the MeshTools::Generation mesh generator to create a uniform
    // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
    // to build a mesh of 8x8 Quad9 elements.  Building these
    // higher-order elements allows us to use higher-order
    // approximation, as in example 3.
    std::string input_mesh = data("input", "NONE");
    if("NONE" == input_mesh)
    {
        std::cout << "Creating mesh" << std::endl;
        double maxx = data("maxx", 1.0);
        double maxy = data("maxy", 1.0);
        int nelx = data("nelx", 0);
        int nely = data("nely", 0);
        MeshTools::Generation::build_square (mesh,
                                           nelx, nely,
                                           0., maxx,
                                           0., maxy,
                                           TRI3);
    }
    else
    {
        std::cout << "Reading mesh" << std::endl;
        int n_refinements = data("nrefs", 0);
        BeatIt::serial_mesh_partition(init.comm(), input_mesh, &mesh, n_refinements);
    }

    std::cout << "Scaling mesh" << std::endl;
    double scale = data("scale", 1.0);
    MeshTools::Modification::scale(mesh, scale, scale, scale);
    // Print information about the mesh to the screen.
    mesh.print_info();
    auto dim = mesh.mesh_dimension();

    // Create an equation systems object.
    std::cout << "Creating equation systems" << std::endl;
    EquationSystems equation_systems(mesh);

    // Declare the system and its variables.
    // Create a transient system named "Stokes"
    std::cout << "Creating Stokes system" << std::endl;
    LinearImplicitSystem &system = equation_systems.add_system < LinearImplicitSystem > ("Stokes");

    //Get parameters for stokes
    // viscosity
    mu = data("mu", 1.0);
    // Stabilization parameter
    tau = data("tau", -1.0);

    // Add the variables "u" & "v" to "Stokes".  They
    // will be approximated using second-order approximation.
    elOrder = FIRST;
    if(tau <= 0.0)
    {
        elOrder = SECOND;
        mesh.all_second_order();
    }
    unsigned int ux = system.add_variable("ux", elOrder);
    unsigned int uy = system.add_variable("uy", elOrder);
    unsigned int uz = 100;
    if(dim == 3) uz = system.add_variable("uz", elOrder);

    // Add the variable "p" to "Stokes". This will
    // be approximated with a first-order basis,
    // providing an LBB-stable pressure-velocity pair.
    unsigned int p = system.add_variable("p", FIRST);

    // Give the system a pointer to the matrix assembly
    // function.
    system.attach_assemble_function(assemble_stokes);

    std::cout << "Adding Dirichlet BC for parabolic velocity profile" << std::endl;
    // Read pressure boundary conditions
    std::string p_bc_id = data("p_bc_id", "NONE");
    std::string p_bc = data("p_bc", "NONE");
    std::vector<libMesh::boundary_id_type> p_bc_id_vec;
    BeatIt::readList(p_bc_id, p_bc_id_vec);
    std::vector<double> p_bc_vec;
    BeatIt::readList(p_bc, p_bc_vec);
    if(p_bc_vec.size() != p_bc_id_vec.size() )
    {
        std::cout << "Pressure boundary conditions are not set up correctly in the input file: " << std::endl;
        std::cout << "Specify the list of sidesets IDs in 'p_bc_id', e.g. p_bc_id = '1, 3, 6' " << std::endl;
        std::cout << "Specify the list of corresponding pressures in 'p_bc', e.g. p_bc_id = '10.0, 12.0, 11.0' " << std::endl;
        std::cout << "The lists p_bc_id and p_pc must have the same number of entries. " << std::endl;
        std::cout << "Aborting. " << std::endl;
        throw std::runtime_error("Pressure Boundary conditions: wrong input file");
    }
    else
    {
        for(unsigned int k = 0; k< p_bc_vec.size(); ++k)
        {
            pressure_bc[p_bc_id_vec[k]] = p_bc_vec[k];
        }
    }
    std::set < libMesh::boundary_id_type > dirichlet_stokes;

    std::string bcs = data("Dbc", "NONE");
    std::set < libMesh::boundary_id_type > bcs_vec;
    BeatIt::readList(bcs, dirichlet_stokes);
    std::vector<unsigned int> vars_vel_stokes(3);
    vars_vel_stokes[0] = ux;
    vars_vel_stokes[1] = uy;
    if(dim == 3) vars_vel_stokes[2] = uz;
    AnalyticFunction<> bczv(zfv);
    libMesh::ZeroFunction < Number > zero;
    libMesh::DirichletBoundary dirichlet_bc_stokes(dirichlet_stokes, vars_vel_stokes, &zero);
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc_stokes);

    // Initialize the data structures for the equation system.
    std::cout << "Initialize equation systems" << std::endl;
    equation_systems.init();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    // Assemble & solve the linear system,
    // then write the solution.
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

    // double scale = 1.0e-4;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Stokes");

    // Get a constant reference to the mesh object.
    const MeshBase &mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    LinearImplicitSystem &system = es.get_system < LinearImplicitSystem > ("Stokes");

    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = system.variable_number("ux");
    const unsigned int v_var = system.variable_number("uy");
    unsigned int w_var = 100;
    if(dim == 3) w_var = system.variable_number("uz");
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
    QGauss qrule(dim, libMesh::FIFTH);

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
    const std::vector<std::vector<RealGradient>> &dpsi = fe_pres->get_dphi();

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

    // BC stuff
    libMesh::UniquePtr < libMesh::FEBase > fe_face(libMesh::FEBase::build(dim, fe_vel_type));
    libMesh::QGauss qface(dim - 1, libMesh::FIFTH);
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
        if(dim == 3) dof_map.dof_indices(elem, dof_indices_w, w_var);
        dof_map.dof_indices(elem, dof_indices_p, p_var);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        unsigned int n_w_dofs = 0;
        if(dim == 3) n_w_dofs = dof_indices_w.size();
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

        Kvu.reposition(v_var * n_u_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition(v_var * n_u_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition(v_var * n_u_dofs, p_var * n_v_dofs, n_v_dofs, n_p_dofs);

        if(dim == 3)Kww.reposition(w_var * n_w_dofs, w_var * n_w_dofs, n_w_dofs, n_w_dofs);
        if(dim == 3)Kwp.reposition(w_var * n_w_dofs, p_var * n_w_dofs, n_w_dofs, n_p_dofs);

        Kpu.reposition(p_var * n_u_dofs, u_var * n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition(p_var * n_u_dofs, v_var * n_u_dofs, n_p_dofs, n_v_dofs);
        if(dim == 3)Kpw.reposition(p_var * n_u_dofs, w_var * n_u_dofs, n_p_dofs, n_w_dofs);
        Kpp.reposition(p_var * n_u_dofs, p_var * n_u_dofs, n_p_dofs, n_p_dofs);

        Fu.reposition(u_var * n_u_dofs, n_u_dofs);
        Fv.reposition(v_var * n_u_dofs, n_v_dofs);
        if(dim == 3) Fw.reposition(w_var * n_u_dofs, n_w_dofs);
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
                    Kuu(i, j) += JxW[qp] * mu * (dphi[i][qp] * dphi[j][qp]);
                    Kvv(i, j) += JxW[qp] * mu * (dphi[i][qp] * dphi[j][qp]);
                    if(dim == 3) Kww(i, j) += JxW[qp] * mu * (dphi[i][qp] * dphi[j][qp]);
                }
            }

            // up coupling
            for (unsigned int i = 0; i < n_u_dofs; i++)
            {
                for (unsigned int j = 0; j < n_p_dofs; j++)
                {
                    Kup(i, j) -= JxW[qp] * psi[j][qp] * dphi[i][qp](0);
                    Kvp(i, j) -= JxW[qp] * psi[j][qp] * dphi[i][qp](1);
                    if(dim == 3) Kwp(i, j) -= JxW[qp] * psi[j][qp] * dphi[i][qp](2);
                }
            }

            // Assemble the pressure row
            // pu coupling
            for (unsigned int i = 0; i < n_p_dofs; i++)
            {
                for (unsigned int j = 0; j < n_u_dofs; j++)
                {
                    Kpu(i, j) -= JxW[qp] * psi[i][qp] * dphi[j][qp](0);
                    Kpv(i, j) -= JxW[qp] * psi[i][qp] * dphi[j][qp](1);
                    if(dim == 3) Kpw(i, j) -= JxW[qp] * psi[i][qp] * dphi[j][qp](2);
                }
            }

            // pp coupling
            if(FIRST == elOrder)
            {
                double h = elem->hmax();

                for (unsigned int i = 0; i < n_p_dofs; i++)
                {
                    for (unsigned int j = 0; j < n_p_dofs; j++)
                    {
                        Kpp(i, j) -= tau * 0.5 * h * h * JxW[qp] * dpsi[i][qp] * dpsi[j][qp];
                    }
                }
            }
            else
            {
                for (unsigned int i = 0; i < n_p_dofs; i++)
                {
                    for (unsigned int j = 0; j < n_p_dofs; j++)
                    {
                        Kpp(i, j) -= 1e-6 * JxW[qp] * psi[i][qp] * psi[j][qp];
                    }
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
                const unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, side);

                // Pressure boundary condition
                if(pressure_bc.find(boundary_id) != pressure_bc.end() )
                {
                    double pressure = pressure_bc.find(boundary_id)->second;
                    fe_face->reinit(elem, side);
                    // Loop over face qp
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        const double xq = qface_point[qp](0);
                        const double yq = qface_point[qp](1);
                        const double zq = qface_point[qp](2);
                        // Assemble Matrix BC
                        for (unsigned int l = 0; l < phi_face.size(); l++)
                        {
                            for (unsigned int k = 0; k < phi_face.size(); k++)
                            {
                                Kuu(l, k) -= JxW_face[qp] * mu * (dphi_face[k][qp] * normals[qp]) * phi_face[l][qp];
                                Kvv(l, k) -= JxW_face[qp] * mu * (dphi_face[k][qp] * normals[qp]) * phi_face[l][qp];
                                if(dim == 3) Kww(l, k) += JxW_face[qp] * mu * (dphi_face[l][qp] * normals[qp]) * phi_face[k][qp];
                            }
                        } // Assemble Matrix BC
                        // Assemble
                        for (unsigned int l = 0; l < phi_face.size(); l++)
                        {
                            Fu(l) -= JxW_face[qp] * pressure * normals[qp](0) * phi_face[l][qp];
                            Fv(l) -= JxW_face[qp] * pressure * normals[qp](1) * phi_face[l][qp];
                            if(dim == 3) Fw(l) -= JxW_face[qp] * pressure * normals[qp](2) * phi_face[l][qp];
                        } // Assemble RHS BC

                    }

                }
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
        system.matrix->add_matrix(Ke, dof_indices, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    } // end of element loop
    std::cout << "End Stokes Assembly" << std::endl;
}


