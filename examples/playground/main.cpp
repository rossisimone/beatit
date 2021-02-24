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

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
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

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

// The definition of a geometric element
#include "libmesh/elem.h"
#include "libmesh/getpot.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_poisson(EquationSystems & es, const std::string & system_name);
// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_stokes(EquationSystems & es, const std::string & system_name);

void zf (DenseVector<Number> & output,
        const Point & p,
        const Real)
{
    output(0) = 0.0;
}

void of (DenseVector<Number> & output,
        const Point & p,
        const Real)
{
    output(0) = 1.0;
}


void inflow (DenseVector<Number> & output,
        const Point & p,
        const Real)
{
    double r = std::sqrt(p(0)*p(0)+p(1)*p(1));
    output(0) = (1.0-r/9e-3*r/9e-3);
    std::cout << "of:  " << output << std::endl;
}

// The main program.
int main(int argc, char ** argv)
{
    // Initialize libMesh.
    LibMeshInit init(argc, argv);

    // across the default MPI communicator.
    Mesh mesh(init.comm());

    MeshTools::Generation::build_square (mesh,
                                       2, 1,
                                       -1., 1.,
                                       -1., 1.,
                                       QUAD4);


    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
         auto p = elem->centroid();
         if(p(0) < 0) elem->subdomain_id() = 1;
         else elem->subdomain_id() = 2;
    }
    mesh.prepare_for_use();
    mesh.print_info();

    EquationSystems equation_systems(mesh);
    LinearImplicitSystem & system_profile = equation_systems.add_system<LinearImplicitSystem>("Poisson");
    std::set< libMesh::subdomain_id_type > left_subdomain;
    std::set< libMesh::subdomain_id_type > right_subdomain;
    left_subdomain.insert(1);
    right_subdomain.insert(2);
    system_profile.add_variable("L", FIRST, LAGRANGE, &left_subdomain);
    system_profile.add_variable("R", FIRST, LAGRANGE, &right_subdomain);


    LinearImplicitSystem & system_profile2 = equation_systems.add_system<LinearImplicitSystem>("Poisson2");
    system_profile2.add_variable("U", FIRST, LAGRANGE);

    std::cout << "Init System" << std::endl;
    equation_systems.init();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    DofMap & dof_map = system_profile.get_dof_map();
    DofMap & dof_map2 = system_profile2.get_dof_map();
    std::vector<libMesh::dof_id_type> l_dofs;
    std::vector<libMesh::dof_id_type> r_dofs;
    std::vector<libMesh::dof_id_type> u_dofs;
    std::vector<libMesh::dof_id_type> dofs;
    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices(elem, dofs);  
        dof_map.dof_indices(elem, l_dofs, 0);  
        dof_map.dof_indices(elem, r_dofs, 1) ; 
        dof_map.dof_indices(elem, u_dofs) ;

       auto id = elem->id();

       std::cout << "On elem " << id << "\n"
                 << "#   dofs: " << dofs.size() << "\n"  
                 << "# l_dofs: " << l_dofs.size() << "\n"  
                 << "# r_dofs: " << r_dofs.size() << "\n"  
                 << "# u_dofs: " << u_dofs.size() << std::endl;  
    }


    libMesh::ExodusII_IO(mesh).write_equation_systems ("poisson.e",
            equation_systems);



    // All done.
    return 0;
}

void assemble_poisson(EquationSystems & es, const std::string & libmesh_dbg_var(system_name))
{
    std::cout << "Start Poisson Assembly" << std::endl;
    // It is a good idea to make sure we are assembling
    // the proper system.

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Poisson");

    // Numeric ids corresponding to each variable in the system
    const unsigned int u_var = system.variable_number("phi");

    // Get the Finite Element type for "u".  Note this will be
    // the same as the type for "v".
    FEType fe_vel_type = system.variable_type(u_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    std::unique_ptr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));

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
    const std::vector<Real> & JxW = fe_vel->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> & dphi = fe_vel->get_dphi();
    const std::vector<std::vector<Number>> & phi = fe_vel->get_phi();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix < Number > Ke;
    DenseVector < Number > Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the active_elem_iterator.
    for (const auto & elem : mesh.active_local_element_ptr_range())
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
        dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The NumericMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    } // end of element loop
    std::cout << "End Poisson Assembly" << std::endl;


}


void assemble_stokes(EquationSystems & es, const std::string & libmesh_dbg_var(system_name))
{
    std::cout << "Start Assembly" << std::endl;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Stokes");

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Stokes");
    LinearImplicitSystem & system_phi = es.get_system<LinearImplicitSystem>("Poisson");

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
    std::unique_ptr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));

    // Build a Finite Element object of the specified type for
    // the pressure variables.
    std::unique_ptr<FEBase> fe_pres(FEBase::build(dim, fe_pres_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    //QGauss qrule(dim, fe_vel_type.default_quadrature_order());
    QGauss qrule(dim, SECOND);

    // Tell the finite element objects to use our quadrature rule.
    fe_vel->attach_quadrature_rule(&qrule);
    fe_pres->attach_quadrature_rule(&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe_vel->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> & dphi = fe_vel->get_dphi();

    // The element shape functions for the pressure variable
    // evaluated at the quadrature points.
    const std::vector<std::vector<Real>> & psi = fe_pres->get_phi();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    const DofMap & dof_map_phi = system_phi.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix < Number > Ke;
    DenseVector < Number > Fe;

    DenseSubMatrix<Number> Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke), Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke), Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke), Kpu(Ke), Kpv(
            Ke), Kpw(Ke), Kpp(Ke);

    DenseSubVector<Number> Fu(Fe), Fv(Fe), Fw(Fe), Fp(Fe);

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_w;
    std::vector<dof_id_type> dof_indices_p;

    std::vector<dof_id_type> dof_indices_phi;
    std::vector<double> profile;
    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  In case users later
    // modify this program to include refinement, we will be safe and
    // will only consider the active elements; hence we use a variant of
    // the active_elem_iterator.
    system_phi.update();
    std::unique_ptr<libMesh::FEBase> fe_face(libMesh::FEBase::build(dim, fe_pres_type));
    libMesh::QGauss qface(dim - 1, libMesh::FIRST);
    fe_face->attach_quadrature_rule(&qface);


    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        double h = elem->hmax();
        dof_map_phi.dof_indices(elem, dof_indices_phi);
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_u, u_var);
        dof_map.dof_indices(elem, dof_indices_v, v_var);
        dof_map.dof_indices(elem, dof_indices_w, w_var);
        dof_map.dof_indices(elem, dof_indices_p, p_var);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        const unsigned int n_w_dofs = dof_indices_w.size();
        const unsigned int n_p_dofs = dof_indices_p.size();

        profile.resize(dof_indices_phi.size());
        system_phi.current_local_solution->get(dof_indices_phi, profile);

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
        Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
        Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);

        Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);

        Kww.reposition (w_var*n_w_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);
        Kwp.reposition (w_var*n_w_dofs, p_var*n_w_dofs, n_w_dofs, n_p_dofs);

        Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
        Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
        Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

        Fu.reposition (u_var*n_u_dofs, n_u_dofs);
        Fv.reposition (v_var*n_u_dofs, n_v_dofs);
        Fv.reposition (w_var*n_u_dofs, n_w_dofs);
        Fp.reposition (p_var*n_u_dofs, n_p_dofs);

        // Now we will build the element matrix.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            // Assemble the u-velocity row
            // uu coupling
            for (unsigned int i=0; i<n_u_dofs; i++)
              for (unsigned int j=0; j<n_u_dofs; j++)
                Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

            // up coupling
            for (unsigned int i=0; i<n_u_dofs; i++)
              for (unsigned int j=0; j<n_p_dofs; j++)
                Kup(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](0);


            // Assemble the v-velocity row
            // vv coupling
            for (unsigned int i=0; i<n_v_dofs; i++)
              for (unsigned int j=0; j<n_v_dofs; j++)
                Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

            // vp coupling
            for (unsigned int i=0; i<n_v_dofs; i++)
              for (unsigned int j=0; j<n_p_dofs; j++)
                Kvp(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](1);

            // Assemble the v-velocity row
            // vv coupling
            for (unsigned int i=0; i<n_v_dofs; i++)
              for (unsigned int j=0; j<n_v_dofs; j++)
                Kww(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

            // vp coupling
            for (unsigned int i=0; i<n_v_dofs; i++)
              for (unsigned int j=0; j<n_p_dofs; j++)
                Kwp(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](2);

            // Assemble the pressure row
            // pu coupling
            for (unsigned int i=0; i<n_p_dofs; i++)
              for (unsigned int j=0; j<n_u_dofs; j++)
                Kpu(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](0);

            // pv coupling
            for (unsigned int i=0; i<n_p_dofs; i++)
              for (unsigned int j=0; j<n_v_dofs; j++)
                Kpv(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](1);

            // pv coupling
            for (unsigned int i=0; i<n_p_dofs; i++)
              for (unsigned int j=0; j<n_w_dofs; j++)
                Kpw(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](2);

            // pv coupling
            // div u = 0
            // q * div u + q * div u' = 0
            // q * div u - nabla q * u' = 0
            // u' = - tau * Res
            // Res = mu Delta u - nabla p
            // u' = tau * nabla p
            // - q * div u - tau * nabla q * nabla p = 0
            double tau = 100*0.5*h*h;
            for (unsigned int i = 0; i < n_p_dofs; i++)
            {
                for (unsigned int j = 0; j < n_p_dofs; j++)
                {
                    Kpp(i, j) += -tau * JxW[qp] * dphi[i][qp] * dphi[j][qp];
                }
            }

        } // end of the quadrature point qp-loop


        {
            for (unsigned int side = 0; side < elem->n_sides(); side++)
            {
                if (elem->neighbor_ptr(side) == libmesh_nullptr)
                {
                    const unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, side);

                    if(boundary_id != 3)
                    {

                        double bc_value_u = 0.0;
                        double bc_value_v = 0.0;
                        double bc_value_w = 1000.0;
                        if(boundary_id == 7 ) bc_value_u = 0.0;


                     const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();
                     const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();
                     int n_phi = phi_face.size();

                     const std::vector<libMesh::Point> & qface_point = fe_face->get_xyz();
                     const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_face->get_dphi();
                     const std::vector<libMesh::Point>& normals = fe_face->get_normals();
                     fe_face->reinit(elem, side);

                     for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                     {
                         double phi = 0.0;
                         for(unsigned int i = 0; i < n_phi; ++i)
                         {
                             phi += profile[i] * phi_face[i][qp];
                         }
                         bc_value_u *= phi;

                         double penalty = 1e6;
                         for (unsigned int i = 0; i < n_phi; i++)
                         {
                             Fu(i) += penalty * JxW_face[qp] * bc_value_u * phi_face[i][qp];
                             Fv(i) += penalty * JxW_face[qp] * bc_value_v * phi_face[i][qp];
                             Fw(i) += penalty * JxW_face[qp] * bc_value_w * phi_face[i][qp];
                         }

                         for (unsigned int i = 0; i < n_phi; i++)
                         {
                             for (unsigned int j = 0; j < n_phi; j++)
                             {
                                 Kuu(i,j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];
                                 Kvv(i,j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];
                                 Kww(i,j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];
                             }
                         }

                     }

                    }
                }
            }

        }

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations.
        //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        //dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The NumericMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    } // end of element loop
    std::cout << "End Assembly" << std::endl;

}
