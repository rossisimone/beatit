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
#include "libmesh/transient_system.h"
#include "Util/IO/io.hpp"
#include "libmesh/getpot.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include "petscmat.h"
#include "petscksp.h"
#include <petsc/private/kspimpl.h>
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
// Bring in everything from the libMesh namespace
using namespace libMesh;

std::map<boundary_id_type, double> pressure_bc;
double mu;
double tau;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_poisson(EquationSystems & es, const std::string & system_name);
 // this fun is called in 110 in main, so dont delete!

void zf(DenseVector<Number> & output, const Point & p, const Real)
{
    output(0) = 0.0;
}
// Function prototype.  This function will assemble the system
// matrix and right-hand-side.
void assemble_stokes(EquationSystems & es, const std::string & system_name); // I think I can delete this??

// The main program.
int main(int argc, char ** argv)
{
    // Initialize libMesh.
    LibMeshInit init(argc, argv);
    GetPot input = BeatIt::readInputFile(argc, argv);
    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());

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
    std::string filename = input("MESH_FILE", "default_meshname");
    mesh.read(filename);
    double scale = input("SCALE", 1.0e-4); // or could set to one
    MeshTools::Modification::scale(mesh, scale, scale, scale);
    //mesh.all_second_order();

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    LinearImplicitSystem & system_profile = equation_systems.add_system<LinearImplicitSystem>("Poisson");
    system_profile.add_variable("phi", FIRST);
    system_profile.attach_assemble_function(assemble_poisson);
    std::set<libMesh::boundary_id_type> dirichlet_poisson;
    int dirichletID= input("dirID", 0); // second arg is default
    dirichlet_poisson.insert(dirichletID); // was for ID 7
    std::vector<unsigned int> vars_p(1);
    vars_p[0] = 0;
    AnalyticFunction<> bcz(zf);
    libMesh::DirichletBoundary dirichlet_bc_poisson(dirichlet_poisson, vars_p, bcz);
    system_profile.get_dof_map().add_dirichlet_boundary(dirichlet_bc_poisson); // add sidesets and un comment this

    // Declare the system and its variables.
    // Create a transient system named "Stokes"
    TransientLinearImplicitSystem & system = equation_systems.add_system<TransientLinearImplicitSystem>("Stokes");

    // Add the variables "u" & "v" to "Stokes".  They
    // will be approximated using first order ( instead of second) approximation. For p1 p1
    system.add_variable("ux", FIRST);
    system.add_variable("uy", FIRST);
    system.add_variable("uz", FIRST);

    // still keeping ux, uy, uz since pressure doesnt vary spatially

    // Add the variable "p" to "Stokes". This will
    // be approximated with a first-order basis,
    // providing an LBB-stable pressure-velocity pair.
    system.add_variable("p", FIRST);

    // Give the system a pointer to the matrix assembly
    // function.
    system.attach_assemble_function(assemble_stokes);

    // adding pressure_BC stuff here

// chanign defaults of the below to test if it's reading it from the input
    std::string p_bc_id = input("PRESS_BC_IDS", "8,9"); // default was 1,3 testing that it's actually reading input
    //std::string p_bc_id = "1,3"; //data("p_bc_id", "NONE");
    std::string p_bc_amp = input("PRESS_BC_AMP", "8, 10"); // also testing, want to read 0,2 from input this will just be the amp that they do //data("p_bc_amp", "NONE");
    std::vector<libMesh::boundary_id_type> p_bc_id_vec;
    BeatIt::readList(p_bc_id, p_bc_id_vec);
    std::vector<double> p_bc_vec;
    BeatIt::readList(p_bc_amp, p_bc_vec);
    for(unsigned int k = 0; k< p_bc_vec.size(); ++k)
    {
        pressure_bc[p_bc_id_vec[k]] = p_bc_vec[k];
    }


    ExplicitSystem & system_mu = equation_systems.add_system<ExplicitSystem>("MU"); // Mu viscosity
      ExplicitSystem & system_tau = equation_systems.add_system<ExplicitSystem>("TAU"); // Tau shear stress
    system_mu.add_variable("mu", CONSTANT, MONOMIAL);
    system_tau.add_variable("tau", CONSTANT, MONOMIAL);
    // Initialize the data structures for the equation system.
    equation_systems.init();
    system_mu.solution->close();
    system_tau.solution->close();

    equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 250;// may edit with input file
    equation_systems.parameters.set<Real>("linear solver tolerance") = TOLERANCE;
    const Real dt = input("DT", 0.015);
    system.time   = 0.0; // current time
    auto n_timesteps = input("N_TS", 10); // may want to change
    equation_systems.parameters.set<Real> ("dt")   = dt;
    equation_systems.parameters.set<Real> ("rho")   = input("RHO", 1.06); // may want to get density from input
    equation_systems.parameters.set<Real> ("mu") = input("MU", 0.0037*1e3); // also mu
    equation_systems.parameters.set<Real> ("tau") = 0; // also mu


    equation_systems.parameters.set<Real> ("penalty") = input("PENALTY", 1.e10);
    equation_systems.parameters.set<Real> ("radius") = input("RADIUS", 50);
    equation_systems.parameters.set<Real> ("scale")   = scale;

    equation_systems.parameters.set<Real> ("v_wall_id") = input("V_WALL_IDS", 2);
    equation_systems.parameters.set<Real> ("inflow_id") = input("INFLOW_ID", 3);
    equation_systems.parameters.set<Real> ("outflow_id") = input("OUTFLOW_ID", 4);


    // Prints information about the system to the screen.
    equation_systems.print_info();


    // Assemble & solve the linear system,
    // then write the solution.
    equation_systems.get_system("Poisson").solve(); // to get BCs // do I still need Poisson??

    ExodusII_IO exporter(mesh); // maybe scaled

    for (unsigned int t_step=1; t_step<=n_timesteps; ++t_step)
    {
        // Increment the time counter, set the time step size as
        // a parameter in the EquationSystem.
        system.time += dt;

        // A pretty update message
        libMesh::out << "\n\n*** Solving time step "
                     << t_step
                     << ", time = "
                     << system.time
                     << " ***"
                    << std::endl;

        // Now we need to update the solution vector from the
        // previous time step.  This is done directly through
        // the reference to the Stokes system.
        *system.old_local_solution = *system.current_local_solution; // updating

        // At the beginning of each solve, reset the linear solver tolerance
        // to a "reasonable" starting value.
        const Real initial_linear_solver_tol = input("LIN_SOLV_TOL", 1.e-6);
        equation_systems.get_system("Stokes").solve();

        const unsigned int write_interval = input("WRITE_ITER", 2); // writing for every iteration
        std::string output_filename = input("Output_FILENAME", "carreau_def.e");
        if ((t_step+1)%write_interval == 0)
          {
            exporter.write_timestep(output_filename,
                                  equation_systems,
                                  t_step+1, // we're off by one since we wrote the IC and the Exodus numbering is 1-based.
                                  system.time);
          }

    }

//#ifdef LIBMESH_HAVE_EXODUS_API
//    ExodusII_IO(mesh).write_equation_systems ("out.e",
//            equation_systems);
//#endif // #ifdef LIBMESH_HAVE_EXODUS_API

    // All done.
    return 0;
}

void assemble_stokes(EquationSystems & es, const std::string & libmesh_dbg_var(system_name))
{
    double scale = es.parameters.get<Real> ("scale");
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Stokes");

    double inflow_director_x = 0.5312; // where did these values come from?
    double inflow_director_z =-0.2540;
    double inflow_director_y = 0.8083;

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system<TransientLinearImplicitSystem>("Stokes");
    LinearImplicitSystem & system_phi = es.get_system<LinearImplicitSystem>("Poisson");
    ExplicitSystem & system_mu = es.get_system<ExplicitSystem>("MU");
    ExplicitSystem & system_tau = es.get_system<ExplicitSystem>("TAU");

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
    QGauss qrule(dim, fe_vel_type.default_quadrature_order());

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
    const std::vector<std::vector<Real>> & phi = fe_vel->get_phi();
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
    const DofMap & dof_map_mu = system_mu.get_dof_map();
    const DofMap & dof_map_tau = system_tau.get_dof_map();

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
    std::vector<dof_id_type> dof_indices_mu;
    std::vector<dof_id_type> dof_indices_tau;

    // Adding BC stuff from stokes example
    libMesh::UniquePtr < libMesh::FEBase > fe_face(libMesh::FEBase::build(dim, fe_vel_type));
    libMesh::QGauss qface(dim - 1, libMesh::FIFTH);
    fe_face->attach_quadrature_rule(&qface);
    const std::vector<libMesh::Real> &JxW_face = fe_face->get_JxW();
    const std::vector<std::vector<libMesh::Real> > &phi_face = fe_face->get_phi();
    int n_phi = phi_face.size();
    const std::vector<libMesh::Point> &qface_point = fe_face->get_xyz();
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_face = fe_face->get_dphi();
    const std::vector<libMesh::Point> &normals = fe_face->get_normals();


    double mu_inf = es.parameters.get<Real> ("mu");
    double tau_not = es.parameters.get<Real> ("tau");

    double rho = es.parameters.get<Real> ("rho");
    double dt = es.parameters.get<Real> ("dt");
    double penalty = es.parameters.get<Real> ("penalty");
    double radius = es.parameters.get<Real> ("radius");
    double v_wall_id = es.parameters.get<Real> ("v_wall_id");
    double inflow_id = es.parameters.get<Real> ("inflow_id");
    double outflow_id = es.parameters.get<Real> ("outflow_id");



    libMesh::out << " ***** Mu:"
                  << mu_inf
                  <<  std::endl;

    double time = system.time;
    double sint = std::sin(3*3.1415*(time+0.025));

    std::vector<double> un;
    std::vector<double> vn;
    std::vector<double> wn;
    RealGradient veln;
    libMesh::TensorValue<libMesh::Number> grad_vel;
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
        dof_map.dof_indices(elem, dof_indices_u, u_var);
        dof_map.dof_indices(elem, dof_indices_v, v_var);
        dof_map.dof_indices(elem, dof_indices_w, w_var);
        dof_map.dof_indices(elem, dof_indices_p, p_var);
        dof_map_mu.dof_indices(elem, dof_indices_mu);
        dof_map_tau.dof_indices(elem, dof_indices_tau);

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
        un.resize(n_dofs);
        vn.resize(n_dofs);
        wn.resize(n_dofs);

        system.current_local_solution->get(dof_indices_u, un);
        system.current_local_solution->get(dof_indices_v, vn);
        system.current_local_solution->get(dof_indices_w, wn);

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
        Kuw.reposition(u_var * n_u_dofs, w_var * n_u_dofs, n_u_dofs, n_w_dofs);
        Kup.reposition(u_var * n_u_dofs, p_var * n_u_dofs, n_u_dofs, n_p_dofs);

        Kvu.reposition(v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition(v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        Kvw.reposition(v_var * n_v_dofs, w_var * n_v_dofs, n_v_dofs, n_w_dofs);
        Kvp.reposition(v_var * n_v_dofs, p_var * n_v_dofs, n_v_dofs, n_p_dofs);

        Kwu.reposition(w_var * n_w_dofs, u_var * n_w_dofs, n_w_dofs, n_u_dofs);
        Kwv.reposition(w_var * n_w_dofs, v_var * n_w_dofs, n_w_dofs, n_v_dofs);
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
        double mu_average = 0.0;
         unsigned int nqp = qrule.n_points();
        // Now we will build the element matrix.
        for (unsigned int qp = 0; qp < nqp; qp++)
        {
            veln *= 0.0;
            grad_vel *= 0.0;
            for (unsigned int l = 0; l < n_u_dofs; ++l)
            {
                for (int jdim = 0; jdim < dim; jdim++)
                {
                    grad_vel(0, jdim) += dphi[l][qp](jdim) * un[l];
                    grad_vel(1, jdim) += dphi[l][qp](jdim) * vn[l];
                    grad_vel(2, jdim) += dphi[l][qp](jdim) * wn[l];
                }
            }

            for (unsigned int i = 0; i < n_u_dofs; i++)
            {
                veln(0) += un[i] * phi[i][qp];
                veln(1) += vn[i] * phi[i][qp];
                veln(2) += wn[i] * phi[i][qp];
            }

            double mu0 = 0.0074*1e3;
            double lambda = 0.033; // both from rheology
            auto D = grad_vel + grad_vel.transpose();
            auto g2 = 0.5*D.contract(D);
            mu = mu_inf + (mu0-mu_inf)/std::sqrt(1+lambda*lambda*g2);
            mu_average += mu;

            // Assemble the u-velocity row
            // uu coupling
            for (unsigned int i = 0; i < n_u_dofs; i++)
            {

                Fu(i) += rho / dt * veln(0) * JxW[qp] * phi[i][qp];
                Fv(i) += rho / dt * veln(1) * JxW[qp] * phi[i][qp];
                Fw(i) += rho / dt * veln(2) * JxW[qp] * phi[i][qp];

                for (unsigned int j = 0; j < n_u_dofs; j++)
                {
                    Kuu(i, j) += mu * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    Kuu(i, j) += rho / dt * JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    Kuu(i, j) += rho * JxW[qp] * phi[i][qp]  * (veln * dphi[j][qp]);

                    Kvv(i, j) += mu * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    Kvv(i, j) += rho / dt * JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    Kvv(i, j) += rho * JxW[qp] * phi[i][qp]  * (veln * dphi[j][qp]);

                    Kww(i, j) += mu * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    Kww(i, j) += rho / dt * JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    Kww(i, j) += rho * JxW[qp] * phi[i][qp]  * (veln * dphi[j][qp]);

                }
            }

            // up coupling
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kup(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](0);

            // vp coupling
            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kvp(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](1);

            // vp coupling
            for (unsigned int i = 0; i < n_w_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kwp(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](2);

            // Assemble the pressure row
            // pu coupling
            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++)
                    Kpu(i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp](0);

            // pv coupling
            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++)
                    Kpv(i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp](1);

            // Assemble the pressure row
            // pu coupling
            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_w_dofs; j++)
                    Kpw(i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp](2);

            // pp coupling
            double h = elem->hmax();

            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kpp(i, j) += -0.1 * h * h * JxW[qp] * dphi[i][qp] * dphi[j][qp];

        } // end of the quadrature point qp-loop
        mu_average /= nqp;
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
            for (auto s : elem->side_index_range()) // boundary condition section
            {
                if (elem->neighbor_ptr(s) == libmesh_nullptr)
                {
                    std::unique_ptr<const Elem> side(elem->build_side_ptr(s)); //
                    const unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, s);

                // force u=0 on vessel sides
              for (auto ns : side->node_index_range())
              {

                    if ( boundary_id == v_wall_id)
                    {
                      const Real penalty =1.e10;
                      const Real u_f=0.0;
                      const Real v_f= 0.0;
                      const Real w_f= 0.0;

                    for( auto n : elem -> node_index_range())
                      if ( elem -> node_id(n) ==side -> node_id(ns))
                      {
                      // Matrix contribution.
                        Kuu(n,n) += penalty;
                        Kvv(n,n) += penalty;
                        Kww(n,n) += penalty;

                        // Right-hand-side contribution.
                        Fu(n) += penalty*u_f;
                        Fv(n) += penalty*v_f;
                        Fw(n) += penalty*w_f;
                      }
// adding this in for tau maybe??
                      for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                      {
                          // const double xq = qface_point[qp](0);
                          // const double yq = qface_point[qp](1);
                          // const double zq = qface_point[qp](2);

                          // begin mu stuff
                          veln *= 0.0;
                          grad_vel *= 0.0;
                          for (unsigned int l = 0; l < n_u_dofs; ++l)
                          {
                              for (int jdim = 0; jdim < dim; jdim++)
                              {
                                  grad_vel(0, jdim) += dphi[l][qp](jdim) * un[l];
                                  grad_vel(1, jdim) += dphi[l][qp](jdim) * vn[l];
                                  grad_vel(2, jdim) += dphi[l][qp](jdim) * wn[l];
                              }
                          }

                          // for (unsigned int i = 0; i < n_u_dofs; i++)
                          // {
                          //     veln(0) += un[i] * phi[i][qp];
                          //     veln(1) += vn[i] * phi[i][qp];
                          //     veln(2) += wn[i] * phi[i][qp];
                          // }

                          double mu0 = 0.0074*1e3;
                          double lambda = 0.033; // both from rheology
                          auto D = grad_vel + grad_vel.transpose();
                          auto g2 = 0.5*D.contract(D);
                          mu = mu_inf + (mu0-mu_inf)/std::sqrt(1+lambda*lambda*g2);
                          auto tau= mu*D*normals[qp];
                  }
                }


                    if(pressure_bc.find(boundary_id) != pressure_bc.end() )
                    {
                        double pressure_amp = pressure_bc.find(boundary_id)->second;
                        double pressure= pressure_amp;
                         // grabbing pressure constant as second arg
                        if (boundary_id == inflow_id)
                        { double inflow_freq = 1;
                          pressure = pressure_amp*std::sin(2*3.1415*(time/inflow_freq))*std::sin(2*3.1415*(time/inflow_freq));
                        }
                        else if (boundary_id == outflow_id)
                        {
                          double inflow_freq = 4;
                          pressure = pressure_amp*std::sin(2*3.1415*(time/inflow_freq))*std::sin(2*3.1415*(time/inflow_freq));

                        }

                        fe_face->reinit(elem, s);
                        // Loop over face qp
                        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                        {
                            const double xq = qface_point[qp](0);
                            const double yq = qface_point[qp](1);
                            const double zq = qface_point[qp](2);

                            // begin mu stuff
                            veln *= 0.0;
                            grad_vel *= 0.0;
                            for (unsigned int l = 0; l < n_u_dofs; ++l)
                            {
                                for (int jdim = 0; jdim < dim; jdim++)
                                {
                                    grad_vel(0, jdim) += dphi[l][qp](jdim) * un[l];
                                    grad_vel(1, jdim) += dphi[l][qp](jdim) * vn[l];
                                    grad_vel(2, jdim) += dphi[l][qp](jdim) * wn[l];
                                }
                            }

                            for (unsigned int i = 0; i < n_u_dofs; i++)
                            {
                                veln(0) += un[i] * phi[i][qp];
                                veln(1) += vn[i] * phi[i][qp];
                                veln(2) += wn[i] * phi[i][qp];
                            }

                            double mu0 = 0.0074*1e3;
                            double lambda = 0.033; // both from rheology
                            auto D = grad_vel + grad_vel.transpose();
                            auto g2 = 0.5*D.contract(D);
                            mu = mu_inf + (mu0-mu_inf)/std::sqrt(1+lambda*lambda*g2);
                            auto tau= mu*D*normals[qp];
                            // end mu stuff
                            // Assemble Matrix BC
                            for (unsigned int l = 0; l < phi_face.size(); l++)
                            {
                              // put mu stuff here like extra for loop for grad etc.


                                for (unsigned int k = 0; k < phi_face.size(); k++)
                                {
                                    Kuu(l, k) += JxW_face[qp] * mu * (dphi_face[k][qp] * normals[qp]) * phi_face[l][qp];
                                    Kvv(l, k) += JxW_face[qp] * mu * (dphi_face[k][qp] * normals[qp]) * phi_face[l][qp];
                                    Kww(l, k) += JxW_face[qp] * mu * (dphi_face[k][qp] * normals[qp]) * phi_face[l][qp];
                                }
                            } // Assemble Matrix BC since mu dynamic, for now making 1
                            for (unsigned int l = 0; l < phi_face.size(); l++)
                            {
                                Fu(l) -= JxW_face[qp] * pressure * normals[qp](0) * phi_face[l][qp];
                                Fv(l) -= JxW_face[qp] * pressure * normals[qp](1) * phi_face[l][qp];
                                Fw(l) -= JxW_face[qp] * pressure * normals[qp](2) * phi_face[l][qp];
                            } // Assemble RHS BC

//// Can delete everything below ///////
                          }

                    } // end face node loop

                    } // end if (elem->neighbor(side) == libmesh_nullptr)
                } // end boundary condition section
            }

            // If this assembly program were to be used on an adaptive mesh,
            // we would have to apply any hanging node constraint equations.
            //dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

            // The element matrix and right-hand-side are now built
            // for this element.  Add them to the global matrix and
            // right-hand-side vector.  The NumericMatrix::add_matrix()
            // and NumericVector::add_vector() members do this for us.
            system.matrix->add_matrix(Ke, dof_indices);
            system.rhs->add_vector(Fe, dof_indices);
            system_mu.solution->set(dof_indices_mu[0],mu_average);
            system_tau.solution->set(dof_indices_tau[0],mu_average);

        } // end of element loop

    // SET UP IS for FIELDSPLIT: we should do this only once,
    //                           but since  the matrix is reassmebled at every timestep
    //                             I will put it here for now
    {
        std::cout << "* Assigning field split information ... " << std::flush;

        IS is_v_local;
        IS is_p_local;
        std::vector<libMesh::dof_id_type> v_indices;
        std::vector<libMesh::dof_id_type> vx_indices;
        std::vector<libMesh::dof_id_type> vy_indices;
        std::vector<libMesh::dof_id_type> vz_indices;
        std::vector<libMesh::dof_id_type> p_indices;
        int var_num = 0;
        dof_map.local_variable_indices(vx_indices, mesh, var_num);
        dof_map.local_variable_indices(vy_indices, mesh, var_num+1);
        dof_map.local_variable_indices(vz_indices, mesh, var_num+2);
        dof_map.local_variable_indices(p_indices, mesh, var_num+3);
        v_indices.insert(v_indices.end(), vx_indices.begin(), vx_indices.end());
        v_indices.insert(v_indices.end(), vy_indices.begin(), vy_indices.end());
        v_indices.insert(v_indices.end(), vz_indices.begin(), vz_indices.end());

        ISCreateGeneral(PETSC_COMM_SELF, v_indices.size(), reinterpret_cast<int*>(&v_indices[0]),PETSC_COPY_VALUES,&is_v_local);
        ISCreateGeneral(PETSC_COMM_SELF, p_indices.size(), reinterpret_cast<int*>(&p_indices[0]),PETSC_COPY_VALUES,&is_p_local);

        typedef libMesh::PetscMatrix<libMesh::Number> PetscMatrix;
        typedef libMesh::PetscLinearSolver<libMesh::Number> PetscLinearSolver;
        auto linear_solver = dynamic_cast<PetscLinearSolver *>(system.get_linear_solver());
        linear_solver->init(dynamic_cast<PetscMatrix *>(system.matrix), "ns_");
        KSPAppendOptionsPrefix(linear_solver->ksp(),"ns_");
        KSPSetFromOptions(linear_solver->ksp());
        PCFieldSplitSetIS(linear_solver->ksp()->pc,"v",is_v_local);
        PCFieldSplitSetIS(linear_solver->ksp()->pc,"p",is_p_local);
        KSPSetFromOptions(linear_solver->ksp());
       std::cout << "  done! " << std::flush;

       int size;
       ISGetSize(is_v_local, &size);
       std::cout << " is_v_local size: " << size << std::flush;
       ISGetSize(is_p_local, &size);
       std::cout << ", is_p_local size: " << size << std::flush;
    }

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
