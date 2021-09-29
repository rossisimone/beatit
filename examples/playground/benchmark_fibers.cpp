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

void solve_poisson(EquationSystems & es);
void eval_fibers(EquationSystems & es);

namespace BCFunction
{
// Define a wrapper for exact_solution that will be needed below
void one (DenseVector<Number> & output,
                             const Point & p,
                             const Real)
{
  output(0) = 1.0;
}

// Define a wrapper for exact_solution that will be needed below
void zero (DenseVector<Number> & output,
                             const Point & p,
                             const Real)
{
  output(0) = 0.0;
}

}

// The main program.
int main(int argc, char ** argv)
{
    // Read input file
    GetPot command_line(argc, argv);
    std::string input_file_name = command_line.follow("input_file", 2, "-i", "--input");
    GetPot input_data(input_file_name);

    // Initialize libMesh.
    LibMeshInit init(argc, argv);

    // Create mesh
    Mesh mesh(init.comm());
    std::string mesh_filename = input_data("mesh", "NONE");
    mesh.read(mesh_filename);
    mesh.print_info();

    // Create equation system and systems
    EquationSystems equation_systems(mesh);
    LinearImplicitSystem & poisson_system = equation_systems.add_system<LinearImplicitSystem>("Poisson");
    poisson_system.add_variable("t", FIRST, LAGRANGE);
    ExplicitSystem & fibers_system = equation_systems.add_system<ExplicitSystem>("fibers");
    fibers_system.add_variable("fibersx", CONSTANT, MONOMIAL);
    fibers_system.add_variable("fibersy", CONSTANT, MONOMIAL);
    fibers_system.add_variable("fibersz", CONSTANT, MONOMIAL);
    ExplicitSystem & xfibers_system = equation_systems.add_system<ExplicitSystem>("xfibers");
    xfibers_system.add_variable("xfibersx", CONSTANT, MONOMIAL);
    xfibers_system.add_variable("xfibersy", CONSTANT, MONOMIAL);
    xfibers_system.add_variable("xfibersz", CONSTANT, MONOMIAL);
    ExplicitSystem & sheets_system = equation_systems.add_system<ExplicitSystem>("sheets");
    sheets_system.add_variable("sheetsx", CONSTANT, MONOMIAL);
    sheets_system.add_variable("sheetsy", CONSTANT, MONOMIAL);
    sheets_system.add_variable("sheetsz", CONSTANT, MONOMIAL);


    // Setup Dirichlet boundary conditions for Poisson problem
    std::vector<unsigned int> variables(1, 0);

    int bc_zero_id = input_data("bc_zero", 0);
    std::set<boundary_id_type> zero_boundary_id;
    zero_boundary_id.insert(bc_zero_id);
    AnalyticFunction<> zero_function(BCFunction::zero);
    DirichletBoundary zero_bc(zero_boundary_id, variables, zero_function);


    int bc_one_id = input_data("bc_one", 1);
    std::set<boundary_id_type> one_boundary_id;
    one_boundary_id.insert(bc_one_id);
    AnalyticFunction<> one_function(BCFunction::one);
    DirichletBoundary one_bc(one_boundary_id, variables, one_function);


    poisson_system.get_dof_map().add_dirichlet_boundary(zero_bc);
    poisson_system.get_dof_map().add_dirichlet_boundary(one_bc);

    // Initialize System
    std::cout << "Init System" << std::endl;
    equation_systems.init();


    // Prints information about the system to the screen.
    equation_systems.print_info();

    solve_poisson(equation_systems);

    eval_fibers(equation_systems);

    libMesh::ExodusII_IO exporter(mesh); \
    exporter.write_equation_systems ("poisson.e", equation_systems);
    exporter.write_element_data (equation_systems);
    return 0;
}

void solve_poisson(EquationSystems & es)
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
    const unsigned int t_var = system.variable_number("t");

    FEType fe_type = system.variable_type(t_var);

    // Build a Finite Element object of the specified type for
    // the velocity variables.
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    //QGauss qrule(dim, fe_vel_type.default_quadrature_order());
    QGauss qrule(dim, FIRST);

    // Tell the finite element objects to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
    const std::vector<std::vector<Number>> & phi = fe->get_phi();

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
        fe->reinit(elem);

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
            }
        } // end of the quadrature point qp-loop

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
    system.matrix->close();
    system.rhs->close();
    std::cout << "End Poisson Assembly" << std::endl;

    // Do not reassemble the matrix when calling solve()
    system.zero_out_matrix_and_rhs = false;
    system.assemble_before_solve = false;

    // Solve system
    system.solve();
}

namespace Parameter
{
    constexpr double degrees_to_radiant = M_PI / 180.0;
    constexpr double alpha_endo = 60.0 * degrees_to_radiant;
    constexpr double alpha_epi = -60.0 * degrees_to_radiant;

    constexpr double r_short_endo = 2.5e-2;
    constexpr double r_long_endo  = 9e-2;

    constexpr double r_short_epi = 3.5e-2;
    constexpr double r_long_epi  = 9.7e-2;
}


void eval_fibers(EquationSystems & es)
{
    std::cout << "Start Fibers Evaluation" << std::endl;
    // It is a good idea to make sure we are assembling
    // the proper system.

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    LinearImplicitSystem & poisson_system = es.get_system<LinearImplicitSystem>("Poisson");
    const unsigned int t_var = poisson_system.variable_number("t");

    ExplicitSystem & fibers_system = es.get_system<ExplicitSystem>("fibers");
    ExplicitSystem & xfibers_system = es.get_system<ExplicitSystem>("xfibers");
    ExplicitSystem & sheets_system = es.get_system<ExplicitSystem>("sheets");

    DofMap& fibers_dof_map = fibers_system.get_dof_map();
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> poisson_dof_indices;
    std::vector<dof_id_type> fibers_dof_indices;

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
        fibers_dof_map.dof_indices(elem, fibers_dof_indices);

        auto centroid = elem->centroid();
        double x = centroid(0);
        double y = centroid(1);
        double z = centroid(2);

        double t = poisson_system.point_value(t_var, centroid, elem);

        double alpha = ( Parameter::alpha_endo   + (Parameter::alpha_epi   - Parameter::alpha_endo  ) * t );
        double rl    = ( Parameter::r_long_endo  + (Parameter::r_long_epi  - Parameter::r_long_endo ) * t);
        double rs    = ( Parameter::r_short_endo + (Parameter::r_short_epi - Parameter::r_short_endo) * t);

        double mu    = -std::acos( x / rl );
        double theta = std::atan( z / y );
        // e_t = d x / d t
        //     = d/dt ( rl cos(mu), rs sin(mu) cos(tehta), rs sin(mu) sin(theta) )
        //     = ( drl/dt cos(mu), drs/dt sin(mu) cos(tehta), drs/dt sin(mu) sin(theta) )
        double drl_dt = Parameter::r_long_epi  - Parameter::r_long_endo;
        double drs_dt = Parameter::r_short_epi - Parameter::r_short_endo;
        Point et ( drl_dt * std::cos(mu),
                   drs_dt * std::sin(mu) * std::cos(theta),
                   drs_dt * std::sin(mu) * std::sin(theta));

        Point emu(-rl * std::sin(mu),
                   rs * std::cos(mu) * std::cos(theta),
                   rs * std::cos(mu) * std::sin(theta));


        Point etheta( 0.0,
                     -rs * std::sin(mu) * std::sin(theta),
                      rs * std::sin(mu) * std::cos(theta));
        if(y > 0)
        {
            et(1) *= -1.0;
            et(2) *= -1.0;

            emu(0) *= -1.0;
        }

        et = et.unit();
        emu = emu.unit();
        etheta = etheta.unit();

//        Point f = et;
//        Point n = emu;
//        Point s = etheta;

        Point f = std::sin(alpha) * emu + std::cos(alpha) * etheta ; f = f.unit();
        Point n = etheta.cross(emu); n = n.unit();
        Point s = f.cross(n); s = s.unit();

        for (int k = 0; k < 3; ++k)
        {
            fibers_system.solution->set(fibers_dof_indices[k],f(k));
            xfibers_system.solution->set(fibers_dof_indices[k],n(k));
            sheets_system.solution->set(fibers_dof_indices[k],s(k));
        }

    } // end of element loop
    std::cout << "End Fiber Evaluation" << std::endl;

}

