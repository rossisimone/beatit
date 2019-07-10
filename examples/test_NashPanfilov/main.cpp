/*
 ============================================================================

 .______    _______     ___   .___________.    __  .___________.
 |   _  \  |   ____|   /   \  |           |   |  | |           |
 |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
 |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |
 |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |
 |______/  |_______/__/     \__\  |__|        |__|     |__|

 BeatIt - code for cardiovascular simulations
 Copyright (C) 2016 Simone Rossi

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ============================================================================
 */

/**
 * \file main.cpp
 *
 * \brief Here we test GetPot
 *
 * \author srossi
 *
 * \version 0.0
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 7, 2016
 *
 */

#include "libmesh/getpot.h"
#include "Util/IO/io.hpp"
#include "Util/Timer.hpp"
// Functions to initialize the library.
#include "libmesh/libmesh.h"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"

// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/vector_value.h"
#include "libmesh/gmv_io.h"
//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"

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

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"


#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"

#include "libmesh/vtk_io.h"
// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
struct Poisson
{
public:
static void assemble_poisson(
        libMesh::EquationSystems & es,
        const std::string & system_name);
};
void assemble_r(
        libMesh::EquationSystems & es,
        const std::string & system_name);

void IC(libMesh::EquationSystems & es, const std::string & system_name);
void ICr(libMesh::EquationSystems & es, const std::string & system_name);

int main(int argc, char ** argv)
{
    BeatIt::printBanner(std::cout);
    GetPot data("data.pot");

    // Initialize the library.  This is necessary because the library
    // may depend on a number of other libraries (i.e. MPI and PETSc)
    // that require initialization before use.  When the LibMeshInit
    // object goes out of scope, other libraries and resources are
    // finalized.
    libMesh::LibMeshInit init(argc, argv);

    // Check for proper usage. The program is designed to be run
    // as follows:
    // ./ex1 -d DIM input_mesh_name [-o output_mesh_name]
    // where [output_mesh_name] is an optional parameter giving
    // a filename to write the mesh into.
//	if (argc < 4)
//	{
//	  // This handy function will print the file name, line number,
//	  // specified message, and then throw an exception.
//	  libmesh_error_msg("Usage: " << argv[0] << " -d 2 in.mesh [-o out.mesh]");
//	}

    // Get the dimensionality of the mesh from argv[2]
    const unsigned int dim = data("DIM", 3);

    // Create a mesh, with dimension to be overridden later, on the
    // default MPI communicator.
    libMesh::Mesh mesh(init.comm());

    // We may need XDR support compiled in to read binary .xdr files
    std::string meshfile = data("input_mesh_name", "Pippo.e");

    // Read the input mesh.
    mesh.read(&meshfile[0]);

    // Print information about the mesh to the screen.
    mesh.print_info();


    libMesh::EquationSystems systems(mesh);

    systems.parameters.set<bool>("test") = true;
    typedef libMesh::Real Real;
    systems.parameters.set<Real>("mu1") = 0.12;
    systems.parameters.set<Real>("mu2") = 0.3;
    systems.parameters.set<Real>("k") = 8.0;
    systems.parameters.set<Real>("a") = 0.1;
    systems.parameters.set<Real>("eps") = 0.01;
    systems.parameters.set<Real>("e0") = 1;
    systems.parameters.set<Real>("D") = 1;

    systems.add_system<libMesh::TransientLinearImplicitSystem>("SimpleSystem");
    unsigned int V_var = systems.get_system("SimpleSystem").add_variable("V",
            libMesh::FIRST);

    // Construct a Dirichlet boundary condition object
    // We impose a "clamped" boundary condition on the
    // "left" boundary, i.e. bc_id = 3
//    std::set<libMesh::boundary_id_type> boundary_ids;
//    boundary_ids.insert(1);
//    boundary_ids.insert(2);
//    boundary_ids.insert(3);
//    boundary_ids.insert(4);
//    std::vector<unsigned int> variables(1);
//    variables[0] = V_var;
//
//    // Create a ZeroFunction to initialize dirichlet_bc
//    libMesh::ZeroFunction<> zf;
//    libMesh::DirichletBoundary dirichlet_bc(boundary_ids, variables, &zf);

    // Give the system a pointer to the matrix assembly
    // function.  This will be called when needed by the
    // library.
    systems.get_system("SimpleSystem").attach_assemble_function(
            &Poisson::assemble_poisson);
    systems.get_system("SimpleSystem").attach_init_function (IC);

    systems.add_system<libMesh::TransientExplicitSystem>("ComplexSystem");
    systems.get_system("ComplexSystem").add_variable("r", libMesh::CONSTANT,
            libMesh::MONOMIAL);
    systems.get_system("ComplexSystem").attach_init_function (ICr);
    systems.get_system("ComplexSystem").attach_assemble_function(
            assemble_r);


    // Define the mesh refinement object that takes care of adaptively
    // refining the mesh.
   libMesh:: MeshRefinement mesh_refinement(mesh);
   // These parameters determine the proportion of elements that will
   // be refined and coarsened. Any element within 30% of the maximum
   // error on any element will be refined, and any element within 30%
   // of the minimum error on any element might be coarsened
   mesh_refinement.refine_fraction()  =  data("refine_fraction", 0.95);
   mesh_refinement.coarsen_fraction() =data("coarsen_fraction", 0.05);
   // We won't refine any element more than 5 times in total
   int max_num_mesh_ref = data("max_num_mesh_ref", 0);
   mesh_refinement.max_h_level()      = max_num_mesh_ref;
   // Refinement parameters
   const unsigned int max_r_steps = max_num_mesh_ref; // Refine the mesh 5 times

    // We must add the Dirichlet boundary condition _before_
    // we call equation_systems.init()
//    systems.get_system("SimpleSystem").get_dof_map().add_dirichlet_boundary(
//            dirichlet_bc);

    // Initialize the data structures for the equation system.
    systems.init();


    // Solve the system "Convection-Diffusion".  This will be done by
    // looping over the specified time interval and calling the
    // solve() member at each time step.  This will assemble the
    // system and call the linear solver.
    Real dt = data("dt", 0.1);

    systems.get_system("SimpleSystem").time = 0.;


    // Prints information about the system to the screen.
    systems.print_info();
    // Write the system.
    std::string output_system = data("output_systems", "");
    if ("" != output_system)
        systems.write(&output_system[0]);
    // Write the output mesh if the user specified an
    // output file name.
    std::string output_file = data("output_mesh", "");
    if ("" != output_file)
        mesh.write(&output_file[0]);

    // TIME LOOP
    std::string exodus_filename = data("exodus_output_file", "transient.e");
    libMesh::ExodusII_IO exo(mesh);
    exo.write_equation_systems (exodus_filename, systems);
    exo.append(true);
    // then append element-based discontinuous plots of the stresses
   // exo.write_element_data(systems);
    std::string gmv_output_file = data("gmv_output_file", "out.gmv");
    libMesh::GMVIO(mesh).write_equation_systems (gmv_output_file+".00000", systems);


    libMesh::TransientLinearImplicitSystem & system = systems.get_system<
            libMesh::TransientLinearImplicitSystem>("SimpleSystem");
    libMesh::TransientExplicitSystem & system_r = systems.get_system<
            libMesh::TransientExplicitSystem>("ComplexSystem");

    unsigned int t_step = 0;

    BeatIt::Timer timer;
    int numSteps = data("numSteps", 100);
    int saveSteps = data("saveSteps", 100);
    int AMR_steps = data("AMR_steps", 1);
    for (; t_step < numSteps; t_step++)
    {
        timer.start();


//        std::ostringstream file_name;
//
//        file_name << "out.pvtu."
//                  << std::setw(3)
//                  << std::setfill('0')
//                  << std::right
//                  << t_step+1;
//        libMesh::VTKIO(mesh).write_equation_systems (file_name.str(), systems);

        // Incremenet the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        systems.get_system("SimpleSystem").time += dt;
        systems.get_system("ComplexSystem").time += dt;

        systems.parameters.set<libMesh::Real> ("time") = systems.get_system("SimpleSystem").time;
        systems.parameters.set<libMesh::Real> ("dt")   = dt;

        // A pretty update message
        libMesh::out << " Solving time step ";
        // Do fancy zero-padded formatting of the current time.
        {
          std::ostringstream out;

          out << std::setw(2)
              << std::right
              << t_step
              << ", time="
              << std::fixed
              << std::setw(6)
              << std::setprecision(3)
              << std::setfill('0')
              << std::left
              << systems.get_system("SimpleSystem").time
              <<  "...";

          libMesh::out << out.str() << std::endl;
        }

        // At this point we need to update the old
        // solution vector.  The old solution vector
        // will be the current solution vector from the
        // previous time step.  We will do this by extracting the
        // system from the EquationSystems object and using
        // vector assignment.  Since only TransientSystems
        // (and systems derived from them) contain old solutions
        // we need to specify the system type when we ask for it.

        *system.old_local_solution = *systems.get_system("SimpleSystem").current_local_solution;
        *system_r.old_local_solution = *systems.get_system("ComplexSystem").current_local_solution;

        // Define the refinement loop
        if( 0 == t_step%AMR_steps )
        {
          for (unsigned int r_step=0; r_step<=max_r_steps; r_step++)
          {
              // Assemble & solve the linear system
              system_r.solve();

              systems.get_system("SimpleSystem").solve();

              // We need to ensure that the mesh is not refined on the last iteration
              // of this loop, since we do not want to refine the mesh unless we are
              // going to solve the equation system for that refined mesh.
              if (r_step != max_r_steps)
              {
                  // Error estimation objects, see Adaptivity Example 2 for details
                  libMesh::ErrorVector error;
                  libMesh::KellyErrorEstimator error_estimator;
                  //libMesh::LaplacianErrorEstimator error_estimator;

                  // Compute the error for each active element
                  error_estimator.estimate_error(system, error);

                  // Output error estimate magnitude
                  libMesh::out << "Error estimate\nl2 norm = "
                               << error.l2_norm()
                               << "\nmaximum = "
                               << error.maximum()
                               << std::endl;

                  // Flag elements to be refined and coarsened
                  mesh_refinement.flag_elements_by_error_fraction (error);

                  // Perform refinement and coarsening
                  mesh_refinement.flag_elements_by_error_fraction (error);
                  mesh_refinement.refine_and_coarsen_elements();
                  // Reinitialize the equation_systems object for the newly refined
                  // mesh. One of the steps in this is project the solution onto the
                  // new mesh
                  systems.reinit();
                }
          }


          }

          timer.stop();

			  if( 0 == t_step%saveSteps )
			  {
				  std::cout << "Exporting to EXODUS" << std::endl;
				  exo.write_timestep (exodus_filename, systems, t_step+1, systems.get_system("SimpleSystem").time);
				  //exo.write_element_data(systems);
			  }
			  if( 0 == t_step%saveSteps )
			  {
				          std::ostringstream file_name;

				          file_name << gmv_output_file + "."
				                    << std::setw(5)
				                    << std::setfill('0')
				                    << std::right
				                    << t_step+1;
						  std::cout << "Exporting to GMV" << std::endl;
				          libMesh::GMVIO(mesh).write_equation_systems (file_name.str(), systems);
			  }

        }

    timer.print(std::cout);


//	  libMesh::VTKIO (mesh).write_equation_systems ("out.pvtu", systems);
    // LibMeshInit object was created first, its destruction occurs
    // last, and it's destructor finalizes any external libraries and
    // checks for leaked memory.
    return 0;
}

// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void
Poisson::assemble_poisson(
        libMesh::EquationSystems & es,
        const std::string & system_name)
{
    using libMesh::UniquePtr;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "SimpleSystem");

    // Get a constant reference to the mesh object.
    const libMesh::MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientLinearImplicitSystem & system = es.get_system<
            libMesh::TransientLinearImplicitSystem>("SimpleSystem");

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientExplicitSystem & system_r = es.get_system<
            libMesh::TransientExplicitSystem>("ComplexSystem");

    const unsigned int V_Var = system.variable_number("V");

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const libMesh::DofMap & dof_map = system.get_dof_map();
    const libMesh::DofMap & dof_map_r = system_r.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    libMesh::FEType fe_type = dof_map.variable_type(0);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a UniquePtr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  UniquePtr's in the context of
    // quadrature rules.
    UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule(dim, libMesh::FIRST);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Declare a special finite element object for
    // boundary integration.
    UniquePtr<libMesh::FEBase> fe_face(libMesh::FEBase::build(dim, fe_type));

    // Boundary integration requires one quadraure rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    libMesh::QGauss qface(dim - 1, libMesh::FIRST);

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

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
    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseVector<libMesh::Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_r(1);

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
        dof_map_r.dof_indices(elem, dof_indices_r);

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
        Ke.resize(dof_indices.size(), dof_indices.size());

        Fe.resize(dof_indices.size());
        const libMesh::Real dt = es.parameters.get<libMesh::Real> ("dt");
        const libMesh::Real mu1 = es.parameters.get<libMesh::Real> ("mu1");
        const libMesh::Real mu2 = es.parameters.get<libMesh::Real> ("mu2");
        const libMesh::Real eps = es.parameters.get<libMesh::Real> ("eps");
        const libMesh::Real k = es.parameters.get<libMesh::Real> ("k");
        const libMesh::Real a = es.parameters.get<libMesh::Real> ("a");
        const libMesh::Real D = es.parameters.get<libMesh::Real> ("D");
        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        double r = (*system_r.current_local_solution)(dof_indices_r[0]);

        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

            libMesh::Number V = 0;
            for(unsigned int i = 0; i < phi.size(); ++i )
            {
                V += phi[i][qp] *  system.old_solution(dof_indices[i]);
            }

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test funcions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i < phi.size(); i++)
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    // Mass term
                    Ke(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]) / dt;
                    // stiffness term
                    Ke(i, j) += JxW[qp] * D * (dphi[i][qp] * dphi[j][qp]);
                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            {
                const libMesh::Real x = q_point[qp](0);
                const libMesh::Real y = q_point[qp](1);
                const libMesh::Real eps = 1.e-3;

                // "fxy" is the forcing function for the Poisson equation.
                // In this case we set fxy to be a finite difference
                // Laplacian approximation to the (known) exact solution.
                //
                // We will use the second-order accurate FD Laplacian
                // approximation, which in 2D is
                //
                // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
                //                u(i-1,j) + u(i+1,j) +
                //                -4*u(i,j))/h^2
                //
                // Since the value of the forcing function depends only
                // on the location of the quadrature point (q_point[qp])
                // we will compute it here, outside of the i-loop
                const libMesh::Real R = - k * V * (V - 1.0) * (V - a) - r * V;
                //const libMesh::Real R = V * (V - 1.0);

                for (unsigned int i = 0; i < phi.size(); i++)
                {
                    // Mass term
                    Fe(i) += JxW[qp] * V * phi[i][qp] / dt;
                    // Reaction term
                    Fe(i) += JxW[qp] * R * phi[i][qp];
                }
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
        // interpolatory bases like the standard Lagrange polynomials,
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
            for (unsigned int side = 0; side < elem->n_sides(); side++)
                if (elem->neighbor_ptr(side) == libmesh_nullptr)
                {
                    // The value of the shape functions at the quadrature
                    // points.
                    const std::vector<std::vector<libMesh::Real> > & phi_face =
                            fe_face->get_phi();

                    // The Jacobian * Quadrature Weight at the quadrature
                    // points on the face.
                    const std::vector<libMesh::Real> & JxW_face =
                            fe_face->get_JxW();

                    // The XYZ locations (in physical space) of the
                    // quadrature points on the face.  This is where
                    // we will interpolate the boundary value function.
                    const std::vector<libMesh::Point> & qface_point =
                            fe_face->get_xyz();

                    // Compute the shape function values on the element
                    // face.
                    fe_face->reinit(elem, side);

                    // Loop over the face quadrature points for integration.
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        // The location on the boundary of the current
                        // face quadrature point.
                        const libMesh::Real xf = qface_point[qp](0);
                        const libMesh::Real yf = qface_point[qp](1);

                        // The penalty value.  \frac{1}{\epsilon}
                        // in the discussion above.
                        const libMesh::Real penalty = 1.e10;

                        // The boundary value.
                        const libMesh::Real value = 0.0;

//                  // Matrix contribution of the L2 projection.
//                  for (unsigned int i=0; i<phi_face.size(); i++)
//                    for (unsigned int j=0; j<phi_face.size(); j++)
//                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
//
//                  // Right-hand-side contribution of the L2
//                  // projection.
//                  for (unsigned int i=0; i<phi_face.size(); i++)
//                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                    }
                }
        }

        // We have now finished the quadrature point loop,
        // and have therefore applied all the boundary conditions.

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The  SparseMatrix::add_matrix()
        // and  NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }

    // All done!
}

libMesh::Number initial_value(
                                const libMesh::Point & p,
                                const libMesh::Parameters & parameters,
                                const std::string &,
                                const std::string &)
{
    if( 40 < p(0) ) return 1.0;
    else return 0.0;
}

void IC(libMesh::EquationSystems & es, const std::string & system_name)
{

    // Get a reference to the Convection-Diffusion system object.
    libMesh::TransientLinearImplicitSystem & system = es.get_system<
            libMesh::TransientLinearImplicitSystem>(system_name);

    // Project initial conditions at time 0
    es.parameters.set<libMesh::Real>("time") = system.time = 0;

    system.project_solution(initial_value, libmesh_nullptr, es.parameters);
}


libMesh::Number initial_value_r(
                                const libMesh::Point & p,
                                const libMesh::Parameters & parameters,
                                const std::string &,
                                const std::string &)
{
   return 0.0;
}


void assemble_r(
        libMesh::EquationSystems & es,
        const std::string & system_name)
{
    using libMesh::UniquePtr;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "ComplexSystem");

    // Get a constant reference to the mesh object.
    const libMesh::MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientExplicitSystem & system_r = es.get_system<
            libMesh::TransientExplicitSystem>("ComplexSystem");


    libMesh::TransientLinearImplicitSystem & system = es.get_system<
            libMesh::TransientLinearImplicitSystem>("SimpleSystem");
    const unsigned int V_Var = system.variable_number("V");
    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const libMesh::DofMap & dof_map_r = system_r.get_dof_map();
    const libMesh::DofMap & dof_map = system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    libMesh::FEType fe_type = dof_map.variable_type(V_Var);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a UniquePtr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  UniquePtr's in the context of
    // quadrature rules.
    UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule(dim, libMesh::FIRST);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Declare a special finite element object for
    // boundary integration.
    UniquePtr<libMesh::FEBase> fe_face(libMesh::FEBase::build(dim, fe_type));

    // Boundary integration requires one quadraure rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    libMesh::QGauss qface(dim - 1, libMesh::CONSTANT);

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

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    libMesh::DenseVector<libMesh::Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices_r(1);
    std::vector<libMesh::dof_id_type> dof_indices;
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
    int elct = 0;
    const libMesh::Real dt = es.parameters.get<libMesh::Real> ("dt");
    const libMesh::Real mu1 = es.parameters.get<libMesh::Real> ("mu1");
    const libMesh::Real mu2 = es.parameters.get<libMesh::Real> ("mu2");
    const libMesh::Real eps = es.parameters.get<libMesh::Real> ("eps");
    const libMesh::Real k = es.parameters.get<libMesh::Real> ("k");
    const libMesh::Real a = es.parameters.get<libMesh::Real> ("a");

    for (; el != end_el; ++el)
    {
        // Store a pointer to the element we are currently
        // working on.  This allows for nicer syntax later.
        const libMesh::Elem * elem = *el;
    	 dof_map.dof_indices(elem, dof_indices);
         // Compute the element-specific data for the current
         // element.  This involves computing the location of the
         // quadrature points (q_point) and the shape functions
         // (phi, dphi) for the current element.
         fe->reinit(elem);
         double V = 0.0;
    	 for(unsigned int i = 0; i < phi.size(); ++i )
        {
            V +=  phi[i][0]*system.old_solution(dof_indices[i]);
        }

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map_r.dof_indices(elem, dof_indices_r);
        elct++;
        double r_old = (*system_r.current_local_solution)(dof_indices_r[0]);
        double fr = (eps + mu1 * r_old / (mu2+V)) * (-r_old-k*V*(V-a-1.));
        double r_new =  r_old+dt*fr;
        system_r.solution->set( dof_indices_r[0],r_new);
	}

   system_r.solution->close();
    // All done!
}





void ICr(libMesh::EquationSystems & es, const std::string & system_name)
{

    // Get a reference to the Convection-Diffusion system object.
    libMesh::TransientExplicitSystem& system = es.get_system<
            libMesh::TransientExplicitSystem>("ComplexSystem");

    // Project initial conditions at time 0
    es.parameters.set<libMesh::Real>("time") = system.time = 0;

    system.project_solution(initial_value_r, libmesh_nullptr, es.parameters);
}


