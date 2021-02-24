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
#include "libmesh/mesh_generation.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"
#include "libmesh/mesh_refinement.h"


#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_linear_solver.h"

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void advance(libMesh::EquationSystems & es,
                      const std::string & system_name);


namespace Parameters
{
    double E = 250.0;
    double nu = 0.3;
    double mu = E / (2* (1+nu));
    double lambda = E*nu / (1+nu) / (1-2*nu);
    double dt = 0.0;
    double time = 0.0;
    double rho = 1.0;
    double alpha = 0;
    double gamma = 0;
    double t0 = 0.0;
    double t1 = 6.26;

    int bottom = 0;
    int right = 1;
    int top = 2;
    int left = 3;
}
    int main (int argc, char ** argv)
{
        using namespace libMesh;
	BeatIt::printBanner(std::cout);
	GetPot data("input.beat");

	// Initialize the library.  This is necessary because the library
	// may depend on a number of other libraries (i.e. MPI and PETSc)
	// that require initialization before use.  When the LibMeshInit
	// object goes out of scope, other libraries and resources are
	// finalized.
	libMesh::LibMeshInit init (argc, argv);

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
	const unsigned int dim = data("DIM", 2);

	// Create a mesh, with dimension to be overridden later, on the
	// default MPI communicator.
	libMesh::Mesh mesh(init.comm());

	// We may need XDR support compiled in to read binary .xdr files
	int nx = data("nx", 2);
	int ny = data("ny",1);
	double maxx = data("maxx", 2.0);
	double maxy = data("maxy", 1.0);

	// Read the input mesh.
	if(ny > 0)
	{
	  MeshTools::Generation::build_square (mesh,
	                                       nx, ny,
	                                       0., maxx,
	                                       0., maxy,
	                                       TRI3);
	}
	else
	{
	    std::string meshfile = data("mesh","NONE");
	    mesh.read(meshfile);
	    MeshRefinement refinement(mesh);
	    int nrefs = data("nrefs",0);
	    refinement.uniformly_refine(nrefs);
	}
	// Print information about the mesh to the screen.
	mesh.print_info();

	libMesh::EquationSystems systems(mesh);

	libMesh::TransientLinearImplicitSystem& primal_momentum_system = systems.add_system<libMesh::TransientLinearImplicitSystem>("primal_momentum");
	unsigned int primal_v1_var = primal_momentum_system.add_variable("v1", libMesh::FIRST);
	unsigned int primal_v2_var = primal_momentum_system.add_variable("v2", libMesh::FIRST);
	primal_momentum_system.add_vector("force");
	libMesh::TransientLinearImplicitSystem& primal_disp_system = systems.add_system<libMesh::TransientLinearImplicitSystem>("primal_disp");
	unsigned int primal_u1_var = primal_disp_system.add_variable("u1", libMesh::FIRST);
	unsigned int primal_u2_var = primal_disp_system.add_variable("u2", libMesh::FIRST);

    libMesh::TransientLinearImplicitSystem& mixed_momentum_system = systems.add_system<libMesh::TransientLinearImplicitSystem>("mixed_momentum");
    unsigned int mixed_v1_var = mixed_momentum_system.add_variable("mv1", libMesh::FIRST);
    unsigned int mixed_v2_var = mixed_momentum_system.add_variable("mv2", libMesh::FIRST);
    mixed_momentum_system.add_vector("force");
    libMesh::TransientLinearImplicitSystem& pressure_system = systems.add_system<libMesh::TransientLinearImplicitSystem>("pressure");
    unsigned int p_var = pressure_system.add_variable("p", libMesh::FIRST);
    libMesh::TransientLinearImplicitSystem& mixed_disp_system = systems.add_system<libMesh::TransientLinearImplicitSystem>("mixed_disp");
    unsigned int mixed_u1_var = mixed_disp_system.add_variable("mu1", libMesh::FIRST);
    unsigned int mixed_u2_var = mixed_disp_system.add_variable("mu2", libMesh::FIRST);

    // Initialize the data structures for the equation system.
    systems.init();
    // Prints information about the system to the screen.
    systems.print_info();

    Parameters::dt = data("dt", 0.01);
    Parameters::time = 0.0;
    int num_iter = data("num_iter",10);
    int save_iter = data("save_iter", 1);
    int save_iter2 = data("save_iter", 1);

    //BC
    Parameters::bottom = data("bottom", 0);
    Parameters::right = data("right", 1);
    Parameters::top = data("top", 2);
    Parameters::left = data("left", 3);

    std::cout << "BCs: " << Parameters::bottom << ", " << Parameters::right << ", "
                         << Parameters::top    << ", " << Parameters::left << std::endl;
    std::string output_file = data("output_file", "out.exo");
    libMesh::ExodusII_IO exo(mesh);
    libMesh::ExodusII_IO exo2(mesh);
    exo.write_equation_systems ("nonsense.exo", systems);
    exo2.write_equation_systems (output_file, systems);

    int iter = 0;
    exo.write_timestep("nonsense.exo", systems, save_iter, Parameters::time);
    int export_iter = 1;
    exo2.write_timestep(output_file, systems, export_iter, Parameters::time);

    std::cout << "TIME LOOP" << std::endl;
    for( ; iter < num_iter; )
    {
        Parameters::time += Parameters::dt;
        iter++;

        //advance
        advance(systems,"");

        if(iter % save_iter == 0)
        {
            save_iter ++;
            exo.write_timestep("nonsense.exo", systems, save_iter, Parameters::time);
        }
        if(iter % save_iter2 == 0)
        {
            std::cout << "TIME: " << Parameters::time << ", iter: " << iter << std::endl;
            std::cout << "saving" << std::endl;
            export_iter ++;
            exo2.write_timestep(output_file, systems, export_iter, Parameters::time);
        }
    }

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
void advance( libMesh::EquationSystems & es,
                      const std::string & system_name)
{
    using namespace libMesh;
	using std::unique_ptr;
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "SimpleSystem");


  // Get a constant reference to the mesh object.
  const  libMesh::MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  libMesh::TransientLinearImplicitSystem& primal_momentum_system = es.get_system<libMesh::TransientLinearImplicitSystem> ("primal_momentum");
  libMesh::TransientLinearImplicitSystem& primal_disp_system = es.get_system<libMesh::TransientLinearImplicitSystem>("primal_disp");
  libMesh::TransientLinearImplicitSystem& mixed_momentum_system = es.add_system<libMesh::TransientLinearImplicitSystem>("mixed_momentum");
  libMesh::TransientLinearImplicitSystem& mixed_disp_system = es.get_system<libMesh::TransientLinearImplicitSystem>("mixed_disp");
  libMesh::TransientLinearImplicitSystem& pressure_system = es.get_system<libMesh::TransientLinearImplicitSystem>("pressure");


  const unsigned int u1_var = primal_momentum_system.variable_number ("v1");
  const unsigned int u2_var = primal_momentum_system.variable_number ("v2");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const libMesh::DofMap & primal_dof_map = primal_momentum_system.get_dof_map();
  const libMesh::DofMap & mixed_dof_map = mixed_momentum_system.get_dof_map();
  const libMesh::DofMap & disp_dof_map = primal_disp_system.get_dof_map();
  const libMesh::DofMap & pressure_dof_map = pressure_system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  libMesh::FEType fe_type = primal_dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.  Introduction Example 4
  // describes some advantages of  std::unique_ptr's in the context of
  // quadrature rules.
  std::unique_ptr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  libMesh::QGauss qrule (dim, libMesh::THIRD);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  std::unique_ptr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));

  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  libMesh::QGauss qface(dim-1, libMesh::THIRD);
  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector< libMesh::Real> & JxW = fe->get_JxW();
  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector< libMesh::Point> & q_point = fe->get_xyz();
  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector< libMesh::Real> > & phi = fe->get_phi();
  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector< libMesh::RealGradient> > & dphi = fe->get_dphi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  libMesh::DenseMatrix< libMesh::Number> Me;
  libMesh::DenseVector< libMesh::Number> Fe; // primal
  libMesh::DenseVector< libMesh::Number> Fme;// mixed
  libMesh::DenseVector< libMesh::Number> Fp;//  pressure
  libMesh::DenseMatrix< libMesh::Number> Mpp;

  libMesh::DenseSubVector<Number>
    Fu1(Fe),
    Fu2(Fe);
  libMesh::DenseSubVector<Number>
    Fmu1(Fme),
    Fmu2(Fme);
  libMesh::DenseSubMatrix<Number> Mu1u1(Me), Mu2u2(Me);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector< libMesh::dof_id_type> primal_dof_indices;
  std::vector< libMesh::dof_id_type> mixed_dof_indices;
  std::vector< libMesh::dof_id_type> disp_dof_indices;
  std::vector< libMesh::dof_id_type> pressure_dof_indices;

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
  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const  libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  // Loop over the elements.  Note that  ++el is preferred to
  // el++ since the latter requires an unnecessary temporary
  // object.
  for ( ; el != end_el ; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const  libMesh::Elem * elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      primal_dof_map.dof_indices (elem, primal_dof_indices);
      pressure_dof_map.dof_indices (elem, pressure_dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).

      // The  DenseMatrix::resize() and the  DenseVector::resize()
      // members will automatically zero out the matrix  and vector.

      auto n_primal_dofs = primal_dof_indices.size();
      auto n_pressure_dofs = pressure_dof_indices.size();
      Me.resize (n_primal_dofs,n_primal_dofs);
      Mpp.resize (n_pressure_dofs,n_pressure_dofs);
      Fe.resize (n_primal_dofs);
      Fme.resize(n_primal_dofs);
      Fp.resize (n_pressure_dofs);

      Mu1u1.reposition (u1_var*n_pressure_dofs, u1_var*n_pressure_dofs, n_pressure_dofs, n_pressure_dofs);
      Mu2u2.reposition (u2_var*n_pressure_dofs, u2_var*n_pressure_dofs, n_pressure_dofs, n_pressure_dofs);
      Fu1.reposition (u1_var*n_pressure_dofs, n_pressure_dofs);
      Fu2.reposition (u2_var*n_pressure_dofs, n_pressure_dofs);
      Fmu1.reposition (u1_var*n_pressure_dofs, n_pressure_dofs);
      Fmu2.reposition (u2_var*n_pressure_dofs, n_pressure_dofs);


      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test funcions (i) against
          // the trial functions (j).
          for (unsigned int i=0; i<phi.size(); i++)
          {
            for (unsigned int j=0; j<phi.size(); j++)
              {
                Mu1u1(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
                Mu2u2(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
                Mpp(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
              }
          }
          // This is the end of the matrix summation loop
          // Now we build the element right-hand-side contribution.
          // This involves a single loop in which we integrate the
          // "forcing function" in the PDE against the test functions.

          // stress component
          {
            libMesh::Gradient grad_u1 = primal_momentum_system.point_gradient(0, q_point[qp], elem);
            libMesh::Gradient grad_u2 = primal_momentum_system.point_gradient(1, q_point[qp], elem);
            libMesh::Tensor E, grad_u, sigma, eye;
            eye(0,0) = 1; eye(1,1) = 1; eye(2,2) = 1;
            grad_u(0,0) = grad_u1(0); grad_u(0,1) = grad_u1(1);
            grad_u(1,0) = grad_u2(0); grad_u(1,1) = grad_u2(1);
            E = 0.5 * ( grad_u + grad_u.transpose() );
            sigma = 2.0 * Parameters::mu * E + Parameters::lambda * E.tr() * eye;

            libMesh::Gradient grad_mu1= mixed_momentum_system.point_gradient(0, q_point[qp], elem);
            libMesh::Gradient grad_mu2= mixed_momentum_system.point_gradient(1, q_point[qp], elem);
            double pressure = pressure_system.point_value(0, q_point[qp], elem);

            libMesh::Tensor Em, grad_mu, sigma_m;
            grad_mu(0,0) = grad_mu1(0); grad_mu(0,1) = grad_mu1(1);
            grad_mu(1,0) = grad_mu2(0); grad_mu(1,1) = grad_mu2(1);
            Em = 0.5 * ( grad_mu + grad_mu.transpose() );
            sigma_m = 2 * Parameters::mu * ( Em - 1 / 3.0 * Em.tr() * eye) + pressure * eye;
            double p_rhs = ( 2.0 / 3.0 * Parameters::mu + Parameters::lambda ) * Em.tr();
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

            for (unsigned int i=0; i<phi.size(); i++)
            {
              auto F_qp = - sigma * dphi[i][qp] * JxW[qp];
              auto Fm_qp= - sigma_m* dphi[i][qp] * JxW[qp];
              Fu1(i) += F_qp(0);
              Fu2(i) += F_qp(1);
              Fmu1(i)+= Fm_qp(0);
              Fmu2(i)+= Fm_qp(1);
              Fp(i)  += JxW[qp]*phi[i][qp] * p_rhs;
            }
          }

          // FORCE
          {
            const  libMesh::Real x = q_point[qp](0);
            const  libMesh::Real y = q_point[qp](1);
            const  libMesh::Real eps = 1.e-3;

            double fx =  2.0 * Parameters::alpha * x;
            double fy =  2.0 * Parameters::gamma * x;

            for (unsigned int i=0; i<phi.size(); i++)
            {
              Fu1(i) += JxW[qp] * fx * phi[i][qp];
              Fu2(i) += JxW[qp] * fy * phi[i][qp];
              Fmu1(i) += JxW[qp] * fx * phi[i][qp];
              Fmu2(i) += JxW[qp] * fy * phi[i][qp];
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
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              auto boundary_id = mesh.get_boundary_info().boundary_id(elem, side);
              //  0 bottom
              //  1 right
              //  2 top
              //  3 left
              libMesh::Gradient traction;
              bool dirichlet = true;
              if(boundary_id == Parameters::bottom)
              {
                  traction(0) = 0.0;
                  traction(1) = 0.0;
                  dirichlet = false;
              }
              else if (boundary_id == Parameters::right)
              {
                  traction(0) = Parameters::t0;
                  traction(1) = Parameters::t1;
                  dirichlet = false;
              }
              else if (boundary_id == Parameters::top)
              {
                  traction(0) = 0.0;
                  traction(1) = 0.0;
                  dirichlet = false;
              }
              // The value of the shape functions at the quadrature
              // points.
              const std::vector<std::vector< libMesh::Real> > & phi_face = fe_face->get_phi();

              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector< libMesh::Real> & JxW_face = fe_face->get_JxW();

              // The XYZ locations (in physical space) of the
              // quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector< libMesh::Point> & qface_point = fe_face->get_xyz();

              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, side);

              // Loop over the face quadrature points for integration.
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // The location on the boundary of the current
                  // face quadrature point.
                  const  libMesh::Real xf = qface_point[qp](0);
                  const  libMesh::Real yf = qface_point[qp](1);

                  // The penalty value.  \frac{1}{\epsilon}
                  // in the discussion above.
                  if(dirichlet)
                  {
                      const   libMesh::Real penalty = 1.e10;
                      const  libMesh::Real value = 0.0;
                      for (unsigned int i=0; i<phi_face.size(); i++)
                      {
                        for (unsigned int j=0; j<phi_face.size(); j++)
                        {
                              Mu1u1(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
                              Mu2u2(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
                        }
                      }
                  }
                  else
                  {
                      for (unsigned int i=0; i<phi.size(); i++)
                      {
                        auto F_qp = traction  * phi_face[i][qp] * JxW_face[qp];
                        auto Fm_qp= traction* phi_face[i][qp] * JxW_face[qp];
                        Fu1(i) += F_qp(0);
                        Fu2(i) += F_qp(1);
                        Fmu1(i)+= Fm_qp(0);
                        Fmu2(i)+= Fm_qp(1);
                      }
                  }
            }
        }
      }
      primal_momentum_system.matrix->add_matrix (Me, primal_dof_indices);
      mixed_momentum_system.matrix->add_matrix (Me, primal_dof_indices);
      primal_momentum_system.rhs->add_vector    (Fe, primal_dof_indices);
      mixed_momentum_system.rhs->add_vector    (Fme, primal_dof_indices);
      pressure_system.matrix->add_matrix    (Mpp, pressure_dof_indices);
      pressure_system.rhs->add_vector    (Fp, pressure_dof_indices);
    }

  primal_momentum_system.matrix->close();
  mixed_momentum_system.matrix->close();
  pressure_system.matrix->close();
  primal_momentum_system.rhs->close();
  mixed_momentum_system.rhs->close();
  pressure_system.rhs->close();

  PetscBool rtol_set;
  double runtime_rtol;
  int ierr;
  double tol = 1e-6;
  int max_its = 100;
  ierr = PetscOptionsGetReal(nullptr, "", "-ksp_rtol", &runtime_rtol, &rtol_set);
  PetscBool max_it_set;
  int runtime_max_it;
  ierr = PetscOptionsGetInt(nullptr, "", "-ksp_max_it", &runtime_max_it, &max_it_set);

//primal solve
  {
      using namespace libMesh;
      PetscLinearSolver<double> * solver = dynamic_cast< PetscLinearSolver<double> * >(primal_momentum_system.get_linear_solver());
      PetscMatrix<double> * M_mat = dynamic_cast<PetscMatrix<double>*>(primal_momentum_system.matrix);
      ierr = KSPSetFromOptions( solver->ksp());
      solver->solve(*M_mat, *M_mat, primal_momentum_system.get_vector("force"), *primal_momentum_system.rhs, rtol_set ? runtime_rtol : tol, max_it_set ? runtime_max_it : max_its);
  }
  //mixed solve
  {
      using namespace libMesh;
      PetscLinearSolver<double> * solver = dynamic_cast< PetscLinearSolver<double> * >(mixed_momentum_system.get_linear_solver());
      PetscMatrix<double> * M_mat = dynamic_cast<PetscMatrix<double>*>(mixed_momentum_system.matrix);
      ierr = PetscOptionsGetInt(nullptr, "", "-ksp_max_it", &runtime_max_it, &max_it_set);
      ierr = KSPSetFromOptions( solver->ksp());
      solver->solve(*M_mat, *M_mat, mixed_momentum_system.get_vector("force"), *mixed_momentum_system.rhs, rtol_set ? runtime_rtol : tol, max_it_set ? runtime_max_it : max_its);
  }
  //pressure solve
  {
      using namespace libMesh;
      PetscLinearSolver<double> * solver = dynamic_cast< PetscLinearSolver<double> * >(pressure_system.get_linear_solver());
      PetscMatrix<double> * M_mat = dynamic_cast<PetscMatrix<double>*>(pressure_system.matrix);
      ierr = PetscOptionsGetInt(nullptr, "", "-ksp_max_it", &runtime_max_it, &max_it_set);
      ierr = KSPSetFromOptions( solver->ksp());
      solver->solve(*M_mat, *M_mat, *pressure_system.solution, *pressure_system.rhs, rtol_set ? runtime_rtol : tol, max_it_set ? runtime_max_it : max_its);
  }

  // close solution
  primal_momentum_system.solution->close();
  mixed_momentum_system.solution->close();
  primal_disp_system.solution->close();
  mixed_disp_system.solution->close();
  // update disp
  primal_momentum_system.solution->add(Parameters::dt/Parameters::rho,  primal_momentum_system.get_vector("force"));
  mixed_momentum_system.solution->add(Parameters::dt/Parameters::rho,  mixed_momentum_system.get_vector("force"));
  primal_disp_system.solution->add(Parameters::dt, *primal_momentum_system.solution);
  mixed_disp_system.solution->add(Parameters::dt, *mixed_momentum_system.solution);
  // All done!
}





