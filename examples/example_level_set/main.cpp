#include <iostream>

#include "Util/IO/io.hpp"


#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "Util/SpiritFunction.hpp"

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
void assemble_poisson(libMesh::EquationSystems & es,
                      const std::string & system_name);

int main(int argc, char ** argv)
{
    BeatIt::printBanner(std::cout);

    using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);
      // allow us to use higher-order approximation.
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh(init.comm());

      int elX = 10;
      int elY = 10;

        MeshTools::Generation::build_square (mesh,
                                       elX, elY,
                                       0., 1.,
                                       0., 1.,
                                       TRI3);

  libMesh::EquationSystems es(mesh);

  libMesh::ExplicitSystem& ls_sys = es.add_system<libMesh::ExplicitSystem>("levet_set");
  ls_sys.add_variable("phi",FIRST,LAGRANGE);

  LinearImplicitSystem&	p_sys = es.add_system<libMesh::LinearImplicitSystem>("s1");
  unsigned int p_var = p_sys.add_variable("u");
//  p_sys.add_variable("a");
	 std::set<libMesh::boundary_id_type> boundary_ids;
	 boundary_ids.insert(1);
	 boundary_ids.insert(2);
	 boundary_ids.insert(3);
	 boundary_ids.insert(4);
	 std::vector<unsigned int> variables(1);
	 variables[0] = p_var;
	     // Create a ZeroFunction to initialize dirichlet_bc
	 libMesh::ZeroFunction< > zf;
	 libMesh::DirichletBoundary dirichlet_bc(boundary_ids, variables, &zf);
	 p_sys.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
  p_sys.attach_assemble_function(assemble_poisson);

	  es.init();
	  es.print_info();

  BeatIt::SpiritFunction init_data;
  std::string func = "0.5*(x-0.25)+0.5*(y-0.5)";
  init_data.read(func);
  ls_sys.project_solution(&init_data);

  es.get_system("s1").solve();

    libMesh::ExodusII_IO exo(mesh);
    exo.write_timestep("ls.exo",es,1,1.0);


    return 0;
}




// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_poisson( libMesh::EquationSystems & es,
                      const std::string & system_name)
{
	using libMesh::UniquePtr;
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "SimpleSystem");


  // Get a constant reference to the mesh object.
  const  libMesh::MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  libMesh::LinearImplicitSystem & system = es.get_system<libMesh::LinearImplicitSystem> ("s1");
//  libMesh::LinearImplicitSystem & ls_sys = es.get_system<libMesh::LinearImplicitSystem> ("level_set");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const libMesh::DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  libMesh::FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.  Introduction Example 4
  // describes some advantages of  UniquePtr's in the context of
  // quadrature rules.
  UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  libMesh::QGauss qrule (dim, libMesh::FIRST);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));

  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  libMesh::QGauss qface(dim-1, libMesh::FIRST);

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
  libMesh::DenseMatrix< libMesh::Number> Ke;
  libMesh::DenseVector< libMesh::Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector< libMesh::dof_id_type> dof_indices;

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
      dof_map.dof_indices (elem, dof_indices);

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
      int num_dofs_u = dof_indices.size();
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test funcions (i) against
          // the trial functions (j).
          for (unsigned int i=0; i<phi.size(); i++)
            for (unsigned int j=0; j<phi.size(); j++)
              {
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
//                Ke(i+num_dofs_u,j+num_dofs_u)  = 1.0;
              }

          // This is the end of the matrix summation loop
          // Now we build the element right-hand-side contribution.
          // This involves a single loop in which we integrate the
          // "forcing function" in the PDE against the test functions.
          {
            const  libMesh::Real x = q_point[qp](0);
            const  libMesh::Real y = q_point[qp](1);
            const  libMesh::Real eps = 1.e-3;


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
            const  libMesh::Real fxy = 1.0;

            for (unsigned int i=0; i<phi.size(); i++)
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
          if (elem->neighbor(side) == libmesh_nullptr)
            {
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
                  const   libMesh::Real penalty = 1.e10;

                  // The boundary value.
                  const  libMesh::Real value = 0.0;

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

