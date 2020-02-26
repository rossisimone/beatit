#include "libmesh/getpot.h"
// Functions to initialize the library.
#include "libmesh/libmesh.h"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/patch_recovery_error_estimator.h"
// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"
#include "libmesh/linear_solver.h"

#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"

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
#include "libmesh/perf_log.h"
// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void assemble_amr(
        libMesh::EquationSystems & es,
        const std::string & system_name);

void assemble_matrix(
        libMesh::EquationSystems & es,
        const std::string & system_name);
void assemble_rhs(
        libMesh::EquationSystems & es,
        const std::string & system_name);

// Initial conditions
void IC(libMesh::EquationSystems & es, const std::string & system_name);

int main(int argc, char ** argv)
{
    GetPot data("data.txt");
      int numElementsX = data("elX", 10);
      int numElementsY = data("elY",  10);
      int numElementsZ = data("elZ",  10);
      int num_steps = data("num_steps", 50);
	  double dt = data("dt", 0.01);
	  int sstep= data("sstep", 10);
      double D = data("D", 1e-2);

      int max_h_level = data("levels", 0);
      bool amr = ( max_h_level > 0) ? true : false;
      double ref_fraction = data("ref", 0.7);
      double coar_fraction = data("coar", 0.3);
      double tol = data("tol", 1e-3);

    libMesh::LibMeshInit init (argc, argv, MPI_COMM_WORLD);


    // Create mesh: UNIT CUBE
    libMesh::Mesh mesh(init.comm());
    libMesh::ElemType elType = libMesh::TET4;
    if(numElementsZ < 1) elType = libMesh::TRI3;
     libMesh::MeshTools::Generation::build_cube (mesh,
    		  	  	  	  	  	  	  	  numElementsX, numElementsY, numElementsZ,
                                         0., 1.0,
                                         0., 1.0,
                                         0., 1.0,
                                         elType);
    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create Systems
    libMesh::EquationSystems es(mesh);
    libMesh::TransientLinearImplicitSystem& system = es.add_system<libMesh::TransientLinearImplicitSystem>("bistable");
    unsigned int var = system.add_variable("V", libMesh::FIRST);
    system.attach_assemble_function(assemble_amr);
    system.attach_init_function (IC);
    es.init();

    system.time = 0.;
    es.parameters.set<double>("dt") = dt;
    es.parameters.set<double>("D") = D;

    // TIME LOOP
    std::string exodus_filename = "output.e";
	libMesh::ExodusII_IO exodus(mesh);
    exodus.write_equation_systems (exodus_filename, es);
	exodus.append(true);
	libMesh::GMVIO gmv(mesh);

    unsigned int t_step = 0;
	unsigned int save_step = 0;

	// Define the mesh refinement object that takes care of adaptively
    // refining the mesh.
   libMesh:: MeshRefinement mesh_refinement(mesh);
   // These parameters determine the proportion of elements that will
   // be refined and coarsened. Any element within 30% of the maximum
   // error on any element will be refined, and any element within 30%
   // of the minimum error on any element might be coarsened
   mesh_refinement.refine_fraction()  = ref_fraction;
   mesh_refinement.coarsen_fraction() = coar_fraction;
   // We won't refine any element more than 5 times in total
   mesh_refinement.max_h_level()      =  max_h_level;
   // Refinement parameters
   unsigned int max_r_steps = 2*max_h_level; // Refine the mesh 5 times

   //Create linear solver
   typedef  libMesh::LinearSolver<libMesh::Number> LinearSolver;
   std::unique_ptr<LinearSolver> linear_solver(LinearSolver::build( es.comm() ) );
   linear_solver->init();

  // Starting the solver
   libMesh::PerfLog perf_log("perf_log");
   perf_log.push("time loop");


//   if(!amr)
//    {
// 	     perf_log.push("no amr assembly");
         assemble_matrix(es,"bistable");
//         perf_log.pop("no amr assembly");
//    }

    for (; t_step < num_steps; t_step++)
    {
        // Output evey 10 timesteps to file.
    	std::cout << "Step: " << t_step << std::endl;
     system.update();

		*system.old_local_solution = *system.current_local_solution;

    if(amr)
    {
    	bool stop_amr = false;
        // Define the refinement loop
          for (unsigned int r_step=0; r_step<=max_r_steps; r_step++)
          {
//       	     perf_log.push("amr solve");
              // Assemble & solve the linear system
        	  //std::cout << "Solving on refinement: "  << r_step << std::endl;
			   assemble_rhs(es,"bistable");

//         perf_log.pop("no amr rhs assembly");
//		 perf_log.push("no amr linear solver");

   std::unique_ptr<LinearSolver> linear_solve_amr(LinearSolver::build( es.comm() ) );
   linear_solve_amr->init();
			linear_solve_amr->solve (*system.matrix, nullptr,
															*system.solution,
															*system.rhs, 1e-12, 2000);
			 //system.solve();
//       	     perf_log.pop("amr solve");

              // We need to ensure that the mesh is not refined on the last iteration
              // of this loop, since we do not want to refine the mesh unless we are
              // going to solve the equation system for that refined mesh.
//			 perf_log.push("amr");

              if (r_step != max_r_steps)
                {
                  // Error estimation objects, see Adaptivity Example 2 for details
                  libMesh::ErrorVector error;
                  libMesh::KellyErrorEstimator error_estimator;
//                  libMesh::PatchRecoveryErrorEstimator error_estimator;
                  //libMesh::LaplacianErrorEstimator error_estimator;

                  // Compute the error for each active element
                  error_estimator.estimate_error(system, error);

                  // Output error estimate magnitude
                  libMesh::out << "refinement step:" << r_step
                		       << ", error estimate: l2 norm = "
                               << error.l2_norm()
                               << ", maximum = "
                               << error.maximum()
                               << std::endl;

                  double maxerr = error.maximum();
                  if(maxerr < tol)
				  {
                	  stop_amr = true;
                	  break;
				  }
                  std::cout << "Refining: " << std::endl;
                  // Flag elements to be refined and coarsened
                  mesh_refinement.flag_elements_by_error_fraction (error);

                  // Perform refinement and coarsening
                  mesh_refinement.flag_elements_by_error_fraction (error);
                  mesh_refinement.refine_and_coarsen_elements();
                  // Reinitialize the equation_systems object for the newly refined
                  // mesh. One of the steps in this is project the solution onto the
                  // new mesh
                  es.reinit();

				  assemble_matrix(es,"bistable");

                }

              if(stop_amr) break;
//    		 perf_log.pop("amr");

          }
    }
    else
    {
//		 perf_log.push("no amr rhs  assembly");
         assemble_rhs(es,"bistable");
//         perf_log.pop("no amr rhs assembly");
//		 perf_log.push("no amr linear solver");

         linear_solver->solve (*system.matrix, nullptr,
														*system.solution,
														*system.rhs, 1e-12, 2000);
//		 perf_log.pop("no amr linear solver");
    }

  	  system.time += dt;
  	  if(t_step%sstep == 0)
  	  {
 		  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>"  << std::endl;
 		  std::cout << "Time: " << system.time << std::endl;
		  std::cout << "Exporting solution"   << std::endl;
  		  exodus.write_timestep (exodus_filename, es, ++save_step, system.time);
  		   gmv.write_equation_systems ( "output.gmv."+std::to_string(save_step),  es);
 		  std::cout << ">>>>>>>>>>>>>>>>>>>>>>>"  << std::endl;

  	  }
//  	  max_r_steps = 2;
  	  max_r_steps = 1;
	}

   perf_log.pop("time loop");

    return 0;
}

// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_amr(
        libMesh::EquationSystems & es,
        const std::string & system_name)
{
    using std::unique_ptr;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "bistable");

    double dt = es.parameters.get<double>("dt");
    double D = es.parameters.get<double>("D");

    // Get a constant reference to the mesh object.
    const libMesh::MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientLinearImplicitSystem & system = es.get_system<
            libMesh::TransientLinearImplicitSystem>("bistable");


    const unsigned int V_Var = system.variable_number("V");

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
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule(dim, libMesh::SECOND);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
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

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
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

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

            libMesh::Number V = 0;
            for(unsigned int i = 0; i < phi.size(); ++i )
            {
                V += phi[i][qp] * ( *system.old_local_solution)(dof_indices[i]);
            }

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test funcions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i < phi.size(); i++)
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    // Mass term
                    Ke(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    // stiffness term
                    // Diffusion
                    Ke(i, j) += dt * D * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            {
                const libMesh::Real R = 8 * V * (V - 1.0) * (0.1 - V);

                for (unsigned int i = 0; i < phi.size(); i++)
                {
                    // Mass term
                    Fe(i) += JxW[qp] * V * phi[i][qp];
                    // Reaction term
                    Fe(i) += dt * JxW[qp] * R * phi[i][qp];
                }
            }
        }
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
    if( p(0) <= 0.25 &&  p(1) <= 0.25 ) return 1.0;
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



// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_matrix(
        libMesh::EquationSystems & es,
        const std::string & system_name)
{
    using std::unique_ptr;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "bistable");

    double dt = es.parameters.get<double>("dt");
    double D = es.parameters.get<double>("D");

    // Get a constant reference to the mesh object.
    const libMesh::MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientLinearImplicitSystem & system = es.get_system<
            libMesh::TransientLinearImplicitSystem>("bistable");


    const unsigned int V_Var = system.variable_number("V");

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
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule(dim, libMesh::SECOND);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
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

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
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

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test funcions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i < phi.size(); i++)
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    // Mass term
                    Ke(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    // stiffness term
                    // Diffusion
                    Ke(i, j) += dt * D * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                }

        }
         // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The  SparseMatrix::add_matrix()
        // and  NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
    }
    // All done!
}


// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_rhs(
        libMesh::EquationSystems & es,
        const std::string & system_name)
{


    using std::unique_ptr;
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "bistable");

    double dt = es.parameters.get<double>("dt");

    // Get a constant reference to the mesh object.
    const libMesh::MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientLinearImplicitSystem & system = es.get_system<
            libMesh::TransientLinearImplicitSystem>("bistable");
	 system.rhs->zero();


    const unsigned int V_Var = system.variable_number("V");

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
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule(dim, libMesh::SECOND);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
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
    libMesh::DenseVector<libMesh::Number> Fe;


    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
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
        Fe.resize(dof_indices.size());

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

            libMesh::Number V = 0;
            for(unsigned int i = 0; i < phi.size(); ++i )
            {
                V += phi[i][qp] * ( *system.old_local_solution)(dof_indices[i]);
            }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            {
                const libMesh::Real R = 8 * V * (V - 1.0) * (0.1 - V);

                for (unsigned int i = 0; i < phi.size(); i++)
                {
                    // Mass term
                    Fe(i) += JxW[qp] * V * phi[i][qp];
                    // Reaction term
                    Fe(i) += dt * JxW[qp] * R * phi[i][qp];
                }
            }
        }
         // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The  SparseMatrix::add_matrix()
        // and  NumericVector::add_vector() members do this for us.
        system.rhs->add_vector(Fe, dof_indices);
    }
    // All done!
}

