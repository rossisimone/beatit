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
#include "libmesh/vector_value.h"

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
#include "libmesh/linear_solver.h"
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

#include <fstream>



double exact_solution(double mu, double alpha, double x,  double time);
double exact_derivative(double mu, double alpha, double x, double time);
void set_exact_ic(double mu, double alpha, double time, libMesh::EquationSystems & es);
void eval_l2_err(double& l2_V_err,
		         double& l2_Q_err,
				 double dt, double mu, double alpha, double time,
				 libMesh::EquationSystems & es);
void set_at_exact_solution(double mu, double alpha, double time, libMesh::EquationSystems & es);

void solve(
        libMesh::EquationSystems & es);
void solveHeun(
        libMesh::EquationSystems & es);


int main(int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace

    using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);


      // Use the MeshTools::Generation mesh generator to create a uniform
      // 3D grid
      // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
      // the number of elements on each side is read from the input file
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh(init.comm());

      GetPot commandLine ( argc, argv );
      std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
      GetPot data(datafile_name);
      // allow us to use higher-order approximation.
      int numElementsX = data("mesh/elX", 1600);
      int numElementsY = data("mesh/elY",   5);
      int numElementsZ = data("mesh/elZ",   4);
      double maxX= data("mesh/maxX", 2.0);
      double maxY = data("mesh/maxY", 0.7);
      double maxZ = data("mesh/maxZ", 0.3);

      MeshTools::Generation::build_line ( mesh, numElementsX, -25.0, 25.0 );

    // Print information about the mesh to the screen.
    mesh.print_info();


    libMesh::EquationSystems systems(mesh);

    double mu = data("mu", 0.2);
    double alpha = data("alpha", 0.1);
    int ms = data("ms", 10);
    double dt = data("dt", 0.003125);
    typedef libMesh::Real Real;
    systems.parameters.set<Real>("mu") = mu;
    systems.parameters.set<Real>("alpha") = alpha;
    systems.parameters.set<Real>("dt") = dt;
    systems.parameters.set<Real>("dummy") = 42.0;
    systems.parameters.set<Real>("nobody") = 0.0;

    typedef libMesh::TransientLinearImplicitSystem Sys;
    Sys& Qsys = systems.add_system<Sys>("SimpleSystem");
    unsigned int Q_var = systems.get_system("SimpleSystem").add_variable("Q",
            libMesh::FIRST);

    Qsys.add_matrix("lumped_mass");
    Qsys.add_matrix("mass");
    Qsys.add_matrix("stiffness");
    Qsys.add_vector("In");
    Qsys.add_vector("I2");
    Qsys.add_vector("Jn");
    Qsys.add_vector("J2");
    Qsys.add_vector("Q2");
    Qsys.add_vector("V2");
    Qsys.add_vector("FnQ");
    Qsys.add_vector("F2Q");
    Qsys.add_vector("MFnQ");
    Qsys.add_vector("MF2Q");
    Qsys.add_vector("dIdVn");
    Qsys.add_vector("dIdV2");
    Qsys.add_vector("KVn");
    Qsys.add_vector("KQn");
    Qsys.add_vector("MQn");
    Qsys.add_vector("KV2");
    Qsys.add_vector("MQ2");

    Sys& Vsys = systems.add_system<Sys>("WaveSystem");
    unsigned int V_var = systems.get_system("WaveSystem").add_variable("V",
            libMesh::FIRST);

    Sys& ATsys = systems.add_system<Sys>("ATSystem");
    unsigned int at_var = systems.get_system("ATSystem").add_variable("at",
            libMesh::FIRST);


//    systems.get_system("SimpleSystem").attach_assemble_function(
//            assemble_poisson);
//    systems.get_system("SimpleSystem").attach_init_function (IC);



    // Initialize the data structures for the equation system.
    systems.init();
    std::cout << "Systems initialized" << std::endl;
    *ATsys.solution = -1;

    systems.get_system("SimpleSystem").time = 0.;

	double time = systems.get_system("SimpleSystem").time;

    double l2_V_err = 0.0;
    double l2_Q_err = 0.0;
   set_exact_ic( mu, alpha, time, systems);
//   Vsys.old_local_solution->print();

   set_at_exact_solution( mu, alpha, time, systems);
    std::string err_output = "error_" + std::to_string(numElementsX)
                           + "_mu_" + std::to_string(mu)
    	  	  	  	  	  	 + "_alpha_" + std::to_string(alpha) + ".txt";
    std::ofstream errfile;
    errfile.open (err_output);
    errfile << "Time  err_V  err_Q\n";

    // Solve the system "Convection-Diffusion".  This will be done by
    // looping over the specified time interval and calling the
    // solve() member at each time step.  This will assemble the
    // system and call the linear solver.


    // Prints information about the system to the screen.
    systems.print_info();



    // TIME LOOP
    std::string exodus_filename = "transient_ex1.e";
    libMesh::ExodusII_IO exo(mesh);

    exo.write_equation_systems (exodus_filename, systems);
    exo.append(true);
//    libMesh::VTKIO(mesh).write_equation_systems ("out.vtu.000", systems);

    Vsys.update();
    Qsys.update();

    unsigned int t_step=0;
    exo.write_timestep (exodus_filename, systems, t_step+1, time);
    t_step++;
    for (;  systems.get_system("SimpleSystem").time < 1.0; t_step++)
    {
    	std::cout << "step: " << t_step << std::endl;
        // Output evey 10 timesteps to file.


        // Incremenet the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        systems.get_system("SimpleSystem").time += dt;

    	time = systems.get_system("SimpleSystem").time;


        solveHeun(systems);

        Vsys.update();
        Qsys.update();

        set_at_exact_solution( mu, alpha, time, systems);
  	    eval_l2_err(l2_V_err, l2_Q_err, dt,
  			      mu, alpha, time,
				  systems);
        errfile << std::fixed << std::setprecision(7) << time << " " << l2_V_err << " "<<  l2_Q_err << std::endl;

        // Define the refinement loop


        auto first_local_index = Vsys.solution->first_local_index();
        auto last_local_index = Vsys.solution->last_local_index();
        for (auto i = first_local_index; i < last_local_index; ++i )
        {
        	double V = (*Vsys.solution)(i);
        	double at = (*ATsys.solution)(i);
        	if( V >= 0.9 && at == -1.0)
        	{
        		ATsys.solution->set(i, time);
        	}
        }

        if(t_step > ms)
        {
        	break;
        }

    	if(t_step%1 == 0)
        exo.write_timestep (exodus_filename, systems, t_step+1, time);


        }

    exo.write_timestep (exodus_filename, systems, t_step+1, time);


    return 0;
}

// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void solve(
        libMesh::EquationSystems & es)
{
    using std::unique_ptr;
    // It is a good idea to make sure we are assembling
    // the proper system.

    // Get a constant reference to the mesh object.
    const libMesh::MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    double D = 1;
    double mu = es.parameters.get<double>("mu");
    double dt = es.parameters.get<double>("dt");
    double time = es.get_system("SimpleSystem").time;
    double alpha = 0.1;

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientLinearImplicitSystem & Qsystem = es.get_system<
            libMesh::TransientLinearImplicitSystem>("SimpleSystem");
    libMesh::TransientLinearImplicitSystem & Vsystem = es.get_system<
            libMesh::TransientLinearImplicitSystem>("WaveSystem");


    Qsystem.get_matrix("lumped_mass").zero();
    Qsystem.get_matrix("mass").zero();
    Qsystem.get_matrix("stiffness").zero();
    Qsystem.get_vector("In").zero();
    Qsystem.get_vector("I2").zero();
    Qsystem.get_vector("Jn").zero();
    Qsystem.get_vector("J2").zero();
    Qsystem.get_vector("Q2").zero();
    Qsystem.get_vector("V2").zero();
    Qsystem.get_vector("FnQ").zero();
    Qsystem.get_vector("F2Q").zero();
    Qsystem.get_vector("dIdVn").zero();
    Qsystem.get_vector("dIdV2").zero();
    Qsystem.get_vector("KVn").zero();
    Qsystem.get_vector("KQn").zero();
    Qsystem.get_vector("MQn").zero();
    Qsystem.get_vector("KV2").zero();
    Qsystem.get_vector("MQ2").zero();
   Qsystem.get_vector("MFnQ").zero();
   Qsystem.get_vector("MF2Q").zero();

    Qsystem.rhs->zero();
    Vsystem.rhs->zero();
    const unsigned int V_Var = Vsystem.variable_number("V");
    const unsigned int Q_Var = Qsystem.variable_number("Q");

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const libMesh::DofMap & dof_map = Qsystem.get_dof_map();

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
    libMesh::DenseMatrix<libMesh::Number> Mate;

    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseMatrix<libMesh::Number> MLe;
    libMesh::DenseVector<libMesh::Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices;

    // setup MATRICES
    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();

//    std::cout << "Assembling matrices" << std::endl;
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
        Me.resize(dof_indices.size(), dof_indices.size());
        MLe.resize(dof_indices.size(), dof_indices.size());
        Fe.resize(dof_indices.size());

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

//            libMesh::Number Vn = 0;
//            libMesh::Number Qn = 0;
//            libMesh::Number gradVn = 0;
//            libMesh::RealGradient nablaVn;
//            for(unsigned int i = 0; i < phi.size(); ++i )
//            {
//                Vn += phi[i][qp] *  (*Vsystem.old_local_solution)(dof_indices[i]);
//                Qn += phi[i][qp] *  (*Qsystem.old_local_solution)(dof_indices[i]);
//                gradVn += dphi[i][qp](0) *  (*Vsystem.old_local_solution)(dof_indices[i]);
//                nablaVn += dphi[i][qp] *  (*Vsystem.old_local_solution)(dof_indices[i]);
//            }

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test funcions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i < phi.size(); i++)
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    // Mass term
//                    Ke(i, j) += (mu + dt*(1+mu)+ dt*dt ) * JxW[qp] * (phi[i][qp] * phi[j][qp]);
//                    // stiffness term
//                    Ke(i, j) += dt * dt * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
//                    Ke(i, j) += (dt+mu) * JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    // stiffness term
//                    Ke(i, j) += dt*dt * D * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);

//					Fe(i) -= JxW[qp] * Vn * (dphi[j][qp] * dphi[i][qp]);
                    //Ke(i, j) += dt*mu / dt * JxW[qp] * (phi[i][qp] * phi[j][qp]);


                	Me(i, j)  += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                	MLe(i,i) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    Ke(i, j)  += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);

                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
//            {
//                const libMesh::Real x = q_point[qp](0);
//                const libMesh::Real y = q_point[qp](1);
//                const libMesh::Real eps = 1.e-3;
//
//                // "fxy" is the forcing function for the Poisson equation.
//                // In this case we set fxy to be a finite difference
//                // Laplacian approximation to the (known) exact solution.
//                //
//                // We will use the second-order accurate FD Laplacian
//                // approximation, which in 2D is
//                //
//                // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
//                //                u(i-1,j) + u(i+1,j) +
//                //                -4*u(i,j))/h^2
//                //
//                // Since the value of the forcing function depends only
//                // on the location of the quadrature point (q_point[qp])
//                // we will compute it here, outside of the i-loop
////                double Iapp = 0;
////                if(time >= 0.03 && time <= 1.03 )
//////				if(time >= 0.0 && time <= 0.004 )
////                {
////                	double r = 0.5;
////                	double x0 = 25.;
////                	if( x >= x0-r && x <= x0+r )
////					{
////                		Iapp = 0.0;
////                        std::cout << time << std::endl;
////					}
////                }
////                for (unsigned int i = 0; i < phi.size(); i++)
////                {
////
////                	double F = ( Vn > 0.1 ) ? 1.0 : 0.0;
////    				F -= Vn;
//////    		        std::cout << F << std::endl;
////
////					Fe(i) += dt*mu / dt * JxW[qp] * Qn * phi[i][qp];
////                	Fe(i) += dt*JxW[qp] * F * phi[i][qp];
////                	Fe(i) += dt*JxW[qp] * Iapp * phi[i][qp];
////					Fe(i) -= dt*JxW[qp] * D * nablaVn * dphi[i][qp];
////
////
////					Fe(i) -= dt * mu * JxW[qp] * Qn * phi[i][qp];
////                }
//            }
        }


        Qsystem.get_matrix("stiffness").add_matrix(Ke, dof_indices);
        Qsystem.get_matrix("mass").add_matrix(Me, dof_indices);
        Qsystem.get_matrix("lumped_mass").add_matrix(MLe, dof_indices);
//        Qsystem.rhs->add_vector(Fe, dof_indices);
    }
//    std::cout << "Done!" << std::endl;

    // All done!
    Qsystem.get_matrix("stiffness").close();
    Qsystem.get_matrix("mass").close();
    Qsystem.get_matrix("lumped_mass").close();
    Qsystem.rhs->zero();
    Qsystem.rhs->close();
    Qsystem.matrix->zero();
    Qsystem.matrix->close();


//    std::cout << "Setup matrix!" << std::endl;

    // MATRIX:
    // ( mu + dt / 2 ) * ML + dt^2 / 4 * K
    Qsystem.matrix->add(( 0.5 * dt +mu ), Qsystem.get_matrix("lumped_mass") );
    Qsystem.matrix->add( 0.25*dt*dt, Qsystem.get_matrix("stiffness" ) );


//    Qsystem.matrix->print();


//    std::cout << "Computing Fn!" << std::endl;


        // STEP 1: Compute Fn
        //                 FnQ = -In - tau * Jn

    	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
    	const libMesh::MeshBase::const_node_iterator end_node =
    			mesh.local_nodes_end();
    	//Vsystem.current_local_solution->print();
    	for (; node != end_node; ++node)
    	{
    		const libMesh::Node * nn = *node;
    		dof_map.dof_indices(nn, dof_indices, 0);
    		double Vn = (*Vsystem.current_local_solution)(dof_indices[0]);
    		double Qn = (*Qsystem.current_local_solution)(dof_indices[0]);
    		double In = ( Vn > alpha ) ? -1.0 : 0.0;
    		In += Vn;
    		Qsystem.get_vector("In").set(dof_indices[0], In);

    		Qsystem.get_vector("dIdVn").set(dof_indices[0], 1.0);
    		Qsystem.get_vector("FnQ").set(dof_indices[0], -In - mu * Qn );

    	}
  	   Qsystem.get_vector("In").close();
  	   Qsystem.get_vector("FnQ").close();
  	   Qsystem.get_vector("dIdVn").close();

//        std::cout << "Computing RHS Q2 !" << std::endl;

    	Qsystem.get_matrix("stiffness").vector_mult_add( Qsystem.get_vector("KVn"),
                                                 *Vsystem.solution );
    	Qsystem.get_matrix("stiffness").vector_mult_add( Qsystem.get_vector("KQn"),
                                                 *Qsystem.solution );
    	Qsystem.get_matrix("lumped_mass").vector_mult_add( Qsystem.get_vector("MQn"),
                                                 *Qsystem.solution );
    	Qsystem.get_matrix("mass").vector_mult_add( Qsystem.get_vector("MFnQ"),
    												Qsystem.get_vector("FnQ") );



 	   Qsystem.get_vector("KVn").close();
 	   Qsystem.get_vector("KQn").close();
 	   Qsystem.get_vector("MQn").close();
 	   Qsystem.get_vector("MFnQ").close();


    	Qsystem.rhs->add( mu - 0.5*dt, Qsystem.get_vector("MQn") );
    	Qsystem.rhs->add( dt, Qsystem.get_vector("MFnQ") );
    	Qsystem.rhs->add( - dt, Qsystem.get_vector("KVn") );
    	Qsystem.rhs->add( -0.25* dt*dt, Qsystem.get_vector("KQn") );

//        std::cout << "Solving for Q2 !" << std::endl;

        std::unique_ptr<libMesh::LinearSolver<libMesh::Number> > linearSolver;
        linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( es.comm() );
       linearSolver->init();

       double tol = 1e-12;
       double max_iter = 2000;

       std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

       rval = linearSolver->solve (*Qsystem.matrix, nullptr,
                                                            Qsystem.get_vector("Q2"),
                                                            *Qsystem.rhs, tol, max_iter);

//       std::cout << "Updating V2 !" << std::endl;

       Qsystem.get_vector("V2") = *Vsystem.solution;
   	   Qsystem.get_vector("V2").add( 0.5*dt, Qsystem.get_vector("Q2") );
   	   Qsystem.get_vector("V2").add( 0.5*dt, *Qsystem.solution );
	   Qsystem.get_vector("V2").close();
//       std::cout << "Updating Vnp1 !" << std::endl;
//	   Vsystem.old_local_solution->print();
//	   *Vsystem.solution = Qsystem.get_vector("V2");

//       std::cout << "Computing F2Q !" << std::endl;

   	node = mesh.local_nodes_begin();

   	for (; node != end_node; ++node)
   	{
   		const libMesh::Node * nn = *node;
   		dof_map.dof_indices(nn, dof_indices, 0);
   		double V2 = Qsystem.get_vector("V2")(dof_indices[0]);
   		double Q2 = Qsystem.get_vector("Q2")(dof_indices[0]);
   		double I2 = ( V2 > alpha ) ? -1.0 : 0.0;
   		I2 += V2;
   		Qsystem.get_vector("I2").set(dof_indices[0], I2);
   		Qsystem.get_vector("dIdV2").set(dof_indices[0], 1.0);
   		Qsystem.get_vector("F2Q").set(dof_indices[0], -I2 - mu * Q2 );
   	}


	   Qsystem.get_vector("I2").close();
	   Qsystem.get_vector("dIdV2").close();
	   Qsystem.get_vector("F2Q").close();


//    std::cout << "Computing Qn+1 RHS !" << std::endl;

	Qsystem.get_matrix("stiffness").vector_mult_add( Qsystem.get_vector("KV2"),
			                                         Qsystem.get_vector("V2"));
	Qsystem.get_matrix("lumped_mass").vector_mult_add( Qsystem.get_vector("MQ2"),
													   Qsystem.get_vector("Q2") );
	Qsystem.get_matrix("mass").vector_mult_add( Qsystem.get_vector("MF2Q"),
												Qsystem.get_vector("F2Q") );


   Qsystem.get_vector("KV2").close();
   Qsystem.get_vector("MQ2").close();
   Qsystem.get_vector("MF2Q").close();




	Qsystem.rhs->zero();
	Qsystem.rhs->close();

	// mu ML

	Qsystem.get_matrix("lumped_mass").vector_mult_add( *Qsystem.rhs,
													   *Qsystem.solution );
	Qsystem.rhs->add( 0.5*dt/mu, Qsystem.get_vector("MFnQ") );
	Qsystem.rhs->add( 0.5*dt/mu, Qsystem.get_vector("MF2Q") );
	Qsystem.rhs->add(-0.5*dt/mu, Qsystem.get_vector("KVn") );
	Qsystem.rhs->add(-0.5*dt/mu, Qsystem.get_vector("KV2") );
	Qsystem.rhs->add(-0.5*dt/mu, Qsystem.get_vector("MQn") );
	Qsystem.rhs->add(-0.5*dt/mu, Qsystem.get_vector("MQ2") );

//    std::cout << "Solving forg Qn+1 !" << std::endl;

    std::unique_ptr<libMesh::LinearSolver<libMesh::Number> > linearSolver2;
    linearSolver2 =  libMesh::LinearSolver<libMesh::Number>::build( es.comm() );
   linearSolver2->init();

   rval = linearSolver2->solve (Qsystem.get_matrix("lumped_mass"), nullptr,
                                                        *Qsystem.solution,
                                                        *Qsystem.rhs, tol, max_iter);

	   //*Qsystem.solution = Qsystem.get_vector("Q2");

   *Vsystem.solution = Qsystem.get_vector("V2");
//   Qsystem.sol = *Vsystem.solution;
//	   Qsystem.get_vector("V2").add( 0.5*dt, Qsystem.get_vector("Q2") );
//	   Qsystem.get_vector("V2").add( 0.5*dt, *Qsystem.solution );
//   Qsystem.get_vector("V2").close();
////       std::cout << "Updating Vnp1 !" << std::endl;
////	   Vsystem.old_local_solution->print();
//   std::cout << "Advancing !" << std::endl;

   *Qsystem.old_local_solution = *Qsystem.solution;
   *Vsystem.old_local_solution = *Vsystem.solution;


}




void solveHeun(
        libMesh::EquationSystems & es)
{
    using std::unique_ptr;
    // It is a good idea to make sure we are assembling
    // the proper system.

    // Get a constant reference to the mesh object.
    const libMesh::MeshBase & mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    double D = 1;
    double mu = es.parameters.get<double>("mu");
    double dt = es.parameters.get<double>("dt");
    double time = es.get_system("SimpleSystem").time;
    double alpha = 0.1;

    // Get a reference to the LinearImplicitSystem we are solving
    libMesh::TransientLinearImplicitSystem & Qsystem = es.get_system<
            libMesh::TransientLinearImplicitSystem>("SimpleSystem");
    libMesh::TransientLinearImplicitSystem & Vsystem = es.get_system<
            libMesh::TransientLinearImplicitSystem>("WaveSystem");


    Qsystem.get_matrix("lumped_mass").zero();
    Qsystem.get_matrix("mass").zero();
    Qsystem.get_matrix("stiffness").zero();
    Qsystem.get_vector("In").zero();
    Qsystem.get_vector("I2").zero();
    Qsystem.get_vector("Jn").zero();
    Qsystem.get_vector("J2").zero();
    Qsystem.get_vector("Q2").zero();
    Qsystem.get_vector("V2").zero();
    Qsystem.get_vector("FnQ").zero();
    Qsystem.get_vector("F2Q").zero();
    Qsystem.get_vector("dIdVn").zero();
    Qsystem.get_vector("dIdV2").zero();
    Qsystem.get_vector("KVn").zero();
    Qsystem.get_vector("KQn").zero();
    Qsystem.get_vector("MQn").zero();
    Qsystem.get_vector("KV2").zero();
    Qsystem.get_vector("MQ2").zero();
   Qsystem.get_vector("MFnQ").zero();
   Qsystem.get_vector("MF2Q").zero();

    Qsystem.rhs->zero();
    Vsystem.rhs->zero();
    const unsigned int V_Var = Vsystem.variable_number("V");
    const unsigned int Q_Var = Qsystem.variable_number("Q");

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const libMesh::DofMap & dof_map = Qsystem.get_dof_map();

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
    libMesh::DenseMatrix<libMesh::Number> Mate;

    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseMatrix<libMesh::Number> MLe;
    libMesh::DenseVector<libMesh::Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices;

    // setup MATRICES
    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();

//    std::cout << "Assembling matrices" << std::endl;
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
        Me.resize(dof_indices.size(), dof_indices.size());
        MLe.resize(dof_indices.size(), dof_indices.size());
        Fe.resize(dof_indices.size());

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

//            libMesh::Number Vn = 0;
//            libMesh::Number Qn = 0;
//            libMesh::Number gradVn = 0;
//            libMesh::RealGradient nablaVn;
//            for(unsigned int i = 0; i < phi.size(); ++i )
//            {
//                Vn += phi[i][qp] *  (*Vsystem.old_local_solution)(dof_indices[i]);
//                Qn += phi[i][qp] *  (*Qsystem.old_local_solution)(dof_indices[i]);
//                gradVn += dphi[i][qp](0) *  (*Vsystem.old_local_solution)(dof_indices[i]);
//                nablaVn += dphi[i][qp] *  (*Vsystem.old_local_solution)(dof_indices[i]);
//            }

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test funcions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i < phi.size(); i++)
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    // Mass term
//                    Ke(i, j) += (mu + dt*(1+mu)+ dt*dt ) * JxW[qp] * (phi[i][qp] * phi[j][qp]);
//                    // stiffness term
//                    Ke(i, j) += dt * dt * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
//                    Ke(i, j) += (dt+mu) * JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    // stiffness term
//                    Ke(i, j) += dt*dt * D * JxW[qp] * (dphi[i][qp] * dphi[j][qp]);

//					Fe(i) -= JxW[qp] * Vn * (dphi[j][qp] * dphi[i][qp]);
                    //Ke(i, j) += dt*mu / dt * JxW[qp] * (phi[i][qp] * phi[j][qp]);


                	Me(i, j)  += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                	MLe(i,i) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                    Ke(i, j)  += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);

                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
//            {
//                const libMesh::Real x = q_point[qp](0);
//                const libMesh::Real y = q_point[qp](1);
//                const libMesh::Real eps = 1.e-3;
//
//                // "fxy" is the forcing function for the Poisson equation.
//                // In this case we set fxy to be a finite difference
//                // Laplacian approximation to the (known) exact solution.
//                //
//                // We will use the second-order accurate FD Laplacian
//                // approximation, which in 2D is
//                //
//                // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
//                //                u(i-1,j) + u(i+1,j) +
//                //                -4*u(i,j))/h^2
//                //
//                // Since the value of the forcing function depends only
//                // on the location of the quadrature point (q_point[qp])
//                // we will compute it here, outside of the i-loop
////                double Iapp = 0;
////                if(time >= 0.03 && time <= 1.03 )
//////				if(time >= 0.0 && time <= 0.004 )
////                {
////                	double r = 0.5;
////                	double x0 = 25.;
////                	if( x >= x0-r && x <= x0+r )
////					{
////                		Iapp = 0.0;
////                        std::cout << time << std::endl;
////					}
////                }
////                for (unsigned int i = 0; i < phi.size(); i++)
////                {
////
////                	double F = ( Vn > 0.1 ) ? 1.0 : 0.0;
////    				F -= Vn;
//////    		        std::cout << F << std::endl;
////
////					Fe(i) += dt*mu / dt * JxW[qp] * Qn * phi[i][qp];
////                	Fe(i) += dt*JxW[qp] * F * phi[i][qp];
////                	Fe(i) += dt*JxW[qp] * Iapp * phi[i][qp];
////					Fe(i) -= dt*JxW[qp] * D * nablaVn * dphi[i][qp];
////
////
////					Fe(i) -= dt * mu * JxW[qp] * Qn * phi[i][qp];
////                }
//            }
        }


        Qsystem.get_matrix("stiffness").add_matrix(Ke, dof_indices);
        Qsystem.get_matrix("mass").add_matrix(Me, dof_indices);
        Qsystem.get_matrix("lumped_mass").add_matrix(MLe, dof_indices);
//        Qsystem.rhs->add_vector(Fe, dof_indices);
    }
//    std::cout << "Done!" << std::endl;

    // All done!
    Qsystem.get_matrix("stiffness").close();
    Qsystem.get_matrix("mass").close();
    Qsystem.get_matrix("lumped_mass").close();
    Qsystem.rhs->zero();
    Qsystem.rhs->close();
    Qsystem.matrix->zero();
    Qsystem.matrix->close();


        // STEP 1: Compute Fn
        //                 FnQ = -In - tau * Jn

    	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
    	const libMesh::MeshBase::const_node_iterator end_node =
    			mesh.local_nodes_end();
    	//Vsystem.current_local_solution->print();
    	for (; node != end_node; ++node)
    	{
    		const libMesh::Node * nn = *node;
    		dof_map.dof_indices(nn, dof_indices, 0);
    		double Vn = (*Vsystem.current_local_solution)(dof_indices[0]);
    		double Qn = (*Qsystem.current_local_solution)(dof_indices[0]);
    		double In = ( Vn > alpha ) ? -1.0 : 0.0;
    		In += Vn;
    		double Jn = 1.0;
    		double V2 = Vn + dt * Qn;

    		Qsystem.get_vector("In").set(dof_indices[0], In);
    		Qsystem.get_vector("dIdVn").set(dof_indices[0], 1.0);
    		Qsystem.get_vector("FnQ").set(dof_indices[0], -In - mu * Jn * Qn );
    		Qsystem.get_vector("In").set(dof_indices[0], In);
    		Qsystem.get_vector("V2").set(dof_indices[0], V2);


    	}
  	   Qsystem.get_vector("In").close();
  	   Qsystem.get_vector("FnQ").close();
  	   Qsystem.get_vector("dIdVn").close();
  	   Qsystem.get_vector("V2").close();

//        std::cout << "Computing RHS Q2 !" << std::endl;

    	Qsystem.get_matrix("stiffness").vector_mult_add( Qsystem.get_vector("KVn"),
                                                 *Vsystem.solution );
    	Qsystem.get_matrix("lumped_mass").vector_mult_add( Qsystem.get_vector("MQn"),
                                                 *Qsystem.solution );
    	Qsystem.get_matrix("mass").vector_mult_add( Qsystem.get_vector("MFnQ"),
    												Qsystem.get_vector("FnQ") );

 	   Qsystem.get_vector("KVn").close();
 	   Qsystem.get_vector("MQn").close();
 	   Qsystem.get_vector("MFnQ").close();


    	Qsystem.rhs->add( -dt/mu, Qsystem.get_vector("MQn") );
    	Qsystem.rhs->add(  dt/mu, Qsystem.get_vector("MFnQ") );
    	Qsystem.rhs->add( -dt/mu, Qsystem.get_vector("KVn") );
    	Qsystem.rhs->add(    1.0, Qsystem.get_vector("MQn") );

//        std::cout << "Solving for Q2 !" << std::endl;

        std::unique_ptr<libMesh::LinearSolver<libMesh::Number> > linearSolver;
        linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( es.comm() );
       linearSolver->init();

       double tol = 1e-12;
       double max_iter = 2000;

       std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

       rval = linearSolver->solve (Qsystem.get_matrix("lumped_mass"), nullptr,
                                                            Qsystem.get_vector("Q2"),
                                                            *Qsystem.rhs, tol, max_iter);

       Vsystem.solution->add(0.5*dt, *Qsystem.solution);
       Vsystem.solution->add(0.5*dt, Qsystem.get_vector("Q2"));

   	node = mesh.local_nodes_begin();

   	for (; node != end_node; ++node)
   	{
   		const libMesh::Node * nn = *node;
   		dof_map.dof_indices(nn, dof_indices, 0);
   		double V2 = Qsystem.get_vector("V2")(dof_indices[0]);
   		double Q2 = Qsystem.get_vector("Q2")(dof_indices[0]);
   		double I2 = ( V2 > alpha ) ? -1.0 : 0.0;
   		I2 += V2;
		double J2 = 1.0;
   		Qsystem.get_vector("I2").set(dof_indices[0], I2);
   		Qsystem.get_vector("dIdV2").set(dof_indices[0], 1.0);
   		Qsystem.get_vector("F2Q").set(dof_indices[0], -I2 - mu * J2 * Q2 );
   	}


	   Qsystem.get_vector("I2").close();
	   Qsystem.get_vector("dIdV2").close();
	   Qsystem.get_vector("F2Q").close();


	    Qsystem.rhs->zero();
	    Qsystem.rhs->close();

//    std::cout << "Computing Qn+1 RHS !" << std::endl;

	Qsystem.get_matrix("stiffness").vector_mult_add( Qsystem.get_vector("KV2"),
			                                         Qsystem.get_vector("V2"));
	Qsystem.get_matrix("lumped_mass").vector_mult_add( Qsystem.get_vector("MQ2"),
													   Qsystem.get_vector("Q2") );
	Qsystem.get_matrix("mass").vector_mult_add( Qsystem.get_vector("MF2Q"),
												Qsystem.get_vector("F2Q") );


   Qsystem.get_vector("KV2").close();
   Qsystem.get_vector("MQ2").close();
   Qsystem.get_vector("MF2Q").close();




	Qsystem.rhs->zero();
	Qsystem.rhs->close();

	Qsystem.rhs->add( -0.5*dt/mu, Qsystem.get_vector("MQn") );
	Qsystem.rhs->add(  0.5*dt/mu, Qsystem.get_vector("MFnQ") );
	Qsystem.rhs->add( -0.5*dt/mu, Qsystem.get_vector("KVn") );
	Qsystem.rhs->add( -0.5*dt/mu, Qsystem.get_vector("MQ2") );
	Qsystem.rhs->add(  0.5*dt/mu, Qsystem.get_vector("MF2Q") );
	Qsystem.rhs->add( -0.5*dt/mu, Qsystem.get_vector("KV2") );
	Qsystem.rhs->add(    1.0, Qsystem.get_vector("MQn") );
	// mu ML

    std::unique_ptr<libMesh::LinearSolver<libMesh::Number> > linearSolver2;
    linearSolver2 =  libMesh::LinearSolver<libMesh::Number>::build( es.comm() );
   linearSolver2->init();


   rval = linearSolver2->solve (Qsystem.get_matrix("lumped_mass"), nullptr,
                                                        *Qsystem.solution,
                                                        *Qsystem.rhs, tol, max_iter);

	   //*Qsystem.solution = Qsystem.get_vector("Q2");

//   Qsystem.sol = *Vsystem.solution;
//	   Qsystem.get_vector("V2").add( 0.5*dt, Qsystem.get_vector("Q2") );
//	   Qsystem.get_vector("V2").add( 0.5*dt, *Qsystem.solution );
//   Qsystem.get_vector("V2").close();
////       std::cout << "Updating Vnp1 !" << std::endl;
////	   Vsystem.old_local_solution->print();
//   std::cout << "Advancing !" << std::endl;

   *Qsystem.old_local_solution = *Qsystem.solution;
   *Vsystem.old_local_solution = *Vsystem.solution;


}




double speed(double mu, double alpha)
{
	return (1-2*alpha) / std::sqrt(mu+(alpha-alpha*alpha)*(mu-1)*(mu-1) );
}

double exact_solution(double mu, double alpha, double x, double time)
{
	double exact_solution;
	double c = speed(mu, alpha);
	double gamma = c * c * mu - 1.0;
	double beta = - c * ( 1 + mu );
	double Delta = beta*beta-4*gamma;
	double lp = ( -beta + std::sqrt(Delta) ) / (2*gamma);
	double lm = ( -beta - std::sqrt(Delta) ) / (2*gamma);

	double z0 = 0;
	double z = x - c *time;

	if(z >= z0)
	{
		exact_solution = alpha * std::exp(lp*z);
	}
	else
	{
		exact_solution = 1 + (alpha-1) * std::exp(lm*z);
	}
	return exact_solution;

}


double exact_derivative(double mu, double alpha, double x,  double time)
{
	double exact_der;
	double c = speed(mu, alpha);
	double gamma = c * c * mu - 1.0;
	double beta = - c * ( 1 + mu );
	double Delta = beta*beta-4*gamma;
	double lp = ( -beta + std::sqrt(Delta) ) / (2*gamma);
	double lm = ( -beta - std::sqrt(Delta) ) / (2*gamma);

	double z0 = 0;
	double z = x - c*time;

	if(z >= z0)
	{
		exact_der = - c * lp * alpha * std::exp(lp*z);
	}
	else
	{
		exact_der = -c * lm * (alpha-1) * std::exp(lm*z);
	}
	return exact_der;

}


void set_exact_ic(double mu, double alpha, double time, libMesh::EquationSystems & es)
{
	const libMesh::MeshBase& mesh = es.get_mesh();
	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	typedef libMesh::TransientLinearImplicitSystem     MonodomainSystem;
	typedef libMesh::TransientExplicitSystem           IonicModelSystem;
	typedef libMesh::ExplicitSystem                     ParameterSystem;

	MonodomainSystem& monodomain_system  =  es.get_system<MonodomainSystem>("SimpleSystem");//Q
	IonicModelSystem& wave_system =  es.get_system<IonicModelSystem>("WaveSystem");//V

	const libMesh::DofMap& dof_map = monodomain_system.get_dof_map();

	std::vector < libMesh::dof_id_type > dof_indices;
	std::cout << "Setting IC" << std::endl;
	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map.dof_indices(nn, dof_indices, 0);
		libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));

		double exsol = exact_solution(mu, alpha, p(0), time);
		double exder = exact_derivative(mu, alpha, p(0), time);
		monodomain_system.solution->set(dof_indices[0], exder);//Q
		wave_system.solution->set(dof_indices[0], exsol);//V
	}
	std::cout << "done" << std::endl;

	monodomain_system.solution->close();
	*monodomain_system.old_local_solution = *monodomain_system.solution;
	monodomain_system.update();
	wave_system.solution->close();
	*wave_system.old_local_solution = *wave_system.solution;
	wave_system.update();
	wave_system.old_local_solution->close();
	std::cout << "IC set" << std::endl;

}

void set_at_exact_solution(double mu, double alpha, double time, libMesh::EquationSystems & es)
{
	const libMesh::MeshBase& mesh = es.get_mesh();
	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	typedef libMesh::ExplicitSystem                     ParameterSystem;

	ParameterSystem& activation_times_system  = es.add_system<ParameterSystem>("ATSystem");


	const libMesh::DofMap& dof_map = activation_times_system.get_dof_map();

	std::vector < libMesh::dof_id_type > dof_indices;

	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map.dof_indices(nn, dof_indices, 0);
		libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));

		double exsol = exact_solution(mu, alpha, p(0), time);
		double exder = exact_derivative(mu, alpha, p(0), time);
//		monodomain_system.solution->set(dof_indices[0], exder);//Q
		activation_times_system.solution->set(dof_indices[0], exder);//V
	}
	activation_times_system.solution->close();
}



void eval_l2_err(double& l2_V_err,
		         double& l2_Q_err,
				 double dt, double mu, double alpha, double time,
				 libMesh::EquationSystems & es)
{
	typedef libMesh::TransientLinearImplicitSystem     MonodomainSystem;
	typedef libMesh::TransientExplicitSystem           IonicModelSystem;
	typedef libMesh::ExplicitSystem                     ParameterSystem;

    const libMesh::MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

	MonodomainSystem& monodomain_system  =  es.get_system<MonodomainSystem>("SimpleSystem");//Q
	monodomain_system.update();
	IonicModelSystem& wave_system =  es.get_system<IonicModelSystem>("WaveSystem");//V
	wave_system.update();
    const libMesh::DofMap & dof_map_monodomain = monodomain_system.get_dof_map();
    libMesh::FEType fe_type = dof_map_monodomain.variable_type(0);
    std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule(dim, libMesh::FOURTH);
    fe->attach_quadrature_rule(&qrule);
    const std::vector<libMesh::Real> & JxW = fe->get_JxW();
    const std::vector<libMesh::Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<libMesh::Real> > & phi = fe->get_phi();

    std::vector<libMesh::dof_id_type> dof_indices;

    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();


    double err_l2_V2 = 0.0;
    double err_l2_Q2 = 0.0;

    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        const unsigned int elem_id = elem->id();
        dof_map_monodomain.dof_indices(elem, dof_indices);
        fe->reinit(elem);

        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

            libMesh::Number Vn = 0;
            libMesh::Number Qn = 0;
            for(unsigned int i = 0; i < phi.size(); ++i )
            {
                Qn += phi[i][qp] *  (*monodomain_system.current_local_solution)(dof_indices[i]);//Q
                Vn += phi[i][qp] *  (*wave_system.current_local_solution)(dof_indices[i]);//V
            }

            double x = q_point[qp](0);
            double exact_V = exact_solution(mu, alpha, x, time);
            double exact_Q = exact_derivative(mu, alpha, x, time);

            err_l2_V2 +=  JxW[qp] * (exact_V - Vn) * (exact_V - Vn);
            err_l2_Q2 +=  JxW[qp] * (exact_Q - Qn) * (exact_Q - Qn);


        }
    }

    double err_l2_V = std::sqrt(err_l2_V2);
    double err_l2_Q = std::sqrt(err_l2_Q2);


    l2_V_err += dt * err_l2_V;
    l2_Q_err += dt * err_l2_Q;
}

