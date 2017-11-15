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
 * \class main
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 11, 2016
 *
 */


// Basic include files needed for the mesh functionality.
#include "Electromechanics/Electromechanics.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"

// Define the Finite Element object.
#include "libmesh/fe.h"
#include "libmesh/dof_map.h"
// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "Util/GenerateFibers.hpp"

#include <iomanip>
#include "Util/Timer.hpp"
#include <fstream>

double exact_solution(double mu, double alpha, double x,  double time);
double exact_derivative(double mu, double alpha, double x, double time);
void set_exact_ic(double mu, double alpha, double time, libMesh::EquationSystems & es);
void eval_l2_err(double& l2_V_err,
		         double& l2_Q_err,
				 double dt, double mu, double alpha, double time,
				 libMesh::EquationSystems & es);
void set_at_exact_solution(double mu, double alpha, double time, libMesh::EquationSystems & es);


int main (int argc, char ** argv)
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
      int numElementsX = data("mesh/elX", 15);
      int numElementsY = data("mesh/elY",   5);
      int numElementsZ = data("mesh/elZ",   4);
      double maxX= data("mesh/maxX", 2.0);
      double maxY = data("mesh/maxY", 0.7);
      double maxZ = data("mesh/maxZ", 0.3);
      double minX= data("mesh/minX", 0.0);

      bool compute_error = data("error", false);

      MeshTools::Generation::build_line ( mesh, numElementsX, minX, maxX );
      // Print information about the mesh to the screen.
      mesh.print_info();


      std::cout << "Mesh done!" << std::endl;

      // Constructor
      std::cout << "Create monodomain ..." << std::endl;
      libMesh::EquationSystems es(mesh);
      BeatIt::Monowave monodomain(es);

      std::cout << "Setup monodomain ..." << std::endl;
      monodomain.setup(data, "monodomain");
      // Setup the equation systems
      monodomain.init(0.0);
      std::cout << "Assembling monodomain ..." << std::endl;
      monodomain.assemble_matrices();



      BeatIt::TimeData datatime;
      datatime.setup(data, "monodomain");
      datatime.print();

//      std::string system_mass = data("monodomain/system_mass", "mass");
//      std::string iion_mass = data("monodomain/iion_mass", "mass");
            std::string system_mass = data("monodomain/diffusion_mass", "mass");
      std::string iion_mass = data("monodomain/reaction_mass", "lumped_mass");
      bool useMidpointMethod = false;
      int step0 = 0;
      int step1 = 1;

      std::cout << "Assembling monodomain ..." << std::endl;
      monodomain.form_system_matrix(datatime.M_dt,false, system_mass);

      std::cout << "Time loop ..." << std::endl;
      BeatIt::Timer timer;
      timer.start();

      double mu = 0.0;
      double alpha = 0.1;
      double l2_V_err = 0.0;
      double l2_Q_err = 0.0;


      std::ofstream errfile;
      if(compute_error)
      {
          mu= data("monodomain/alpha", 0.1);
          alpha= data("monodomain/BistablePiecewiseLinear/a", 0.1);
    	  set_exact_ic( mu, alpha, datatime.M_time, es);
    	  set_at_exact_solution( mu, alpha, datatime.M_time, es);
          std::string err_output = "error_" + std::to_string(numElementsX)
                                 + "_mu_" + std::to_string(mu)
          	  	  	  	  	  	 + "_alpha_" + std::to_string(alpha) + ".txt";
          std::cout << "Writing errors in " << err_output << std::endl;
          errfile.open (err_output);
          errfile << "Time  err_V  err_Q\n";
      }
      int save_iter = 0;
      std::cout << "Initializing output monodomain ..." << std::endl;

		  monodomain.init_exo_output();
		  save_iter++;

      for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
      {
          datatime.advance();

          monodomain.advance();

          //monodomain.update_pacing(datatime.M_time);
          std::cout << "Reaction" << std::endl;
          monodomain.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, iion_mass);
          std::cout << "Diffusion" << std::endl;
          monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
          std::cout << "save" << std::endl;

          if(compute_error) set_at_exact_solution( mu, alpha, datatime.M_time, es);
          else monodomain.update_activation_time(datatime.M_time, 0.9);
          if(compute_error)
          {
        	  eval_l2_err(l2_V_err, l2_Q_err, datatime.M_dt,
        			      mu, alpha, datatime.M_time,
						  es);
              errfile << datatime.M_time << " " << l2_V_err << " "<<  l2_Q_err << std::endl;
          }


		  if( 0 == datatime.M_iter%datatime.M_saveIter )
          {
             std::cout << "* Test EM: Time: " << datatime.M_time << std::endl;
             //monodomain.save_potential(save_iter+1, datatime.M_time);
             save_iter++;
             monodomain.save_exo_timestep(save_iter, datatime.M_time);
          }

      }
      timer.stop();
      timer.print(std::cout);
      if(compute_error)
      {
    	  errfile.close();
      }
      std::cout << "Saving monodomain parameters ..." << std::endl;
      monodomain.save_parameters();
      //monodomain.save_exo(1, datatime.M_time);
      //double last_activation_time = monodomain.last_activation_time();
      //double potential_norm = monodomain.potential_norm();
      //std::cout << std::setprecision(25) << "pot norm = " << potential_norm << std::endl;
      return 0;

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

	typedef libMesh::TransientLinearImplicitSystem     ElectroSystem;
	typedef libMesh::TransientExplicitSystem           IonicModelSystem;
	typedef libMesh::ExplicitSystem                     ParameterSystem;

	ElectroSystem& monodomain_system  =  es.get_system<ElectroSystem>("monowave");//Q
	ElectroSystem& wave_system =  es.get_system<ElectroSystem>("wave");//V

	const libMesh::DofMap& dof_map = monodomain_system.get_dof_map();

	std::vector < libMesh::dof_id_type > dof_indices;

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
	monodomain_system.solution->close();
	monodomain_system.update();
	wave_system.solution->close();
	wave_system.update();
}

void set_at_exact_solution(double mu, double alpha, double time, libMesh::EquationSystems & es)
{
	const libMesh::MeshBase& mesh = es.get_mesh();
	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	typedef libMesh::ExplicitSystem                     ParameterSystem;

	ParameterSystem& activation_times_system  = es.add_system<ParameterSystem>("activation_times");


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
		activation_times_system.solution->set(dof_indices[0], exsol);//V
	}
	activation_times_system.solution->close();
}



void eval_l2_err(double& l2_V_err,
		         double& l2_Q_err,
				 double dt, double mu, double alpha, double time,
				 libMesh::EquationSystems & es)
{
	typedef libMesh::TransientLinearImplicitSystem     ElectroSystem;
	typedef libMesh::TransientExplicitSystem           IonicModelSystem;
	typedef libMesh::ExplicitSystem                     ParameterSystem;

    const libMesh::MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    ElectroSystem& monodomain_system  =  es.get_system<ElectroSystem>("monowave");//Q
	monodomain_system.update();
	ElectroSystem& wave_system =  es.get_system<ElectroSystem>("wave");//V
	wave_system.update();
    const libMesh::DofMap & dof_map_monodomain = monodomain_system.get_dof_map();
    libMesh::FEType fe_type = dof_map_monodomain.variable_type(0);
    libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule(dim, libMesh::FIRST);
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




