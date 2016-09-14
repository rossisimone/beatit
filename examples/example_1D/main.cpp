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
#include "Electrophysiology/Monodomain/Monodomain.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"


#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"

#include "Util/CTestUtil.hpp"
#include <iomanip>
//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
enum class TestCase
{
	NP,                      // Nash Panfilov model
	NP_AMR,           // Nash Panfilov model using AMR
	ORd,                   // ORd model no AMR
};

double get_reference_value(const BeatIt::Monodomain& monodomain, bool usingAMR);

int main (int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace

    using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);

      // Create a mesh, with dimension to be overridden later, distributed
      // across the default MPI communicator.
      Mesh mesh(init.comm());

      // Use the MeshTools::Generation mesh generator to create a uniform
      // 3D grid
      // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
      // the number of elements on each side is read from the input file
      GetPot commandLine ( argc, argv );
      std::string datafile_name = commandLine.follow ( "nash_panfilov.pot", 2, "-i", "--input" );
      GetPot data(datafile_name);
      // allow us to use higher-order approximation.
      int numElementsX = data("mesh/elX", 15);
      int numElementsY = data("mesh/elY",   5);
      int numElementsZ = data("mesh/elZ",   4);
      double maxX= data("mesh/maxX", 2.0);
      double maxY = data("mesh/maxY", 0.7);
      double maxZ = data("mesh/maxZ", 0.3);

      MeshTools::Generation::build_line ( mesh,
    		  	  	  	  	  	  	  	  	  	  	  	  	  	      numElementsX,
                                                                      0., maxX );
      std::cout << "Mesh done!" << std::endl;

     //IonicModel* M_ionicModelPtr =  BeatIt::IonicModel::IonicModelFactory::Create("NashPanfilov");
      std::string reaction_mass = data("monodomain/reaction_mass", "mass");
      std::string diffusion_mass = data("monodomain/diffusion_mass", "lumped_mass");
      bool useMidpointMethod = data("monodomain/time/use_midpoint", true);
      int step0 = 0;
      int step1 = 1;
      // Create an equation systems object.
      libMesh::EquationSystems es(mesh);
      es.add_system<libMesh::TransientLinearImplicitSystem>("Test");
      es.get_system("Test").add_variable("Pippo",libMesh::FIRST);
      es.init();
      std::string exodus_filename = "pippo.e";
      libMesh::ExodusII_IO(mesh).write_equation_systems (exodus_filename, es);

      // Constructor
      std::cout << "Create monodomain ..." << std::endl;
      BeatIt::Monodomain monodomain(mesh);
      // Setup the equation systems
      std::cout << "Setup monodomain ..." << std::endl;
      monodomain.setup(data, "monodomain");
      std::cout << "Initialize monodomain ..." << std::endl;
      monodomain.init(0.0);

      BeatIt::TimeData datatime;
      datatime.setup(data, "monodomain/");
      datatime.print();

      int save_iter = 0;
      monodomain.init_exo_output();
//      return 0;
      save_iter++;
      monodomain.assemble_matrices();

      libMesh::PerfLog perf_log ("Solving");

      for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
      {

		  datatime.advance();
		  monodomain.advance();
//
		  {
			  if(useMidpointMethod)
			  {
				  perf_log.push("update pacing");
				  monodomain.update_pacing(datatime.M_time-0.75*datatime.M_dt);
				  perf_log.pop("update pacing");
				  perf_log.push("solving reactions");
				  monodomain.solve_reaction_step(datatime.M_dt/2, datatime.M_time,step0, useMidpointMethod, reaction_mass);
				  perf_log.pop("solving reactions");
				  perf_log.push("solving diffusion");
				  monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, diffusion_mass);
				  perf_log.pop("solving diffusion");
				  perf_log.push("update pacing");
				  monodomain.update_pacing(datatime.M_time-0.25*datatime.M_dt);
				  perf_log.pop("update pacing");
				  perf_log.push("solving reactions");
				  monodomain.solve_reaction_step(datatime.M_dt/2, datatime.M_time, step1, useMidpointMethod,  reaction_mass);
				  perf_log.pop("solving reactions");
				  perf_log.push("update activation times");
				  monodomain.update_activation_time(datatime.M_time);
				  perf_log.pop("update activation times");
			  }
			  else
			  {
				  perf_log.push("update pacing");
				  monodomain.update_pacing(datatime.M_time);
				  perf_log.pop("update pacing");
				  perf_log.push("solving reactions");
				  monodomain.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, reaction_mass);
				  perf_log.pop("solving reactions");
				  perf_log.push("solving diffusion");
				  monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, diffusion_mass);
				  perf_log.pop("solving diffusion");
				  perf_log.push("update activation times");
				  monodomain.update_activation_time(datatime.M_time);
				  perf_log.pop("update activation times");
			  }
		  }


          if( 0 == datatime.M_iter%datatime.M_saveIter )
          {
              std::cout << "* Test Monodomain: Time: " << datatime.M_time << std::endl;
//              monodomain.save(save_iter++);
              perf_log.push("saving");

             monodomain.save_exo(save_iter++, datatime.M_time);
              perf_log.pop("saving");
          }
      }

      monodomain.save_parameters();

      double last_activation_time = monodomain.last_activation_time();
      double potential_norm = monodomain.potential_norm();
      std::cout << std::setprecision(25) << "pot norm = " << potential_norm << std::endl;
      // For ctest
//      const double reference_value = get_reference_value(monodomain, usingAMR);

  	//We check only up to 12th
      //We check only up to 12th
//  	return BeatIt::CTest::check_test(potential_norm, reference_value, 1e-12);
      return 0;
}

double get_reference_value(const BeatIt::Monodomain& monodomain, bool usingAMR)
{
	if(usingAMR) return 104.1266875838878434024082;
	else
	{
		std::string im = monodomain.get_ionic_model_name();
		if("ORd" == im) return 129127.463204705185489729;
		else if("TP06" == im) return 120756.656022836803458631;
		else if("Grandi11" == im) return 107227.470842540860758163;
		else return  56.92915449886483258978842;
	}
}



