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
#include "Electrophysiology/Monodomain/Monowave.hpp"
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

//      MeshTools::Generation::build_line ( mesh,
//    		  	  	  	  	  	  	  	  	  	  	  	  	  	      numElementsX,
//                                                                      0., maxX );
      MeshTools::Generation::build_square ( mesh,
    		  	  	  	  	  	  	  	  	  	  	  	  	  	      numElementsX, numElementsY,
                                                                      0., maxX, 0.0, maxY, TRI3 );
      std::cout << "Mesh done!" << std::endl;

     //IonicModel* M_ionicModelPtr =  BeatIt::IonicModel::IonicModelFactory::Create("NashPanfilov");
      std::string reaction_mass = data("monodomain/reaction_mass", "mass");
      std::string diffusion_mass = data("monodomain/diffusion_mass", "lumped_mass");
      bool useMidpointMethod = false;
      int step0 = 0;
      int step1 = 1;

      // Constructor
      std::cout << "Create monodomain ..." << std::endl;
      libMesh::EquationSystems es1(mesh);
      BeatIt::Monowave monodomain(es1);


      std::cout << "Setup monodomain ..." << std::endl;
	  monodomain.setup(data, "monodomain");
      // Setup the equation systems
      monodomain.init(0.0);
      std::cout << "Assembling monodomain ..." << std::endl;
      monodomain.assemble_matrices();
      monodomain.restart();

      int save_iter = 0;
      std::cout << "Initializing output monodomain ..." << std::endl;
      monodomain.init_exo_output();
//      return 0;
      save_iter++;

      BeatIt::TimeData datatime;
      datatime.setup(data, "monodomain/");
      datatime.print();
      libMesh::PerfLog perf_log ("Solving");

      for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
      {

		  datatime.advance();
		  monodomain.advance();

		  monodomain.update_pacing(datatime.M_time);
		  monodomain.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, reaction_mass);
		  monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, diffusion_mass);
		  monodomain.update_activation_time(datatime.M_time);


          if( 0 == datatime.M_iter%datatime.M_saveIter )
          {
              std::cout << "* Test Monowave: Time: " << datatime.M_time << std::endl;
             monodomain.save_exo(save_iter++, datatime.M_time);
          }
      }

      monodomain.save_parameters();
      double last_activation_time = monodomain.last_activation_time();
      double potential_norm = monodomain.potential_norm();
      std::cout << std::setprecision(25) << "pot norm = " << potential_norm << std::endl;
      return 0;
}





