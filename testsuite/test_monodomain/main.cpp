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
#include "Util/IO/io.hpp"
#include <iomanip>

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
      libMesh::PerfLog perf_log ("Monodomain Solver");

      LibMeshInit init (argc, argv, MPI_COMM_WORLD);

      // Create a mesh, with dimension to be overridden later, distributed
      // across the default MPI communicator.
	  perf_log.push("Creating Mesh");

      Mesh mesh(init.comm());

      // Use the MeshTools::Generation mesh generator to create a uniform
      // 3D grid
      // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
      // the number of elements on each side is read from the input file
      GetPot commandLine ( argc, argv );
      std::string datafile_name = commandLine.follow ( "nash_panfilov.beat", 2, "-i", "--input" );
      GetPot data(datafile_name);

      // Reading Meshfile
      std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");
      // Read the input mesh.
      if( meshfile != "Pippo.e" )  mesh.read (&meshfile[0]);
      // or generate the mesh
      else
      {
		  // allow us to use higher-order approximation.
		  int numElementsX = data("mesh/elX", 15);
		  int numElementsY = data("mesh/elY",   5);
		  int numElementsZ = data("mesh/elZ",   4);
		  double maxX= data("mesh/maxX", 2.0);
		  double maxY = data("mesh/maxY", 0.7);
		  double maxZ = data("mesh/maxZ", 0.3);

		  MeshTools::Generation::build_cube (mesh,
											  numElementsX, numElementsY, numElementsZ,
											 0., maxX,
											 0., maxY,
											 0., maxZ,
											 TET4);
      }
   	  perf_log.pop("Creating Mesh");

   	  perf_log.push("Creating Mesh Refinement");
      libMesh:: MeshRefinement mesh_refinement(mesh);

      std::string refine_fractions = data("mesh/refine_fraction", "0.97");
      std::string coarsen_fractions= data("mesh/coarsen_fraction", "0.03");
      std::vector<double> refine_fraction;
      std::vector<double> coarsen_fraction;

      BeatIt::readList(refine_fractions, refine_fraction);
      BeatIt::readList(coarsen_fractions, coarsen_fraction);

      int max_num_mesh_ref = data("mesh/max_num_mesh_ref", 0);
      if(max_num_mesh_ref >0)
      {
          if( refine_fraction.size() != max_num_mesh_ref) max_num_mesh_ref = refine_fraction.size();
      }

      int AMRstep = data("mesh/step", 10);
     std::string error_estimator  =  data("mesh/error_estimator","kelly");
      mesh_refinement.max_h_level()      = max_num_mesh_ref;
      const unsigned int max_r_steps = max_num_mesh_ref;
      bool usingAMR = false;
      if(max_num_mesh_ref > 0) usingAMR = true;
   	  perf_log.pop("Creating Mesh Refinement");
      //IonicModel* M_ionicModelPtr =  BeatIt::IonicModel::IonicModelFactory::Create("NashPanfilov");

      std::string reaction_mass = data("monodomain/reaction_mass", "mass");
      std::string diffusion_mass = data("monodomain/diffusion_mass", "lumped_mass");
      bool useMidpointMethod = data("monodomain/time/use_midpoint", true);
      int step0 = 0;
      int step1 = 1;
      // Create an equation systems object.

//      equation_systems.add_system<libMesh::TransientLinearImplicitSystem>("Test");
//      equation_systems.get_system("Test").add_variable("Pippo",libMesh::FIRST);


   	  perf_log.push("Monodomain setup");

      libMesh::EquationSystems es(mesh);
      BeatIt::Monodomain monodomain(es);
      // Setup the equation systems
      monodomain.setup(data, "monodomain");
      monodomain.init(0.0);

            // Constructor
      std::string poisson_rhs = data("poisson/rhs", "NOPE");
      if( poisson_rhs != "NOPE")
	  {
    	  monodomain.generate_fibers(data, "poisson");
	  }
      std::string V0_boundaries  = data("monodomain/init_boundaries", "NOPE");
      std::vector<int> init_boundaries;
      bool found_init_boundaries = BeatIt::readList(V0_boundaries, init_boundaries);
      if( found_init_boundaries )
      {
    	  for (auto&& b_id : init_boundaries )  monodomain.set_potential_on_boundary(b_id, 1.0);
      }


      BeatIt::TimeData datatime;
      datatime.setup(data, "monodomain");
      datatime.print();

      int save_iter = 0;
      monodomain.init_exo_output();
      save_iter++;
      monodomain.assemble_matrices();
      monodomain.form_system_matrix(datatime.M_dt, useMidpointMethod, diffusion_mass);
   	  perf_log.pop("Monodomain setup");
   	  bool reassemble = usingAMR;

//      libMesh::PerfLog perf_log ("Solving");

      int current_max_r_steps = max_r_steps;
      for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
      {

		  datatime.advance();
		  monodomain.advance();

          if( 0 == (datatime.M_iter-1)%AMRstep && true == usingAMR) current_max_r_steps = max_r_steps;
          else current_max_r_steps = 1;

          for (unsigned int r_step=0; r_step <= current_max_r_steps; r_step++)
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
                    monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, diffusion_mass, reassemble);
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
                    monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, diffusion_mass, reassemble);
                    perf_log.pop("solving diffusion");
                    perf_log.push("update activation times");
                    monodomain.update_activation_time(datatime.M_time);
                    perf_log.pop("update activation times");
              }

              if (r_step != current_max_r_steps && true == usingAMR && 0 == (datatime.M_iter-1)%AMRstep)
              {
                    std::cout << "\n Refinement step: " << r_step << std::endl;
        			perf_log.push("amr");
                    mesh_refinement.refine_fraction()  =  refine_fraction[r_step];
                    mesh_refinement.coarsen_fraction() = coarsen_fraction[r_step];

                    monodomain.amr(mesh_refinement, error_estimator);
                    monodomain.assemble_matrices();
                    monodomain.reinit_linear_solver();
        			perf_log.pop("amr");
              }
          }



          if( 0 == datatime.M_iter%datatime.M_saveIter )
          {
              std::cout << "* Test Monodomain: Time: " << datatime.M_time << std::endl;
              perf_log.push("saving");

              if(usingAMR == true)
			  {
            	  monodomain.save(save_iter);
			  }
              else monodomain.save_exo(save_iter++, datatime.M_time);
              perf_log.pop("saving");
          }
      }

      monodomain.save_parameters();

      double last_activation_time = monodomain.last_activation_time();
      double potential_norm = monodomain.potential_norm();
      std::cout << std::setprecision(25) << "pot norm = " << potential_norm << std::endl;
      // For ctest
      const double reference_value = get_reference_value(monodomain, usingAMR);

  	//We check only up to 12th
      //We check only up to 12th
  	return BeatIt::CTest::check_test(potential_norm, reference_value, 1e-6);

}

double get_reference_value(const BeatIt::Monodomain& monodomain, bool usingAMR)
{
	if(usingAMR) return 116.1200918988604939841025;
	else
	{
		std::string im = monodomain.get_ionic_model_name();
		if("ORd" == im) return 129127.4632047053601127118;
		else if("TP06" == im) return 120756.6560228370362892747;
		else if("Grandi11" == im) return 107207.7281702340114861727;
		else return  56.92915449886483258978842;
	}
}



