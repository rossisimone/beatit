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


#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"

#include "Util/CTestUtil.hpp"
#include "Util/GenerateFibers.hpp"

#include <iomanip>


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
      GetPot commandLine ( argc, argv );
      std::string datafile_name = commandLine.follow ( "em.beat", 2, "-i", "--input" );
      GetPot data(datafile_name);
      // allow us to use higher-order approximation.
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh(init.comm());

      // We may need XDR support compiled in to read binary .xdr files
      //std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");

      // Read the input mesh.
      //mesh.read (&meshfile[0]);

      int numElementsX = data("mesh/elX", 15);
      int numElementsY = data("mesh/elY",  5);
      int numElementsZ = data("mesh/elZ",  4);
      double maxX= data("mesh/maxX", 2.0);
      double maxY = data("mesh/maxY", 0.7);
      double maxZ = data("mesh/maxZ", 0.3);

//      MeshTools::Generation::build_cube (mesh,
//    		  	  	  	  	  	  	  	  numElementsX, numElementsY, numElementsZ,
//                                         0., maxX,
//                                         0., maxY,
//                                         0., maxZ,
//                                         TET4);

      MeshTools::Generation::build_square ( mesh,
                                            numElementsX, numElementsY,
                                            0., maxX,
                                            0., maxY,
                                            TRI3 );

      // Print information about the mesh to the screen.
      mesh.print_info();

      // Constructor
      libMesh::EquationSystems es(mesh);
      BeatIt::Electromechanics em(es, "electromechanics");

      em.setup(data, "monodomain", "elasticity", "activation");

      em.init(0.0);


      int save_iter = 0;
      std::cout << "Initializing output monodomain ..." << std::endl;
      em.M_monowave->init_exo_output();
      std::cout << "Saving monodomain parameters ..." << std::endl;
      em.M_monowave->save_parameters();
//      return 0;

      BeatIt::TimeData datatime;
      datatime.setup(data, "em");
      datatime.print();

      std::string system_mass = data("monodomain/system_mass", "mass");
      std::string iion_mass = data("monodomain/iion_mass", "lumped_mass");
      bool useMidpointMethod = false;
      int step0 = 0;
      int step1 = 1;

      em.M_monowave->assemble_matrices();
      em.M_monowave->form_system_matrix(datatime.M_dt,false, system_mass);

      std::cout << "Assembling monodomain ..." << std::endl;

      std::cout << "Time loop ..." << std::endl;

      int mech_iter = 1.0 / datatime.M_dt;
      for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
      {
          datatime.advance();
          // Electrophysiology

          em.M_monowave->advance();

          em.M_monowave->update_pacing(datatime.M_time);
          em.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, iion_mass);
          em.M_monowave->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);

          em.M_monowave->update_activation_time(datatime.M_time);

          // Activation part
          em.compute_activation(datatime.M_dt);
          // mechanics part
          if( 0 == datatime.M_iter%mech_iter)
          {
              std::cout << "* Test EM: Time: " << datatime.M_time << std::endl;
              em.solve_mechanics();
              std::cout << "Solving mechanics ... " << std::endl;
          }

          if( 0 == datatime.M_iter%datatime.M_saveIter )
          {
             //em.M_monowave->save_potential(save_iter+1, datatime.M_time);
             save_iter++;
             em.save_exo(save_iter, datatime.M_time);
          }

      }

      std::cout << "Saving monodomain parameters ..." << std::endl;
      em.M_monowave->save_parameters();
      em.M_monowave->save_exo(1, datatime.M_time);
      double last_activation_time = em.M_monowave->last_activation_time();
      double potential_norm = em.M_monowave->potential_norm();

      typedef libMesh::TransientExplicitSystem           ActivationSystem;
      std::cout << std::setprecision(25) << "pot norm = " << potential_norm << std::endl;
	  libMesh::LinearImplicitSystem& system  =  em.M_equationSystems.get_system<libMesh::LinearImplicitSystem>(em.M_elasticity->M_myName);
	  libMesh::LinearImplicitSystem& asystem  =  em.M_equationSystems.get_system<libMesh::LinearImplicitSystem>("activation");

	  auto norm = system.solution->linfty_norm ();
	  auto anorm = asystem.solution->linfty_norm ();
    std::cout << "norm is mechanics: " << norm << std::endl;
    std::cout << "activation norm is mechanics: " << anorm << std::endl;

    double reference_potential_norm = 1.180132432165309719351853;
    double reference_norm = 17.90409631663879252982952;
    double reference_anorm = 0.06755865198526839199288929;
    double reference_value = reference_anorm + reference_norm + reference_potential_norm;
    return BeatIt::CTest::check_test(norm+anorm+potential_norm, reference_value, 1e-10);

}





