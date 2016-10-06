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
      std::string datafile_name = commandLine.follow ( "nash_panfilov.pot", 2, "-i", "--input" );
      GetPot data(datafile_name);
      // allow us to use higher-order approximation.
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh(init.comm());

      // We may need XDR support compiled in to read binary .xdr files
      std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");

      // Read the input mesh.
      mesh.read (&meshfile[0]);

//      int numElementsX = data("mesh/elX", 15);
//      int numElementsY = data("mesh/elY",   5);
//      int numElementsZ = data("mesh/elZ",   4);
//      double maxX= data("mesh/maxX", 2.0);
//      double maxY = data("mesh/maxY", 0.7);
//      double maxZ = data("mesh/maxZ", 0.3);
//
//      MeshTools::Generation::build_cube (mesh,
//    		  	  	  	  	  	  	  	  numElementsX, numElementsY, numElementsZ,
//                                         0., maxX,
//                                         0., maxY,
//                                         0., maxZ,
//                                         TET4);


      // Print information about the mesh to the screen.
      mesh.print_info();


      std::string reaction_mass = data("monodomain/reaction_mass", "mass");
      std::string diffusion_mass = data("monodomain/diffusion_mass", "lumped_mass");
      bool useMidpointMethod = data("monodomain/time/use_midpoint", true);
      int step0 = 0;
      int step1 = 1;

      // Constructor
      BeatIt::Monodomain monodomain(mesh);
      // Setup the equation systems
      monodomain.setup(data, "monodomain");
      monodomain.init(0.0);

//      monodomain.set_fibers( BeatIt::Util::generate_gradient_field( mesh, data, "poisson") );
//      monodomain.set_fibers( BeatIt::Util::generate_gradient_field( mesh, data, "poisson") );
      monodomain.generate_fibers(data, "poisson");
      int init_bc = data("monodomain/init_bc", 1);

      monodomain.set_potential_on_boundary(init_bc, 1.0);

      int save_iter = 0;
      monodomain.init_exo_output();
      monodomain.save_parameters();

      double fiber_norm = monodomain.get_fibers()->l1_norm();
      std::cout << std::setprecision(25) << "fiber norm = " << fiber_norm << std::endl;
      // For ctest
      const double reference_value = 2447.322388569795293733478;
      //We check only up to 12th
    return BeatIt::CTest::check_test(fiber_norm, reference_value, 1e-12);
}





