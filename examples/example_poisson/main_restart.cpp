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
#include "PoissonSolver/Poisson.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"

#include "Util/CTestUtil.hpp"
#include <iomanip>

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

      MeshTools::Generation::build_cube (mesh,
    		  	  	  	  	  	  	  	  numElementsX, numElementsY, numElementsZ,
                                         0., maxX,
                                         0., maxY,
                                         0., maxZ,
                                         TET4);

      BeatIt::Poisson poisson(mesh);

      poisson.setup(data, "poisson");
      poisson.read_equation_system();
      std::string output_filename = "poisson_restarted.exo";
      poisson.save_exo( "poisson_restarted.exo");

      return 0;

}



