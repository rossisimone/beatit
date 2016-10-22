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

/*
 * main.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */
// Basic include files needed for the mesh functionality.
#include "Elasticity/Elasticity.hpp"

#include "libmesh/linear_implicit_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
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
      std::string datafile_name = commandLine.follow ( "data.pot", 2, "-i", "--input" );
      GetPot data(datafile_name);
      // allow us to use higher-order approximation.
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh(init.comm());

  double nu = data("nu", 0.3);
  double E = data("E", 5e5);
  int elX = data("elX", 4);
  std::string elTypeName = data("elType", "TRI6");
  std::map<std::string, ElemType> orderMap;
  orderMap["TRI3"] = TRI3;
  orderMap["QUAD4"] = QUAD4;
  orderMap["TRI6"] = TRI6;
  orderMap["QUAD9"] = QUAD9;
  auto elType = orderMap.find(elTypeName)->second;
  auto order = FIRST;
  int elY = elX;
  if(elType == TRI6 || elType == QUAD9 ) order = SECOND;
  auto elType2 = TRI3;
  if( elType == QUAD9 ) elType2 = QUAD4;

  MeshTools::Generation::build_square (mesh,
                                       elX, elY,
                                       0., 1.,
                                       0., 0.2,
                                       elType);

  libMesh::EquationSystems es(mesh);

  BeatIt::Elasticity elas(es, "Elasticity");
  elas.setup(data,"elasticity");
  elas.init_exo_output("elas.exo");
  elas.newton();
  elas.save_exo("elas.exo", 1, 1.0);

	  return 0;
}

