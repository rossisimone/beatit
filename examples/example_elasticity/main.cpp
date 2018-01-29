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
#include "Elasticity/MixedElasticity.hpp"
#include "Elasticity/DynamicElasticity.hpp"
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
#include "Util/TimeData.hpp"
#include <libmesh/mesh_function.h>
#include "BoundaryConditions/BCData.hpp"

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
      std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
      GetPot data1(datafile_name);
      datafile_name = commandLine.follow ( "data.beat", 2, "-j", "--jinput" );
      GetPot data2(datafile_name);
      // allow us to use higher-order approximation.
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh1(init.comm());
      libMesh::Mesh mesh2(init.comm());

      libMesh::EquationSystems es1(mesh1);
      libMesh::EquationSystems es2(mesh2);

      std::string meshfile1 = data1("mesh", "NOMESH");
      mesh1.read (&meshfile1[0]);
      std::string meshfile2 = data2("mesh", "NOMESH");
      mesh2.read (&meshfile2[0]);

      BeatIt::TimeData datatime;
      datatime.setup(data1, "elasticity");
      datatime.print();
      int save_iter = 1;

      BeatIt::DynamicElasticity dynamic(es1, "Elasticity");
      dynamic.setup(data1,"elasticity");
      dynamic.init_exo_output("elas_dynamic.exo");
      dynamic.save_exo("elas_dynamic.exo", save_iter, datatime.M_time);

      libMesh::LinearImplicitSystem& system  =  dynamic.M_equationSystems.get_system<libMesh::LinearImplicitSystem>("displacement");
      std::vector<unsigned int> vars(3,0);
      vars[1] = 1;
      vars[2] = 2;
      libMesh::MeshFunction fe_function(es1, *system.solution, system.get_dof_map(), vars );
      fe_function.init();


      BeatIt::Elasticity quasistatic(es2, "Elasticity");
      quasistatic.setup(data2,"elasticity");
      quasistatic.M_bch.get_bc(5)->M_fe_function = &fe_function;
      quasistatic.init_exo_output("elas_static.exo");
      quasistatic.save_exo("elas_static.exo", save_iter, datatime.M_time);
      libMesh::LinearImplicitSystem& system2  =  quasistatic.M_equationSystems.get_system<libMesh::LinearImplicitSystem>(quasistatic.M_myName);

      for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
      {
          std::cout << "time: " << datatime.M_time << std::endl;

          datatime.advance();
          dynamic.setTime(datatime.M_time);
          dynamic.newton(datatime.M_dt);
          dynamic.advance();
          quasistatic.newton();
          if (0 == datatime.M_iter % datatime.M_saveIter)
          {
              std::cout << "* Test dynamic elasticity: Time: " << datatime.M_time << std::endl;
              save_iter++;
              dynamic.save_exo("elas_dynamic.exo", save_iter, datatime.M_time);
              quasistatic.save_exo("elas_static.exo", save_iter, datatime.M_time);
//                elas.save("dg.vtk", save_iter);
//                elas.save("dg.gmv", save_iter);
          }
          system2.solution->zero();

      }
      return 0;
}

