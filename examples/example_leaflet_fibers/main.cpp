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
#include "PoissonSolver/Poisson.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

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

      // We may need XDR support compiled in to read binary .xdr files
      std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");

      // Read the input mesh.
      mesh.read (&meshfile[0]);
      libMesh::EquationSystems es(mesh);
      BeatIt::Poisson poisson1(es);

      std::cout << "Calling setup: ..." << std::flush;
      std::string pois1 = "poisson1";
      poisson1.setup(data, pois1);
      std::cout << " Done!" << std::endl;
      std::cout << "Calling assemble system: ..." << std::flush;
      poisson1.assemble_system();
      std::cout << " Done!" << std::endl;
      std::cout << "Calling solve system: ..." << std::flush;
      poisson1.solve_system();
      std::cout << " Done!" << std::endl;
      std::cout << "Calling gradient: ..." << std::flush;
      poisson1.compute_elemental_solution_gradient();
      std::cout << " Done!" << std::endl;


      BeatIt::Poisson poisson2(es);
      std::cout << "Calling setup: ..." << std::flush;
      std::string pois2 = "poisson2";
      poisson2.setup(data, pois2);
      std::cout << " Done!" << std::endl;
      std::cout << "Calling assemble system: ..." << std::flush;
      poisson2.assemble_system();
      std::cout << " Done!" << std::endl;
      std::cout << "Calling solve system: ..." << std::flush;
      poisson2.solve_system();
      std::cout << " Done!" << std::endl;
      std::cout << "Calling gradient: ..." << std::flush;
      poisson2.compute_elemental_solution_gradient();
      std::cout << " Done!" << std::endl;

      BeatIt::Poisson poisson3(es);

      std::cout << "Calling setup: ..." << std::flush;
      std::string pois3 = "poisson3";
      poisson3.setup(data, pois3);
      auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").solution;
      auto& grad2 = poisson2. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois2+"_gradient").solution;
      auto& grad3 = poisson3. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois3+"_gradient").solution;

      auto first = grad1->first_local_index();
      auto last = grad1->last_local_index();

     auto normalize = [](double& x, double& y, double& z,
                                           double X, double Y, double Z)
    {
        double norm = std::sqrt( x * x + y * y + z * z);
        if(norm >= 1e-12 )
        {
            x /= norm;
            y /= norm;
            z /= norm;
        }
        else
        {
            x = X;
            y = Y;
            z = Z;
        }
    };

      for(int i = first; i < last; )
      {
    	  int j = i;
    	  double v1_x = (*grad1)(i);
    	  double v2_x = (*grad2)(i);
    	  i++;
    	  double v1_y = (*grad1)(i);
    	  double v2_y = (*grad2)(i);
    	  i++;
    	  double v1_z = (*grad1)(i);
    	  double v2_z = (*grad2)(i);
    	  i++;
    	  normalize(v1_x,v1_y,v1_z,1.0,0.0,0.0);
    	  normalize(v2_x,v2_y,v2_z,0.0,1.0,0.0);
          double cross_x = v1_y *  v2_z -  v1_z * v2_y;
          double cross_y = v1_z *  v2_x -  v1_x * v2_z;
          double cross_z = v1_x *  v2_y -  v1_y * v2_x;
          normalize(cross_x, cross_y,cross_z, 0.0, 0.0, 1.0);

          grad1->set(j,v1_x);
    	  grad2->set(j,v2_x);
          grad3->set(j,cross_x);
    	  j++;
    	  grad1->set(j,v1_y);
    	  grad2->set(j,v2_y);
          grad3->set(j,cross_y);
    	  j++;
    	  grad1->set(j,v1_z);
    	  grad2->set(j,v2_z);
          grad3->set(j,cross_z);
      }

      std::cout << "Calling exporter: ..." << std::flush;
        poisson1.save_exo();
        std::cout << " Done!" << std::endl;
        poisson1.write_equation_system();
//     std::cout << "Calling exporter: ..." << std::flush;
//      poisson2.save_exo();
//      std::cout << " Done!" << std::endl;
//      poisson2.write_equation_system();
//      std::cout << "Calling exporter: ..." << std::flush;
//      poisson3.save_exo();
//      std::cout << " Done!" << std::endl;
//      poisson3.write_equation_system();
//      std::cout << "Calling exporter: ..." << std::flush;
//      poisson4.save_exo();
//      std::cout << " Done!" << std::endl;
//      poisson4.write_equation_system();
    	return 0;
}

