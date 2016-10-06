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
      libMesh::EquationSystems& es1(mesh);
      BeatIt::Poisson poisson1(es1);

      std::cout << "Calling setup: ..." << std::flush;
      poisson1.setup(data, "poisson1");
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



      libMesh::Mesh mesh2(mesh);
      libMesh::EquationSystems& es2(mesh2);

      BeatIt::Poisson poisson2(es2);

      std::cout << "Calling setup: ..." << std::flush;
      poisson2.setup(data, "poisson2");
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


      libMesh::Mesh mesh4(mesh);
      libMesh::EquationSystems& es4(mesh4);
      BeatIt::Poisson poisson4(es4);

      std::cout << "Calling setup: ..." << std::flush;
      poisson4.setup(data, "poisson4");
      std::cout << " Done!" << std::endl;
      std::cout << "Calling assemble system: ..." << std::flush;
      poisson4.assemble_system();
      std::cout << " Done!" << std::endl;
      std::cout << "Calling solve system: ..." << std::flush;
      poisson4.solve_system();
      std::cout << " Done!" << std::endl;
      std::cout << "Calling gradient: ..." << std::flush;
      poisson4.compute_elemental_solution_gradient();
      std::cout << " Done!" << std::endl;


      libMesh::Mesh mesh3(mesh);
      libMesh::EquationSystems& es3(mesh3);
      BeatIt::Poisson poisson3(es3);

      std::cout << "Calling setup: ..." << std::flush;
      poisson3.setup(data, "poisson3");

      auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>("gradient").solution;
      auto& grad2 = poisson2. M_equationSystems.get_system<libMesh::ExplicitSystem>("gradient").solution;
      auto& grad3 = poisson3. M_equationSystems.get_system<libMesh::ExplicitSystem>("gradient").solution;
      auto& grad4 = poisson4. M_equationSystems.get_system<libMesh::ExplicitSystem>("gradient").solution;

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

    double cx = data ("poisson3/centerline_x", 0.0);
    double cy = data ("poisson3/centerline_y", 0.0);
    double cz = data ("poisson3/centerline_z", 1.0);

      for(int i = first; i < last; )
      {
    	  int j = i;
    	  double v1_x = (*grad1)(i);
    	  double v2_x = (*grad2)(i);
    	  double v4_x = (*grad4)(i);
    	  i++;
    	  double v1_y = (*grad1)(i);
    	  double v2_y = (*grad2)(i);
    	  double v4_y = (*grad4)(i);
    	  i++;
    	  double v1_z = (*grad1)(i);
    	  double v2_z = (*grad2)(i);
    	  double v4_z = (*grad4)(i);
    	  i++;
    	  normalize(v1_x,v1_y,v1_z,1.0,0.0,0.0);
    	  normalize(v2_x,v2_y,v2_z,0.0,1.0,0.0);
    	  normalize(v4_x,v4_y,v4_z,1.0,0.0,0.0);
    	  grad1->set(j,v1_x);
    	  grad2->set(j,v2_x);
    	  grad4->set(j,v4_x);
    	  j++;
    	  grad1->set(j,v1_y);
    	  grad2->set(j,v2_y);
    	  grad4->set(j,v4_y);
    	  j++;
    	  grad1->set(j,v1_z);
    	  grad2->set(j,v2_z);
    	  grad4->set(j,v4_z);

      }



   libMesh::MeshBase::const_element_iterator el =
             mesh.active_local_elements_begin();
     const libMesh::MeshBase::const_element_iterator end_el =
             mesh.active_local_elements_end();
  const libMesh::DofMap & dof_map = poisson1.M_equationSystems.get_system<libMesh::ExplicitSystem>("gradient").get_dof_map();
       std::vector<libMesh::dof_id_type> dof_indices;

     for (; el != end_el; ++el)
     {
    	 const libMesh::Elem * elem = *el;
         dof_map.dof_indices(elem, dof_indices);
    	 const auto blockID = elem->subdomain_id();
    	 if( blockID == 3)
    	 {
    	  double sx = (*grad2)(dof_indices[0]);
    	  double sy = (*grad2)(dof_indices[1]);
    	  double sz = (*grad2)(dof_indices[2]);

    	  double  cdot = cx * sx + cy * sy + cz * sz;
    	  double  xfx = cx - cdot * sx;
          double xfy = cy - cdot * sy;
          double xfz = cz - cdot * sz;
          normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

        double fx = sy * xfz - sz * xfy;
        double fy = sz * xfx - sx * xfz;
        double fz = sx * xfy - sy * xfx;
        normalize(fx, fy, fz, 1.0, 0.0, 0.0);

    	  grad3->set(dof_indices[0],fx);
    	  grad3->set(dof_indices[1],fy);
    	  grad3->set(dof_indices[2],fz);
    	 }
    	 else if( blockID == 7 || blockID == 8 || blockID == 10 || blockID == 11   )
    	 {
    	  double v1_x = (*grad4)(dof_indices[0]);
    	  double v2_x = (*grad2)(dof_indices[0]);
    	  double v1_y = (*grad4)(dof_indices[1]);
    	  double v2_y = (*grad2)(dof_indices[1]);
    	  double v1_z = (*grad4)(dof_indices[2]);
    	  double v2_z = (*grad2)(dof_indices[2]);

    	  double cross_x = v1_y *  v2_z -  v1_z * v2_y;
    	  double cross_y = v1_z *  v2_x -  v1_x * v2_z;
    	  double cross_z = v1_x *  v2_y -  v1_y * v2_x;
    	  normalize(cross_x, cross_y,cross_z,0.0,0.0,1.0);

    	  grad3->set(dof_indices[0],cross_x);
    	  grad3->set(dof_indices[1],cross_y);
    	  grad3->set(dof_indices[2],cross_z);

    	 }
    	 else
    	 {
    	  double v1_x = (*grad1)(dof_indices[0]);
    	  double v2_x = (*grad2)(dof_indices[0]);
    	  double v1_y = (*grad1)(dof_indices[1]);
    	  double v2_y = (*grad2)(dof_indices[1]);
    	  double v1_z = (*grad1)(dof_indices[2]);
    	  double v2_z = (*grad2)(dof_indices[2]);

    	  double cross_x = v1_y *  v2_z -  v1_z * v2_y;
    	  double cross_y = v1_z *  v2_x -  v1_x * v2_z;
    	  double cross_z = v1_x *  v2_y -  v1_y * v2_x;
    	  normalize(cross_x, cross_y,cross_z,0.0,0.0,1.0);

    	  grad3->set(dof_indices[0],cross_x);
    	  grad3->set(dof_indices[1],cross_y);
    	  grad3->set(dof_indices[2],cross_z);
    	 }
     }
	  std::cout << "Calling exporter: ..." << std::flush;
      poisson1.save_exo();
      std::cout << " Done!" << std::endl;
      poisson1.write_equation_system();
     std::cout << "Calling exporter: ..." << std::flush;
      poisson2.save_exo();
      std::cout << " Done!" << std::endl;
      poisson2.write_equation_system();
      std::cout << "Calling exporter: ..." << std::flush;
      poisson3.save_exo();
      std::cout << " Done!" << std::endl;
      poisson3.write_equation_system();
      std::cout << "Calling exporter: ..." << std::flush;
      poisson4.save_exo();
      std::cout << " Done!" << std::endl;
      poisson4.write_equation_system();
    	return 0;
}

