// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1>Subdomains Example 2 - Subdomain-Restricted Variables</h1>
// \author Benjamin S. Kirk
// \date 2011
//
// This example builds on the fourth example program by showing how
// to restrict solution fields to a subdomain (or union of
// subdomains).


// C++ include files that we need
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
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
#include "Util/Timer.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;


// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries, like in example 2.
  LibMeshInit init (argc, argv);


  // Use the MeshTools::Generation mesh generator to create a uniform
  // 3D grid
  // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
  // the number of elements on each side is read from the input file
  GetPot commandLine(argc, argv);
  std::string datafile_name = commandLine.follow("data.beat", 2, "-i",
          "--input");
  GetPot data(datafile_name);

  BeatIt::TimeData datatime;
  datatime.setup(data, "bidomain");
  datatime.print();

  // Create a mesh with user-defined dimension on the default MPI
  // communicator.
  int dim = 2;
  Mesh mesh (init.comm(), dim);


  // Read FE order from command line
  std::string order = "FIRST";

  // Read FE Family from command line
  std::string family = "LAGRANGE";

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of 8x8 Quad9 elements in 2D, or Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.

  Real halfwidth = dim > 1 ? 1. : 0.;
  Real halfheight = dim > 2 ? 1. : 0.;

  double center_x = 0.0;
  if ((family == "LAGRANGE") && (order == "FIRST"))
  {
      // allow us to use higher-order approximation.
        int elX = data("mesh/elX", 15);
        int elY = data("mesh/elY", 5);
        int elZ = data("mesh/elZ", 4);
        double maxX = data("mesh/maxX", 2.0);
        center_x = 0.5*maxX;
        double maxY = data("mesh/maxY", 0.7);
        double maxZ = data("mesh/maxZ", 0.3);
      // No reason to use high-order geometric elements if we are
      // solving with low-order finite elements.
      MeshTools::Generation::build_cube (mesh,
              elX, elY, elZ,
                                         0., maxX,
                                         0., maxY,
                                         0., maxZ,
                                         TRI3);
  }
  std::cout << "Mesh done!" << std::endl;

  std::cout << "Creating subdomains done!" << std::endl;

  {
    double z_interface = data("mesh/z_interface", 0.0);

    MeshBase::element_iterator       el     = mesh.elements_begin();
    const MeshBase::element_iterator end_el = mesh.elements_end();

    for ( ; el != end_el; ++el)
    {
        Elem * elem = *el;
        const Point cent = elem->centroid();
        // BATH
        if ( cent(2) > z_interface )
        {
              elem->subdomain_id() = 1;
        }
        // TISSUE
        else
        {
            elem->subdomain_id() = 0;
        }

//        for (unsigned int side=0; side<elem->n_sides(); side++)
//        {
//            if ( elem->neighbor_ptr(side) )
//            {
//                auto * nn = elem->side(side)->get_node(0);
//                double x = std::abs( (*nn)(0) );
//                auto * nn2 = elem->side(side)->get_node(1);
//                double x2 = std::abs( (*nn2)(0) );
//                if ( x < 1e-12 && x2 < 1e-12 )
//                {
//                    std::cout << "x: " << x << std::endl;
//                    mesh.get_boundary_info().add_side(elem, side, 10);
//                }
//            }
//        }
    }
  }

//  mesh.get_boundary_info().print_info();
  mesh.get_boundary_info().regenerate_id_sets();

  // Print information about the mesh to the screen.
  mesh.print_info();


  // Create an equation systems object.
  std::cout << "Create equation system ..." << std::endl;
  EquationSystems es (mesh);
  // Constructor
  std::cout << "Create bidomain with bath..." << std::endl;
  BeatIt::BidomainWithBath bidomain(es);
  std::cout << "Calling setup..." << std::endl;
  bidomain.setup(data, "bidomain");
  std::cout << "Calling init ..." << std::endl;
  bidomain.init(0.0);
  std::cout << "Assembling matrices" << std::endl;
  bidomain.assemble_matrices(datatime.M_dt);

  // Declare the system and its variables.
  // Create a system named "Poisson"

  int save_iter = 1;
std::cout << "Init Output" << std::endl;
  bidomain.init_exo_output();
  bidomain.save_exo_timestep(save_iter++, datatime.M_time);


  std::string system_mass = data("bidomain/diffusion_mass", "mass");
  std::string iion_mass = data("bidomain/reaction_mass", "lumped_mass");
  bool useMidpointMethod = false;
  int step0 = 0;
  int step1 = 1;

    std::cout << "Time loop starts:" << std::endl;
  for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
  {
      datatime.advance();
      std::cout << "Time:" << datatime.M_time << std::endl;
      bidomain.advance();
      //std::cout << "Reaction:" << datatime.M_time << std::endl;
      bidomain.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, iion_mass);

      //std::cout << "Diffusion:" << datatime.M_time << std::endl;

        bidomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        //std::cout << "at:" << datatime.M_time << std::endl;
        bidomain.update_activation_time(datatime.M_time);
        //std::cout << "at done:" << datatime.M_time << std::endl;

        if( 0 == datatime.M_iter%datatime.M_saveIter )
        {
            save_iter++;
          bidomain.save_potential(save_iter, datatime.M_time);
          bidomain.save_exo_timestep(save_iter, datatime.M_time);

        }

  }
  return 0;
}



