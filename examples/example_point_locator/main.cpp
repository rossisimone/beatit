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
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/elem.h"
#include <libmesh/node.h>
#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"

#include "Util/CTestUtil.hpp"
#include <iomanip>
//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "Util/Timer.hpp"

#include <libmesh/point_locator_tree.h>

int main(int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace

    using namespace libMesh;
    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    libMesh::PerfLog perf_log("Timing");

   ///////////////////////////
   //  __  __ ___ ___ _  _  //
   // |  \/  | __/ __| || | //
   // | |\/| | _|\__ \ __ | //
   // |_|  |_|___|___/_||_| //
   ///////////////////////////

    // Timer for the mesh:
    perf_log.push("mesh");

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    ParallelMesh mesh(init.comm());

    // Use the MeshTools::Generation mesh generator to create a uniform
    // 3D grid
    // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
    // the number of elements on each side is read from the input file
    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data.beat", 2, "-i", "--input");
    GetPot data(datafile_name);

    /////////////
    // RESTART //
    /////////////
    // allow us to use higher-order approximation.
    int numElementsX = data("mesh/elX", 15);
    int numElementsY = data("mesh/elY", 5);
    int numElementsZ = data("mesh/elZ", 4);
    double maxX = data("mesh/maxX", 2.0);
    double maxY = data("mesh/maxY", 0.7);
    double maxZ = data("mesh/maxZ", 0.3);
    double rotation = data("mesh/rotation", 0.0);
    double x_translation = data("mesh/x_translation", 0.0);
    double y_translation = data("mesh/y_translation", 0.0);
    double z_translation = data("mesh/z_translation", 0.0);

    std::map < std::string, ElemType > orderMap;
    orderMap["TRI3"] = TRI3;
    orderMap["QUAD4"] = QUAD4;
    orderMap["TRI6"] = TRI6;
    orderMap["QUAD9"] = QUAD9;
    std::string mesh_type = data("mesh/type", "TRI3");
    auto elType = orderMap.find(mesh_type)->second;
    if (numElementsZ > 0)
        elType = TET4;
    else
        elType = TRI3;

    MeshTools::Generation::build_cube(mesh, numElementsX, numElementsY, numElementsZ, 0., maxX, 0.0, maxY, 0.0, maxZ, elType);
    std::cout << "Mesh done!" << std::endl;
    // End timer Mesh
    perf_log.pop("mesh");
    ///////////////////////////////////////////////////
    //    __  __ ___ ___ _  _   ___ _  _ ___  ___    //
    //   |  \/  | __/ __| || | | __| \| |   \/ __|   //
    //   | |\/| | _|\__ \ __ | | _|| .` | |) \__ \   //
    //   |_|  |_|___|___/_||_| |___|_|\_|___/|___/   //
    ///////////////////////////////////////////////////


    //   BLOCK IDS

    {
      double z_interface = 0.5*maxZ;

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

      }
    }
    mesh.get_boundary_info().regenerate_id_sets();



    //    ___ ___ ___ ___ _  _ ___ __  __ ___ _  _ _____ ___
    //   | _ \ __| __|_ _| \| | __|  \/  | __| \| |_   _/ __|
    //   |   / _|| _| | || .` | _|| |\/| | _|| .` | | | \__ \
    //   |_|_\___|_| |___|_|\_|___|_|  |_|___|_|\_| |_| |___/

    perf_log.push("mesh refinement");
    int n_refinements = data("mesh/n_ref", 0);
    for (int k = 0; k < n_refinements; ++k)
    {
        std::cout << "refinement: " << k << std::endl;
        MeshRefinement refinement(mesh);
        refinement.uniformly_refine();
    }

    auto nel = mesh.n_active_elem();
    std::cout << "Mesh has " << nel << " active elements" << std::endl;
    perf_log.pop("mesh refinement");

    //    ___ ___ ___ ___ _  _ ___ __  __ ___ _  _ _____ ___   ___ _  _ ___
    //   | _ \ __| __|_ _| \| | __|  \/  | __| \| |_   _/ __| | __| \| |   \
    //   |   / _|| _| | || .` | _|| |\/| | _|| .` | | | \__ \ | _|| .` | |) |
    //   |_|_\___|_| |___|_|\_|___|_|  |_|___|_|\_| |_| |___/ |___|_|\_|___/

    double x = data("x", 0.0);
    double y = data("y", 0.0);
    double z = data("z", 0.0);
    libMesh::Point point(x,y,z);

    std::set<unsigned short> subdomains;
    subdomains.insert(0);

    libMesh::PointLocatorTree locator (mesh);
    locator.enable_out_of_mesh_mode();
    //locator.init(Trees::NODES);
    double tolerance = data("tol", 1e-8);
    auto * node = locator.locate_node(point,&subdomains,tolerance);
    node->print_info();

    libMesh::ExodusII_IO(mesh).write("test.exo");
    return 0;
}

