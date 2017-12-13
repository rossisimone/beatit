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

#include "libmesh/linear_implicit_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
 #include "libmesh/mesh_modification.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
#include <iomanip>


enum Formulation { Primal,  Mixed, Incompressible };
enum CL { Linear,  NH };

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

    int elX = data("elX", 10);
    int elY = data("elY", 2);
    int elZ = data("elZ", 0);
    std::string elTypeName = data("elType", "TRI6");
    //std::map<std::string, ElemType> orderMap;
    //orderMap["TRI6"] = TRI6;
    //auto elType = orderMap.find(elTypeName)->second;
    auto elType = TRI3;
    if(elTypeName == "TRI6") elType = TRI6;
    if(elZ > 0) elType = TET10;
//    auto order = FIRST;
//    if(elType == TRI6 || elType == QUAD9 ) order = SECOND;
//    auto elType2 = TRI3;
//    if( elType == QUAD9 ) elType2 = QUAD4;

    double maxX = data("maxX", 1.);
    double maxY = data("maxY", 1.);
    double maxZ = data("maxZ", 1.);
    std::cout << "Element type chosen:  " << elType << std::endl;
    std::cout << "Max x:  " << maxX << std::endl;
    std::cout << "Max y:  " << maxY << std::endl;

    // Build grid mesh
    MeshTools::Generation::build_cube(mesh,
            elX,
            elY,
            elZ,
              0.0,maxX,
              0.0,maxY,
              0.0,maxZ,
              elType );

    //Map unit square onto cook's membrane
//    int idx = 0;


    bool cm = data("cm_mesh", true);


    if(cm)
    {

      //Map unit square onto cook's membrane
      MeshBase::const_node_iterator nd = mesh.nodes_begin();
      const MeshBase::const_node_iterator end_nd = mesh.nodes_end();
      for (; nd != end_nd; ++nd)
      {
        libMesh::Point s = **nd;
        (**nd)(0) =  4.8 * s(0) + 2.6;
        (**nd)(1) = -2.8 * s(0) * s(1) + 4.4 * s(0) + 4.4 * s(1) + 2.0;
        if(elZ > 0 ) (**nd)(2) = s(2) + 4.5;

      }

      std::string units = data("units", "cm");
      if(units == "mm")
      {
          libMesh::MeshTools::Modification::scale (mesh, 10, 10, 10);
      }
    }

    libMesh::EquationSystems es(mesh);

    Formulation f;
    std::string formulation = data("elasticity/formulation", "NO_FORMULATION");
    std::cout << "Creating " << formulation << std::endl;

    CL cl;
    std::string cl_name = data("elasticity/materials", "NO_MATERIAL");
    std::cout << "Material " << cl_name << std::endl;

    BeatIt::Elasticity* elas = nullptr;
    if(formulation == "primal")
    {
        elas = new BeatIt::Elasticity(es, "Elasticity");
        f = Primal;
        if(cl_name == "linear") cl = Linear;
        else cl = NH;

    }
    else
    {
        elas = new BeatIt::MixedElasticity(es, "Elasticity");
        double nu = 0.0;
        if(cl_name == "linear")
        {
            cl = Linear;
            nu = data("elasticity/materials/linear/nu", -1.0);
        }
        else
        {
            cl = NH;
            nu = data("elasticity/materials/"+cl_name+"/nu", -1.0);
        }
        if(nu < -0.5) return EXIT_FAILURE;
        else if(nu < 0.5) f = Mixed;
        else f = Incompressible;

    }
    std::cout << "Setting up ... " << std::endl;
    elas->setup(data,"elasticity");
    std::cout << "Initializing output ... " << std::endl;
    elas->init_exo_output("solution.exo");
    int save_iter = 1;
    double ramp_dt = data("ramp/dt", 0.1);
    double ramp_end_time = data("ramp/end_time", 0.9);
    double time = 0.0;
    while(time < ramp_end_time)
    {
        std::cout << "Solving ... " << std::endl;
        time+=ramp_dt;
        elas->setTime(time);
        std::cout << "Time: " << time << std::endl;
        elas->newton();
        elas->save_exo("solution.exo", save_iter, time);
        save_iter++;
    }

    delete elas;
    return 0;
}


