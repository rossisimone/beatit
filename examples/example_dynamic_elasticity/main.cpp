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
#include "Elasticity/MixedDynamicElasticity.hpp"

#include "libmesh/linear_implicit_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
#include <iomanip>
#include "Util/TimeData.hpp"

int main(int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace

    using namespace libMesh;
    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // Use the MeshTools::Generation mesh generator to create a uniform
    // 3D grid
    // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
    // the number of elements on each side is read from the input file
    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data.pot", 2, "-i", "--input");
    GetPot data(datafile_name);
    libMesh::ParallelMesh mesh(init.comm());

    double nu = data("nu", 0.3);
    double E = data("E", 5e5);
    int elX = data("elX", 10);
    int elY = data("elY", 2);
    int elZ = data("elZ", -1);
    std::string elTypeName = data("elType", "TRI6");
    std::map<std::string, ElemType> orderMap;
    orderMap["TRI3"] = TRI3;
    orderMap["QUAD4"] = QUAD4;
    orderMap["TRI6"] = TRI6;
    orderMap["QUAD9"] = QUAD9;
    auto elType = orderMap.find(elTypeName)->second;
    auto order = FIRST;
    if (elType == TRI6 || elType == QUAD9)
        order = SECOND;
    auto elType2 = TRI3;
    if (elType == QUAD9)
        elType2 = QUAD4;

    std::string meshfile = data("mesh", "NOMESH");
    std::cout << "Meshfile: " << meshfile << std::endl;
    if(meshfile == "NOMESH")
    {
        if(elZ > 0)
        {
            MeshTools::Generation::build_cube(mesh, elX, elY, elZ,  0., 10., 0., 10., 0.0, 10.0, TET4);
        }
        else
        {
            MeshTools::Generation::build_square(mesh, elX, elY, 0., 10., 0., 10., elType);
        }
    }
    else mesh.read (&meshfile[0]);

    libMesh::EquationSystems es(mesh);

    std::string formulation = data("elasticity/formulation", "primal");
    std::cout << "Creating " << formulation << std::endl;
    int save_iter = 0;

    BeatIt::Elasticity * elas_ptr;
    std::string out_name;
    if (formulation == "primal")
    {
        std::cout << "Setup primal" << std::endl;
        elas_ptr = new BeatIt::DynamicElasticity(es, "Elasticity");
        out_name = "elas_primal.e";
    }
    else if(formulation == "mixed")
    {
        std::cout << "Setup mixed" << std::endl;
        elas_ptr = new BeatIt::MixedDynamicElasticity(es, "Elasticity");
        out_name = "elas_mixed.e";
    }



    {
        BeatIt::Elasticity& elas = (*elas_ptr);

        elas.setup(data, "elasticity");
        std::cout << "Initializing output " << std::endl;
        elas.init_exo_output(out_name);

        elas.save_exo(out_name, ++save_iter, 0.0);
        std::cout << "Initializing data time " << std::endl;
        BeatIt::TimeData datatime;
        datatime.setup(data, "elasticity");
        datatime.print();
        std::cout << std::setprecision(25) << std::endl;

        for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
        {
            std::cout << "time: " << datatime.M_time << std::endl;

            datatime.advance();
            elas.setTime(datatime.M_time);
            elas.newton(datatime.M_dt);
            elas.advance();
            if (0 == datatime.M_iter % datatime.M_saveIter)
            {
                std::cout << "* Test dynamic elasticity: Time: " << datatime.M_time << std::endl;
                save_iter++;
                elas.save_exo(out_name, save_iter, datatime.M_time);
//                elas.save("dg.vtk", save_iter);
//                elas.save("dg.gmv", save_iter);
            }

        }
    }

    delete elas_ptr;
    return 0;
}

