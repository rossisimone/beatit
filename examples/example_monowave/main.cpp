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

enum class TestCase
{
    NP,                      // Nash Panfilov model
    NP_AMR,           // Nash Panfilov model using AMR
    ORd,                   // ORd model no AMR
};

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
    std::string datafile_name = commandLine.follow("nash_panfilov.pot", 2, "-i", "--input");
    GetPot data(datafile_name);

    /////////////
    // RESTART //
    /////////////
    bool do_restart = data("monodomain/restart/restart", false);
    ExodusII_IO importer(mesh);
    if (do_restart)
    {
        std::string restart_file = data("monodomain/restart/restart_file", "NONE");
        if (restart_file != "NONE")
        {
            importer.read(restart_file);
            mesh.prepare_for_use();
        }
        else
        {
            do_restart = false;
        }

    }

    if (do_restart == false)
    {
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
        if(rotation != 0.0) MeshTools::Modification::rotate(mesh, rotation);
        if( x_translation != 0.0 || y_translation != 0.0 || z_translation != 0.0 ) MeshTools::Modification::translate(mesh, x_translation, y_translation, z_translation);
    }
    std::cout << "Mesh done!" << std::endl;
    // End timer Mesh
    perf_log.pop("mesh");
    ///////////////////////////////////////////////////
    //    __  __ ___ ___ _  _   ___ _  _ ___  ___    //
    //   |  \/  | __/ __| || | | __| \| |   \/ __|   //
    //   | |\/| | _|\__ \ __ | | _|| .` | |) \__ \   //
    //   |_|  |_|___|___/_||_| |___|_|\_|___/|___/   //
    ///////////////////////////////////////////////////



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



    // Constructor
    perf_log.push("create es");
    std::cout << "Create monodomain ..." << std::endl;
    libMesh::EquationSystems es1(mesh);
    perf_log.pop("create es");


    perf_log.push("create monodomain");
    BeatIt::Monowave monodomain(es1);
    perf_log.pop("create monodomain");

    perf_log.push("setup");
    std::cout << "Setup monodomain ..." << std::endl;
    monodomain.setup(data, "monodomain");
    perf_log.pop("setup");

    // Setup the equation systems
    perf_log.push("init");
    monodomain.init(0.0);
    std::cout << "Assembling monodomain ..." << std::endl;
    perf_log.pop("init");

    perf_log.push("assemble matrix");
    monodomain.assemble_matrices();
    perf_log.pop("assemble matrix");


    perf_log.push("restart");
    if (do_restart)
    {
        int restart_step = data("monodomain/restart/step", 2);
        monodomain.restart(importer, restart_step);
    }
    perf_log.pop("rester");

    perf_log.push("init output");
    std::cout << "Initializing output monodomain ..." << std::endl;
    monodomain.init_exo_output();
    perf_log.pop("init output");


    perf_log.push("Data time");
    BeatIt::TimeData datatime;
    datatime.setup(data, "monodomain");
    datatime.print();
    perf_log.pop("Data time");

    perf_log.push("options");
    std::string system_mass = data("monodomain/diffusion_mass", "mass");
    std::string iion_mass = data("monodomain/reaction_mass", "lumped_mass");
    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;

    int save_iter = 0;
    save_iter++;
    perf_log.pop("options");

    perf_log.push("form system matrix");
    monodomain.form_system_matrix(datatime.M_dt, useMidpointMethod, system_mass);
    perf_log.pop("form system matrix");

    perf_log.push("cut");
    bool cut = data("monodomain/cut", false);
    double cut_time = -5.0;
    if (cut)
        cut_time = data("monodomain/c/cut_time", -1.0);
    std::string cut_function;
    if (cut)
        cut_function = data("monodomain/c/function", "NO_FUNCTION");
    std::cout << "cut_time: " << cut_time << ", function: " << cut_function << std::endl;
    perf_log.pop("cut");






    perf_log.push("time loop");

    BeatIt::Timer timer;
    timer.start();

    unsigned int  bID = data("monodomain/bID", 333);
    std::cout << "bID: " << bID << std::endl;
    monodomain.set_potential_on_boundary(bID);


    for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
    {


        perf_log.push("advancing");
        datatime.advance();
        monodomain.advance();
        perf_log.pop("advancing");

        //perf_log.push("pacing");
        //monodomain.update_pacing(datatime.M_time);
        //perf_log.pop("pacing");
        perf_log.push("reaction");
        monodomain.solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        perf_log.pop("reaction");

//          if( 0 == datatime.M_iter%datatime.M_saveIter )
//          {
//              std::cout << "* Test Monowave: Time: " << datatime.M_time << std::endl;
//             monodomain.save_potential(save_iter++, datatime.M_time-0.5*datatime.M_dt);
//          }
        perf_log.push("diffusion");
        monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        perf_log.pop("diffusion");

//          if( 0 == datatime.M_iter%datatime.M_saveIter )
//          {
//              std::cout << "* Test Monowave: Time: " << datatime.M_time << std::endl;
//             monodomain.save_potential(save_iter++, datatime.M_time);
//          }

        perf_log.push("at");
        monodomain.update_activation_time(datatime.M_time);
        perf_log.pop("at");

        perf_log.push("at");

        if (0 == datatime.M_iter % 8)
        {
            std::cout << "Time: " << datatime.M_time << std::endl;
        }
        if (0 == datatime.M_iter % datatime.M_saveIter)
        {
            std::cout << "* Test Monowave: Time: " << datatime.M_time << std::endl;
            save_iter++;
            perf_log.push("output");

            monodomain.save_potential(save_iter, datatime.M_time);
            //monodomain.save(save_iter);
            perf_log.pop("output");
        }
        if (cut && datatime.M_time >= cut_time && datatime.M_time - datatime.M_dt <= cut_time)
        {
            monodomain.cut(datatime.M_time, cut_function);
        }

    }
    timer.stop();
    timer.print(std::cout);
    perf_log.pop("time loop");

    perf_log.push("export parameters");
    monodomain.save_parameters();
    perf_log.pop("export parameters");
    perf_log.push("export solution");
    save_iter++;
    monodomain.save_exo_timestep(save_iter, datatime.M_time);
    monodomain.save_potential(save_iter, datatime.M_time);
    perf_log.pop("export solution");
    monodomain.save_activation_times(1);
//      double last_activation_time = monodomain.last_activation_time();
//      double potential_norm = monodomain.potential_norm();
//      std::cout << std::setprecision(25) << "pot norm = " << potential_norm << std::endl;
    return 0;
}

