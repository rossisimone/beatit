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
#include "Electrophysiology/Bidomain/Bidomain.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"
#include "Electrophysiology/Monodomain/Monowave.hpp"

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
#include "Util/GenerateFibers.hpp"
#include <iomanip>
//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "Util/Timer.hpp"
//#include "libmesh/exodusII_io_helper.h"


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

    // We may need XDR support compiled in to read binary .xdr files
    //std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");

    // Empty mesh
    std::string solver = data("solver", "monodomain");

    // Importer
    libMesh::ExodusII_IO importer(mesh);
    std::string restart_file = data(solver+"/restart/restart_file", "NONE");
    importer.read(restart_file);
    mesh.prepare_for_use();

    libMesh::MeshTools::Modification::scale(mesh, 0.1, 0.1, 0.1);

    libMesh::EquationSystems es(mesh);

    //
    BeatIt::TimeData datatime;

    if(solver == "bidomain")
    {
        BeatIt::Bidomain bidomain(es);
        bidomain.setup(data, "bidomain");

        bidomain.restart(importer, 1, true);

        // DATA TIME
        datatime.setup(data, "bidomain");
        datatime.print();

        std::string system_mass = data("bidomain/diffusion_mass", "mass");
        std::string iion_mass = data("bidomain/reaction_mass", "lumped_mass");

        std::cout << "Assembling matrices ... " << std::flush;

        bidomain.assemble_matrices(datatime.M_dt);
        std::cout << "Done" << std::endl;


        int save_iter = 1;
        std::cout << "Init Output" << std::endl;
        bidomain.init_exo_output();
        bidomain.M_EXOExporter->write_element_data(es);

        bidomain.save_exo_timestep(save_iter++, datatime.M_time);

        bool useMidpointMethod = false;
        int step0 = 0;
        int step1 = 1;

        std::cout << "Time loop starts:" << std::endl;
        for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
        {
            datatime.advance();
            std::cout << "Time:" << datatime.M_time << std::endl;
            bidomain.advance();
            bidomain.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, iion_mass);
            bidomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
            bidomain.update_activation_time(datatime.M_time);
            if( 0 == datatime.M_iter%datatime.M_saveIter )
            {
                save_iter++;
                bidomain.save_potential(save_iter, datatime.M_time);
                bidomain.save_exo_timestep(save_iter, datatime.M_time);
            }
        }
    }
    else
    {
        BeatIt::Monowave monodomain(es);
        monodomain.setup(data, "monodomain");
        monodomain.init(0.0);

        monodomain.restart(importer, 1, true);

        // DATA TIME
        datatime.setup(data, "monodomain");
        datatime.print();

        std::string system_mass = data("monodomain/diffusion_mass", "mass");
        std::string iion_mass = data("monodomain/reaction_mass", "lumped_mass");

        std::cout << "Assembling matrices ... " << std::flush;

        monodomain.assemble_matrices(datatime.M_dt);
        monodomain.form_system_matrix(datatime.M_dt, false);
        std::cout << "Done" << std::endl;


        int save_iter = 1;
        std::cout << "Init Output" << std::endl;
        monodomain.init_exo_output();
        monodomain.M_EXOExporter->write_element_data(es);
        monodomain.save_parameters();
        monodomain.save_exo_timestep(save_iter++, datatime.M_time);

        bool useMidpointMethod = false;
        int step0 = 0;
        int step1 = 1;

        std::cout << "Time loop starts:" << std::endl;
        for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
        {
            datatime.advance();
            std::cout << "Time:" << datatime.M_time << std::endl;
            monodomain.advance();
            monodomain.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, iion_mass);
            monodomain.solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
            monodomain.update_activation_time(datatime.M_time);
            if( 0 == datatime.M_iter%datatime.M_saveIter )
            {
                save_iter++;
                monodomain.save_potential(save_iter, datatime.M_time);
                monodomain.save_exo_timestep(save_iter, datatime.M_time);
            }
        }
        monodomain.save_activation_times(1);

        //monodomain.save_parameters();
    }
//    save_iter++;
//    monodomain.save_exo(save_iter, datatime.M_time);
//    libMesh::ExodusII_IO at_exporter(mesh);
//    std::vector<std::string> output;
//    output.push_back("activation_times");
//    output.push_back("fibersx");
//    output.push_back("fibersy");
//    output.push_back("fibersz");
//    at_exporter.set_output_variables(output);
//    at_exporter.write_equation_systems ("activation_times.exo", es);
//    at_exporter.write_timestep(  "activation_times.exo", es, 1, datatime.M_time);
//    at_exporter.write_element_data(es);
    return 0;
}

