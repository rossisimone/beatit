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
    std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");
    mesh.read (&meshfile[0]);

    int scale = 0.1;
    //MeshTools::Modification::scale(mesh, scale, scale, scale);
    perf_log.pop("mesh");

    libMesh::EquationSystems es(mesh);
    //es.init();



    std::string pois1 = "poisson1";
    BeatIt::Poisson poisson1(es, pois1);

    std::cout << "Calling setup: ..." << std::flush;
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

    std::string pois2 = "poisson2";
    BeatIt::Poisson poisson2(es, pois2);

    std::cout << "Calling setup: ..." << std::flush;
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


    std::string pois3 = "poisson3";
    BeatIt::Poisson poisson3(es, pois3);

    std::cout << "Calling setup: ..." << std::flush;
    poisson3.setup(data, pois3);
    std::cout << " Done!" << std::endl;
    std::cout << "Calling assemble system: ..." << std::flush;
    poisson3.assemble_system();
    std::cout << " Done!" << std::endl;
    std::cout << "Calling solve system: ..." << std::flush;
    poisson3.solve_system();
    std::cout << " Done!" << std::endl;
    std::cout << "Calling gradient: ..." << std::flush;
    poisson3.compute_elemental_solution_gradient();
    std::cout << " Done!" << std::endl;


    // Create Fiber fields
    //  long direction
    auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").solution;
    // through thickness
    auto& grad2 = poisson2. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois2+"_gradient").solution;
    // short direction
    auto& grad3 = poisson3. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois3+"_gradient").solution;

    BeatIt::Util::normalize(*grad1, 1.0, 0.0, 0.0);
    BeatIt::Util::normalize(*grad2, 0.0, 1.0, 0.0);
    BeatIt::Util::normalize(*grad3, 0.0, 0.0, 1.0);

    poisson1.save_exo("poisson1.exo");

    //
    BeatIt::Bidomain bidomain(es);
    bidomain.setup(data, "bidomain");
    bidomain.init(0.0);

    // Fibers
    auto& fibers = bidomain.M_equationSystems.get_system<libMesh::ExplicitSystem>("fibers").solution;
    // Sheets
    auto& sheets = bidomain.M_equationSystems.get_system<libMesh::ExplicitSystem>("sheets").solution;
    // XFibers
    auto& xfibers = bidomain.M_equationSystems.get_system<libMesh::ExplicitSystem>("xfibers").solution;

    // Set fibers;
    *fibers = *grad1;
    *sheets = *grad2;
    *xfibers = *grad3;

    poisson1.deleteSystems();
    poisson2.deleteSystems();
    poisson3.deleteSystems();

    // DATA TIME
    BeatIt::TimeData datatime;
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
    bidomain.save_exo(save_iter++, datatime.M_time);

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
        bidomain.save_exo(save_iter, datatime.M_time);
      }
  }

    return 0;
}

