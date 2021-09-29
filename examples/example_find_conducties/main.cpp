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
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/linear_implicit_system.h"

#include "Util/IO/io.hpp"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

#include <random>


int main(int argc, char ** argv)
{
    // Read input file
    GetPot data = BeatIt::readInputFile(argc, argv);

    // Use libMesh stuff without writing libMesh:: everytime
    using namespace libMesh;
    // Initialize libMesh
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // Create empty mesh
    ParallelMesh mesh(init.comm());
    mesh.allow_renumbering(false);

    // If we want to import the fiber from the meshfile we need to create
    // I/O object we will use to import the data
    // Before reading the fibers we need to create the correct
    // FE spaces to store them
    // This will be done after initializing the Electrophysiology solver
    libMesh::ExodusII_IO importer(mesh);
    //libMesh::Nemesis_IO nemesis_importer(mesh);

    // Create mesh:
    // we run on a cube
//    {
        // number of elements in the x,y and z direction
        int nelx = data("nelx", 10);
        int nely = data("nely", 0);
        // if nelz = 0 we run in 2D
        int nelz = data("nelz", 10);
        // the cube dimensions are defined as [minx, maxx] x [miny, maxy] x [minz, maxz]
        double maxx = data("maxx", 1.0);
        double maxy = data("maxy", 1.0);
        double maxz = data("maxz", 1.0);

        double minx = data("minx", 0.0);
        double miny = data("miny", 0.0);
        double minz = data("minz", 0.0);

        // Create a tetrahedral mesh
        auto elType = TET4;
        // If we are in 2D, create a triangle mesh
        if (nely == 0)
        {
            elType = EDGE2;
        }
        else if (nelz == 0)
        {
             elType = TRI3;
        }

        std::cout << "Creating the cube [" << minx << ", " << maxx << "] x ["
                                           << miny << ", " << maxy << "] x ["
                                           << minx << ", " << maxx << "] " << std::endl;
        std::cout << "Using " << nelx << " x " << nely << " x " << nelz << " elements." << std::endl;
        if(TET4 == elType) std::cout << "Element type TET4" << std::endl;
        else if(TRI3 == elType) std::cout << "Element type TRI3" << std::endl;
        else if(EDGE2 == elType) std::cout << "Element type EDGE2" << std::endl;
        else std::cout << "NO ELEMENT TYPE!!!" << std::endl;

        // Create mesh
        MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, minx, maxx, miny, maxy, minz, maxz, elType);
        // Usually we run using cm as dimensions
        // use this part to scale the mesh if it exported in mm (or um)
//    }

    // output the details about the mesh
    mesh.print_info();

    // Information about the time, such as
    // current iteration, timestep, final time etc
    // can be conveniently read from the input file
    // and stored in a TimeData abject
    // In the input file specify
    //
    //    [time]
    //        # Timestep
    //        dt = 0.125      # Default: 1.0
    //        # Simulation initial time
    //        init_time = 0.0 # Default: 0.0
    //        # Simulation end time
    //        final_time = 10  # Default: 1.0
    //        # Maximum number of timesteps
    //        max_iter = 200000000   # Default: 99999999
    //        # Export the solution every save_iter iterations
    //        save_iter = 8   # Default: 1
    //    [../]
    //
    // Create the TimeData object
    BeatIt::TimeData datatime;
    // Set it up using the input file
    datatime.setup(data, "");
    // Output on screen the stored variables
    datatime.print();

    // Create libMesh Equations systems
    // This will hold the mesh and create the corresponding
    // finite element spaces
    libMesh::EquationSystems es(mesh);


    // Create the Electrophyiosiology model based on the input file
    //
    //   model = monowave
    //   # The input parameters of the model are written
    //   # in a section with the same name as 'model'
    //   [monowave]
    //       # PARAMETERS HERE
    //   [../]
    //
    // Read the model from the input file
    std::string model = data("model", "monowave");
    std::cout << "Create Electrophysiology model ..." << std::endl;
    // Create the model
    BeatIt::ElectroSolver* solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, es);
    // Setup the EP model using the input file
    std::cout << "Calling setup..." << std::endl;
    solver->setup(data, model);
    // Initialize systems
    std::cout << "Calling init ..." << std::endl;
    es.print_info();
    // Set up initial conditions at time
    solver->init(datatime.M_startTime);
    // Now read the fibers if wanted


    // output file counter
    int save_iter = 0;

    // Old parameters that define the method, not to be changed
    // TODO: clean this part
    std::string system_mass = data(model + "/diffusion_mass", "mass");
    std::string iion_mass = data(model + "/reaction_mass", "lumped_mass");
    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;

    double cv_target = data("target_cv", 60.0);
    double cv_tol_percentage = data("cv_tolerance", 5.0);
    double current_cv = cv_target /2 ;
    double res = std::abs( std::abs(cv_target - current_cv) / cv_target - 1.0 ) * 100;


    libMesh::ExplicitSystem& activation_times_system = es.get_system < libMesh::ExplicitSystem > ("activation_times");
    libMesh::ExplicitSystem& conductivity_system = es.get_system<libMesh::ExplicitSystem>("conductivity");
    libMesh::ExplicitSystem& CV_system = es.get_system < libMesh::ExplicitSystem > ("CV");
    libMesh::TransientLinearImplicitSystem& monodomain_system = es.get_system<libMesh::TransientLinearImplicitSystem>("monowave");
    libMesh::TransientLinearImplicitSystem& wave_system = es.add_system<libMesh::TransientLinearImplicitSystem>("wave");


    int iter = 0;
    int max_iter = data("max_iter",10);
    double threshold = data("threshold", 0.9);

    std::cout << "cv_target: " <<  cv_target << std::endl;
    std::cout << "current_cv: " <<  current_cv << std::endl;
    std::cout << "res: " <<  res <<  ", tol: " <<  cv_tol_percentage << std::endl;
    double factor = 1.0;

    while( res > cv_tol_percentage)
    {
        if(iter > max_iter) break;
        iter++;

        // Set up initial conditions at time
        solver->init(datatime.M_startTime);

        conductivity_system.solution->close();
        (*conductivity_system.solution) *= factor;

        std::cout << "Assembling matrices, ITER: " << iter << std::endl;
        solver->assemble_matrices(datatime.M_dt);

        // Start loop in time
        std::cout << "Time loop starts:" << std::endl;
        // Control the time loop using the TimeData object
        for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
        {
        // We are doing a new iteration
        // let's update first the information in the TimeData object
        datatime.advance();
        // Ouput to screen the current time and the current iteration
        //std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter << std::endl;
        // Advance the solution in time: u_n <- u_n+1
        solver->advance();
        // Solve ionic model and evaluate ionic currents
        solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        // Solve monodomain model
        solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        // Update the activation times
        solver->update_activation_time(datatime.M_time, threshold);
        }



        double min_activation_time = activation_times_system.solution->min();
        double max_activation_time = activation_times_system.solution->max();
        std::cout << "Iter:" << iter << ", min activation time: " << min_activation_time << std::endl;
        std::cout << "Iter:" << iter << ", max activation time: " << max_activation_time << std::endl;
        double sigma = conductivity_system.solution->max();
        double old_sigma = sigma;
        if(min_activation_time > 0)
        {
            double cv_at = (maxx - minx) / (max_activation_time - min_activation_time) * 1000;
            std::cout << "Iter:" << iter << ", average cv: " << cv_at << std::endl;
            solver->evaluate_conduction_velocity();
            double cv_min = CV_system.solution->min();
            double cv_max = CV_system.solution->max();
            double cv_mid = CV_system.solution->el( static_cast<int>( CV_system.solution->size()/6) );
            current_cv = cv_at;

            sigma *=  cv_target * cv_target / current_cv / current_cv;
            factor *=  cv_target * cv_target / current_cv / current_cv;

            std::cout << "Iter: " << iter << ", cv_min: " << cv_min << ", cv_max: " << cv_max << ", cv_mid: " << cv_mid << ", factor: " << factor << std::endl;
            res = std::abs( std::abs(cv_target - current_cv) / cv_target ) * 100;
        }
        else
        {
            current_cv = 0.0;
            res = 100.0;
            factor *= 1.5;
            sigma *= 1.5;
        }
        std::cout << "Iter: " << iter << ", CV: " << current_cv << ". old_sigma: " << old_sigma << ", sigma: " << sigma << std::endl;
        // reset
        activation_times_system.solution->close();
        (*activation_times_system.solution) = -1.0;
        activation_times_system.update();
        std::cout << "Resetting: " << std::endl;

        monodomain_system.solution->zero();
        ( *monodomain_system.old_local_solution ) = ( *monodomain_system.solution );
        monodomain_system.update();

        wave_system.solution->zero();
        ( *wave_system.old_local_solution ) = ( *wave_system.solution );
        wave_system.update();

        datatime.reset();
        std::cout << "res: " <<  res <<  ", tol: " <<  cv_tol_percentage << std::endl;
    }

    // delete solver before ending the simulation
    // avoiding memory leaks
    delete solver;

    // The end
    std::cout << "Good luck with your simulation :P" << std::endl;
    return 0;
}
