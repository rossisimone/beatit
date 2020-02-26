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
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"

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
    Mesh mesh(init.comm());

    // Read mesh filename from input file
    std::string mesh_name = data("mesh", "NONE");

    // Create mesh:
    // If we passed a specific filename for the mesh let'd read that
    if ("NONE" != mesh_name)
    {
        // read mesh
        mesh.read(&mesh_name[0]);

        // We may want to refine this mesh n times
        // if the original mesh is too coarse
        // read the number of times we will refine the mesh, 0 by default
        int num_refs = data("refs", 0);
        std::cout << "Refining the mesh " << num_refs << " times. " << std::endl;
        // Refine the mesh
        MeshRefinement(mesh).uniformly_refine(num_refs);
        // Finalize the mesh to run simulations
        mesh.prepare_for_use();
    }
    // If no mesh file has been specified
    // we run on a cube
    else
    {
        // number of elements in the x,y and z direction
        int nelx = data("nelx", 10);
        int nely = data("nely", 10);
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
        if (nelz == 0)
            elType = TRI3;

        std::cout << "Creating the cube [" << minx << ", " << maxx << "] x ["
                                           << miny << ", " << maxy << "] x ["
                                           << minx << ", " << maxx << "] " << std::endl;
        std::cout << "Using " << nelx << " x " << nely << " x " << nelz << " elements." << std::endl;
        if(TET4 == elType) std::cout << "Element type TET4" << std::endl;
        else if(TRI3 == elType) std::cout << "Element type TRI3" << std::endl;
        else std::cout << "NO ELEMENT TYPE!!!" << std::endl;

        // Create mesh
        MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, minx, maxx, miny, maxy, minz, maxz, elType);
        // Usually we run using cm as dimensions
        // use this part to scale the mesh if it exported in mm (or um)
    }
    double scale = data("scale", 1.0);
    MeshTools::Modification::scale(mesh, scale, scale, scale);

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
    // Export simulation parameters
    solver->save_parameters();
    // Assemble matrices
    std::cout << "Assembling matrices" << std::endl;
    solver->assemble_matrices(datatime.M_dt);
    // output file counter
    int save_iter = 0;
    // Export initial condition at time
    solver->save_potential(save_iter, datatime.M_startTime);

    // Parameters to save the activation times
    // A node is activated if the transmembrane potential > threshold
    double threshold = data("threshold", -10.0);
    // We export the activation times every at_save_iter iterations
    int at_save_iter = data("at_save_iter", 25);

    // Old parameters that define the method, not to be changed
    // TODO: clean this part
    std::string system_mass = data(model + "/diffusion_mass", "mass");
    std::string iion_mass = data(model + "/reaction_mass", "lumped_mass");
    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;


    // Start loop in time
    std::cout << "Time loop starts:" << std::endl;
    // Control the time loop using the TimeData object
    for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
    {
        // We are doing a new iteration
        // let's update first the information in the TimeData object
        datatime.advance();
        // Ouput to screen the current time and the current iteration
        std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter << std::endl;
        // Advance the solution in time: u_n <- u_n+1
        solver->advance();
        // Solve ionic model and evaluate ionic currents
        solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        // Solve monodomain model
        solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        // Update the activation times
        solver->update_activation_time(datatime.M_time, threshold);
        //Export the solution if at the right timestep
        if (0 == datatime.M_iter % datatime.M_saveIter)
        {
            // update output file counter
            save_iter++;
            // export current solution
            solver->save_potential(save_iter, datatime.M_time);
        }
        // export the activation times if at the corresponding timestep
        if (0 == datatime.M_iter % (at_save_iter * datatime.M_saveIter))
        {
            // export activation times
            solver->save_activation_times(save_iter);
        }

    }

    // // export activation times again
    solver->save_activation_times(save_iter);

    // delete solver before ending the simulation
    // avoiding memory leaks
    delete solver;

    // The end
    std::cout << "Good luck with your simulation :P" << std::endl;
    return 0;
}
