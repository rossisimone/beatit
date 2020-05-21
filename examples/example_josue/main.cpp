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

#include "Electromechanics/Electromechanics.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"


#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"

#include "Util/IO/io.hpp"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

#include <random>
#include "libmesh/vtk_io.h"
#include <iomanip>      // std::setw

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
    mesh.allow_renumbering(false);

    // Read mesh filename from input file
    std::string mesh_name = data("mesh", "NONE");

    // If we want to import the fiber from the meshfile we need to create
    // I/O object we will use to import the data
    // Before reading the fibers we need to create the correct
    // FE spaces to store them
    // This will be done after initializing the Electrophysiology solver
    libMesh::ExodusII_IO importer(mesh);

    // Create mesh:
    // If we passed a specific filename for the mesh let'd read that
    if ("NONE" != mesh_name)
    {
        // Use the Importer to read
        importer.read(mesh_name);
        // If we did not use the importer we would have done
        //mesh.read(mesh_name);

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
    std::cout << "Create Electrophysiology model ..." << std::endl;
    // Create the model
    BeatIt::Electromechanics em(es, "electromechanics");
    // Setup the EP model using the input file
    std::cout << "Calling setup..." << std::endl;
    em.setup(data, "monowave", "elasticity", "activation");
    // Initialize systems
    std::cout << "Calling init ..." << std::endl;
    es.print_info();
    // Set up initial conditions at time
    std::cout << "\nInitializing em: ..." << std::flush;
    em.init(0.0);
    // Now read the fibers if wanted
    bool read_fibers = data("read_fibers", false);
    if(read_fibers)
    {
        // First show the elemental vcariables that can be imported
        auto elemental_variables = importer.get_elem_var_names();
        for (auto && var : elemental_variables) std::cout << var << std::endl;
        em.M_monowave->read_fibers(importer,1);
    }
    // Export simulation parameters
    // This will also export the fiber field
    int save_iter = 0;
    int save_iter_exo = 0;
    std::cout << "Initializing output monodomain ..." << std::endl;
//    exporter.append(true);
    //em.M_monowave->init_exo_output();
    std::cout << "Saving monodomain parameters ..." << std::endl;
    em.M_monowave->save_parameters();

    // Assemble matrices
    std::cout << "Assembling matrices" << std::endl;
    em.M_monowave->assemble_matrices(datatime.M_dt);
    // output file counter
    // Export initial condition at time
//    solver->save_exo_timestep(save_iter, datatime.M_time);
    em.M_monowave->save_potential(save_iter, datatime.M_startTime);
//    solver->save_parameters();

    // Parameters to save the activation times
    // A node is activated if the transmembrane potential > threshold
    double threshold = data("threshold", -50.0);
    // We export the activation times every at_save_iter iterations
    int at_save_iter = data("at_save_iter", 25);
    std::string output_folder = data("output_folder", "Output");

    // Old parameters that define the method, not to be changed
    // TODO: clean this part
    std::string system_mass = "lumped_mass";
    std::string iion_mass = "mass";
    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;
    double emdt = data("time/emdt", 1.0);
    int em_iter = static_cast<int>(emdt / datatime.M_dt);

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
        em.M_monowave->advance();
        // Solve ionic model and evaluate ionic currents
        std::cout << "solve_reaction_step ... " << std::endl;
        em.M_monowave->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        // Solve monodomain model
        std::cout << "solve_diffusion_step ... " << std::endl;
        em.M_monowave->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        // Update the activation times
        std::cout << "update_activation_time ... " << std::endl;
        em.M_monowave->update_activation_time(datatime.M_time, threshold);
        std::cout << "compute_activation ... " << std::endl;
        em.compute_activation(datatime.M_dt);
        // mechanics part
        std::cout << "Done ... " << std::endl;
        if (0 == datatime.M_iter % em_iter)
        {
            std::cout << "Solving mechanics ... " << std::endl;
            em.M_elasticity->setTime(datatime.M_time);
            em.solve_mechanics();
            std::cout << "* Test EM: Time: " << datatime.M_time << std::endl;
        }

        //Export the solution if at the right timestep
        if (0 == datatime.M_iter % datatime.M_saveIter)
        {
            std::cout << "Exporting  ... " << std::endl;
            // update output file counter
            save_iter++;
            // export current solution
            //em.M_monowave->save_potential(save_iter, datatime.M_time);
            std::cout << "Exporting EXODUS  ... " << std::endl;
            std::set<std::string> exported_variables;
            exported_variables.insert("monowave");
            exported_variables.insert("wave");
            exported_variables.insert("istim");
            exported_variables.insert(em.M_elasticity->M_myName);
            exported_variables.insert("I4f");
            exported_variables.insert("activation");

            std::ostringstream ss;
            ss << std::setw(4) << std::setfill('0') << save_iter;
            std::string step_str = ss.str();
            libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/em" + step_str + ".pvtu", es, &exported_variables);

            //exporter.write_equation_systems("em_" + step_str + ".pvtu", es, &exported_variables);
            std::cout << "done " << std::endl;
            //exporter.write_timestep("em.exo", es,  save_iter,  datatime.M_time );
        }
        // export the activation times if at the corresponding timestep
        if (0 == datatime.M_iter % (at_save_iter * datatime.M_saveIter))
        {
            std::cout << "Exporting  Activation Times... " << std::endl;
            // export activation times
            em.M_monowave->save_activation_times(save_iter);
        }

    }

    // // export activation times again
    em.M_monowave->save_activation_times(save_iter);


    // The end
    std::cout << "Good luck with your simulation :P" << std::endl;
    return 0;
}
