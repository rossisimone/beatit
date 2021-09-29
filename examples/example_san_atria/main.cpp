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
#include "PoissonSolver/Poisson.hpp"

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
#include "libmesh/dof_map.h"

#include "libmesh/numeric_vector.h"
#include <libmesh/boundary_mesh.h>
#include <libmesh/serial_mesh.h>

#include "libmesh/analytic_function.h"
#include "Util/CTestUtil.hpp"
#include "Util/GenerateFibers.hpp"
#include <iomanip>
//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "Util/Timer.hpp"
#include "Util/IO/io.hpp"
#include "Util/SetNumericVectorValues.hpp"
#include "Util/Noise.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Setup a mesh like this
/*
 *
 *     ---------------------------         y_max
 *     |                          |
 *     |                          |
 *     |             4            |
 *     |                          |
 *     ---------------------------         y_interface
 *     |          |     |         |
 *     |    1     |  3  |    2    |
 *     |          |     |         |
 *     |          |     |         |
 *     ---------------------------        y_min
 *     x_min                     x_max
 *           x_san_interface
 *                       x_atria_interface
 */
void create_mesh(Mesh& mesh, GetPot& data)
{
    std::cout << "Calling build_square: " << std::endl;
    double x_min = data("x_min", 0.0);
    double x_max = data("x_max", 1.0);

    double y_min = data("y_min", 0.0);
    double y_max = data("y_max", 1.0);

    int elX = data("elX", 20);
    int elY = data("elY", 20);

    MeshTools::Generation::build_square(mesh, elX, elY, x_min, x_max, y_min, y_max, QUAD4);
}


void setup_mesh_subdomains(Mesh& mesh, GetPot& data)
{
    std::cout << "setup_mesh_subdomains: " << std::endl;
    double x_san_interface = data("x_san_interface", 0.45);
    double x_atria_interface = data("x_atria_interface", 0.55);
    double y_interface = data("y_interface", 0.5);

    for(auto & elem : mesh.active_element_ptr_range() )
    {
        Point c = elem->centroid();
        double x = c(0);
        double y = c(1);
        if( y > y_interface ) elem->subdomain_id() = 4;
        else
        {
            if(x < x_san_interface) elem->subdomain_id() = 1;
            else if(x >= x_atria_interface) elem->subdomain_id() = 2;
            else elem->subdomain_id() = 3;
        }

    }
}

void setup_fibrosis_subdomains(EquationSystems& es, GetPot& data)
{
    std::cout << "setup_fibrosis_subdomains: " << std::endl;

    auto& noise_system = es.get_system<libMesh::ExplicitSystem>("Noise");
    auto& dof_map = noise_system.get_dof_map();
    std::vector<libMesh::dof_id_type> dof_indices;
    double threshold = data("threshold", 0.5);

    for(auto & elem : es.get_mesh().active_element_ptr_range() )
    {
        dof_map.dof_indices(elem, dof_indices);
        double noise = (*noise_system.solution)(dof_indices[0]);
        if(noise > threshold && elem->subdomain_id() == 1)
        {
            elem->subdomain_id() = 3;
        }

    }
}


// Begin the main program.
int main(int argc, char ** argv)
{
    // Initialize libMesh and any dependent libaries, like in example 2.
    LibMeshInit init(argc, argv);

    // Use the MeshTools::Generation mesh generator to create a uniform
    // 3D grid
    // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
    // the number of elements on each side is read from the input file
    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data.beat", 2, "-i", "--input");
    GetPot data(datafile_name);

    bool export_data = data("export", false);
    bool compare_data = data("compare", false);

    // Time Data
    BeatIt::TimeData datatime;
    datatime.setup(data);
    datatime.print();

    // Create a mesh with user-defined dimension on the default MPI
    // communicator.
    Mesh mesh(init.comm());

    std::string meshname = data("mesh/input_mesh_name", "NONE");
    bool generate_fibers = true;


    create_mesh(mesh, data);

    std::cout << "Mesh done!" << std::endl;
    mesh.print_info();
    // Create an equation systems object.
    std::cout << "Create equation system ..." << std::endl;
    EquationSystems es(mesh);
    // setup fibrosis
    Noise noise;
    noise.setup(data, "");
    noise.generate_noise_field(es);

    setup_mesh_subdomains(mesh, data);
    setup_fibrosis_subdomains(es, data);

      std::string output_file = data("output","out.e");
      ExodusII_IO exporter_lattice(mesh);
      exporter_lattice.write_equation_systems(output_file,
              es);
      exporter_lattice.write_element_data(es);

    ///////////////////
     ///////////////////
     ///////////////////
     ///////////////////
     ///////////////////
     // Constructor
     std::cout << "Create bidomain with bath..." << std::endl;
     std::string model = data("model", "monowave");
     BeatIt::ElectroSolver* solver = nullptr;

     std::cout << "CREATING VOLUMETRIC SOLVER" << std::endl;
     solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, es);

     std::string section = data("section", "monowave");
     std::cout << "Calling setup..." << std::endl;
     solver->setup(data, model);
     std::cout << "Calling init ..." << std::endl;
     // Set up initial conditions at time
     solver->init(datatime.M_startTime);

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
 //            solver->save_exo_timestep(save_iter, datatime.M_time);
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

    return 0;
}

