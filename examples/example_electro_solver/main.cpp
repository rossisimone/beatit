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
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/elem.h"

#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"

#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"

#include "Util/CTestUtil.hpp"
#include <iomanip>
//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "Util/Timer.hpp"
#include "Util/IO/io.hpp"
// Bring in everything from the libMesh namespace
using namespace libMesh;

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

    BeatIt::TimeData datatime;
    datatime.setup(data);
    datatime.print();

    // Create a mesh with user-defined dimension on the default MPI
    // communicator.
    ParallelMesh mesh(init.comm());

    // We may need XDR support compiled in to read binary .xdr files
    std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");
    int n_refinements = data("mesh/n_ref", 0);
    std::cout << "n_refs: " << n_refinements << std::endl;
    BeatIt::serial_mesh_partition(init.comm(), meshfile, &mesh, n_refinements);


    double xscale = data("mesh/x_scale", 1);
    double yscale = data("mesh/y_scale", 1);
    double zscale = data("mesh/z_scale", 1);
    if(xscale != 1 || yscale != 1 || zscale != 1) libMesh::MeshTools::Modification::scale(mesh, xscale, yscale, zscale);

//    libMesh::ExodusII_IO importer(mesh);
//
//    double center_x = 0.0;
//    //if ((family == "LAGRANGE") && (order == "FIRST"))
//    {
//        std::string meshname = data("mesh/input_mesh_name", "NONE");
//        if (meshname == "NONE")
//        {
//            // allow us to use higher-order approximation.
//            int elX = data("mesh/elX", 10);
//            int elY = data("mesh/elY", 10);
//            int elZ = data("mesh/elZ", 4);
//            double maxX = data("mesh/maxX", 1.0);
//            center_x = 0.5 * maxX;
//            double maxY = data("mesh/maxY", 1.0);
//            double maxZ = data("mesh/maxZ", 0.08);
//            // No reason to use high-order geometric elements if we are
//            // solving with low-order finite elements.
//            MeshTools::Generation::build_cube(mesh, elX, elY, elZ, 0., maxX, 0., maxY, 0., maxZ, TET4);
//
//            std::cout << "Creating subdomains!" << std::endl;
//
//            {
//                double z_interface = data("mesh/z_interface", 0.0);
//
//                MeshBase::element_iterator el = mesh.elements_begin();
//                const MeshBase::element_iterator end_el = mesh.elements_end();
//
//                for (; el != end_el; ++el)
//                {
//                    Elem * elem = *el;
//                    const Point cent = elem->centroid();
//                    // BATH
//                    if (cent(2) > z_interface)
//                    {
//                        elem->subdomain_id() = 2;
//                    }
//                    // TISSUE
//                    else
//                    {
//                        elem->subdomain_id() = 1;
//                    }
//
//                }
//            }
//            mesh.get_boundary_info().regenerate_id_sets();
//            std::cout << "Creating subdomains done!" << std::endl;
//
//        }
//        else
//        {
//            importer.read(meshname);
//            mesh.prepare_for_use();
//        }
//    }
    std::cout << "Mesh done!" << std::endl;

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    std::cout << "Create equation system ..." << std::endl;
    EquationSystems es(mesh);
    // Constructor
    std::cout << "Create bidomain with bath..." << std::endl;
    std::string model = data("model" , "monowave");
    BeatIt::ElectroSolver* solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, es);
    std::string section = data("section", "monowave");
    std::cout << "Calling setup..." << std::endl;
    solver->setup(data, model);
    std::cout << "Calling init ..." << std::endl;
    solver->init(0.0);
    std::cout << "Import fibers ... " << std::endl;
    //bidomain.restart(importer, 1);
    std::cout << "Assembling matrices" << std::endl;
    solver->assemble_matrices(datatime.M_dt);

    // Declare the system and its variables.
    // Create a system named "Poisson"
    std::cout << "Create Boundary Structure" << std::endl;
    std::set<boundary_id_type> endoIDs;
    endoIDs.insert(1);
    std::set<unsigned short> subdomains;
    subdomains.insert(1);
    //bidomain.init_endocardial_ve(endoIDs, subdomains);
    int save_iter = 0;
    int save_iter_ve = 1;
    //bidomain.save_ve_timestep(save_iter_ve, datatime.M_time);
    std::cout << "Init Output" << std::endl;
    //solver->init_exo_output();
    //solver->save_exo_timestep(save_iter++, datatime.M_time);
    solver->save_potential(save_iter, 0.0);
    solver->save_parameters();
    std::string system_mass = data(model+"/diffusion_mass", "mass");
    std::string iion_mass = data(model+"/reaction_mass", "lumped_mass");


    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;

    auto & V = *es.get_system<libMesh::TransientLinearImplicitSystem>("wave").solution;
    auto & Vn = *es.get_system<libMesh::TransientLinearImplicitSystem>("wave").old_local_solution;
    auto & Vnm1 = *es.get_system<libMesh::TransientLinearImplicitSystem>("wave").older_local_solution;
    auto & Iionn = *es.get_system < libMesh::TransientExplicitSystem > ("iion").solution;
    auto & Iionnm1 = *es.get_system < libMesh::TransientExplicitSystem > ("iion").old_local_solution;
    std::cout << "Time loop starts:" << std::endl;

    //export initial condition

    for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
    {
        datatime.advance();
        std::cout << "Time:" << datatime.M_time << std::endl;
        solver->advance();
        //std::cout << "Reaction:" << datatime.M_time << std::endl;
        solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        //solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, "lumpedmass");

        //std::cout << "Diffusion:" << datatime.M_time << std::endl;
        solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        //solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, "mass");

        //std::cout << "at:" << datatime.M_time << std::endl;
        solver->update_activation_time(datatime.M_time);
        //std::cout << "at done:" << datatime.M_time << std::endl;

        //++save_iter_ve;
        //bidomain.save_ve_timestep(save_iter_ve, datatime.M_time);

        if (0 == datatime.M_iter % datatime.M_saveIter)
        {
            save_iter++;
            solver->save_potential(save_iter,datatime.M_time);
            solver->save_activation_times(save_iter);

            //solver->save_potential(save_iter, datatime.M_time);
            //solver->save_exo_timestep(save_iter, datatime.M_time);

        }

    }
    return 0;
}

