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

// <h1>Systems Example 1 - Stokes Equations</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example shows how a simple, linear system of equations
// can be solved in parallel.  The system of equations are the familiar
// Stokes equations for low-speed incompressible fluid flow.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <sys/stat.h>
// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/cell_tet4.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

// The definition of a geometric element
#include "libmesh/elem.h"
#include "libmesh/getpot.h"
#include "libmesh/metis_partitioner.h"
#include "libmesh/checkpoint_io.h"

#include "PoissonSolver/Poisson.hpp"
#include "Util/IO/io.hpp"
#include "Util/GenerateFibers.hpp"
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// The main program.
int main(int argc, char **argv)
{
    // Initialize libMesh.
    LibMeshInit init(argc, argv);

    int np = init.comm().size();
    int rank = init.comm().rank();

    libMesh::Parallel::Communicator rank_comm;
    std::cout << "Split comm" << std::endl;
    init.comm().split(rank, rank, rank_comm);

    std::string exodus_mesh_output = "test_mesh.e";
    std::string mesh_output_name = "split_mesh.cpr";
    std::string mesh_output_name_with_cpu = "split_mesh.cpr/" + std::to_string(np);
    std::string es_output_name = "test_system.xdr";
    bool skip_renumber_nodes_and_elements = true;

    struct stat buffer;
    bool mesh_already_exists = stat(mesh_output_name_with_cpu.c_str(), &buffer) == 0;
    bool es_already_exists = stat(mesh_output_name.c_str(), &buffer) == 0;
    if (0 == rank && (!mesh_already_exists || !es_already_exists))
    {
        libMesh::Mesh mesh(rank_comm);
        mesh.set_mesh_dimension(3);
        mesh.set_spatial_dimension(3);
        mesh.reserve_elem(46154317);
        mesh.reserve_nodes(8736713);

        std::cout << "add nodes" << std::endl;
        std::ifstream nodes_list("03meshLA.pts");
        std::string line;

        double x, y, z;
        int n = 0;
        if (nodes_list.is_open())
        {
            // first line is the number of nodes:
            getline(nodes_list, line);
            while (getline(nodes_list, line))
            {
                std::istringstream ss(line);
                //line will have
                ss >> x >> y >> z;
                x *= 1e-4;
                y *= 1e-4;
                z *= 1e-4;
                mesh.add_point(libMesh::Point(x, y, z), n);
                n++;
            }
        }
        nodes_list.close();


        std::cout << "add elements" << std::endl;
        std::ifstream elements_list("03meshLA.elem");
        int n1, n2, n3, n4, bID, n0;
        std::string elType;
        int c = 0;
        if (elements_list.is_open())
        {
            getline(elements_list, line);
            while (getline(elements_list, line))
            {
                std::istringstream ss(line);
                libMesh::Elem *elem = mesh.add_elem(new libMesh::Tet4);
                ss >> elType >> n0 >> n1 >> n2 >> n3 >> bID;
                elem->set_node(0) = mesh.node_ptr(n0);
                elem->set_node(1) = mesh.node_ptr(n1);
                elem->set_node(2) = mesh.node_ptr(n2);
                elem->set_node(3) = mesh.node_ptr(n3);
                elem->subdomain_id() = bID;
                c++;
            }
        }
        elements_list.close();
        std::cout << "prepare for use" << std::endl;
        mesh.prepare_for_use(skip_renumber_nodes_and_elements);
        //std::cout << "write to exo" << std::endl;
        //mesh.write(exodus_mesh_output);

        mesh.print_info();
        if (!es_already_exists)
        {
            EquationSystems equation_systems(mesh);
            std::cout << "Create system" << std::endl;
            ExplicitSystem &system = equation_systems.add_system < ExplicitSystem > ("fibers");
            std::cout << "Add variables" << std::endl;
            system.add_variable("fibersx", CONSTANT, MONOMIAL);
            system.add_variable("fibersy", CONSTANT, MONOMIAL);
            system.add_variable("fibersz", CONSTANT, MONOMIAL);
            std::cout << "Init System" << std::endl;
            equation_systems.init();

            std::cout << "Opening fiber file " << std::endl;
            std::ifstream fibers_file("03meshLA.lon");
            double fx, fy, fz;
            int elID = 0;
            std::vector<unsigned int> indices;
            std::cout << "Reading and assigning fibers .. this will take a while! " << std::endl;
            if (fibers_file.is_open())
            {
                while (getline(fibers_file, line))
                {
                    std::istringstream ss(line);
                    //line will have
                    ss >> fx >> fy >> fz;

                    Elem *el = mesh.elem_ptr(elID);
                    system.get_dof_map().dof_indices(el, indices);
                    system.solution->set(indices[0], fx);
                    system.solution->set(indices[1], fy);
                    system.solution->set(indices[2], fz);
                    elID++;
                }
            }
            fibers_file.close();
            //std::cout << "exporting fibers exo... " << std::endl;
            //ExodusII_IO(mesh).write_equation_systems("test_fibers.e", equation_systems);
            std::cout << "exporting es... " << std::endl;
            equation_systems.write(es_output_name);
        }
        std::cout << "Create partitioner " << std::endl;
        MetisPartitioner partitioner;
        std::cout << "partition " << std::endl;
        partitioner.partition(mesh, np);

        std::cout << "write to file " << std::endl;
        CheckpointIO cpr(mesh);
        cpr.current_processor_ids().clear();
        for (unsigned int i = 0; i < np; i++)
            cpr.current_processor_ids().push_back(i);
        cpr.current_n_processors() = np;
        cpr.parallel() = true;
        cpr.binary() = true;
        std::cout << mesh_output_name << std::endl;
        cpr.write(mesh_output_name);
    }
    init.comm().barrier();

    // across the default MPI communicator.
//    std::string mesh_name = "../example_AFib_fibers/fastl/fastl03.e";
//    ReplicatedMesh * mesh_fibers;
//    ReplicatedMesh * mesh_to_split;
//
//
//    libMesh::Parallel::Communicator cpu_0_comm;
//    std::cout << "Split comm" << std::endl;
//
//    ExodusII_IO * importer_ptr;
//    std::string mesh_output_name = "split_mesh.cpr";
//    std::string es_output_name = "test_system.xdr";
//
//    if(rank == 0)
//    {
//        std::cout << "Create replicated mesh" << std::endl;
//        mesh_fibers = new Mesh(cpu_0_comm);
//        std::cout << "Create importer" << std::endl;
//        //ParallelObject()
//        importer_ptr = new ExodusII_IO(*mesh_fibers);
//        std::cout << "Read mesh" << std::endl;
//        importer_ptr->read(mesh_name);
//        mesh_fibers->print_info();
//
//
//        {
//            std::cout << "Create ES" << std::endl;
//            EquationSystems equation_systems(*mesh_fibers);
//            std::cout << "Create system" << std::endl;
//            ExplicitSystem & system = equation_systems.add_system<ExplicitSystem>("fs");
//            std::cout << "Add variables" << std::endl;
//            system.add_variable("fibersx", CONSTANT, MONOMIAL);
//            system.add_variable("fibersy", CONSTANT, MONOMIAL);
//            system.add_variable("fibersz", CONSTANT, MONOMIAL);
//            std::cout << "Init System" << std::endl;
//            equation_systems.init();
//            // Prints information about the system to the screen.
//            equation_systems.print_info();
//            std::cout << "Copy Solutions" << std::endl;
////            importer_ptr->copy_elemental_solution(system,"fibersx","fibersx");
////            importer_ptr->copy_elemental_solution(system,"fibersy","fibersy");
////            importer_ptr->copy_elemental_solution(system,"fibersz","fibersz");
//            std::cout << "Copied solution" << std::endl;
//
//            std::cout << "Opening file " << std::endl;
//            std::ifstream fibers_file("03meshLA.lon");
//            std::string line;
//            double fx, fy, fz;
//            int elID = 0;
//            std::vector<unsigned int> indices;
//            std::cout << "Reading and assigning fibers .. this will take a while! " << std::endl;
//            if(fibers_file.is_open())
//            {
//                  while ( getline (fibers_file,line) )
//                  {
//                      std::istringstream ss(line);
//                      //line will have
//                      ss >> fx >> fy >> fz;
//
//                      Elem  * el = mesh_fibers->elem_ptr(elID);
//                      system.get_dof_map().dof_indices(el, indices);
//                      system.solution->set(indices[0], fx);
//                      system.solution->set(indices[1], fy);
//                      system.solution->set(indices[2], fz);
//                      elID++;
//                  }
//            }
//            fibers_file.close();
//            std::cout << "exporting ... " << std::endl;
//            //ExodusII_IO(*mesh_fibers).write_equation_systems("test_fibers.e", equation_systems);
//            equation_systems.write(es_output_name);
//
//
//            std::cout << "Create partitioner " << std::endl;
//            MetisPartitioner partitioner;
//            std::cout << "partition " << std::endl;
//            partitioner.partition(*mesh_to_split, np);
//
//            std::cout << "write to file " << std::endl;
//            CheckpointIO cpr(*mesh_to_split);
//            cpr.current_processor_ids().clear();
//            for (unsigned int i = 0; i < np; i++)
//                cpr.current_processor_ids().push_back(i);
//            cpr.current_n_processors() = np;
//            cpr.parallel() = true;
//            cpr.binary() = true;
//            std::cout << mesh_output_name << std::endl;
//            cpr.write(mesh_output_name);
//        }
//    }

//    init.comm().barrier();

//    std::string mesh_output_name = "mesh.cpr";
//    if(rank == 0)
//    {
//        mesh_to_split = new ReplicatedMesh(cpu_0_comm);
//        mesh_to_split->read(mesh_name);
//        MetisPartitioner partitioner;
//        partitioner.partition(*mesh_to_split, np);
//
//        CheckpointIO cpr(*mesh_to_split);
//        cpr.current_processor_ids().clear();
//        for (unsigned int i = 0; i < np; i++)
//            cpr.current_processor_ids().push_back(i);
//        cpr.current_n_processors() = np;
//        cpr.parallel() = true;
//        cpr.binary() = true;
//        std::cout << mesh_output_name << std::endl;
//        cpr.write(mesh_output_name);
//    }

    // Read the model from the input file
    GetPot data_ep = BeatIt::readInputFile(argc, argv);

    std::cout << "Create parallel mesh" << std::endl;
    ParallelMesh parallel_mesh(init.comm());
    parallel_mesh.read(mesh_output_name);
    parallel_mesh.print_info();

    std::cout << "Create parallel ES" << std::endl;
    EquationSystems parallel_equation_systems(parallel_mesh);
    std::cout << "read fiber system" << std::endl;
    parallel_equation_systems.read(es_output_name);

    ExplicitSystem &fiber_system = parallel_equation_systems.get_system < ExplicitSystem > ("fibers");
    double norm = fiber_system.solution->linfty_norm();
    std::cout << "fiber L infty norm: " << norm << std::endl;

    bool export_fibers = data_ep("export_fibers", true);
    if(export_fibers)
    {
        std::cout << "exporting ... " << std::endl;
        ExodusII_IO exporter(parallel_mesh);
        exporter.write_equation_systems("parallel_fibers.e", parallel_equation_systems);
        exporter.write_element_data(parallel_equation_systems);
    }
    //parallel_equation_systems.init();
    norm = fiber_system.solution->linfty_norm();
    std::cout << "fiber L infty norm: " << norm << std::endl;
    // Prints information about the system to the screen.
    parallel_equation_systems.print_info();

    std::cout << "Create EP ..." << std::endl;

    std::string model = data_ep("model", "monowave");
    // Create the TimeData object
    std::cout << "Create Datatime object ..." << std::endl;
    BeatIt::TimeData datatime;
    // Set it up using the input file
    datatime.setup(data_ep, "");
    // Output on screen the stored variables
    datatime.print();
    std::cout << "Create Electrophysiology model ..." << std::endl;
    // Create the model
    BeatIt::ElectroSolver *solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, parallel_equation_systems);
    // Setup the EP model using the input file
    std::cout << "Calling setup..." << std::endl;
    solver->setup(data_ep, model);
    norm = fiber_system.solution->linfty_norm();
    std::cout << "fiber L infty norm: " << norm << std::endl;
    // Initialize systems
    std::cout << "Calling init ..." << std::endl;
    // Set up initial conditions at time
    norm = fiber_system.solution->linfty_norm();
    std::cout << "fiber L infty norm: " << norm << std::endl;
    solver->init(datatime.M_startTime);
    ////////////////////////////////////////////
    // keep going with the EP
    ////////////////////////////////////////////
    norm = fiber_system.solution->linfty_norm();
    std::cout << "fiber L infty norm: " << norm << std::endl;
    ExplicitSystem &sheets_system = parallel_equation_systems.get_system < ExplicitSystem > ("sheets");
    norm = sheets_system.solution->linfty_norm();
    std::cout << "sheets L infty norm: " << norm << std::endl;
    ExplicitSystem &xfiber_system = parallel_equation_systems.get_system < ExplicitSystem > ("xfibers");
    norm = xfiber_system.solution->linfty_norm();
    std::cout << "xfiber L infty norm: " << norm << std::endl;
    parallel_equation_systems.print_info();
    std::cout << "Assembling matrices" << std::endl;
    solver->assemble_matrices(datatime.M_dt);
    // output file counter
    int save_iter = 0;
    // Export initial condition at time
    //    solver->save_exo_timestep(save_iter, datatime.M_time);
    solver->save_potential(save_iter, datatime.M_startTime);
    // Parameters to save the activation times
    // A node is activated if the transmembrane potential > threshold
    double threshold = data_ep("threshold", -10.0);
    // We export the activation times every at_save_iter iterations
    int at_save_iter = data_ep("at_save_iter", 25);

    // Old parameters that define the method, not to be changed
    // TODO: clean this part
    std::string system_mass = data_ep(model + "/diffusion_mass", "mass");
    std::string iion_mass = data_ep(model + "/reaction_mass", "lumped_mass");
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
            solver->save_activation_times_nemesis(save_iter);
        }

    }

    // // export activation times again
    solver->save_activation_times_nemesis(save_iter);
    //solver->save_activation_times(save_iter);
    // delete solver before ending the simulation
    // avoiding memory leaks
    delete solver;

//    delete mesh_fibers;
//    delete mesh_to_split;
    // All done.
}
