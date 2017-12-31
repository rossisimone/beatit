

// C++ include files that we need
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include <libmesh/boundary_mesh.h>
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
#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "Util/Timer.hpp"
#include <libmesh/point_locator_tree.h>


#include "libmesh/dof_map.h"

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

    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data_at.beat", 2, "-j", "--jnput");
    GetPot data_data(datafile_name);

    std::string simulation_datafile_name = commandLine.follow("data.beat", 2, "-i", "--input");
    GetPot data(simulation_datafile_name);



    // IMPORT DATA
    std::cout << "Importing DATA: " << std::endl;
    Mesh mesh_data(init.comm());
    ExodusII_IO data_importer(mesh_data);
    std::string data_file = data_data("data_file", "NONE");
    if ("NONE" != data_file )
     {
        data_importer.read(data_file);
        mesh_data.allow_renumbering (false);
        mesh_data.allow_remote_element_removal(true);
        mesh_data.prepare_for_use();
    }
    std::cout << "Mesh imported! " << std::endl;
    auto& nodal_var_names = data_importer.get_nodal_var_names ();
    std::cout << "available fields:" << std::endl;
    for(auto && var : nodal_var_names) std::cout << var << std::endl;

    libMesh::EquationSystems data_es(mesh_data);
    auto& data_system = data_es.add_system<libMesh::ExplicitSystem>("data");
    data_system.add_variable("Recording", libMesh::FIRST, libMesh::LAGRANGE);
    std::cout << "Data Equation system created! " << std::endl;

    auto& sim_system = data_es.add_system<libMesh::ExplicitSystem>("simulation");
    sim_system.add_variable("Simulation", libMesh::FIRST, libMesh::LAGRANGE);
    std::cout << "Sim Equation system created! " << std::endl;
    data_es.init();
    std::cout << "Systems created! " << std::endl;

    std::cout << "Copying solution: " << std::endl;
    data_importer.copy_nodal_solution(data_system, "Recording", "ActivationTimes");
    std::cout << "Done! " << std::endl;
    std::cout << "n_nodes: " << mesh_data.n_local_nodes() << std::endl;
    std::cout << "parallel_n_nodes : " << mesh_data.parallel_n_nodes () << std::endl;
    mesh_data.renumber_nodes_and_elements();
    std::cout << "n_nodes: " << mesh_data.n_local_nodes() << std::endl;
    std::cout << "parallel_n_nodes : " << mesh_data.parallel_n_nodes () << std::endl;
    libMesh::MeshTools::Modification::scale(mesh_data, 0.1, 0.1, 0.1);


    // IMPORT SIMULATION

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());
    // READ SOLUTION
    ExodusII_IO importer(mesh);
    std::string restart_file = data_data("restart_file", "NONE");
    if ("NONE" != restart_file )
   {
       importer.read(restart_file);
       mesh.prepare_for_use();
   }
   else
   {
       std::cout << "restart file has not been specified: " << restart_file << std::endl;
       std::cout << "ABORTING!" << std::endl;
       return 0;
   }

    // Constructor
    perf_log.push("create es");
    std::cout << "Create equation system ..." << std::endl;
    libMesh::EquationSystems es(mesh);
    perf_log.pop("create es");

    std::cout << "Create a bidomain with bath solver to import the solution " << std::endl;
    BeatIt::BidomainWithBath solver(es);
    solver.setup(data,"bidomain");
    solver.init(0.0);
    int restart_step = data_data("step", 1);
    solver.restart(importer, restart_step);



    // Populate the new Simulation field in the data mesh
    std::cout << "Populating the mesh" << std::endl;
    std::set<unsigned short> subdomains;
    subdomains.insert(1);

    BoundaryMesh endo(init.comm());
    std::set<boundary_id_type> endoIDs;
    endoIDs.insert(1);

    mesh.boundary_info->sync(endoIDs, endo, subdomains);
    std::map< dof_id_type, dof_id_type > node_id_map;
    std::map< dof_id_type, unsigned char > side_id_map;
    mesh.boundary_info->get_side_and_node_maps(endo, node_id_map, side_id_map);
    std::map< dof_id_type, dof_id_type > reverse_node_id_map;
    for(auto it = node_id_map.begin(); it != node_id_map.end(); ++it)
    {
        reverse_node_id_map[it->second] = it->first;
    }
    ExodusII_IO(endo).write("endo.exo");

    int n_data_nodes = mesh_data.n_local_nodes();
    std::map<int, int> node_to_node_map;

    //libMesh::PointLocatorTree locator(mesh, Trees::NODES);
    libMesh::PointLocatorTree locator(mesh, Trees::ELEMENTS);


    libMesh::MeshBase::const_node_iterator node = mesh_data.local_nodes_begin();
    const libMesh::MeshBase::const_node_iterator end_node =
            mesh_data.local_nodes_end();
    double tolerance = data_data("tolerance", 1e-2);

    ExplicitSystem& activation_times_system = es.get_system
            < ExplicitSystem > ("activation_times");

    const libMesh::DofMap & dof_map = sim_system.get_dof_map();
    const libMesh::DofMap & dof_map_at = activation_times_system.get_dof_map();

    std::vector < libMesh::dof_id_type > dof_indices_sim;
    std::vector < libMesh::dof_id_type > dof_indices_at;

    std::cout << "Looping over the nodes" << std::endl;
    for (; node != end_node; ++node)
    {
        const libMesh::Node * nn = *node;
        libMesh::Point p_data((*nn)(0), (*nn)(1), (*nn)(2));

        libMesh::MeshBase::const_node_iterator enode = endo.local_nodes_begin();
        const libMesh::MeshBase::const_node_iterator end_enode =
                endo.local_nodes_end();

        double dist_min = 1e9;
        int ID = -1;
        libMesh::MeshBase::const_node_iterator closest_node = enode;
        for(; enode != end_enode; ++enode)
        {
            const libMesh::Node * enn = *enode;
            libMesh::Point p((*enn)(0), (*enn)(1), (*enn)(2));

            p -= p_data;
            double current_dist = p.norm_sq();
            if(current_dist < dist_min)
            {
                closest_node = enode;
                ID = (*closest_node)->id();
                dist_min = current_dist;
            }
        }

        dof_map.dof_indices(nn, dof_indices_sim, 0);


        const libMesh::Node * enn = mesh.node_ptr(reverse_node_id_map.find(ID)->second);

        dof_map_at.dof_indices(enn, dof_indices_at, 0);
        double value = (*activation_times_system.solution)(dof_indices_at[0]);
        sim_system.solution->set(dof_indices_sim[0], value);

//
//        libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));
//
//        auto * node = locator.locate_node(p,&subdomains,tolerance);
//        //auto * node = locator.locate_node(p,nullptr,tolerance);
//        if(node)
//        {
//            //std::cout << "Found node!" << std::endl;
//            dof_map_at.dof_indices(node, dof_indices_at, 0);
//            double value = (*activation_times_system.solution)(dof_indices_at[0]);
//            sim_system.solution->set(dof_indices_sim[0], value);
//        }
    }

    std::cout << "Exporting " << std::endl;
    ExodusII_IO exporter(mesh_data);
    exporter.write_equation_systems("comparison.exo",data_es);
    std::cout << "Done " << std::endl;









    return 0;
}
