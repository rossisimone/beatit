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
#include "libmesh/parallel_mesh.h"
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

// Bring in everything from the libMesh namespace
using namespace libMesh;

namespace Data
{
static std::string path = "";
}
void ic_patient2_plane_wave(DenseVector<Number> & output, const Point & p, const Real)
{
    double x = p(0);
    double y = p(1);
    double z = p(2);

    Point N(-0.11, 0.95, 0.28);
    N.unit();
    Point X0(-1.86, -4.6, 9.4);

    double ic = -81.2;

    if( N.contract(p-X0) < 0 )
    {
        ic = 10.0;
    }
    output(0) = ic;
}

void ic_patient2_ext(DenseVector<Number> & output, const Point & p, const Real)
{
    double x = p(0);
    double y = p(1);
    double z = p(2);
    double ic = -81.2;
    if ((x +0.4520240907997266) * (x +0.4520240907997266) + (y+6.631625) * (y+6.631625) + (z - 10.13692) * (z - 10.13692) <= 4.0)
    {
        ic = 10.0;
    }
    output(0) = ic;
}

struct EndoRankMap
{
    EndoRankMap() : data_unique_ID(), mesh_ID(), rank(), endo_ID() {}
    std::vector<int>  data_unique_ID;
    std::vector<int>  mesh_ID;
    std::vector<int>  endo_ID;
    std::vector<int>  rank;
    bool check()
    {
        int size1 = data_unique_ID.size();
        int size2 = mesh_ID.size();
        int size3 = rank.size();
        int size4 = endo_ID.size();
        if(size1 == size2 && size2 == size3 && size3 == size4) return true;
        else return false;
    }
};


void init_files( libMesh::MeshBase& endo,
                 EndoRankMap& id_map, std::map<dof_id_type, dof_id_type>& reverse_node_id_map, std::string meshfile, int unipole)
{
    std::cout << "initializing: " << meshfile << std::endl;
    ReplicatedMesh mesh_unip1(endo.comm());
    mesh_unip1.read(meshfile, nullptr, true, true);
    mesh_unip1.print_info(std::cout);
    libMesh::MeshBase::const_node_iterator node_data = mesh_unip1.nodes_begin();
    const libMesh::MeshBase::const_node_iterator end_node_data = mesh_unip1.nodes_end();
    unsigned int my_rank = endo.comm().rank();

    //std::cout << "endo num nodes: " << endo.n_nodes() << std::endl;
    int n_unip_nodes = mesh_unip1.n_nodes();
    std::vector<int> data_IDs;
    std::vector<int> mesh_IDs;
    std::vector<int> ranks;
    std::vector<int> endo_IDs;
    int c = 0;
    for (; node_data != end_node_data; ++node_data)
    {
        c++;
        const libMesh::Node * nn = *node_data;
        auto dataID = nn->id();

        libMesh::Point p_data((*nn)(0), (*nn)(1), (*nn)(2));
        //std::cout << p_data(0) << ", " << p_data(1) << ", " << p_data(2) << std::endl;
        libMesh::MeshBase::const_node_iterator enode = endo.pid_nodes_begin(my_rank);
        const libMesh::MeshBase::const_node_iterator end_enode = endo.pid_nodes_end(my_rank);
        double dist_min = 1000000.0;
        int ID = -1;
        //libMesh::MeshBase::const_node_iterator closest_node = enode;
        libMesh::Node * cs_node;
        for (; enode != end_enode; ++enode)
        {
            libMesh::Node * enn = *enode;
            libMesh::Point p((*enn)(0), (*enn)(1), (*enn)(2));
            //std::cout << "p_data = " << p_data(0) << ", " << p_data(1) << ", " << p_data(2) << std::endl;
            //std::cout << "p = " << p(0) << ", " << p(1) << ", " << p(2) << std::endl;
            p -= p_data;
            double current_dist = p.norm_sq();
            //std::cout << "p_diff = " << p(0) << ", " << p(1) << ", " <<  p(2) << ", norm: " << current_dist << std::endl;
            //std::cout << current_dist << std::endl;
            if (current_dist < dist_min)
            {
                cs_node = enn;
                ID = cs_node->id();
                dist_min = current_dist;
            }
        }

        unsigned int rank = endo.comm().rank();

        //std::cout << my_rank << ": cn" << &closest_node << ", ID: " << ID << ", dist: " <<  dist_min << std::endl;
        //std::cout << my_rank << ": " << rank << ", min: " << dist_min << std::endl;
        //endo.comm().maxloc(dist_min, rank);
        double global_dist_min = dist_min;
        //endo.comm().min(global_dist_min);

        struct {
            double val;
            int rank;
        } mp, mp_global;

        mp.val=global_dist_min;
        mp.rank=my_rank;

        MPI_Allreduce(&mp, &mp_global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        global_dist_min = mp_global.val;
        rank = mp_global.rank;
        //std::cout << my_rank << ": minloc " << mp_global.rank << ", min: " << mp_global.val << ", ID: " << ID << std::endl;
        //std::cout << my_rank << ": " << rank << std::endl;
        //std::cout <<  "rank " << rank << " writes the file" << std::endl;

        //if (rank == my_rank)

        if(dist_min == global_dist_min && rank == my_rank)
        {
            int mesh_ID = reverse_node_id_map.find(ID)->second;

            data_IDs.push_back(dataID);
            mesh_IDs.push_back(mesh_ID);
            ranks.push_back(rank);
            endo_IDs.push_back(ID);


            std::string filename;
            if(unipole == 1)
            {
                filename = Data::path + "/uni1_node_" + std::to_string(dataID) + ".txt";
            }
            else if(unipole == 2)
            {
                filename = Data::path + "/uni2_node_" + std::to_string(dataID) + ".txt";
            }
            else
            {
                filename = Data::path + "/at_node_" + std::to_string(dataID) + ".txt";
            }

            std::ofstream file(filename);

            std::streamsize ss = file.precision();
            file.precision(10);
            file << "Data Point: " << dataID << std::endl;
            file << (*nn)(0) << " " << (*nn)(1) << " " << (*nn)(2) << std::endl;
            file << "Point: " << mesh_ID << std::endl;
            //const libMesh::Node * cnn = *closest_node;
            file << (*cs_node)(0) << " " << (*cs_node)(1) << " " << (*cs_node)(2) << std::endl;
            file.precision(ss);

            file << "time Ve V Q at" << std::endl;
            file.close();
            ID = mesh_ID;

        }
        endo.comm().broadcast(ID, rank);
        //std::cout <<  my_rank <<": rank " << rank << " ID: " << ID << ", dataID: " << dataID << ", dist min: " << global_dist_min << std::endl;
    }
    endo.comm().barrier();
    endo.comm().allgather(ranks);
    endo.comm().allgather(data_IDs);
    endo.comm().allgather(mesh_IDs);
    endo.comm().allgather(endo_IDs);

    int size1 = ranks.size();
    int endo_n_nodes = endo.n_nodes();
    if(size1 != n_unip_nodes)
    {
        throw std::runtime_error("Wrong num nodes");
    }

    int size2 = data_IDs.size();
    std::cout << "data_IDs: " << size2 << std::endl;
    if(size2 != n_unip_nodes) throw std::runtime_error("Wrong num nodes");
    int size3 = mesh_IDs.size();
    std::cout << "mesh_IDs: " << size3 << std::endl;
    if(size3 != n_unip_nodes) throw std::runtime_error("Wrong num nodes");

    std::cout << "swapping vectors" << std::endl;
    id_map.data_unique_ID.swap(data_IDs);
    id_map.mesh_ID.swap(mesh_IDs);
    id_map.rank.swap(ranks);
    id_map.endo_ID.swap(endo_IDs);
    std::cout << "done" << std::endl;

    if(id_map.check())
    {
        std::cout << "size id map: " << id_map.data_unique_ID.size() << std::endl;
    }
    else
    {
        throw std::runtime_error("WWRONG!!!");
    }
    std::cout << "initializing: " << meshfile << ". done."<< std::endl;
}


void find_point(libMesh::MeshBase& endo,
                libMesh::Point& point,
                int data_ID,
                libMesh::Point& cs_point,
                EndoRankMap& id_map,
                std::map<dof_id_type, dof_id_type>& reverse_node_id_map,
                std::string file_name = "" )
{
    std::vector<int> mesh_IDs;
    std::vector<int> ranks;
    std::vector<int> endo_IDs;

    unsigned int my_rank = endo.comm().rank();
    libMesh::MeshBase::const_node_iterator enode = endo.pid_nodes_begin(my_rank);
    const libMesh::MeshBase::const_node_iterator end_enode = endo.pid_nodes_end(my_rank);

    // Define minimum distance and ID of the point to find on the endo mesh
    double dist_min = 1000000.0;
    int ID = -1;
    // Fine closest node
    libMesh::Node * cs_node;
    for (; enode != end_enode; ++enode)
    {
        libMesh::Node * enn = *enode;
        libMesh::Point p((*enn)(0), (*enn)(1), (*enn)(2));
        p -= point;
        double current_dist = p.norm_sq();
        if (current_dist < dist_min)
        {
            cs_node = enn;
            ID = cs_node->id();
            dist_min = current_dist;
        }
    }

    // On each processor we have found the closest node
    // Now we need to take closest over all processors
    unsigned int rank = endo.comm().rank();
    double global_dist_min = dist_min;
    // We define a local class that will contain the local and global values
    struct {
        double val;
        int rank;
    } mp, mp_global;
    // Initialize the local values
    mp.val=global_dist_min;
    mp.rank=my_rank;
    // Find the minimum using AllReduce
    MPI_Allreduce(&mp, &mp_global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    // Give the values just found to everyone using local variables
    global_dist_min = mp_global.val;
    rank = mp_global.rank;

    double x = (*cs_node)(0);
    double y = (*cs_node)(1);
    double z = (*cs_node)(2);

    endo.comm().broadcast(x, rank);
    endo.comm().broadcast(y, rank);
    endo.comm().broadcast(z, rank);

    cs_point(0) = x;
    cs_point(1) = y;
    cs_point(2) = z;
    // Populate EndoRank map
    // Do this only if we are the processor that holds the closest point
    if(dist_min == global_dist_min && rank == my_rank)
    {
        // Given the ID of the surface mesh, find the ID on the volumetric mesh
        int mesh_ID = reverse_node_id_map.find(ID)->second;
        // store the data ID
        id_map.data_unique_ID.push_back(data_ID);
        // store the volumetric mesh ID
        id_map.mesh_ID.push_back(mesh_ID);
        // store rank where this point is
        id_map.rank.push_back(rank);
        // store the surface mesh ID
        id_map.endo_ID.push_back(ID);

        // If the filename is not empty,
        // initialize the file for output
        if(file_name != "")
        {
            std::ostringstream data_ss;
            data_ss << std::setw(4) << std::setfill('0') << data_ID;
            std::string data_ID_str = data_ss.str();

            std::string filename;
            filename = Data::path + "/" + file_name + data_ID_str + ".txt";
            std::ofstream file(filename);
            std::streamsize ss = file.precision();
            file.precision(10);
            // Write the coordinates of the original point
            file << "Data Point: " << data_ID << std::endl;
            file << point(0) << " " << point(1) << " " << point(2) << std::endl;
            // Write the coordinates of the mesh point
            file << "Point: " << mesh_ID << std::endl;
            file << (*cs_node)(0) << " " << (*cs_node)(1) << " " << (*cs_node)(2) << std::endl;
            file.precision(ss);
            file << "time Ve V Q at" << std::endl;
            file.close();
        }
    }
}

void project_unipolar_data( libMesh::MeshBase& endo,
                             EndoRankMap& unipolar1_id_map,
                             EndoRankMap& unipolar2_id_map,
                             EndoRankMap& bipolar_id_map,
                             std::map<dof_id_type, dof_id_type>& reverse_node_id_map)
{
    std::string unip1_data = "uni1_unprojected.csv";
    std::string unip2_data = "uni2_unprojected.csv";
    std::string bip_data = "bipolar_unprojected.csv";
    std::string line_uni1;
    std::string line_uni2;
    std::string line_bip;

    std::ifstream unip1_file (unip1_data);
    std::ifstream unip2_file (unip2_data);
    std::ifstream bip_file (bip_data);

    libMesh::Point bip_point;
    libMesh::Point unip1_point;
    libMesh::Point unip2_point;

    // Read files
    int data_ID = 0;

    EndoRankMap tmp_bipolar_map;
    EndoRankMap tmp_unipolar1_map;
    EndoRankMap tmp_unipolar2_map;

    std::ofstream bipolar_projected_file("bipolar_projected_points.csv");
    std::ofstream unipolar1_projected_file("unipolar1_projected_points.csv");
    std::ofstream unipolar2_projected_file("unipolar2_projected_points.csv");

    while( std::getline(bip_file,line_bip) )
    {
        data_ID++;

        std::getline(unip1_file,line_uni1);
        std::getline(unip2_file,line_uni2);

        std::istringstream bip_str(line_bip);
        std::istringstream unip1_str(line_uni1);
        std::istringstream unip2_str(line_uni2);
        std::string coordinate;
        std::stringstream value;

        // populate points
        for(int i = 0; i < 3; ++i)
        {
            std::getline(bip_str,coordinate,',');
            bip_point(i) = std::stod( coordinate );

            std::getline(unip1_str,coordinate,',');
            unip1_point(i) = std::stod( coordinate );

            std::getline(unip2_str,coordinate,',');
            unip2_point(i) = std::stod( coordinate );
        } // populate points

        std::cout << "\nPoint ID: " << data_ID << std::endl;
        bip_point.print();
        std::cout << std::endl;
        // Find the closest point to the bipolar position, and push back data to tmp_bipolar_map
        libMesh::Point bipolar_closest_point;
        find_point(endo, bip_point, data_ID, bipolar_closest_point, tmp_bipolar_map, reverse_node_id_map, "at_node_");
        bipolar_projected_file << bipolar_closest_point(0) << ", ";
        bipolar_projected_file << bipolar_closest_point(1) << ", ";
        bipolar_projected_file << bipolar_closest_point(2)  << ", ";
        bipolar_projected_file << data_ID << std::endl;
        bipolar_closest_point.print();

        // Compute displacement vecotr between bipolar and unipolar positions
        libMesh::Point unipolar1_dislacement(unip1_point);
        unipolar1_dislacement -= bip_point;
        // Displace the projected bipolar point with the displacement vector
        libMesh::Point displaced_unipolar1_point(bipolar_closest_point);
        displaced_unipolar1_point += unipolar1_dislacement;
        // Find point closest to the displaced bipolar point
        libMesh::Point unipolar1_closest_point;
        find_point(endo, displaced_unipolar1_point, data_ID, unipolar1_closest_point, tmp_unipolar1_map, reverse_node_id_map, "uni1_node_");
        unipolar1_projected_file << unipolar1_closest_point(0) << ", ";
        unipolar1_projected_file << unipolar1_closest_point(1) << ", ";
        unipolar1_projected_file << unipolar1_closest_point(2) << ", ";
        unipolar1_projected_file << data_ID << std::endl;

        // Repeat for unipolar position 2:
        // Compute displacement vecotr between bipolar and unipolar positions
        libMesh::Point unipolar2_dislacement(unip2_point);
        unipolar2_dislacement -= bip_point;
        // Displace the projected bipolar point with the displacement vector
        libMesh::Point displaced_unipolar2_point(bipolar_closest_point);
        displaced_unipolar2_point += unipolar2_dislacement;
        // Find point closest to the displaced bipolar point
        libMesh::Point unipolar2_closest_point;
        find_point(endo, displaced_unipolar2_point, data_ID, unipolar2_closest_point, tmp_unipolar2_map, reverse_node_id_map, "uni2_node_");
        unipolar2_projected_file << unipolar2_closest_point(0) << ", ";
        unipolar2_projected_file << unipolar2_closest_point(1) << ", ";
        unipolar2_projected_file << unipolar2_closest_point(2) << ", ";
        unipolar2_projected_file << data_ID << std::endl;

        //if (data_ID == 2) exit(-1);
    } // Read file line

    bip_file.close();
    unip1_file.close();
    unip2_file.close();

    bipolar_projected_file.close();
    unipolar1_projected_file.close();
    unipolar2_projected_file.close();
    // Gather the the maps:
    endo.comm().barrier();
    endo.comm().allgather(tmp_bipolar_map.rank);
    endo.comm().allgather(tmp_unipolar1_map.rank);
    endo.comm().allgather(tmp_unipolar2_map.rank);

    endo.comm().allgather(tmp_bipolar_map.data_unique_ID);
    endo.comm().allgather(tmp_unipolar1_map.data_unique_ID);
    endo.comm().allgather(tmp_unipolar2_map.data_unique_ID);

    endo.comm().allgather(tmp_bipolar_map.mesh_ID);
    endo.comm().allgather(tmp_unipolar1_map.mesh_ID);
    endo.comm().allgather(tmp_unipolar2_map.mesh_ID);

    endo.comm().allgather(tmp_bipolar_map.endo_ID);
    endo.comm().allgather(tmp_unipolar1_map.endo_ID);
    endo.comm().allgather(tmp_unipolar2_map.endo_ID);

    // swap the gathered vector with the actual ones
    bipolar_id_map.data_unique_ID.swap( tmp_bipolar_map.data_unique_ID);
    bipolar_id_map.rank.swap( tmp_bipolar_map.rank);
    bipolar_id_map.mesh_ID.swap( tmp_bipolar_map.mesh_ID);
    bipolar_id_map.endo_ID.swap( tmp_bipolar_map.endo_ID);

    unipolar1_id_map.data_unique_ID.swap( tmp_unipolar1_map.data_unique_ID);
    unipolar1_id_map.rank.swap( tmp_unipolar1_map.rank);
    unipolar1_id_map.mesh_ID.swap( tmp_unipolar1_map.mesh_ID);
    unipolar1_id_map.endo_ID.swap( tmp_unipolar1_map.endo_ID);

    unipolar2_id_map.data_unique_ID.swap( tmp_unipolar2_map.data_unique_ID);
    unipolar2_id_map.rank.swap( tmp_unipolar2_map.rank);
    unipolar2_id_map.mesh_ID.swap( tmp_unipolar2_map.mesh_ID);
    unipolar2_id_map.endo_ID.swap( tmp_unipolar2_map.endo_ID);

}

void write_to_file(libMesh::MeshBase& mesh, const EndoRankMap& id_map, double time, libMesh::TransientLinearImplicitSystem& V_sys,
        libMesh::TransientLinearImplicitSystem& bid_sys, libMesh::ExplicitSystem& at_sys, int unipole)
{
    std::cout << "write_to_file" << std::endl;
    int n = id_map.data_unique_ID.size();
    int my_rank = mesh.comm().rank();
//    for (auto && elm : rank_point_map)
    std::cout << my_rank << ": my rank size: " << n << std::endl;
    for(int i = 0; i < n; ++i)
    {
        //std::cout << i << std::endl;

        if (id_map.rank[i] == my_rank)
        {

            int ID = id_map.mesh_ID[i];
            int dataID = id_map.data_unique_ID[i];
            //std::cout << my_rank << ": " << id_map.rank[i] << ", " << ID << ", " << id_map.endo_ID[i] << ", " << dataID << std::endl;
            const libMesh::Node * nn = mesh.query_node_ptr(ID);
//            // get V
            libMesh::DofMap * dof_map = &V_sys.get_dof_map();
            std::vector<libMesh::dof_id_type> dof_indices;
            dof_map->dof_indices(nn, dof_indices, 0);
            double V = (*V_sys.current_local_solution)(dof_indices[0]);
//
//            // get at
            dof_indices.clear();
            dof_map = &at_sys.get_dof_map();
            dof_map->dof_indices(nn, dof_indices, 0);
            double at = (*at_sys.current_local_solution)(dof_indices[0]);

            // get Q & Ve
            dof_indices.clear();
            dof_map = &bid_sys.get_dof_map();
            dof_map->dof_indices(nn, dof_indices);
            double Q = (*bid_sys.current_local_solution)(dof_indices[0]);
            double Ve = (*bid_sys.current_local_solution)(dof_indices[1]);

            std::ostringstream ss;
            ss << std::setw(4) << std::setfill('0') << dataID;
            std::string data_ID_str = ss.str();
            if (unipole == 1)
            {
                std::ofstream file(Data::path + "/uni1_node_" + data_ID_str + ".txt", std::ios::app);
                //file << "time Ve V Q at" << std::endl;
                file << time << " " << Ve << " " << V << " " << Q << " " << at << std::endl;
                file.close();
            }
            else if (unipole == 2)
            {
                std::ofstream file(Data::path + "/uni2_node_" + data_ID_str + ".txt", std::ios::app);
                //file << "time Ve V Q at" << std::endl;
                file << time << " " << Ve << " " << V << " " << Q << " " << at << std::endl;
                file.close();
            }
            else
            {
                std::ofstream file(Data::path + "/at_node_" + data_ID_str + ".txt", std::ios::app);
                //file << "time Ve V Q at" << std::endl;
                file << time << " " << Ve << " " << V << " " << Q << " " << at << std::endl;
                file.close();
            }
        }
    }
    std::cout << "write_to_file done" << std::endl;

}


void extract_endocardial_mesh(libMesh::BoundaryMesh& endo, libMesh::MeshBase& mesh, int bdID, int subdomainID)
{
//    endo.set_n_partitions() = mesh.n_partitions();
//    std::map<dof_id_type, dof_id_type> node_id_map;
//    std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> side_id_map;
//
//    MeshBase::const_element_iterator el = mesh.active_subdomain_elements_begin();
//    const MeshBase::const_element_iterator end_el = mesh.elements_end();
//
//    this->_find_id_maps(requested_boundary_ids, 0, &node_id_map, 0, &side_id_map, subdomains_relative_to);

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

    BeatIt::TimeData datatime;
    datatime.setup(data);
    datatime.print();

    // Create a mesh with user-defined dimension on the default MPI
    // communicator.
    ParallelMesh mesh(init.comm());

    std::string meshname = data("mesh/input_mesh_name", "NONE");
    if (meshname == "NONE")
    {
        // allow us to use higher-order approximation.
        int elX = data("mesh/elX", 15);
        int elY = data("mesh/elY", 5);
        int elZ = data("mesh/elZ", 4);
        double maxX = data("mesh/maxX", 2.0);
        double maxY = data("mesh/maxY", 0.7);
        double maxZ = data("mesh/maxZ", 0.3);
        // No reason to use high-order geometric elements if we are
        // solving with low-order finite elements.
        if (elZ > 0)
            MeshTools::Generation::build_cube(mesh, elX, elY, elZ, 0., maxX, 0., maxY, 0., maxZ, TET4);
        else
            MeshTools::Generation::build_cube(mesh, elX, elY, elZ, 0., maxX, 0., maxY, 0., maxZ, TRI3);

        std::cout << "Creating subdomains!" << std::endl;

        {
            double z_interface = data("mesh/z_interface", 10000.0);
            double r_interface = data("mesh/r_interface", 0.0);
            double y_interface = data("mesh/y_interface", 10000.0);
            std::cout << "z_interface: " << z_interface << std::endl;

            MeshBase::element_iterator el = mesh.elements_begin();
            const MeshBase::element_iterator end_el = mesh.elements_end();

            for (; el != end_el; ++el)
            {
                Elem * elem = *el;
                const Point cent = elem->centroid();
                // BATH
                if (cent(2) > z_interface || cent(1) > y_interface)
                {
                    elem->subdomain_id() = 2;
                    double h = elem->hmax();
                    if (cent(2) < r_interface)
                        elem->set_refinement_flag(libMesh::Elem::REFINE);
                }
                // TISSUE
                else
                {
                    elem->subdomain_id() = 1;
                    elem->set_refinement_flag(libMesh::Elem::REFINE);
                }

            }
        }
        mesh.get_boundary_info().regenerate_id_sets();
        std::cout << "Creating subdomains done!" << std::endl;

    }
    else
    {
        // We may need XDR support compiled in to read binary .xdr files
        int n_refinements = data("mesh/n_ref", 0);
        std::cout << "n_refs: " << n_refinements << std::endl;
        BeatIt::serial_mesh_partition(init.comm(), meshname, &mesh, n_refinements);

        double xscale = data("mesh/x_scale", 1);
        double yscale = data("mesh/y_scale", 1);
        double zscale = data("mesh/z_scale", 1);
        if (xscale != 1 || yscale != 1 || zscale != 1)
            libMesh::MeshTools::Modification::scale(mesh, xscale, yscale, zscale);
    }

    MeshRefinement refinement(mesh);
    refinement.refine_elements();
    std::cout << "Mesh done!" << std::endl;
    mesh.print_info();

//    auto* node = mesh.query_node_ptr(1);
//    std::cout << mesh.comm().rank() << ": node 1 proc id: " << node->processor_id() << std::endl;
    //mesh.get_boundary_info().print_info(std::cout);
    // Create an equation systems object.
    std::cout << "Create equation system ..." << std::endl;
    EquationSystems es(mesh);




    //////////////////
    // CREATING FIBERS
    //////////////////
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


    auto& grad1 = poisson1.M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1 + "_gradient").solution;
    // through thickness
    auto& grad2 = poisson2.M_equationSystems.get_system<libMesh::ExplicitSystem>(pois2 + "_gradient").solution;
    // short direction
    auto& grad3 = poisson3.M_equationSystems.get_system<libMesh::ExplicitSystem>(pois3 + "_gradient").solution;

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    const libMesh::DofMap & dof_map = es.get_system<libMesh::ExplicitSystem>(pois1 + "_gradient").get_dof_map();

    std::vector<libMesh::dof_id_type> dof_indices;

    for (; el != end_el; ++el)
    {
        libMesh::Elem * elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        const auto blockID = elem->subdomain_id();

        if (blockID == 1)
        {

            double sx = (*grad2)(dof_indices[0]);
            double sy = (*grad2)(dof_indices[1]);
            double sz = (*grad2)(dof_indices[2]);

            double xfx = (*grad3)(dof_indices[0]);
            double xfy = (*grad3)(dof_indices[1]);
            double xfz = (*grad3)(dof_indices[2]);

            double fx = xfy * sz - xfz * sy;
            double fy = xfz * sx - xfx * sz;
            double fz = xfx * sy - xfy * sx;

            grad1->set(dof_indices[0], fx);
            grad1->set(dof_indices[1], fy);
            grad1->set(dof_indices[2], fz);
        }
    }
    auto& phi = poisson2.get_P0_solution();

    BeatIt::Util::normalize(*grad1, 1.0, 0.0, 0.0);
    BeatIt::Util::normalize(*grad2, 0.0, 1.0, 0.0);
    BeatIt::Util::normalize(*grad3, 0.0, 0.0, 1.0);

    poisson1.save_exo("poisson1.exo");

    ///////////////////
    ///////////////////
    ///////////////////
    ///////////////////
    ///////////////////
    // Constructor
    std::cout << "Create bidomain with bath..." << std::endl;
    std::string model = data("model", "monowave");
    BeatIt::ElectroSolver* solver = nullptr;

    bool solve_on_endo = data("endo_solve", false);
    BoundaryMesh *  endo_mesh = nullptr;
    EquationSystems * endo_es = nullptr;
    if(solve_on_endo)
    {
        std::cout << "CREATING ENDO MESH" << std::endl;
        endo_mesh = new BoundaryMesh(init.comm());
        std::set<boundary_id_type> endoIDs;
        endoIDs.insert(1);
        std::set<unsigned short> subdomains;
        subdomains.insert(1);
        //mesh.boundary_info->print_info(std::cout);
        //mesh.boundary_info->print_summary(std::cout);
        std::cout << "Sync endo" << std::endl;
        mesh.boundary_info->sync(endoIDs, *endo_mesh, subdomains);
        //endo.prepare_for_use();

        endo_mesh->print_info(std::cout);
        std::cout << "CREATING NEW EQUATION SYSTEM" << std::endl;
        endo_es = new EquationSystems(*endo_mesh);

        std::cout << "CREATING ENDOCARDIAL SOLVER" << std::endl;
        solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, *endo_es);
    }
    else
    {
        std::cout << "CREATING VOLUMETRIC SOLVER" << std::endl;
        BeatIt::ElectroSolver* solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, es);
    }
    std::string section = data("section", "monowave");
    std::cout << "Calling setup..." << std::endl;
    solver->setup(data, model);
    std::cout << "Calling init ..." << std::endl;
    solver->init(0.0);

    // set fibers
    // Fibers
    auto& fibers = solver->M_equationSystems.get_system<libMesh::ExplicitSystem>("fibers").solution;
    // Sheets
    auto& sheets = solver->M_equationSystems.get_system<libMesh::ExplicitSystem>("sheets").solution;
    // XFibers
    auto& xfibers = solver->M_equationSystems.get_system<libMesh::ExplicitSystem>("xfibers").solution;

    // Set fibers;
    if(solve_on_endo)
    {
        std::cout << "CREATING ENDOCARDIAL FIBERS" << std::endl;

        libMesh::MeshBase::const_element_iterator el_endo = endo_mesh->local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el_endo = endo_mesh->local_elements_end();

        const auto& endo_fibers_dof_map = solver->M_equationSystems.get_system<libMesh::ExplicitSystem>("fibers").get_dof_map();
        const auto& grad1_dof_map = poisson1.M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1 + "_gradient").get_dof_map();
        std::vector<unsigned int> fibers_dofs;
        std::vector<unsigned int> grad1_dofs;
        for( ; el_endo != end_el_endo; el_endo++)
        {
            const libMesh::Elem * elem = *el_endo;
            endo_fibers_dof_map.dof_indices(elem, fibers_dofs);

            const  libMesh::Elem * interior_parent = elem->interior_parent();
            grad1_dof_map.dof_indices(interior_parent, grad1_dofs);

            double fx = (*grad1)(grad1_dofs[0]);
            double fy = (*grad1)(grad1_dofs[1]);
            double fz = (*grad1)(grad1_dofs[2]);
            fibers->set(fibers_dofs[0], fx);
            fibers->set(fibers_dofs[1], fy);
            fibers->set(fibers_dofs[2], fz);

            double xfx = (*grad3)(grad1_dofs[0]);
            double xfy = (*grad3)(grad1_dofs[1]);
            double xfz = (*grad3)(grad1_dofs[2]);
            xfibers->set(fibers_dofs[0], xfx);
            xfibers->set(fibers_dofs[1], xfy);
            xfibers->set(fibers_dofs[2], xfz);

            double sx = (*grad2)(grad1_dofs[0]);
            double sy = (*grad2)(grad1_dofs[1]);
            double sz = (*grad2)(grad1_dofs[2]);
            sheets->set(fibers_dofs[0], sx);
            sheets->set(fibers_dofs[1], sy);
            sheets->set(fibers_dofs[2], sz);
        }
        std::cout << "DONE!" << std::endl;
        std::cout << "n sys: " << solver->M_equationSystems.n_systems() << std::endl;
        std::cout << "n sys: " << endo_es->n_systems() << std::endl;
        solver->M_equationSystems.print_info(std::cout);

    }
    else
    {
        *fibers = *grad1;
        *sheets = *grad2;
        *xfibers = *grad3;
    }

    int boundary_ic = data(model + "/boundary_ic", -1);
    if(boundary_ic > 0)
    {
        double boundary_value = data(model + "/boundary_value", 30.0);
        if(solve_on_endo)
        {
            std::map< dof_id_type, dof_id_type > node_id_map;
            std::map< dof_id_type, unsigned char > side_id_map;
            mesh.boundary_info->get_side_and_node_maps  (  *endo_mesh, node_id_map, side_id_map);

            libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
            libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

            libMesh::MeshBase::const_element_iterator el_endo = endo_mesh->active_local_elements_begin();
            const libMesh::MeshBase::const_element_iterator end_el_endo = endo_mesh->active_local_elements_end();

            auto & wave_system = solver->M_equationSystems.get_system<TransientLinearImplicitSystem>("wave");
            const libMesh::DofMap & dof_map = wave_system.get_dof_map();
            std::vector < libMesh::dof_id_type > v_dof_indices;

            for ( ; el != end_el; el++)
            {
                const libMesh::Elem * interior_parent = *el;
                for (unsigned int side = 0; side < interior_parent->n_sides(); side++)
                {
                    const unsigned int boundary_id = mesh.boundary_info->boundary_id(interior_parent, side);
                    if(boundary_id == boundary_ic)
                    {
                        std::unique_ptr<const Elem> side_el = interior_parent->build_side_ptr(side);
                        for(int nn = 0; nn < side_el->n_nodes() ; ++nn )
                        {
                            unsigned int nodeID = side_el->node_id(nn);
                            auto it_n = node_id_map.find(nodeID);
                            if(it_n != node_id_map.end())
                            {
                                unsigned int endo_nodeID = it_n->second;
                                Node * endo_node = endo_mesh->node_ptr(endo_nodeID);
                                dof_map.dof_indices(endo_node, v_dof_indices);
                                wave_system.solution->set(v_dof_indices[0], boundary_value);

                            }
                        }
                    }
                }

            }
            wave_system.solution->close();


        }
        else
        {
            solver->set_potential_on_boundary(boundary_ic, boundary_value, 1);
        }
    }


    poisson1.deleteSystems();
    poisson2.deleteSystems();
    poisson3.deleteSystems();

    Data::path = solver->M_outputFolder;
    int patient = data("patient", -1);
    std::map<int, int> rank_uni1_map;
    std::map<int, int> rank_uni2_map;
    std::map<int, int> rank_at_map;
    EndoRankMap uni1_id_map;
    EndoRankMap uni2_id_map;
    EndoRankMap at_id_map;

    if (patient > 0)
    {
        // Create an AnalyticFunction object that we use to project the BC
        // This function just calls the function exact_solution via exact_solution_wrapper
        if (patient == 2)
        {
            std::cout << "Found Patient 2!" << std::endl;

            //libMesh::AnalyticFunction<> ic_function(ic_patient2_ext);
            libMesh::AnalyticFunction<> ic_function(ic_patient2_plane_wave);
            solver->setup_ic(ic_function);
            if(!endo_mesh)
            {
                endo_mesh = new BoundaryMesh(init.comm());
                std::set<boundary_id_type> endoIDs;
                endoIDs.insert(1);
                std::set<unsigned short> subdomains;
                subdomains.insert(1);
                //mesh.boundary_info->print_info(std::cout);
                //mesh.boundary_info->print_summary(std::cout);
                std::cout << "Sync endo" << std::endl;
                mesh.boundary_info->sync(endoIDs, *endo_mesh, subdomains);
                //endo.prepare_for_use();
                endo_mesh->print_info(std::cout);
            }
            // Create map to access point from boundary mesh to point of the whole mesh
            std::map<dof_id_type, dof_id_type> reverse_node_id_map;
            std::cout << "Creating reverse map" <<std::endl;

            libMesh::MeshBase::const_node_iterator node_endo = endo_mesh->local_nodes_begin();
            const libMesh::MeshBase::const_node_iterator end_node_endo = endo_mesh->local_nodes_end();
            bool stop = false;
            for(; node_endo != end_node_endo; ++node_endo)
            {
                bool found = false;
                libMesh::MeshBase::const_node_iterator node_interior = mesh.nodes_begin();
                const libMesh::MeshBase::const_node_iterator end_node_interior = mesh.nodes_end();
                for(; node_interior != end_node_interior; ++node_interior)
                {
                    libMesh::Node * enode = *node_endo;
                    libMesh::Node * inode = *node_interior;
                    libMesh::Node nodo(*enode);
                    nodo -= *inode;
                    if(nodo.norm()< 1e-6)
                    {
                        found = true;
                        reverse_node_id_map[enode->id()] = inode->id();
                        break;
                    }
                }
                if(found == false)
                {
                    libMesh::Node * enode = *node_endo;
                    stop = true;
                    std::cout << endo_mesh->comm().rank() << ": Node: " << enode->id() << " was not found in the interior mesh" << std::endl;
                    throw std::runtime_error("error");
                }


            }
            if(stop) exit(-1);

            std::cout << "Creating reverse map done: "  << reverse_node_id_map.size() <<std::endl;
            //init_files(endo, uni1_id_map, reverse_node_id_map, "unip1.e", 1);
            //init_files(endo, uni2_id_map, reverse_node_id_map, "unip2.e", 2);
            //init_files(endo, at_id_map, reverse_node_id_map, "at_data.e", 0);
            project_unipolar_data(*endo_mesh, uni1_id_map, uni2_id_map, at_id_map, reverse_node_id_map);

        }
    }

    //return 0 ;
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
    //return 0;
    std::string system_mass = data(model + "/diffusion_mass", "mass");
    std::string iion_mass = data(model + "/reaction_mass", "lumped_mass");
    //bidomain.restart(importer, 1);
    std::cout << "Assembling matrices" << std::endl;
    solver->assemble_matrices(datatime.M_dt);


    if (export_data)
    {
        solver->save_parameters();
        solver->save_potential(save_iter, 0.0);
    }


    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;

    auto & V_sys = solver->M_equationSystems.get_system<libMesh::TransientLinearImplicitSystem>("wave");
    auto & bid_sys = solver->M_equationSystems.get_system<libMesh::TransientLinearImplicitSystem>(solver->model());
    auto & at_sys = solver->M_equationSystems.get_system<libMesh::ExplicitSystem>("activation_times");

    if (patient > 0)
    {
        init.comm().barrier();
        write_to_file(mesh, uni1_id_map, datatime.M_time, V_sys, bid_sys, at_sys, 1);
        write_to_file(mesh, uni2_id_map, datatime.M_time, V_sys, bid_sys, at_sys, 2);
        write_to_file(mesh, at_id_map, datatime.M_time, V_sys, bid_sys, at_sys, 0);
        init.comm().barrier();
    }

    std::cout << "Time loop starts:" << std::endl;
    //return 0;
    //export initial condition
    for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
    {
        datatime.advance();
        std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter << std::endl;
        solver->advance();
        //std::cout << "Reaction:" << datatime.M_time << std::endl;
        solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        //solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, "lumpedmass");

        //std::cout << "Diffusion:" << datatime.M_time << std::endl;
        solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        //solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, "mass");

        //std::cout << "at:" << datatime.M_time << std::endl;
        solver->update_activation_time(datatime.M_time, -5.0);
        //std::cout << "at done:" << datatime.M_time << std::endl;

        //++save_iter_ve;
        //bidomain.save_ve_timestep(save_iter_ve, datatime.M_time);

        if ( 0 == datatime.M_iter % static_cast<int>(1.0/datatime.M_dt) )
        {
            if (patient > 0)
            {
                init.comm().barrier();
                write_to_file(mesh, uni1_id_map, datatime.M_time, V_sys, bid_sys, at_sys, 1);
                write_to_file(mesh, uni2_id_map, datatime.M_time, V_sys, bid_sys, at_sys, 2);
                write_to_file(mesh, at_id_map, datatime.M_time, V_sys, bid_sys, at_sys, 0);
                init.comm().barrier();
            }
        }

        if (0 == datatime.M_iter % datatime.M_saveIter && export_data)
        {
            save_iter++;
            solver->save_potential(save_iter, datatime.M_time);

            //solver->save_potential(save_iter, datatime.M_time);
            //solver->save_exo_timestep(save_iter, datatime.M_time);

        }

    }
//    if (export_data)
        solver->save_activation_times(save_iter);
    delete solver;
    if(endo_es) delete endo_es;
    if(endo_mesh) delete endo_mesh;
    return 0;
}

