/*
 * main_remap_regions.cpp
 *
 *  Created on: Jul 18, 2019
 *      Author: srossi
 */

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include <libmesh/point_locator_tree.h>
#include "libmesh/elem.h"
#include "libmesh/getpot.h"
#include <algorithm>
#include <tuple>
#include "Util/IO/io.hpp"
#include <libmesh/boundary_mesh.h>

int main(int argc, char ** argv)
{
    // Initialize libmesh
    using namespace libMesh;
    // Import input file
    GetPot data = BeatIt::readInputFile(argc, argv);

    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // Create mesh:
    // one for the conforming mesh
    // one with the regions to be mapped
    libMesh::Mesh mesh(init.comm());

    std::string mesh_file = data("input_mesh", "NONE");
    if ("NONE" != mesh_file)
    {
        std::cout << "Reading mesh " << mesh_file << std::endl;
        mesh.read(mesh_file);
    }
    else
    {
        std::cout << "Set input_mesh in input file!" << std::endl;
        throw std::runtime_error("No mesh given");
    }

    double tolerance = data("tolerance", 1e-3);
    int default_id = data("default_id", 666);

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // set bc flags
    //bool set_bc_flags = data("bc", false);
    if(true)
    {
        std::string bc_meshes = data("bcs", "NONE");
        std::string bcs_ids = data("bcs_ids", "");
        std::vector < std::string > bc_files_vec;
        BeatIt::readList(bc_meshes, bc_files_vec);
        for(auto && m : bc_files_vec) std::cout << m << std::endl;
        std::vector < unsigned int > bc_ids_vec;
        BeatIt::readList(bcs_ids, bc_ids_vec);
        for(auto && m : bc_ids_vec) std::cout << m << std::endl;
        int num_bc_subregions = bc_files_vec.size();
        int num_ids = bc_ids_vec.size();
        std::cout << "Num bc meshes: " << num_bc_subregions << ", num ids: " << num_ids << std::endl;
        double bc_tolerance = data("bc_tolerance", 1e-2);

        std::cout << "Reading bc meshes and creating point locator ...  " << std::endl;
        std::vector < std::unique_ptr<libMesh::Mesh> > bc_region_meshes_ptr(num_bc_subregions);
        //std::vector < std::unique_ptr<libMesh::PointLocatorTree> > bc_locator_ptr(num_bc_subregions);
        for (int i = 0; i < num_bc_subregions; ++i)
        {
            std::cout << "Reading bc mesh " << bc_files_vec[i] << ", " << i << std::endl;
            bc_region_meshes_ptr[i].reset(new libMesh::Mesh(init.comm()));
            bc_region_meshes_ptr[i]->read(bc_files_vec[i]);
            bc_region_meshes_ptr[i]->print_info();
            //bc_locator_ptr[i].reset(new libMesh::PointLocatorTree(*bc_region_meshes_ptr[i]));
            //bc_locator_ptr[i]->set_close_to_point_tol(bc_tolerance);
        }

        int bc_missing_elements = 0;
        std::set<int> bc_missing_elem;

        std::cout << "Looping over elements: " << std::endl;

        int c = 0;
        for (; el != end_el; ++el)
        {
            c++;
            if( c % 100000 == 0) std::cout << "Elem: " << c << std::endl;
            Elem * elem = *el;
            for (int s = 0; s < elem->n_sides(); ++s)
            {
                if (nullptr == elem->neighbor_ptr(s))
                {
                    auto side_elem_ptr  = elem->build_side_ptr(s);
                    auto node0 = side_elem_ptr->node_ptr(0);
                    auto node1 = side_elem_ptr->node_ptr(1);
                    auto node2 = side_elem_ptr->node_ptr(2);

                    bool found_elem = false;
                    int region = -1;
                    // Let's loop over the number of subregion meshes and look for the centroid
                    // As soon as we find it, we can pass to the next element
                    for (int i = 0; i < num_bc_subregions; ++i)
                    {
                        //std::cout << "i = " << i << " " << std::flush;
                        bool found_node = false;
                        libMesh::MeshBase::const_node_iterator node =  bc_region_meshes_ptr[i]->active_nodes_begin();
                        const libMesh::MeshBase::const_node_iterator end_nodes =  bc_region_meshes_ptr[i]->active_nodes_end();
                        for (; node != end_nodes; ++node)
                        {
                            libMesh::Node * bc_node = *node;
                            libMesh::Point v(*bc_node - *node0);
                            double distance = v.norm();
                            if(distance < tolerance)
                            {
                                found_node = true;
                                region = i;
                                //std::cout << "\n Found on mesh: " << i << std::endl;
                                break;
                            }
                        }
                        if(found_node)
                        {
                            libMesh::MeshBase::const_element_iterator bc_el = bc_region_meshes_ptr[i]->active_local_elements_begin();
                            const libMesh::MeshBase::const_element_iterator end_bc_el = bc_region_meshes_ptr[i]->active_local_elements_end();
                            for (; bc_el != end_bc_el; ++bc_el)
                            {
                                Elem * bc_elem = *bc_el;
                                auto bc_node0 = bc_elem->node_ptr(0);
                                auto bc_node1 = bc_elem->node_ptr(1);
                                auto bc_node2 = bc_elem->node_ptr(2);

                                libMesh::Point v0(*bc_node0 - *node0);
                                double distance0 = v0.norm();
                                libMesh::Point v1(*bc_node1 - *node0);
                                double distance1 = v1.norm();
                                libMesh::Point v2(*bc_node2 - *node0);
                                double distance2 = v2.norm();
                                if(distance0 < tolerance || distance1 < tolerance || distance2 < tolerance)
                                {
                                    v0  = *bc_node0 - *node1;
                                    distance0 = v0.norm();
                                    v1 = *bc_node1 - *node1;
                                    distance1 = v1.norm();
                                    v2 = *bc_node2 - *node1;
                                    distance2 = v2.norm();
                                    if(distance0 < tolerance || distance1 < tolerance || distance2 < tolerance)
                                    {
                                         v0  = *bc_node0 - *node2;
                                         distance0 = v0.norm();
                                         v1 = *bc_node1 - *node2;
                                         distance1 = v1.norm();
                                         v2 = *bc_node2 - *node2;
                                         distance2 = v2.norm();
                                         if(distance0 < tolerance || distance1 < tolerance || distance2 < tolerance)
                                         {
                                             found_elem = true;
                                             unsigned int sideset_id = i+1;
                                             if(num_ids == num_bc_subregions) sideset_id = bc_ids_vec[i];
                                             mesh.boundary_info->add_side (elem, s, sideset_id);
                                             break;
                                         }
                                    }
                                } // end distance check
                            } // end loop over elem bcs
                        } // end if found node
                        if(found_elem) break;
                    } // end loop over mesh subregions

//                        const Elem* subregion_element = (*bc_locator_ptr[i])(centroid);
//                        if (subregion_element)
//                        {
//                            unsigned int sideset_id = i+1;
//                            if(num_ids == num_bc_subregions) sideset_id = bc_ids_vec[i];
//                            mesh.boundary_info->add_side (elem, s, sideset_id);
//                            found = true;
//                            break;
//                        }
//                    }

                    if (!found_elem)
                    {
                        mesh.boundary_info->add_side (elem, s, default_id);
                        int elID = elem->id();
                        bc_missing_elements++;
                        bc_missing_elem.insert(elID);
                        std::cout << "\n Missed elem: " << bc_missing_elements << std::endl;
                    }

                }
            }
        }
        std::cout << "Missed " << bc_missing_elem.size() << " boundary elements." << std::endl;

//        // Fix boundaries
//        BoundaryMesh * boundary_mesh = new BoundaryMesh(init.comm());
//        mesh.boundary_info->sync(*boundary_mesh);
//        boundary_mesh->print_info(std::cout);
//
//        num_missed_els = 10;
//        loop = 0;
//        while (num_missed_els >= 1)
//        {
//            loop++;
//            el = boundary_mesh.active_local_elements_begin();
//            const libMesh::MeshBase::const_element_iterator b_end_el = boundary_mesh.active_local_elements_end();
//            num_missed_els = 0;
//            for (; el != b_end_el; ++el)
//            {
//                Elem * elem = *el;
//                auto subID = elem->subdomain_id();
//                if (default_id == subID)
//                {
//                    std::vector<int> ids;
//                    for (unsigned int side = 0; side < elem->n_sides(); side++)
//                    {
//                        Elem * neighbor = elem->neighbor_ptr(side);
//                        if (neighbor)
//                            ids.push_back(neighbor->subdomain_id());
//                    }
//                    int min = *std::min_element(ids.begin(), ids.end());
//
//                    elem->subdomain_id() = min;
//                    if (default_id == elem->subdomain_id() || 0 == elem->subdomain_id())
//                        num_missed_els++;
//                }
//            }
//            std::cout << "Missed: " << c << " elements, still missing " << num_missed_els << std::endl;
//            if (loop >= 10)
//                break;
//        }

        mesh.prepare_for_use();
    }
    std::cout << "Export to new mesh" << std::endl;
    std::string output_file = data("output", "remapped_mesh.e");
    ExodusII_IO exporter(mesh);
    exporter.write(output_file);

    std::cout << "Good luck!" << std::endl;
    return 0;
}

