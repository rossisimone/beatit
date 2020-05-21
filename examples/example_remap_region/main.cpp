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

    std::string vtk_mesh_file = data("input_mesh", "NONE");
    if ("NONE" != vtk_mesh_file)
    {
        std::cout << "Reading mesh " << vtk_mesh_file << std::endl;
        mesh.read(vtk_mesh_file);
    }
    else
    {
        std::cout << "Set input_mesh in input file!" << std::endl;
        throw std::runtime_error("No mesh given");
    }

    std::string subregion_meshes = data("subregions", "NONE");
    bool remap_volume = true;
    if("NONE" == subregion_meshes) remap_volume = false;
    std::vector < std::string > mesh_files_vec;
    BeatIt::readList(subregion_meshes, mesh_files_vec);

    std::cout << subregion_meshes << std::endl;
    int num_subregions = mesh_files_vec.size();
    std::cout << "num subregions: " << num_subregions << std::endl;
    for (auto && m : mesh_files_vec)
        std::cout << m << std::endl;
    double tolerance = data("tolerance", 1e-2);

    int default_id = data("default_id", 666);

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    if(remap_volume)
    {
        std::cout << "Reading meshes and creating point locator ...  " << std::endl;
        std::vector < std::unique_ptr<libMesh::Mesh> > region_meshes_ptr(num_subregions);
        std::vector < std::unique_ptr<libMesh::PointLocatorTree> > locator_ptr(num_subregions);
        for (int i = 0; i < num_subregions; ++i)
        {
            std::cout << "Reading mesh " << mesh_files_vec[i] << ", " << i << std::endl;
            region_meshes_ptr[i].reset(new libMesh::Mesh(init.comm()));
            region_meshes_ptr[i]->read(mesh_files_vec[i]);
            region_meshes_ptr[i]->print_info();
            locator_ptr[i].reset(new libMesh::PointLocatorTree(*region_meshes_ptr[i]));
            locator_ptr[i]->set_close_to_point_tol(tolerance);
        }

        int missing_elements = 0;
        std::set<int> missing_elem;

        // Loop over each element of the mesh to remap
        // and look for the centroid of the element in one
        // of the subregion meshes
        // By default, the blockID is set to:  default_id (666)
        // In the meantime, count the number of missed elements
        // storing the corresponding ID in: missing_elem
        std::cout << "Looping over elements: " << std::endl;
        int c = 0;
        for (; el != end_el; ++el)
        {
            c++;
            Elem *elem = *el;

            auto centroid = elem->centroid();
            // Set element ID to the number of the beast
            elem->subdomain_id() = default_id;
            // Let's see if we can find the centroid in one of the meshes
            // Otherwise the ID will stay the number of the beast
            bool found = false;

            // Let's loop over the number of subregion meshes and look for the centroid
            // As soon as we find it, we can pass to the next element
            for (int i = 0; i < num_subregions; ++i)
            {

                const Elem *subregion_element = (*locator_ptr[i])(centroid);
                if (subregion_element)
                {
                    elem->subdomain_id() = i + 1;
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                int elID = elem->id();
                missing_elements++;
                missing_elem.insert(elID);
            }

        }

        std::cout << "Missed " << missing_elem.size() << " elements." << std::endl;

        std::cout << "Export first remapped mesh" << std::endl;
        ExodusII_IO(mesh).write("remapped_step_1.e");

        int num_missed_els = 10;
        int loop = 0;
        while (num_missed_els >= 1)
        {
            loop++;
            el = mesh.active_local_elements_begin();

            num_missed_els = 0;
            for (; el != end_el; ++el)
            {
                Elem *elem = *el;
                auto subID = elem->subdomain_id();
                if (default_id == subID)
                {
                    std::vector<int> ids;
                    for (unsigned int side = 0; side < elem->n_sides(); side++)
                    {
                        Elem *neighbor = elem->neighbor_ptr(side);
                        if (neighbor)
                            ids.push_back(neighbor->subdomain_id());
                    }
                    int min = *std::min_element(ids.begin(), ids.end());

                    elem->subdomain_id() = min;
                    if (default_id == elem->subdomain_id() || 0 == elem->subdomain_id())
                        num_missed_els++;
                }
            }
            std::cout << "Missed: " << c << " elements, still missing " << num_missed_els << std::endl;
            if (loop >= 10)
                break;
        }

        std::cout << "Export first remapped mesh" << std::endl;
        ExodusII_IO(mesh).write("remapped_step_2.e");

        // Fix boundary mesh for embryo meshes
        bool fix_boundary = data("fix", false);
        if (fix_boundary)
        {
            std::set < libMesh::dof_id_type > boundary_node_set;

            el = mesh.active_local_elements_begin();
            for (; el != end_el; ++el)
            {
                Elem *elem = *el;
                auto subID = elem->subdomain_id();
                if (1 == subID)
                {
                    for (int s = 0; s < elem->n_sides(); ++s)
                    {
                        if (nullptr == elem->neighbor_ptr(s))
                        {
                            auto side_elem_ptr = elem->build_side_ptr(s);
                            for (int n = 0; n < side_elem_ptr->n_nodes(); ++n)
                            {
                                boundary_node_set.insert(side_elem_ptr->node_ptr(n)->id());
                            }
                        }
                    }
                }
            }
            std::cout << "Boundary nodes on 1: " << boundary_node_set.size() << std::endl;

            el = mesh.active_local_elements_begin();
            for (; el != end_el; ++el)
            {
                Elem *elem = *el;
                auto subID = elem->subdomain_id();
                if (2 == subID)
                {
                    for (int n = 0; n < elem->n_nodes(); ++n)
                    {
                        auto it = boundary_node_set.find(elem->node_ptr(n)->id());
                        if (it != boundary_node_set.end())
                        {
                            std::cout << "changing element id on boundary" << std::endl;
                            elem->subdomain_id() = 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    // set bc flags
    bool set_bc_flags = data("bc", false);
    if(set_bc_flags)
    {
        std::string bc_meshes = data("bcs", "NONE");
        std::string bcs_ids = data("bcs_ids", "");
        std::vector < std::string > bc_files_vec;
        BeatIt::readList(bc_meshes, bc_files_vec);
        std::vector < unsigned int > bc_ids_vec;
        BeatIt::readList(bcs_ids, bc_ids_vec);
        int num_bc_subregions = bc_files_vec.size();
        int num_ids = bc_ids_vec.size();
        double bc_tolerance = data("bc_tolerance", 1e-2);

        std::cout << "Reading bc meshes and creating point locator ...  " << std::endl;
        std::vector < std::unique_ptr<libMesh::Mesh> > bc_region_meshes_ptr(num_bc_subregions);
        std::vector < std::unique_ptr<libMesh::PointLocatorTree> > bc_locator_ptr(num_bc_subregions);
        for (int i = 0; i < num_bc_subregions; ++i)
        {
            std::cout << "Reading bc mesh " << bc_files_vec[i] << ", " << i << std::endl;
            bc_region_meshes_ptr[i].reset(new libMesh::Mesh(init.comm()));
            bc_region_meshes_ptr[i]->read(bc_files_vec[i]);
            bc_region_meshes_ptr[i]->print_info();
            bc_locator_ptr[i].reset(new libMesh::PointLocatorTree(*bc_region_meshes_ptr[i]));
            bc_locator_ptr[i]->set_close_to_point_tol(bc_tolerance);
        }

        el = mesh.active_local_elements_begin();
        int bc_missing_elements = 0;
        std::set<int> bc_missing_elem;

        std::cout << "Looping over elements: " << std::endl;

        int c = 0;
        for (; el != end_el; ++el)
        {
            Elem * elem = *el;
            for (int s = 0; s < elem->n_sides(); ++s)
            {
                if (nullptr == elem->neighbor_ptr(s))
                {
                    auto centroid = elem->side_ptr(s)->centroid();
                    bool found = false;

                    // Let's loop over the number of subregion meshes and look for the centroid
                    // As soon as we find it, we can pass to the next element
                    for (int i = 0; i < num_bc_subregions; ++i)
                    {

                        const Elem* subregion_element = (*bc_locator_ptr[i])(centroid);
                        if (subregion_element)
                        {
                            unsigned int sideset_id = i+1;
                            if(num_ids == num_bc_subregions) sideset_id = bc_ids_vec[i];
                            mesh.boundary_info->add_side (elem, s, sideset_id);
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                    {
                        mesh.boundary_info->add_side (elem, s, default_id);
                        int elID = elem->id();
                        bc_missing_elements++;
                        bc_missing_elem.insert(elID);
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

