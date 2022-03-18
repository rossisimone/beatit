/*
 * SetNumericVectorValues.cpp
 *
 *  Created on: May 14, 2018
 *      Author: srossi
 */
#include "Util/SetNumericVectorValues.hpp"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/system.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
 #include "libmesh/boundary_info.h"

namespace BeatIt
{

namespace SetValues
{

    void set_scalar_solution_on_boundary(libMesh::EquationSystems &es, libMesh::System & system, unsigned int boundID, double value, int subdomain)
    {
        std::cout << "* SetValues::  : WARNING:  set_potential_on_boundary works only for TET4" << std::flush;
        libMesh::MeshBase & mesh = es.get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        const libMesh::DofMap & dof_map = system.get_dof_map();
        std::vector < libMesh::dof_id_type > dof_indices;

        libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

        if(subdomain >= 0 )
        {
            el = mesh.active_local_subdomain_elements_begin(subdomain);
            end_el = mesh.active_local_subdomain_elements_end(subdomain);
        }

        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;

            {
                dof_map.dof_indices(elem, dof_indices);
                for (unsigned int side = 0; side < elem->n_sides(); side++)
                {
                    if (elem->neighbor_ptr(side) == libmesh_nullptr)
                    {
                        unsigned int n_boundary_ids=mesh.boundary_info->n_boundary_ids(elem,side);
                        std::vector<short int> boundary_ids_vec(n_boundary_ids);
                        mesh.boundary_info->boundary_ids(elem,side, boundary_ids_vec);
                        const unsigned int boundary_id = boundary_ids_vec[0];
                        if (boundary_id == boundID)
                        {
                            system.solution->set(dof_indices[0], value);
                            system.solution->set(dof_indices[1], value);
                            system.solution->set(dof_indices[2], value);
                        }
                    }
                }
            }
        }
        std::cout << "done" << std::endl;
        system.solution->close();
        std::cout << "done" << std::endl;
    }


} // SetValues

}// BeatIt



