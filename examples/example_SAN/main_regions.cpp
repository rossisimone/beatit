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

int main(int argc, char ** argv)
{
    // Initialize libmesh
    using namespace libMesh;
    LibMeshInit init (argc, argv, MPI_COMM_WORLD);

    // Read input file
    GetPot commandLine ( argc, argv );
    std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
    GetPot data(datafile_name);

    // Create mesh:
    // one for the conforming mesh
    // one with the regions to be mapped
    libMesh::Mesh SAN_mesh(init.comm());
    libMesh::Mesh regions_mesh(init.comm());

    // Read meshes
    // first take name from the input file
    std::string SAN_meshfile = data("SAN_mesh", "NONE");
    std::cout << "Heart mesh file: " << SAN_meshfile << std::endl;
    std::string regions_meshfile = data("regions_mesh", "NONE");
    std::cout << "Regions mesh file: " << regions_meshfile << std::endl;

    // Read the input meshes.
    std::cout << "Reading Heart Mesh" << std::endl;
    SAN_mesh.read (&SAN_meshfile[0]);
    std::cout << "Reading Regions Mesh" << std::endl;
    regions_mesh.read (&regions_meshfile[0]);

    // For each element in the heart mesh
    // we look for the centroid
    // and we try to locate it in the regions mesh

	// REGIONS MAP
	//  1 - SAN
	//  2 - atrium

    // We only need to split the atria and the ventricles!!!!
    // Not touching the valves and the vessels.
    // We are going to do three loop
    // 1) Ventricles
    // 2) LA
    // 3) RA
    // Whish me luck

    // We need to create a point locator
    std::cout << "Creating Point Locator" << std::endl;
    libMesh::PointLocatorTree locator (regions_mesh);
    //locator.enable_out_of_mesh_mode ();
    double tolerance = data("tolerance", 1e-2);
    locator.set_close_to_point_tol(tolerance);

    libMesh::MeshBase::const_element_iterator el = SAN_mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = SAN_mesh.active_local_elements_end();
    std::cout << "Looping over heart mesh elements ... " << std::endl;
    int missing_elements = 0;
    std::set<int > missing_elem;

    for( ; el != end_el; ++el)
    {
        Elem * SAN_elem = *el;

        auto centroid = SAN_elem->centroid();

        int blockID = SAN_elem->subdomain_id();
        std::unique_ptr< std::set<subdomain_id_type> > subdomains;

        bool remap = true;

        if(remap)
        {
			const Elem* regions_element = locator(centroid,subdomains.get());
			if(!regions_element)
			{
				int elID = SAN_elem->id();
				std::cout << blockID << ":" << elID << ", " << std::flush;
				std::cout << "Centroid not found! " << std::endl;
				centroid.print();
				std::cout << std::endl;
				missing_elements++;
				missing_elem.insert(elID);
			}
			else
			{
				//std::cout << regions_element->subdomain_id() << ":" << regions_element->id() << ", " << std::flush;
				SAN_elem->subdomain_id() = regions_element->subdomain_id();
			}
	}
    }

    std::cout << "Fixing missed elements" << std::endl;
    int loops = 0;
    while(missing_elem.size() > 0)
    {
    	std::cout << "Loop: " << loops++ << ", over " << missing_elem.size() << " elements"<< std::endl;
		for(auto&& missedID : missing_elem)
		{
			Elem * missed = SAN_mesh.elem_ptr(missedID);
			bool unmissed = false;
			for (unsigned int side=0; side< missed->n_sides(); side++)
			{
				Elem * neighbor = missed->neighbor_ptr(side);
				bool found = false;
				if(neighbor)
				{
					int subdID = neighbor->subdomain_id();
					//if(subdID > 100)
					{
						missed->subdomain_id() = subdID;
						found = true;
					}
				}
				if(found)
				{
					unmissed = true;
					break;
				}
			}
			if(unmissed) missing_elem.erase(missedID);
		}
		if(loops > 10) break;
    }

    std::cout << "Missed " << missing_elements << " elements." << std::endl;
    std::cout << "Export to new mesh" << std::endl;
    ExodusII_IO exporter(SAN_mesh);
    exporter.write("regioned_SAN_mesh.e");

    std::cout << "Good luck!" << std::endl;
    return 0;
}





