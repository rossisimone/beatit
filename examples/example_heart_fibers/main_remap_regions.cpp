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
    libMesh::Mesh heart_mesh(init.comm());
    libMesh::Mesh regions_mesh(init.comm());

    // Read meshes
    // first take name from the input file
    std::string heart_meshfile = data("heart_mesh", "NONE");
    std::cout << "Heart mesh file: " << heart_meshfile << std::endl;
    std::string regions_meshfile = data("regions_mesh", "NONE");
    std::cout << "Regions mesh file: " << regions_meshfile << std::endl;

    // Read the input meshes.
    std::cout << "Reading Heart Mesh" << std::endl;
    heart_mesh.read (&heart_meshfile[0]);
    std::cout << "Reading Regions Mesh" << std::endl;
    regions_mesh.read (&regions_meshfile[0]);

    // For each element in the heart mesh
    // we look for the centroid
    // and we try to locate it in the regions mesh

	// REGIONS MAP
	//  1 - MV (Mitral Valve)
	//  2 - TV (Tricuspid Valve)
	//  3, 4, 15 - PV (Pulmonic Valve)
	//  5 - PA (Pulmonary Arteries)
	//  6, 7, 8 - AV (Aortic Valve)
	//  9 - A (Aorta)
	// 10, 11, 12 - RPM (Right Papillary Muscles)
	// 13, 14 - LPM (Left Papillary Muscles)
	// 16 - LV (Left Ventricle)
	// 17 - RV (Right Ventricle)
	// 18, 19 - VC (Vena Cava)
	// 20 - TVR (Tricuscpid Valve Ring)
	// 21 - RAA (Right Atrial Appendage)
	// 22 - SAN (Sino-Atrial Node)
	// 23 - CT (Crista Terminalis)
	// 24 - RAPW (Right Atrial Posterior Wall)
	// 25 - RAAW (Right Atrial Anterior Wall)
	// 26 - RAR (Right Atrial Roof)
	// 27, 28, 31, 37 - PVs (Pulmonary Veins)
	// 29 - MVR (Mitral Valve Ring)
	// 30 - LAA (Left Atrial Appendage)
	// 32 - LA (Left Atrium)
	// 33, 35 - Left and Right Carina
	// 34, 36 - Left and Right Antra
    int laIDs[] = {27, 28, 29, 30, 31,32, 33, 34, 35, 36, 37};
    int raIDs[] = {18, 19, 20, 21, 22, 23, 24, 25, 26};
    int vIDs[]  = {16, 17};
    // HEART MAP
    // 1 - LA
    // 2 - RA
    // 3 - LV+RV
    // 40 - A
    // 41, 42, 43 - AV
    // 50 - PA
    // 51, 52, 53 - PV
    // 60 - MV
    // 61, 62 - LPM
    // 70 - TV
    // 71, 72, 73 - RPM

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

    libMesh::MeshBase::const_element_iterator el = heart_mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = heart_mesh.active_local_elements_end();
    std::cout << "Looping over heart mesh elements ... " << std::endl;
    int missing_elements = 0;
    std::set<int > missing_elem;

    for( ; el != end_el; ++el)
    {
        Elem * heart_elem = *el;

        auto centroid = heart_elem->centroid();

        int blockID = heart_elem->subdomain_id();
        std::unique_ptr< std::set<subdomain_id_type> > subdomains;

        bool remap = true;
        switch(blockID)
        {
        	case 1:
        	{
        		subdomains.reset(new  std::set<subdomain_id_type>(laIDs, laIDs+11));
        		break;
        	}
        	case 2:
        	{
        		subdomains.reset(new std::set<subdomain_id_type>(raIDs, raIDs+9));
        		break;
        	}
        	case 3:
        	{
        		subdomains.reset(new std::set<subdomain_id_type>(vIDs, vIDs+2));
        		break;
        	}
        	default:
        	{
        		remap = false;
        		break;
        	}

        }

        if(remap)
        {
			const Elem* regions_element = locator(centroid,subdomains.get());
			if(!regions_element)
			{
				int elID = heart_elem->id();
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
				heart_elem->subdomain_id() = regions_element->subdomain_id()+100;
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
			Elem * missed = heart_mesh.elem_ptr(missedID);
			bool unmissed = false;
			for (unsigned int side=0; side< missed->n_sides(); side++)
			{
				Elem * neighbor = missed->neighbor_ptr(side);
				bool found = false;
				if(neighbor)
				{
					int subdID = neighbor->subdomain_id();
					if(subdID > 100)
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
    ExodusII_IO exporter(heart_mesh);
    exporter.write("regioned_heart.e");

    std::cout << "Good luck!" << std::endl;
    return 0;
}





