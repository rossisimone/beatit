/*
 * main_remap_regions.cpp
 *
 *  Created on: Jul 18, 2019
 *      Author: srossi
 */


#include "PoissonSolver/Poisson.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
#include <iomanip>
#include "Util/GenerateFibers.hpp"
#include "libmesh/exodusII_io.h"
#include <libmesh/point_locator_tree.h>
#include <fstream>
#include <sstream>
#include "libmesh/xdr_io.h"
#include "libmesh/enum_xdr_mode.h"

using namespace libMesh;

const int N = 9;
void evaluate(double f[], double s[], double n[],
              double phi[], double dphi[][3], int blockID);

void cross(double v1[], double v2[], double cross_product[])
{
	cross_product[0] = v1[1] * v2[2] - v1[2] * v2[1];
	cross_product[1] = v1[2] * v2[0] - v1[0] * v2[2];
	cross_product[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return;
}


void cross(double v1[], int n, double v2[], int m, double dphi[][3], double cross_product[])
{
    auto normalize = [](double& x, double& y, double& z,
                        double X, double Y, double Z)
    {
        double norm = std::sqrt( x * x + y * y + z * z);
        if(norm >= 1e-12 )
        {
            x /= norm;
            y /= norm;
            z /= norm;
        }
        else
        {
            x = X;
            y = Y;
            z = Z;
        }
    };

	v1[0] = dphi[n][0];
	v1[1] = dphi[n][1];
	v1[2] = dphi[n][2];
	normalize(v1[0], v1[1], v1[2], 1.0, 0.0, 0.0);

	v2[0] = dphi[m][0];
	v2[1] = dphi[m][1];
	v2[2] = dphi[m][2];
	normalize(v2[0], v2[1], v2[2], 0.0, 1.0, 0.0);
	// compute f as the cross product
	cross(v1, v2, cross_product);
	normalize(cross_product[0], cross_product[1], cross_product[2], 0.0, 0.0, 1.0);
	return;
}


int main(int argc, char ** argv)
{
    // Initialize libmesh
    using namespace libMesh;
    typedef libMesh::ExodusII_IO EXOExporter;

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
    libMesh::Mesh final_mesh(init.comm());

    // Read meshes
    // first take name from the input file
    std::string heart_meshfile = data("heart_mesh", "NONE");
    std::string regions_meshfile = data("regions_mesh", "NONE");
    std::string final_meshfile = data("final_mesh", "NONE");

    std::cout << "Heart mesh file: " << heart_meshfile << std::endl;
    // Read the input meshes.
    std::cout << "Reading Heart Mesh" << std::endl;
    heart_mesh.read (&heart_meshfile[0]);
    std::cout << "Reading Region Mesh" << std::endl;
    regions_mesh.read (&regions_meshfile[0]);
    std::cout << "Reading Final Mesh" << std::endl;


    int num_refs = data("valves_refs", 0);
    if(num_refs > 0)
    {
        std::cout << "Refinement" << std::endl;
        libMesh::MeshRefinement mesh_refinement(heart_mesh);
        std::cout << "Refining Valves" << std::endl;
        libMesh::MeshBase::const_element_iterator el = heart_mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el = heart_mesh.active_local_elements_end();
        libMesh::PointLocatorTree point_locator(regions_mesh);
        std::set<libMesh::subdomain_id_type> allowed_subdomains;
        regions_mesh.subdomain_ids(allowed_subdomains);
        std::cout << "Setting Flags Valves" << std::endl;
        int c = 0;
        for (; el != end_el; ++el)
        {
            c++;
            libMesh::Elem * elem = *el;
            auto elID = elem->id();
            const Point centroid = elem->centroid();
            const Elem * subregion_element = point_locator(centroid);
            unsigned int blockID = 0;
            if(subregion_element) blockID = subregion_element->subdomain_id();
            else
            {
                //std::cout << "Subregion Element NOT FOUND! " << elID << " " << c << " performing linear search: " << std::endl;
                const Elem * subregion_element_linear_search = point_locator.perform_linear_search(centroid, &allowed_subdomains, true, 5e-1);
                if(subregion_element_linear_search) blockID = subregion_element_linear_search->subdomain_id();
                else
                {
                    //std::cout << "Subregion Element NOT FOUND! " << elID << " " << c << "performing fuzzy linear search:" <<   std::endl;
                    std::set<const libMesh::Elem*> subregion_element_fuzzy_linear_search = point_locator.perform_fuzzy_linear_search(centroid, &allowed_subdomains, 5e-1);
                    if(subregion_element_fuzzy_linear_search.size() > 0)
                    {
                        double distance = 1e9;
                        for(auto && found_elem : subregion_element_fuzzy_linear_search)
                        {
                            Point ptest(centroid);
                            ptest -= found_elem->centroid();
                            double new_distance = ptest.norm();
                            if( new_distance < distance )
                                blockID = found_elem->subdomain_id();
                        }
                    }
                    else
                    {
                        std::cout << "Subregion Element NOT FOUND! " << elID << " " << c << "NOT FOUND!" <<   std::endl;
                        throw std::runtime_error("subregion element not found");
                    }
                }
            }

            switch(blockID)
            {
                // AV
                case 41:
                case 42:
                case 43:
                // PV
                case 51:
                case 52:
                case 53:
                // MV
                case 60:
                // TV
                case 70:
                {
                    //std::cout << "Setting Refinement Flag in " << blockID << std::endl;
                    elem->set_refinement_flag(libMesh::Elem::REFINE);
                    break;
                }
                default:
                {
                    elem->set_refinement_flag(libMesh::Elem::DO_NOTHING);
                    break;
                }
            }

        }
        // refine mesh
        mesh_refinement.refine_elements();
        std::cout << "Prepare" << std::endl;
        heart_mesh.prepare_for_use(false);
        std::cout << "output" << std::endl;
        EXOExporter(heart_mesh).write("refined_valves.e");
        return 0;
    }

//    final_mesh.read (&final_meshfile[0]);
    //heart_mesh.read("full_heart_mesh.e");

    std::cout << "Changing sidesets names" << std::endl;
    // set boundary_ids names
    heart_mesh.get_boundary_info().sideset_name(103) = "right_ventricle_septum";
    heart_mesh.get_boundary_info().sideset_name(104) = "right_ventricle_septum";
    heart_mesh.get_boundary_info().sideset_name(105) = "aorta_top_surface";
    heart_mesh.get_boundary_info().sideset_name(106) = "pulmonary_artery_left_branch";
    heart_mesh.get_boundary_info().sideset_name(107) = "pulmonary_artery_right_branch";
    heart_mesh.get_boundary_info().sideset_name(440) = "aortic_valve_side";
    heart_mesh.get_boundary_info().sideset_name(441) = "aortic_valve_bottom";
    heart_mesh.get_boundary_info().sideset_name(442) = "aortic_valve_top";
    heart_mesh.get_boundary_info().sideset_name(540) = "pulmonic_valve_side";
    heart_mesh.get_boundary_info().sideset_name(541) = "pulmonic_valve_bottom";
    heart_mesh.get_boundary_info().sideset_name(542) = "pulmonic_valve_top";
    heart_mesh.get_boundary_info().sideset_name(640) = "mitral_valve_side";
    heart_mesh.get_boundary_info().sideset_name(642) = "mitral_valve_bottom";
    heart_mesh.get_boundary_info().sideset_name(740) = "tricuspid_valve_side";
    heart_mesh.get_boundary_info().sideset_name(742) = "tricuspid_valve_bottom";
    heart_mesh.get_boundary_info().sideset_name(743) = "tricuspid_valve_side";

    bool export_new_mesh_only = data("export_new_mesh_only", false);
    if(export_new_mesh_only)
    {
        heart_mesh.get_boundary_info().remove_id(440);
        heart_mesh.get_boundary_info().remove_id(540);
        heart_mesh.get_boundary_info().remove_id(640);
        heart_mesh.get_boundary_info().remove_id(743);
        heart_mesh.write ("heart_mesh_v5.xdr");
        return 0;
    }
    libMesh::EquationSystems es(heart_mesh);

    std::string pois[N];
    for(int i = 0; i < N; i++) pois[i] = "poisson" + std::to_string(i);

    typedef BeatIt::Poisson Poisson;
    typedef std::shared_ptr<Poisson> PoissonPtr;
    Poisson* pv[N];
    // 1: LEFT TRANSMURAL
    // 2: RIGHT TRANSMURAL
    // 3: PAPILLARY MUSCLES
    // 4: TRANS-INFLOWS
    // 5: TRANS-OUTFLOWS

    for(int i = 0; i < N; i++)
    {
		std::cout << "Solving Poisson " << i << std::endl;
		std::string pois1 = pois[i];
		BeatIt::Poisson poisson(es, pois1);
		std::cout << "Calling setup: ..." << std::flush;
		poisson.setup(data, pois1);
		std::cout << " Done!" << std::endl;
		std::cout << "Calling assemble system: ..." << std::flush;
		poisson.assemble_system();
		std::cout << " Done!" << std::endl;
		std::cout << "Calling solve system: ..." << std::flush;
		poisson.solve_system();
		std::cout << " Done!" << std::endl;
		std::cout << "Calling gradient: ..." << std::flush;
		poisson.compute_elemental_solution_gradient();
		std::cout << " Done!" << std::endl;
		pv[i] = &poisson;
    }

    // Normalize gradients
	std::cout << "Normalize: ..." << std::flush;
	for (int k = 0; k < N; ++k)
	{
		auto& grad = es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution;
		BeatIt::Util::normalize(*grad, 1.0, 0.0, 0.0);
	}
	std::cout << " Done!" << std::endl;

	// DEFINE FIBER SYSTEMS:
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& f_sys = es.add_system<FiberSystem>("fibers");
    f_sys.add_variable( "fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable( "fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable( "fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.init();
    auto& f_v = f_sys.solution;

    FiberSystem& s_sys = es.add_system<FiberSystem>("sheets");
    s_sys.add_variable( "sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.add_variable( "sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.add_variable( "sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.init();
    auto& s_v = s_sys.solution;

    FiberSystem& n_sys = es.add_system<FiberSystem>("xfibers");
    n_sys.add_variable( "xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.add_variable( "xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.add_variable( "xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.init();
    auto& n_v = n_sys.solution;

    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& l_sys = es.add_system<FiberSystem>("lambda");
    l_sys.add_variable( "lambda", libMesh::FIRST, libMesh::LAGRANGE);
    l_sys.init();
    auto& l_v = l_sys.solution;


    // Define auxiliary vectors used in the computations of the fibers
	double fx = 0,fy = 0, fz = 0;
	double f[3];
	double sx = 0, sy = 0, sz = 0;
    double s[3];
	double nx = 0, ny = 0, nz = 0;
    double n[3];
    double lambda;
	double phi[N];
	double dphi[N][3];


	// Loop over each element and based on the blockID
	// assign the fibers in some way
    libMesh::MeshBase::const_element_iterator el = heart_mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = heart_mesh.active_local_elements_end();
    std::cout << "\nGetting dof map:  ... " << std::flush;
    const libMesh::DofMap & dof_map = f_sys.get_dof_map();
    const libMesh::DofMap & dof_map_p = es.get_system<libMesh::ExplicitSystem>(pois[0]+"_P0").get_dof_map();
    const libMesh::DofMap & dof_map_phi = es.get_system<libMesh::ExplicitSystem>(pois[0]).get_dof_map();
	std::cout << " done \n " << std::flush;
    std::vector < libMesh::dof_id_type > dof_indices;
    std::vector < libMesh::dof_id_type > dof_indices_p;
    std::vector < libMesh::dof_id_type > dof_indices_phi;
    std::vector < libMesh::Real > phi_values;
    std::cout << "\nLooping over elements: ... \n" << std::flush;


    std::set<libMesh::subdomain_id_type> allowed_subdomains;
    regions_mesh.subdomain_ids(allowed_subdomains);
    libMesh::PointLocatorTree point_locator(regions_mesh);
    int c = 0;
    for (; el != end_el; ++el)
    {
        c++;
        libMesh::Elem * elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        dof_map_p.dof_indices(elem, dof_indices_p);
        dof_map_phi.dof_indices(elem, dof_indices_phi);
        auto elID = elem->id();
        const Point centroid = elem->centroid();
        const Elem * subregion_element = point_locator(centroid);
        unsigned int blockID = 0;
        if(subregion_element) blockID = subregion_element->subdomain_id();
        else
        {
            //std::cout << "Subregion Element NOT FOUND! " << elID << " " << c << " performing linear search: " << std::endl;
            const Elem * subregion_element_linear_search = point_locator.perform_linear_search(centroid, &allowed_subdomains, true, 5e-1);
            if(subregion_element_linear_search) blockID = subregion_element_linear_search->subdomain_id();
            else
            {
                //std::cout << "Subregion Element NOT FOUND! " << elID << " " << c << "performing fuzzy linear search:" <<   std::endl;
                std::set<const libMesh::Elem*> subregion_element_fuzzy_linear_search = point_locator.perform_fuzzy_linear_search(centroid, &allowed_subdomains, 5);
                if(subregion_element_fuzzy_linear_search.size() > 0)
                {
                    double distance = 1e9;
                    for(auto && found_elem : subregion_element_fuzzy_linear_search)
                    {
                        Point ptest(centroid);
                        ptest -= found_elem->centroid();
                        double new_distance = ptest.norm();
                        if( new_distance < distance )
                            blockID = found_elem->subdomain_id();
                    }
                }
                else
                {
                    std::cout << "Subregion Element NOT FOUND! " << elID << " " << c << "NOT FOUND!" <<   std::endl;
                    throw std::runtime_error("subregion element not found");
                }
            }
        }

    	for (int k = 0; k < N; ++k)
		{
            phi[k] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_P0").solution)(dof_indices_p[0]);
    		dphi[k][0] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[0]);
    		dphi[k][1] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[1]);
    		dphi[k][2] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[2]);
		}

    	evaluate(f, s, n, phi, dphi, blockID);

        // simplify subdomains
        bool left_side = true;
        std::string subdomain_name = "";
    	{

    	    unsigned int new_block_ID = 666;
    	    switch(blockID)
    	    {
                // LV
                case 103:
                {
                    subdomain_name = "left_ventricle";
                    new_block_ID = 1;
                    break;
                }
                // LV Pap
                case 63:
                case 64:
                {
                    subdomain_name = "left_ventricle_papillary_muscles";
                    new_block_ID = 10;
                    break;
                }
                // RV
                case 114:
                {
                    subdomain_name = "right_ventricle";
                    left_side = false;
                    new_block_ID = 2;
                    break;
                }
                // RV Pap
                case 74:
                case 75:
                case 76:
                {
                    subdomain_name = "right_ventricle_papillary_muscles";
                    left_side = false;
                    new_block_ID = 20;
                    break;
                }
                // LA
                case 202:
                case 212:
                case 216:
                {
                    subdomain_name = "left_atrium";
                    new_block_ID = 3;
                    break;
                }
                // PVs
                case 213:
                case 218:
                case 219:
                case 220:
                {
                    subdomain_name = "pulmonary_veins";
                    double threshold = 0.2;
                    if(blockID == 220 )
                        threshold = 0.15;
                    double potential = phi[5];
                    if(potential <= threshold)
                        new_block_ID = 30;
                    else
                        new_block_ID = 3;
                    break;
                }
                // RA
                case 301:
                case 305:
                case 308:
                case 309:
                case 317:
                case 321:
                {
                    subdomain_name = "right_atrium";
                    left_side = false;
                    new_block_ID = 4;
                    break;
                }
                // IVC + SVC
                case 306:
                {
                    subdomain_name = "vena_cava";
                    left_side = false;
                    double potential = phi[5];
                    // use phi 1 for s and phi 1 for n
                    if(potential >= 0.35)
                    {
                        new_block_ID = 4;
                    }
                    else
                    {
                        new_block_ID = 40;
                    }
                    break;
                }
                case 307:
                {
                    subdomain_name = "vena_cava";
                    left_side = false;
                    double potential = phi[5];
                    // use phi 1 for s and phi 1 for n
                    if(potential >= 0.3)
                    {
                        new_block_ID = 4;
                    }
                    else
                    {
                        new_block_ID = 40;
                    }
                    break;
                }
                // Aorta
                case 40:
                {
                    subdomain_name = "aorta";
                    new_block_ID = 5;
                    break;
                }
                // AV
                case 41:
                case 42:
                case 43:
                {
                    subdomain_name = "aortic_valve";
                    new_block_ID = 50;
                    break;
                }
                // Pulmonary Artery
                case 50:
                {
                    subdomain_name = "pulmonary_artery";
                    left_side = false;
                    new_block_ID = 6;
                    break;
                }
                // PV
                case 51:
                case 52:
                case 53:
                {
                    subdomain_name = "pulmonic_valve";
                    left_side = false;
                    new_block_ID = 60;
                    break;
                }
                // MV
                case 60:
                {
                    subdomain_name = "mitral_valve";
                    new_block_ID = 13;

                    // leaflet close to the AV
                    if (centroid(0)*0.16+centroid(1)*0.8-centroid(2)*0.58 < 111*0.16+101*0.8-74*0.58)
                    {
                        subdomain_name = "anterior_mitral_leaflet";
                        new_block_ID = 13;
                    }
                    else
                    {
                        subdomain_name = "posterior_mitral_leaflet";
                        new_block_ID = 14;
                    }
                    break;
                }
                // TV
                case 70:
                {
                    subdomain_name = "tricuspid_valve";
                    left_side = false;
                    new_block_ID = 24;
                    break;
                }
                // Skeleton
                case 112:
                {
                    subdomain_name = "left_cardiac_skeleton";
                    new_block_ID = 1234;
                    break;
                }
                case 310:
                {
                    subdomain_name = "right_cardiac_skeleton";
                    new_block_ID = 1235;
                    break;
                }
                // Mitral Anterior Marginal Chordae
                case 6504:
                case 6510:
                case 6511:
                case 6518:
                {
                    subdomain_name = "mitral_anterior_marginal_chordae";
                    new_block_ID = 661;
                    break;
                }
                // Mitral Posterior Marginal Chordae
                case 6502:
                case 6503:
                case 6505:
                case 6508:
                case 6509:
                case 6512:
                case 6513:
                case 6514:
                {
                    subdomain_name = "mitral_posterior_marginal_chordae";
                    new_block_ID = 662;
                    break;
                }
                // Mitral Anterior Strut Chordae
                case 6501:
                case 6517:
                {
                    subdomain_name = "mitral_anterior_strut_chordae";
                    new_block_ID = 663;
                    break;
                }
                // Mitral Anterior Strut Chordae
                case 6516:
                {
                    // if true we are on the marginal chordae
                    if( centroid(0)*0.52+centroid(1)*0.28+centroid(2)*0.81 < 136.476*0.52+79.8388*0.28+61.3883*0.81 )
                    {
                        subdomain_name = "mitral_anterior_marginal_chordae";
                        new_block_ID = 661;
                    }
                    // we are on the strut
                    else
                    {
                        subdomain_name = "mitral_anterior_strut_chordae";
                        new_block_ID = 663;
                    }
                    break;
                }
                // Mitral Posterior Basal Chordae
                case 6506:
                case 6507:
                case 6515:
                 {
                     subdomain_name = "mitral_posterior_basal_chordae";
                     new_block_ID = 665;
                     break;
                 }
                 // Mitral Antterior Basal Chordae NOT AVAILABLE
                // chordae
                 // Mitral Posterior Basal Chordae
                 case 7701:
                 case 7702:
                 case 7703:
                 case 7704:
                  {
                      subdomain_name = "septal_tricuspid_chordae";
                      new_block_ID = 667;
                      break;
                  }
                default:
                {
                    subdomain_name = "tricuspid_chordae";
                    break;
                }
    	    }
    	    elem->subdomain_id() = new_block_ID;
    	    heart_mesh.subdomain_name(new_block_ID) = subdomain_name;
    	}

    	for(auto && dof : dof_indices_phi)
    	{
    	    double val = 0.0;
            if( left_side )
                val = (*es.get_system<libMesh::ExplicitSystem>(pois[0]).solution)(dof);
            else
                val = (*es.get_system<libMesh::ExplicitSystem>(pois[1]).solution)(dof);
            l_v->set(dof, val);
    	}


        f_v->set(dof_indices[0], f[0]);
        f_v->set(dof_indices[1], f[1]);
        f_v->set(dof_indices[2], f[2]);

        s_v->set(dof_indices[0], s[0]);
        s_v->set(dof_indices[1], s[1]);
        s_v->set(dof_indices[2], s[2]);

        n_v->set(dof_indices[0], n[0]);
        n_v->set(dof_indices[1], n[1]);
        n_v->set(dof_indices[2], n[2]);
    }
    f_v->close();
    s_v->close();
    n_v->close();


    // Export the solution
    std::cout << "Calling exporter: ..."  << ". " << std::flush;
	EXOExporter exporter(heart_mesh);
	// If we want to export only the fiber field without
	// the solution of the poisson problems
	bool export_only_fibers = data("only_fibers", false);
	if(export_only_fibers)
	{
		std::vector<std::string> varname(10);
		varname[0] = "fibersx";
		varname[1] = "fibersy";
		varname[2] = "fibersz";
		varname[3] = "sheetsx";
		varname[4] = "sheetsy";
		varname[5] = "sheetsz";
		varname[6] = "xfibersx";
		varname[7] = "xfibersy";
		varname[8] = "xfibersz";
        varname[9] = "lambda";
		exporter.set_output_variables (varname);
	}
    heart_mesh.get_boundary_info().remove_id(440);
    heart_mesh.get_boundary_info().remove_id(540);
    heart_mesh.get_boundary_info().remove_id(640);
    heart_mesh.get_boundary_info().remove_id(743);

    std::string output = data("output", "new_heart_with_fibers.e");
	exporter.write_equation_systems(output, es);
	exporter.write_element_data(es);



    EXOExporter exporter2(heart_mesh);
    exporter2.write_equation_systems("complete_solution.e", es);
    exporter2.write_element_data(es);


    for(int i = 0; i < N; i++)
    {
        std::cout << "Deleting Poisson " << i << std::endl;
        const int N = 1;
        std::string pois1 = pois[i];
        es.delete_system(pois1);

    }

    std::cout << "ES2 " << std::endl;
    libMesh::EquationSystems es2(heart_mesh);
    std::cout << "F2 " << std::endl;
    FiberSystem& f_sys2 = es2.add_system<FiberSystem>("fibers");
    f_sys2.add_variable( "fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys2.add_variable( "fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys2.add_variable( "fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    std::cout << "S2 " << std::endl;
    FiberSystem& s_sys2 = es2.add_system<FiberSystem>("sheets");
    s_sys2.add_variable( "sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys2.add_variable( "sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys2.add_variable( "sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    std::cout << "N2 " << std::endl;
    FiberSystem& n_sys2 = es2.add_system<FiberSystem>("xfibers");
    n_sys2.add_variable( "xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys2.add_variable( "xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys2.add_variable( "xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    std::cout << "L2 " << std::endl;
    FiberSystem& l_sys2 = es2.add_system<FiberSystem>("lambda");
    l_sys2.add_variable( "lambda", libMesh::FIRST, libMesh::LAGRANGE);
    std::cout << "init " << std::endl;
    es2.init();
    auto& l_v2 = l_sys2.solution;
    auto& f_v2 = f_sys2.solution;
    auto& n_v2 = n_sys2.solution;
    auto& s_v2 = s_sys2.solution;
    std::cout << "copy " << std::endl;
    *l_sys2.solution =  *l_sys.solution;
    *f_sys2.solution =  *f_sys.solution;
    *s_sys2.solution =  *s_sys.solution;
    *n_sys2.solution =  *n_sys.solution;
    std::cout << "Exporting to XDR " << std::endl;
    std::string output_file = data("xdr_output", "new_heart_with_fibers.xdr");
    XdrMODE xdr_mode =  WRITE;
    const int write_mode = EquationSystems::WRITE_DATA;
    es2.write(output_file, xdr_mode, write_mode, /*partition_agnostic*/ true);
    return 0;
}

void evaluate(double f[], double s[], double n[],
              double phi[], double dphi[][3], int blockID)
{

    auto normalize = [](double& x, double& y, double& z,
                        double X, double Y, double Z)
    {
        double norm = std::sqrt( x * x + y * y + z * z);
        if(norm >= 1e-12 )
        {
            x /= norm;
            y /= norm;
            z /= norm;
        }
        else
        {
            x = X;
            y = Y;
            z = Z;
        }
    };

//    std::cout << "blockID: " << blockID << std::endl;
    switch(blockID)
    {
		// Left Ventricle Checked
		case 103:
		// right ventricle Checked
		case 114:
		{
			double potential = phi[0];
			double sx = dphi[0][0];
			double sy = dphi[0][1];
			double sz = dphi[0][2];
			if ( blockID == 114 || blockID == 74 || blockID == 75 || blockID == 76 )
			{
				potential = phi[1];
				sx = dphi[1][0];
				sy = dphi[1][1];
				sz = dphi[1][2];
			}

			double cx = 0.6599407598247158; //data (section+"/centerline_x", 0.0);
			double cy = -0.6674877611992615; //data (section+"/centerline_y", 0.0);
			double cz = -0.3448742990876162; // (section+"/centerline_z", 1.0);

			double epi_angle = -60;//data (section+"/epi_angle", -60.0);
			double endo_angle = 60;//data (section+"/endo_angle", 60.0);



			double cdot = cx * sx + cy * sy + cz * sz;
			double xfx = cx - cdot * sx;
			double xfy = cy - cdot * sy;
			double xfz = cz - cdot * sz;

			normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

			double fx = sy * xfz - sz * xfy;
			double fy = sz * xfx - sx * xfz;
			double fz = sx * xfy - sy * xfx;
			normalize(fx, fy, fz, 1.0, 0.0, 0.0);

			double teta1 = M_PI * epi_angle / 180.0;
			double teta2 = M_PI * endo_angle / 180.0;
			double m = (teta1 - teta2 );
			double q = teta2;
			double teta =  m * potential + q;

			double sa = std::sin (teta);
			double sa2 = 2.0 * std::sin (0.5 * teta) * std::sin (0.5 * teta);

			double W11 = 0.0; double W12 = -sz; double W13 =  sy;
			double W21 =  sz; double W22 = 0.0; double W23 = -sx;
			double W31 = -sy; double W32 =  sx; double W33 = 0.0;
			//
			double R11 = 1.0 + sa * W11 + sa2 * ( sx * sx - 1.0 );
			double R12 = 0.0 + sa * W12 + sa2 * ( sx * sy );
			double R13 = 0.0 + sa * W13 + sa2 * ( sx * sz );
			double R21 = 0.0 + sa * W21 + sa2 * ( sy * sx );
			double R22 = 1.0 + sa * W22 + sa2 * ( sy * sy - 1.0 );
			double R23 = 0.0 + sa * W23 + sa2 * ( sy * sz );
			double R31 = 0.0 + sa * W31 + sa2 * ( sz * sx );
			double R32 = 0.0 + sa * W32 + sa2 * ( sz * sy );
			double R33 = 1.0 + sa * W33 + sa2 * ( sz * sz - 1.0 );

			double f0x =   R11 * fx + R12 * fy + R13 * fz;
			double f0y =   R21 * fx + R22 * fy + R23 * fz;
			double f0z =   R31 * fx + R32 * fy + R33 * fz;
			normalize(f0x, f0y, f0z, 1.0, 0.0, 0.0);

			xfx = f0y * sz - f0z * sy;
			xfy = f0z * sx - f0x * sz;
			xfz = f0x * sy - f0y * sx;
			normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

			f[0] = f0x; f[1] = f0y; f[2] = f0z;
			s[0] = sx;  s[1] = sy;  s[2] = sz;
			n[0] = xfx; n[1] = xfy; n[2] = xfz;

			break;
		}
        // left atrial appendage Checked
        case 216:
        // right atria appendage Checked
        case 308:
        {
			double cx =0.6742238029343747; //data (section+"/centerline_x", 0.0);
			double cy = -0.7296511576681224; //data (section+"/centerline_y", 0.0);
			double cz =0.11415538388651877; // (section+"/centerline_z", 1.0);

            double sx = dphi[0][0];
            double sy = dphi[0][1];
            double sz = dphi[0][2];

            if(blockID == 308)
            {
                cx = 0.75314410973367;
                cy = -0.49830792973321275;
                cz = 0.42949174280593244;
                sx = dphi[1][0];
                sy = dphi[1][1];
                sz = dphi[1][2];
            }

            double cdot = cx * sx + cy * sy + cz * sz;
            double xfx = cx - cdot * sx;
            double xfy = cy - cdot * sy;
            double xfz = cz - cdot * sz;

            normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

            double fx = sy * xfz - sz * xfy;
            double fy = sz * xfx - sx * xfz;
            double fz = sx * xfy - sy * xfx;
            normalize(fx, fy, fz, 1.0, 0.0, 0.0);

            f[0] = fx;  f[1] = fy;  f[2] = fz;
            s[0] = sx;  s[1] = sy;  s[2] = sz;
            n[0] = xfx; n[1] = xfy; n[2] = xfz;
            break;
        }
        // Left Ventricle Ppapillary muscles Checked
        case 63:
        case 64:
        {
            double sx = dphi[0][0];
            double sy = dphi[0][1];
            double sz = dphi[0][2];

            double fx = -sy;
            double fy = sx;
            double fz = 0.0;
            normalize(fx, fy, fz, 1.0, 0.0, 0.0);

            double xfx = sy * fz - sz * fy;
            double xfy = sz * fx - sx * fz;
            double xfz = sx * fy - sy * fx;
            normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

            f[0] = fx;  f[1] = fy;  f[2] = fz;
            s[0] = sx;  s[1] = sy;  s[2] = sz;
            n[0] = xfx; n[1] = xfy; n[2] = xfz;
            break;
        }
        // Right Ventricle Ppapillary muscles Checked
        case 74:
        case 75:
        case 76:
        {
            double sx = dphi[1][0];
            double sy = dphi[1][1];
            double sz = dphi[1][2];

            double fx = -sy;
            double fy = sx;
            double fz = 0.0;
            normalize(fx, fy, fz, 1.0, 0.0, 0.0);

            double xfx = sy * fz - sz * fy;
            double xfy = sz * fx - sx * fz;
            double xfz = sx * fy - sy * fx;
            normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

            f[0] = fx;  f[1] = fy;  f[2] = fz;
            s[0] = sx;  s[1] = sy;  s[2] = sz;
            n[0] = xfx; n[1] = xfy; n[2] = xfz;
            break;
        }
        // Aortic Valve Checked
        case 41:
        case 42:
        case 43:
		// Mitral Valve Checked
		case 60:
        // Pulmonic Valve Checked
        case 51:
        case 52:
        case 53:
        // Tricuspid Valve Checked
        case 70:
		{
			// use phi2 for s and phi 7 for n
			cross(s, 2, n, 7, dphi, f);
			break;
		}
        // Aorta Checked
		case 40:
		// compute n as the cross product
		{
			// use phi0 for s and phi 3 for f
			cross(s, 0, f, 3, dphi, n);
			break;
		}
		// Pulmonary Artery Checked
		case 50:
		// compute n as the cross product
		{
			// use phi1 for s and phi 3 for f
			cross(s, 1, f, 3, dphi, n);
			break;
		}
		// Mitral Valve Ring Checked
		case 112:
        // Left atrial floor Checked
        case 212:
        // Tricuspid Valve Ring Checked
        case 310:
		{
			// use phi 0 for s and phi 3 for n
			cross(s, 0, n, 8, dphi, f);
			break;
		}
		// ICV Checked
		case 306:
        {
            // use phi 5 for f
            double potential = phi[5];
            // use phi 1 for s and phi 1 for n
            if(potential >= 0.35)
            {
                cross(s, 1, n, 4, dphi, f);
            }
            else
            {
                cross(s, 1, f, 4, dphi, n);
            }
            break;
        }
        // SCV Checked
		case 307:
		{
            // use phi 5 for f
			double potential = phi[5];
            // use phi 1 for s and phi 1 for n
			if(potential >= 0.3)
			{
				cross(s, 1, n, 4, dphi, f);
			}
			else
			{
				cross(s, 1, f, 4, dphi, n);
			}
			break;
		}
		// Papillary Muscles
//		case 61:
//		case 62:
//		case 71:
//		case 72:
//		case 73:
//		{
//			f[0] = dphi[2][0];
//			f[1] = dphi[2][1];
//			f[2] = dphi[2][2];
//			normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);
//
//			s[0] =-dphi[2][1];
//			s[1] = dphi[2][0];
//			s[2] = 0.0;
//			normalize(s[0], s[1], s[2], 0.0, 1.0, 0.0);
//
//			cross(f, s, n);
//			break;
//		}

		//  (LA)
		case 202:
		{
		    double potential = phi[4];
            double low_threshold = 0.4;
            double high_threshold = 0.7;
	        // Left + Right Carina (LA)
            if(potential <= low_threshold || potential >= high_threshold)
            {
                // use phi4 for s and phi 1 for n
                cross(s, 0, n, 6, dphi, f);
            }
            // use phi4 for s and phi 1 for n
            else
            {
                potential = phi[3];
                double threshold = 0.2565;
                if(potential <= threshold)
                {
                    cross(s, 0, n, 8, dphi, f);
                }
                else
                    cross(s, 0, n, 4, dphi, f);
            }
			break;
		}
    	// Pulmonary Veins Checked
		case 213:
		case 218:
		case 219:
		case 220:
		{
			double threshold = 0.2;
            if(blockID == 220 )
                threshold = 0.15;
			double potential = phi[5];
			//std::cout << "Potential: " << potential << std::endl;
            if(potential <= threshold)
            {
              //  std::cout << "cross(s, 0, f, 4, dphi, n): " << std::flush;
                cross(s, 0, f, 4, dphi, n);
//                std::cout << "done" << std::endl;
            }
            else
            {
//                std::cout << "cross(s, 0, f, 4, dphi, n): " << std::flush;
                cross(s, 0, n, 4, dphi, f);
//                std::cout << "done" << std::endl;
            }
			break;
		}
        // Eustachian Valve checked
        case 305:
        {
            // use phi 1 for s and phi 5 for n
            cross(s, 1, n, 5, dphi, f);
            break;
        }
        // RA anterior wall checked
        case 301:
        {
            // use phi 1 for s and phi 5 for n
            cross(s, 1, f, 6, dphi, n);
            break;
        }
        // RA superior wall checked
        case 321:
        {
            // use phi 1 for s and phi 5 for n
            cross(s, 1, n, 4, dphi, f);
            break;
        }
        // right atrium pw checked
        case 317:
		{
		    double potential = phi[1];
            double threshold = 0.5;
			if(potential > 0.5)
			cross(s, 1, n, 3, dphi, f);
			else
            cross(s, 1, f, 5, dphi, n);
			break;
		}
		// right atrium roof checked
		case 309:
		{
			cross(s, 1, f, 4, dphi, n);
			break;
		}
		//SAN
		case 122:
		// Crista Terminalis
        case 123:
		// right atrium anterior wall
		case 125:
		{
			cross(s, 1, n, 2, dphi, f);
			break;
		}
		// CHORDAE
        case 6501:
        case 6502:
        case 6503:
        case 6504:
        case 6505:
        case 6506:
        case 6507:
        case 6508:
        case 6509:
        case 6510:
        case 6511:
        case 6512:
        case 6513:
        case 6514:
        case 6515:
        case 6516:
        case 6517:
		case 6518:
        case 7701:
        case 7702:
        case 7703:
        case 7704:
        case 7705:
        case 7706:
        case 7707:
        case 7708:
        case 7709:
        case 7710:
        case 7711:
        case 7712:
        {
            // use phi 5 for f
            f[0] = dphi[5][0];
            f[1] = dphi[5][1];
            f[2] = dphi[5][2];
            normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);

            // Take any direction for s and n
            // since it does not matter in the chordae
            s[0] = f[1];
            s[1] = -f[0];
            s[2] = 0.0;
            normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);

            n[0] = f[1] * s[2] - f[2] * s[1];
            n[1] = f[2] * s[0] - f[0] * s[2];
            n[2] = f[0] * s[1] - f[1] * s[0];
            normalize(n[0], n[1], n[2], 0.0, 0.0, 1.0);
            break;
        }
        default:
        {
            std::cout << "You should not be here!!!!" << blockID <<  std::endl;
            s[0] = dphi[0][0];
            s[1] = dphi[0][1];
            s[2] = dphi[0][2];
            normalize(s[0], s[1], s[2], 1.0, 0.0, 0.0);

            f[0] =-dphi[0][1];
            f[1] = dphi[0][0];
            f[2] = 0.0;
            normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);

            n[0] = f[1] * s[2] - f[2] * s[1];
            n[1] = f[2] * s[0] - f[0] * s[2];
            n[2] = f[0] * s[1] - f[1] * s[0];
            normalize(n[0], n[1], n[2], 0.0, 0.0, 1.0);
            break;
        }
    }
}
