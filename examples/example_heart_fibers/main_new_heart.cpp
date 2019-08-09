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

#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
#include <iomanip>
#include "Util/GenerateFibers.hpp"
#include "libmesh/exodusII_io.h"

#include <fstream>
#include <sstream>

using namespace libMesh;

const int N = 7;
void evaluate(double f[], double s[], double n[],
              double phi[], double dphi[][3], int blockID, Point& x);

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
    final_mesh.read (&final_meshfile[0]);
    //heart_mesh.read("full_heart_mesh.e");


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
		const int N = 1;
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

    // Define auxiliary vectors used in the computations of the fibers
	double fx = 0,fy = 0, fz = 0;
	double f[3];
	double sx = 0, sy = 0, sz = 0;
    double s[3];
	double nx = 0, ny = 0, nz = 0;
    double n[3];

	double phi[N];
	double dphi[N][3];


	// Loop over each element and based on the blockID
	// assign the fibers in some way
    libMesh::MeshBase::const_element_iterator el = heart_mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = heart_mesh.active_local_elements_end();
    std::cout << "\nGetting dof map:  ... " << std::flush;
    const libMesh::DofMap & dof_map = f_sys.get_dof_map();
    const libMesh::DofMap & dof_map_p = es.get_system<libMesh::ExplicitSystem>(pois[0]+"_P0").get_dof_map();
	std::cout << " done \n " << std::flush;
    std::vector < libMesh::dof_id_type > dof_indices;
    std::vector < libMesh::dof_id_type > dof_indices_p;
    std::cout << "\nLooping over elements: ... \n" << std::flush;

    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        dof_map_p.dof_indices(elem, dof_indices_p);
        auto elID = elem->id();
        Point centroid = elem->centroid();
        const auto blockID = regions_mesh.elem_ptr(elID)->subdomain_id();

    	for (int k = 0; k < N; ++k)
		{
            phi[k] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_P0").solution)(dof_indices_p[0]);
    		dphi[k][0] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[0]);
    		dphi[k][1] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[1]);
    		dphi[k][2] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[2]);
		}

    	evaluate(f, s, n, phi, dphi, blockID, centroid);
    	// change subdomain ID Pulmonary Veins and Vena Cava:
		if(blockID == 127 || blockID == 128 || blockID == 131 || blockID == 137
		|| blockID == 118 || blockID == 119)
		{
			double threshold = 0.75;
			double potential = phi[5];
			if(blockID == 127) threshold = 0.15;
			if(blockID == 128) threshold = 0.2;
			int new_blockID = 10;
			if(blockID == 118 || blockID == 119)
			{
				threshold = 0.6;
				potential = phi[2];
				new_blockID = 20;
			}

			bool change_blockID = false;
			if(blockID == 127 || blockID == 128)
				if(potential < threshold)
					change_blockID = true;

			if(blockID == 131 || blockID == 137)
				if(potential > threshold)
					change_blockID = true;

			if(blockID == 118 || blockID == 119)
				if(potential > threshold)
					change_blockID = true;

			if(change_blockID)
			{
				final_mesh.elem_ptr(elID)->subdomain_id() = new_blockID;
			}
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

    // Export the solution
	typedef libMesh::ExodusII_IO EXOExporter;
    std::cout << "Calling exporter: ..."  << ". " << std::flush;
	EXOExporter exporter(final_mesh);
	// If we want to export only the fiber field without
	// the solution of the poisson problems
	bool export_only_fibers = data("only_fibers", false);
	if(export_only_fibers)
	{
		std::vector<std::string> varname(9);
		varname[0] = "fibersx";
		varname[1] = "fibersy";
		varname[2] = "fibersz";
		varname[3] = "sheetsx";
		varname[4] = "sheetsy";
		varname[5] = "sheetsz";
		varname[6] = "xfibersx";
		varname[7] = "xfibersy";
		varname[8] = "xfibersz";
		exporter.set_output_variables (varname);
	}
	exporter.write_equation_systems("new_heart_with_fibers.e", es);
	exporter.write_element_data(es);


    return 0;
}

void evaluate(double f[], double s[], double n[],
              double phi[], double dphi[][3], int blockID, Point& x)
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

    switch(blockID)
    {
		// Left Ventricle
		case 116:
		// right ventricle
		case 117:
		{
			double potential = phi[0];
			double sx = dphi[0][0];
			double sy = dphi[0][1];
			double sz = dphi[0][2];
			if (blockID == 117)
			{
				potential = phi[1];
				sx = dphi[1][0];
				sy = dphi[1][1];
				sz = dphi[1][2];
			}

			double cx = 0.53393899516939; //data (section+"/centerline_x", 0.0);
			double cy = -0.674105753868089; //data (section+"/centerline_y", 0.0);
			double cz = -0.510382779920559; // (section+"/centerline_z", 1.0);

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
        // left atrial appendage
        case 130:
        // right atria appendage
        case 121:
        {
			double cx =0.34020660984721; //data (section+"/centerline_x", 0.0);
			double cy = -0.823146463920235; //data (section+"/centerline_y", 0.0);
			double cz =0.454631016926783; // (section+"/centerline_z", 1.0);

            double sx = dphi[0][0];
            double sy = dphi[0][1];
            double sz = dphi[0][2];

            if(blockID == 121)
            {
                cx = -0.726767523370148;
                cy = 0.381478658634116;
                cz = -0.571211869608061;
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
        // Aortic Valve
        case 41:
        case 42:
        case 43:
		// Mitral Valve
		case 60:
		{
			// use phi4 for s and phi 0 for n
			cross(s, 4, n, 0, dphi, f);
			break;
		}

        // Pulmonic Valve
		case 51:
		case 52:
		case 53:
		// Tricuspid Valve
		case 70:
		// compute f as the cross product
        {
			// use phi4 for s and phi 1 for n
			cross(s, 4, n, 1, dphi, f);
			break;
        }
        // Aorta
		case 40:
		// compute n as the cross product
		{
			// use phi4 for s and phi 1 for n
			cross(s, 0, f, 3, dphi, n);
			break;
		}
		// Pulmonary Artery
		case 50:
		// compute n as the cross product
		{
			// use phi4 for s and phi 1 for n
			cross(s, 1, f, 3, dphi, n);
			break;
		}
		// Mitral Valve Ring
		case 129:
		{
			// use phi4 for s and phi 1 for n
			cross(s, 0, n, 2, dphi, f);
			break;
		}
		// Tricuspid Valve Ring
		case 120:
		{
			// use phi4 for s and phi 1 for n
			cross(s, 1, n, 2, dphi, f);
			break;
		}

		// ICV + SCV
		case 118:
		case 119:
		{
			double potential = phi[2];
			if(phi[2] <= 0.6)
			{
				cross(s, 1, n, 2, dphi, f);
			}
			else
			{
				cross(s, 1, f, 2, dphi, n);
			}
			break;
		}
		// Papillary Muscles
		case 61:
		case 62:
		case 71:
		case 72:
		case 73:
		{
			f[0] = dphi[2][0];
			f[1] = dphi[2][1];
			f[2] = dphi[2][2];
			normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);

			s[0] =-dphi[2][1];
			s[1] = dphi[2][0];
			s[2] = 0.0;
			normalize(s[0], s[1], s[2], 0.0, 1.0, 0.0);

			cross(f, s, n);
			break;
		}
		// Left + Right Carina (LA)
		case 133:
		case 135:
		{
			// use phi4 for s and phi 1 for n
			cross(s, 0, n, 4, dphi, f);
			break;
		}
		default:
		{
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
		// left atrium
        case 132:
		// Left + Right Antra (LA)
		case 134:
		case 136:
        {
			// use phi4 for s and phi 0 for n
			cross(s, 0, n, 5, dphi, f);
			break;
        }
    	// Pulmonary Veins
		case 127:
		case 128:
		case 131:
		case 137:
		{
			double threshold = 0.75;
			double potential = phi[5];
			if(blockID == 131 || blockID == 137)
			{
				if(phi[5] <= threshold)
				{
					cross(s, 0, n, 5, dphi, f);
				}
				else
				{
					cross(s, 0, f, 5, dphi, n);
				}
				break;
			}

			if(blockID == 127)
			{
				threshold = 0.15;
			}
			if(blockID == 128)
			{
				threshold = 0.2;
			}
			if(phi[5] >= threshold)
			{
				cross(s, 0, n, 5, dphi, f);
			}
			else
			{
				cross(s, 0, f, 5, dphi, n);
			}

			break;
		}
        // right atrium pw
        case 124:
		{
//			if(phi[1] > 0.5)
			cross(s, 1, n, 6, dphi, f);
//			else
//			cross(s, 1, f, 6, dphi, n);
			break;
		}
		// right atrium roof
		case 126:
		{
			cross(s, 1, n, 6, dphi, f);
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

    }
}
