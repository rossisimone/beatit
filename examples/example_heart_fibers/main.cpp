/*
 ============================================================================

 .______    _______     ___   .___________.    __  .___________.
 |   _  \  |   ____|   /   \  |           |   |  | |           |
 |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
 |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |     
 |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |     
 |______/  |_______/__/     \__\  |__|        |__|     |__|     
 
 BeatIt - code for cardiovascular simulations
 Copyright (C) 2016 Simone Rossi

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ============================================================================
 */

/*
 * main.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */
// Basic include files needed for the mesh functionality.
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

const int N = 7;


void evaluate(double f[], double s[], double n[],
              double phi[], double dphi[][3], int blockID);

int main (int argc, char ** argv)
{
     // Bring in everything from the libMesh namespace

      using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);


      // Use the MeshTools::Generation mesh generator to create a uniform
      // 3D grid
      // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
      // the number of elements on each side is read from the input file
      GetPot commandLine ( argc, argv );
      std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
      GetPot data(datafile_name);

      // allow us to use higher-order approximation.
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh(init.comm());

      // We may need XDR support compiled in to read binary .xdr files
      std::string meshfile = data("mesh/input_mesh_name", "duke_heart_v5.e");

      // Read the input mesh.
      mesh.read (&meshfile[0]);


      libMesh::EquationSystems es(mesh);

      std::string pois[N];
      for(int i = 0; i < N; i++) pois[i] = "poisson" + std::to_string(i);

      typedef BeatIt::Poisson Poisson;
      typedef std::shared_ptr<Poisson> PoissonPtr;
      Poisson* pv[N];

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

		std::cout << "Get Gradients: ... " << std::flush;
		//auto& grad1 = es.get_system<libMesh::ExplicitSystem>(pois[1]+"_gradient").solution;
		std::cout << " Done!" << std::endl;


		std::cout << "Normalize: ..." << std::flush;
		for (int k = 0; k < N; ++k)
		{
			auto& grad = es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution;
			BeatIt::Util::normalize(*grad, 1.0, 0.0, 0.0);
		}
		//BeatIt::Util::normalize(*grad1, 0.0, 1.0, 0.0);
		std::cout << " Done!" << std::endl;




		auto& grad = es.get_system<libMesh::ExplicitSystem>(pois[0]+"_gradient").solution;
        //first = grad->first_local_index();
        //last = grad->last_local_index();

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





    	double fx = 0;
    	double fy = 0;
    	double fz = 0;
    	double f[3];


    	double sx = 0;
    	double sy = 0;
    	double sz = 0;
        double s[3];

    	double nx = 0;
    	double ny = 0;
    	double nz = 0;
        double n[3];

    	double phi[N];
    	double dphi[N][3];



        libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
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
            const auto blockID = elem->subdomain_id();

        	for (int k = 0; k < N; ++k)
			{
                phi[k] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_P0").solution)(dof_indices_p[0]);
        		dphi[k][0] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[0]);
        		dphi[k][1] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[1]);
        		dphi[k][2] = (*es.get_system<libMesh::ExplicitSystem>(pois[k]+"_gradient").solution)(dof_indices[2]);
			}

        	evaluate(f, s, n, phi, dphi, blockID);


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


		std::cout << "Calling exporter: ..."  << ". " << std::flush;

		typedef libMesh::ExodusII_IO EXOExporter;
		EXOExporter exporter(mesh);
		exporter.write_equation_systems("HeartPois0/poisson0.exo", es);
		exporter.write_element_data(es);

        std::cout << " Done!" << std::endl;



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

    switch(blockID)
    {
        // Left Ventricle
        case 3:
        // right ventricle
        case 14:
        {
            double potential = phi[0];
            double sx = dphi[0][0];
            double sy = dphi[0][1];
            double sz = dphi[0][2];
            if (blockID == 14)
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
        case 16:
        {
            double cx =-0.523696611217322; //data (section+"/centerline_x", 0.0);
            double cy = 0.848062778126582; //data (section+"/centerline_y", 0.0);
            double cz =-0.0808169769028613; // (section+"/centerline_z", 1.0);

            double sx = dphi[0][0];
            double sy = dphi[0][1];
            double sz = dphi[0][2];

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
        // upper left pulmonary vein
        case 13:
        // mitral valve ring
        case 12:
        // lower left pulmonary vein
        case 18:
        // lower right pulmonary vein
        case 19:
        // upper right pulmonary vein
        case 20:
        {
            double v1_x = dphi[0][0];
            double v2_x = dphi[2][0];
            double v1_y = dphi[0][1];
            double v2_y = dphi[2][1];
            double v1_z = dphi[0][2];
            double v2_z = dphi[2][2];

            double cdot = v1_x * v2_x + v1_y * v2_y + v1_z * v2_z;
            v1_x = v1_x - cdot * v2_x;
            v1_y = v1_y - cdot * v2_y;
            v1_z = v1_z - cdot * v2_z;

            normalize(v1_x, v1_y, v1_z, 0.0, 0.0, 1.0);
            double cross_x = v1_y * v2_z - v1_z * v2_y;
            double cross_y = v1_z * v2_x - v1_x * v2_z;
            double cross_z = v1_x * v2_y - v1_y * v2_x;
            normalize(cross_x, cross_y, cross_z, 0.0, 0.0, 1.0);

            f[0] = v1_x;  f[1] = v1_y;  f[2] = v1_z;
            s[0] = v2_x;  s[1] = v2_y;  s[2] = v2_z;
            n[0] = cross_x; n[1] = cross_y; n[2] = cross_z;
            break;
        }
        // Left Atrium
        case 2:
        {
           double v1_x = dphi[3][0];
           double v2_x = dphi[2][0];
           double v1_y = dphi[3][1];
           double v2_y = dphi[2][1];
           double v1_z = dphi[3][2];
           double v2_z = dphi[2][2];

           double cdot = v1_x * v2_x + v1_y * v2_y + v1_z * v2_z;

           v1_x = v1_x - cdot * v2_x;
           v1_y = v1_y - cdot * v2_y;
           v1_z = v1_z - cdot * v2_z;
           normalize(v1_x, v1_y, v1_z, 0.0, 0.0, 1.0);
           double cross_x = v1_y * v2_z - v1_z * v2_y;
           double cross_y = v1_z * v2_x - v1_x * v2_z;
           double cross_z = v1_x * v2_y - v1_y * v2_x;
           normalize(cross_x, cross_y, cross_z, 0.0, 0.0, 1.0);

           f[0] = v1_x;  f[1] = v1_y;  f[2] = v1_z;
           s[0] = v2_x;  s[1] = v2_y;  s[2] = v2_z;
           n[0] = cross_x; n[1] = cross_y; n[2] = cross_z;
           break;
        }
        // Right Atrium
        case 1:
        // inter atrial band
        case 4:
        // eustachian valve
        case 5:
        // inferior vena cava
        case 6:
        // superior vena cava
        case 7:
        // right atria appendage
        case 8:
        // right atria roof muscle
        case 9:
        // tricuspid valve ring
        case 10:
        // aorta
        case 11:
        // pulmonary artery
        case 15:
        // cresta terminalis
        case 17:
        default:
        {
            s[0] = dphi[0][0];
            s[1] = dphi[0][1];
            s[2] = dphi[0][2];

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


