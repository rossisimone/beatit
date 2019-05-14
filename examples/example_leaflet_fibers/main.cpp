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

#include <fstream>
#include <sstream>


const int N = 4;


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

      std::ifstream right_atrium("right_atrium0.csv");
      int objID, elID, elID2;
      char v, v2;
      std::string line;
      if(right_atrium.is_open())
      {
          getline (right_atrium,line);

            while ( getline (right_atrium,line) )
            {
                std::istringstream ss(line);
                //line will have
                ss >> objID >> v >> elID;
//                    if(fn == 10) elID = elID2;
                //std::cout << elID << std::endl;
//                  if(elID == 16896) trick = true;
//                if(fn == 10 && trick == true) elID = objID;
                auto * elem = mesh.query_elem(elID-1);
                elem->subdomain_id() = 21;
            }
      }
      right_atrium.close();

      libMesh::EquationSystems es(mesh);

      std::string pois[N];
      for(int i = 0; i < N; i++) pois[i] = "poisson" + std::to_string(i);

      typedef BeatIt::Poisson Poisson;
      typedef std::shared_ptr<Poisson> PoissonPtr;
      Poisson* pv[N];





      for(int i = 0; i < N; i++)
      {
		std::cout << "Solving Poisson " << i << std::endl;
		//const int N = 1;
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
//        exporter.set_output_variables (varname);
		exporter.write_equation_systems("pa.exo", es);
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
        // Pulmonary Valve Leaflets
        case 51:
        case 52:
        case 53:
        // Aortic Valve Leaflets
        case 41:
        case 42:
        case 43:
		//tricuspid valve
		case 70:
               // mitral valve 
               case 60:
        {
            // s
            double v1_x = dphi[0][0];
            double v1_y = dphi[0][1];
            double v1_z = dphi[0][2];

            // n
            double v2_x = dphi[1][0];
            double v2_y = dphi[1][1];
            double v2_z = dphi[1][2];

            double cross_x = v1_y * v2_z - v1_z * v2_y;
            double cross_y = v1_z * v2_x - v1_x * v2_z;
            double cross_z = v1_x * v2_y - v1_y * v2_x;
            normalize(cross_x, cross_y, cross_z, 0.0, 0.0, 1.0);

            s[0] = v1_x;  s[1] = v1_y;  s[2] = v1_z;
            n[0] = v2_x;  n[1] = v2_y;  n[2] = v2_z;
            f[0] = cross_x; f[1] = cross_y; f[2] = cross_z;
            break;
        }
        // Pulmonary Trunk
        case 50:
        // Aorta
        case 40:
        {
            s[0] = dphi[3][0];
            s[1] = dphi[3][1];
            s[2] = dphi[3][2];
            normalize(s[0], s[1], s[2], 1.0, 0.0, 0.0);

            f[0] = dphi[2][0];
            f[1] = dphi[2][1];
            f[2] = dphi[2][2];

            normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);

            n[0] = f[1] * s[2] - f[2] * s[1];
            n[1] = f[2] * s[0] - f[0] * s[2];
            n[2] = f[0] * s[1] - f[1] * s[0];
            normalize(n[0], n[1], n[2], 0.0, 0.0, 1.0);
            break;
        }
        // RV papillaary muscles
        case 221:
        case 222:
        case 223:
        // LV papillaary muscles
        case 64:
        case 65:
        case 66:
        case 67:
	// RV Chordae
        case 7301:
        case 7302:
        case 7303:
        case 7304:
        case 7305:
        case 7306:
        case 7307:
        case 7308:
        case 7309:
        case 7310:
        case 7311:
        case 7312:
        // LV Chordae
        case 6301:
        case 6302:
        case 6303:
        case 6304:
        case 6305:
        case 6306:
        case 6307:
        case 6308:
        case 6309:
        case 6310:
        case 6311:
        case 6312:
        case 6313:
        case 6314:
        case 6315:
        case 6316:
        {
            // f
            double f_x = dphi[1][0];
            double f_y = dphi[1][1];
            double f_z = dphi[1][2];

            // s
            double s_x = -dphi[1][1];
            double s_y = dphi[1][0];
            double s_z = 0;

            double cross_x = f_y * s_z - f_z * s_y;
            double cross_y = f_z * s_x - f_x * s_z;
            double cross_z = f_x * s_y - f_y * s_x;
            normalize(cross_x, cross_y, cross_z, 0.0, 0.0, 1.0);

            s[0] = s_x;  s[1] = s_y;  s[2] = s_z;
            f[0] = f_x;  f[1] = f_y;  f[2] = f_z;
            n[0] = cross_x; n[1] = cross_y; n[2] = cross_z;
            break;

        	break;
        }
        default:
        {
            throw std::runtime_error("Wrong!!!");
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


