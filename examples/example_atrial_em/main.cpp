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
#include "Electromechanics/Electromechanics.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "Util/Timer.hpp"
#include "libmesh/dof_map.h"
#include "libmesh/exodusII_io.h"

#include <iomanip>

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
      std::string datafile_name = commandLine.follow ( "em.beat", 2, "-i", "--input" );
      GetPot data(datafile_name);
      GetPot data_fibers("data.beat");

      // allow us to use higher-order approximation.
      // Create a mesh, with dimension to be overridden later, on the
      // default MPI communicator.
      libMesh::Mesh mesh(init.comm());

      // We may need XDR support compiled in to read binary .xdr files
      std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");
      bool buildFibers = data("monodomain/build_fibers", true);
      bool readMesh = data("mesh/read_mesh", false);
      // Read the input mesh.
      bool do_restart = data("monodomain/restart/restart", false);
      ExodusII_IO importer(mesh) ;
      if(do_restart)
      {
        std::string restart_file = data("monodomain/restart/restart_file", "NONE");
        if(restart_file != "NONE")
        {
            importer.read(restart_file);
            mesh.prepare_for_use();
        }
        else
        {
            do_restart = false;
        }

      }

      if(do_restart == false)
      {
          if(meshfile != "Pippo.e")
          {
              std::cout << "EXXXXXXX" << std::endl;
              mesh.read (&meshfile[0]);
              std::cout << "EXXXXXXX" << std::endl;
          }
          else
          {
              std::cout << "Mesh: " << meshfile << std::flush;
              std::cout << "Either you did not specify the mesh in mesh/input_mesh_name" << std::endl;
              std::cout << "Or you want to do a simulation of Pippo,  which (I agree) it would be nice." << std::endl;
              throw std::runtime_error("WRONG MESH!");
          }
      }
      std::cout << "\nConstructing es: ..." << std::flush;
      libMesh::EquationSystems es(mesh);
      std::cout << "\ninit es: ..." << std::flush;

      //es.init();

      std::cout << "\nConstructing em: ..." << std::flush;
      BeatIt::Electromechanics em(es, "electromechanics");
      // Setup the equation systems
      std::cout << "\nSetting it up em: ..." << std::flush;
      em.setup(data, "monodomain", "elasticity", "activation");
      std::cout << "\nInitializing em: ..." << std::flush;
      em.init(0.0);

      if(do_restart)
      {
          int restart_step = data("monodomain/restart/step", 2);
          //monodomain.restart(importer, restart_step);
          std::cout << "ERROR: restart not coded yet" << std::endl;
          throw std::runtime_error("error");
          std::string fibers_file = data("monodomain/restart/fibers_file", "NONE");
          int fibers_step = data("monodomain/restart/fibers_step", 1);
//        if(fibers_file != "NONE")
//        {
//            ExodusII_IO fiber_importer(mesh) ;
//            fiber_importer.read(fibers_file);
//          monodomain.readFibers(fiber_importer, fibers_step);
//          buildFibers = false;
//        }

      }

      BeatIt::Timer timer2;

      timer2.start();

      if(buildFibers)
        {
              std::string pois1 = "poisson1";
              BeatIt::Poisson poisson1(es, pois1);

              std::cout << "Calling setup: ..." << std::flush;
              poisson1.setup(data_fibers, "poisson1");
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
              auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").solution;
              auto first = grad1->first_local_index();
              auto last = grad1->last_local_index();
              std::cout << "First: " << first << ", Last: " << last << std::endl;


              std::string pois2 = "poisson2";
              BeatIt::Poisson poisson2(es, pois2);

              std::cout << "Calling setup: ..." << std::flush;
              poisson2.setup(data_fibers, "poisson2");
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

              std::string pois4 = "poisson4";
              BeatIt::Poisson poisson4(es, pois4);

              std::cout << "Calling setup: ..." << std::flush;
              poisson4.setup(data_fibers, "poisson4");
              std::cout << " Done!" << std::endl;
              std::cout << "Calling assemble system: ..." << std::flush;
              poisson4.assemble_system();
              std::cout << " Done!" << std::endl;
              std::cout << "Calling solve system: ..." << std::flush;
              poisson4.solve_system();
              std::cout << " Done!" << std::endl;
              std::cout << "Calling gradient: ..." << std::flush;
              poisson4.compute_elemental_solution_gradient();
              std::cout << " Done!" << std::endl;


              std::string pois3 = "poisson3";
              BeatIt::Poisson poisson3(es, pois3);

              std::cout << "Calling setup: ..." << std::flush;
              poisson3.setup(data_fibers, "poisson3");

              std::cout << "\n\nCreating Fibers: ... \n" << std::flush;

        //      auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").solution;
              auto& grad2 = poisson2. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois2+"_gradient").solution;
              auto& grad3 = poisson3. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois3+"_gradient").solution;
              auto& grad4 = poisson4. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois4+"_gradient").solution;
              first = grad1->first_local_index();
              last = grad1->last_local_index();
              std::cout << "First: " << first << ", Last: " << last << std::endl;

        //      auto first = grad1->first_local_index();
        //      auto last = grad1->last_local_index();

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

            double cx = data_fibers ("poisson3/centerline_x", 0.0);
            double cy = data_fibers ("poisson3/centerline_y", 0.0);
            double cz = data_fibers ("poisson3/centerline_z", 1.0);

              std::cout << "\n\nNormalizing gradients: ... \n" << std::flush;

              for(int i = first; i < last; )
              {
                  int j = i;
                  double v1_x = (*grad1)(i);
                  double v2_x = (*grad2)(i);
                  double v4_x = (*grad4)(i);
                  i++;
                  double v1_y = (*grad1)(i);
                  double v2_y = (*grad2)(i);
                  double v4_y = (*grad4)(i);
                  i++;
                  double v1_z = (*grad1)(i);
                  double v2_z = (*grad2)(i);
                  double v4_z = (*grad4)(i);
                  i++;
                  normalize(v1_x,v1_y,v1_z,1.0,0.0,0.0);
                  normalize(v2_x,v2_y,v2_z,0.0,1.0,0.0);
                  normalize(v4_x,v4_y,v4_z,1.0,0.0,0.0);
                  grad1->set(j,v1_x);
                  grad2->set(j,v2_x);
                  grad4->set(j,v4_x);
                  j++;
                  grad1->set(j,v1_y);
                  grad2->set(j,v2_y);
                  grad4->set(j,v4_y);
                  j++;
                  grad1->set(j,v1_z);
                  grad2->set(j,v2_z);
                  grad4->set(j,v4_z);

              }



           libMesh::MeshBase::const_element_iterator el =
                     mesh.active_local_elements_begin();
             const libMesh::MeshBase::const_element_iterator end_el =
                     mesh.active_local_elements_end();
        //           auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").solution;
              std::cout << "\nGetting dof map:  ... " << std::flush;

          const libMesh::DofMap & dof_map = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").get_dof_map();
              std::cout << " done \n " << std::flush;
               std::vector<libMesh::dof_id_type> dof_indices;

              std::cout << "\nLooping over elements: ... \n" << std::flush;

             for (; el != end_el; ++el)
             {
                 const libMesh::Elem * elem = *el;
                 dof_map.dof_indices(elem, dof_indices);
                 const auto blockID = elem->subdomain_id();
                 if( blockID == 3)
                 {
                  double sx = (*grad2)(dof_indices[0]);
                  double sy = (*grad2)(dof_indices[1]);
                  double sz = (*grad2)(dof_indices[2]);

                  double  cdot = cx * sx + cy * sy + cz * sz;
                  double  xfx = cx - cdot * sx;
                  double xfy = cy - cdot * sy;
                  double xfz = cz - cdot * sz;
                  normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

                double fx = sy * xfz - sz * xfy;
                double fy = sz * xfx - sx * xfz;
                double fz = sx * xfy - sy * xfx;
                normalize(fx, fy, fz, 1.0, 0.0, 0.0);

                  grad3->set(dof_indices[0],fx);
                  grad3->set(dof_indices[1],fy);
                  grad3->set(dof_indices[2],fz);
                 }
                 else if( blockID == 7 || blockID == 8 || blockID == 10 || blockID == 11   )
                 {
                  double v1_x = (*grad4)(dof_indices[0]);
                  double v2_x = (*grad2)(dof_indices[0]);
                  double v1_y = (*grad4)(dof_indices[1]);
                  double v2_y = (*grad2)(dof_indices[1]);
                  double v1_z = (*grad4)(dof_indices[2]);
                  double v2_z = (*grad2)(dof_indices[2]);

                  double cross_x = v1_y *  v2_z -  v1_z * v2_y;
                  double cross_y = v1_z *  v2_x -  v1_x * v2_z;
                  double cross_z = v1_x *  v2_y -  v1_y * v2_x;
                  normalize(cross_x, cross_y,cross_z,0.0,0.0,1.0);

                  grad3->set(dof_indices[0],cross_x);
                  grad3->set(dof_indices[1],cross_y);
                  grad3->set(dof_indices[2],cross_z);

                 }
                 else
                 {
                  double v1_x = (*grad1)(dof_indices[0]);
                  double v2_x = (*grad2)(dof_indices[0]);
                  double v1_y = (*grad1)(dof_indices[1]);
                  double v2_y = (*grad2)(dof_indices[1]);
                  double v1_z = (*grad1)(dof_indices[2]);
                  double v2_z = (*grad2)(dof_indices[2]);

                  double cross_x = v1_y *  v2_z -  v1_z * v2_y;
                  double cross_y = v1_z *  v2_x -  v1_x * v2_z;
                  double cross_z = v1_x *  v2_y -  v1_y * v2_x;
                  normalize(cross_x, cross_y,cross_z,0.0,0.0,1.0);

                  grad3->set(dof_indices[0],cross_x);
                  grad3->set(dof_indices[1],cross_y);
                  grad3->set(dof_indices[2],cross_z);
                 }
             }






             std::cout << "Calling exporter: ..." << std::flush;

              std::cout << "Calling exporter: ..." << std::flush;
              poisson4.save_exo("poisson4.exo");
              std::cout << " Done!" << std::endl;
              poisson4.write_equation_system();

              std::cout << "\nGetting fibers: ..." << std::flush;
              //auto& fibers = em.M_monowave-> M_equationSystems.get_system<libMesh::ExplicitSystem>("fibers").solution;
              auto& fibers = es.get_system<libMesh::ExplicitSystem>("fibers").solution;

              std::cout << "\nSetting fibers: ..." << std::flush;
              for(int i = first; i < last; ++i)
              {
                  fibers->set(i, (*grad3)(i));
              }
              poisson1.deleteSystems();
              poisson2.deleteSystems();
              poisson3.deleteSystems();
              poisson4.deleteSystems();
        }

      timer2.stop();

      int save_iter = 0;
      int save_iter_exo = 0;
      std::cout << "Initializing output monodomain ..." << std::endl;
      em.M_monowave->init_exo_output();
      std::cout << "Saving monodomain parameters ..." << std::endl;
      em.M_monowave->save_parameters();
//      return 0;

      BeatIt::TimeData datatime;
      datatime.setup(data, "em");
      datatime.print();

      std::string system_mass = data("monodomain/system_mass", "mass");
      std::string iion_mass = data("monodomain/iion_mass", "lumped_mass");
      bool useMidpointMethod = false;
      int step0 = 0;
      int step1 = 1;

      std::cout << "Assembling monodomain ..." << std::endl;
      em.M_monowave->assemble_matrices();
      em.M_monowave->form_system_matrix(datatime.M_dt,false, system_mass);
      std::cout << " done" << std::endl;

      double emdt = data("em/time/emdt", 1.0);
      int em_iter = static_cast<int>(emdt/datatime.M_dt);
      double emdt_exo = data("em/time/save_exo", 10.0);
      int em_exo_iter = static_cast<int>(emdt_exo/datatime.M_dt);

      std::cout << "Time loop ..." << std::endl;
      for( ; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime ; )
      {
          datatime.advance();
          // Electrophysiology

          em.M_monowave->advance();

          em.M_monowave->update_pacing(datatime.M_time);
          em.solve_reaction_step(datatime.M_dt, datatime.M_time,step0, useMidpointMethod, iion_mass);
          em.M_monowave->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);

          em.M_monowave->update_activation_time(datatime.M_time);
          // Activation part
          em.compute_activation(datatime.M_dt);
          //em.compute_activation(datatime.M_dt);
          // mechanics part
          if( 0 == datatime.M_iter%em_iter )
          {
              std::cout << "Solving mechanics ... " << std::endl;
              em.M_elasticity->setTime(datatime.M_time);
              em.solve_mechanics();
             std::cout << "* Test EM: Time: " << datatime.M_time << std::endl;
          }

          if( 0 == datatime.M_iter%datatime.M_saveIter )
          {
             save_iter++;
             em.save_gmv(save_iter, datatime.M_time);
          }
         if( 0 ==  datatime.M_iter%em_exo_iter )
         {
             save_iter_exo++;
             em.save_exo(save_iter_exo, datatime.M_time);
         }

      }

      std::cout << "Saving monodomain parameters ..." << std::endl;
      em.M_monowave->save_parameters();
      em.M_monowave->save_exo(1, datatime.M_time);

      return 0;
}

