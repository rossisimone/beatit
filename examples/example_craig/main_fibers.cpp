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

/**
 * \file main.cpp
 *
 * \class main
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 11, 2016
 *
 */

// Basic include files needed for the mesh functionality.
#include "PoissonSolver/Poisson.hpp"
#include "Electrophysiology/Bidomain/Bidomain.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"
#include "Electrophysiology/Monodomain/Monowave.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"

#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "Util/CTestUtil.hpp"
#include "Util/GenerateFibers.hpp"
#include <iomanip>
//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "Util/Timer.hpp"


void rotate_fibers( libMesh::NumericVector<double>& phi,
                    libMesh::NumericVector<double>& fibers,
                    libMesh::NumericVector<double>& sheets,
                    libMesh::NumericVector<double>& xfibers,
                    double endo_angle = 0.0, double epi_angle = 0.0)
{
    std::cout << "ROTATING FIBERS: " << std::endl;
    double potential = 0.0;

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    double xfx = 0.0;
    double xfy = 0.0;
    double xfz = 0.0;


    double W11 = 0.0;
    double W12 = 0.0;
    double W13 = 0.0;
    double W21 = 0.0;
    double W22 = 0.0;
    double W23 = 0.0;
    double W31 = 0.0;
    double W32 = 0.0;
    double W33 = 0.0;
    //
    double R11 = 0.0;
    double R12 = 0.0;
    double R13 = 0.0;
    double R21 = 0.0;
    double R22 = 0.0;
    double R23 = 0.0;
    double R31 = 0.0;
    double R32 = 0.0;
    double R33 = 0.0;
    double sa  = 0.0;
    double sa2 = 0.0;
    double teta1 = 0.0;
    double teta2 = 0.0;
    double teta  = 0.0;
    double m = 0.0;
    double q = 0.0;
    double f0x = 0.0;
    double f0y = 0.0;
    double f0z = 0.0;
    auto first_local_index_sol = phi.first_local_index();
    auto last_local_index_sol = phi.last_local_index();

    auto j = first_local_index_sol;
    auto first_local_index = fibers.first_local_index();
    auto last_local_index = fibers.last_local_index();

    for(auto i = first_local_index; i < last_local_index; )
    {

        potential = phi(j);
        j++;

        fx = fibers(i);
        fy = fibers(i+1);
        fz = fibers(i+2);

        sx = sheets(i);
        sy = sheets(i+1);
        sz = sheets(i+2);

        xfx = xfibers(i);
        xfy = xfibers(i+1);
        xfz = xfibers(i+2);

        teta1 = M_PI * epi_angle / 180.0;
        teta2 = M_PI * endo_angle / 180.0;
        m = (teta1 - teta2 );
        q = teta2;
        teta =  m * potential + q;


        //*************************************************************//
        // The fiber field F is a rotation of the flat fiber field f
        // F = R f
        // where R is the rotation matrix.
        // To compute R we need the sin(teta) and
        // the sin(teta)^2 and the cross-product matrix W (check
        // rodrigues formula on wikipedia :) )
        //*************************************************************//
        sa = std::sin (teta);
        sa2 = 2.0 * std::sin (0.5 * teta) * std::sin (0.5 * teta);

        W11 = 0.0; W12 = -sz; W13 =  sy;
        W21 =  sz; W22 = 0.0; W23 = -sx;
        W31 = -sy; W32 =  sx; W33 = 0.0;
        //
        R11 = 1.0 + sa * W11 + sa2 * ( sx * sx - 1.0 );
        R12 = 0.0 + sa * W12 + sa2 * ( sx * sy );
        R13 = 0.0 + sa * W13 + sa2 * ( sx * sz );
        R21 = 0.0 + sa * W21 + sa2 * ( sy * sx );
        R22 = 1.0 + sa * W22 + sa2 * ( sy * sy - 1.0 );
        R23 = 0.0 + sa * W23 + sa2 * ( sy * sz );
        R31 = 0.0 + sa * W31 + sa2 * ( sz * sx );
        R32 = 0.0 + sa * W32 + sa2 * ( sz * sy );
        R33 = 1.0 + sa * W33 + sa2 * ( sz * sz - 1.0 );

        f0x =   R11 * fx + R12 * fy + R13 * fz;
        f0y =   R21 * fx + R22 * fy + R23 * fz;
        f0z =   R31 * fx + R32 * fy + R33 * fz;
        BeatIt::Util::normalize(f0x, f0y, f0z, 1.0, 0.0, 0.0);


        xfx = f0y * sz - f0z * sy;
        xfy = f0z * sx - f0x * sz;
        xfz = f0x * sy - f0y * sx;
        BeatIt::Util::normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

        xfibers.set(i,   xfx);
        xfibers.set(i+1, xfy);
        xfibers.set(i+2, xfz);

        fibers.set(i,   f0x);
        fibers.set(i+1, f0y);
        fibers.set(i+2, f0z);

        i += 3;
    }
}


int main(int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace

    using namespace libMesh;
    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    libMesh::PerfLog perf_log("Timing");

   ///////////////////////////
   //  __  __ ___ ___ _  _  //
   // |  \/  | __/ __| || | //
   // | |\/| | _|\__ \ __ | //
   // |_|  |_|___|___/_||_| //
   ///////////////////////////

    // Timer for the mesh:
    perf_log.push("mesh");

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    ParallelMesh mesh(init.comm());

    // Use the MeshTools::Generation mesh generator to create a uniform
    // 3D grid
    // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
    // the number of elements on each side is read from the input file
    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data_fibers.beat", 2, "-i", "--input");
    GetPot data(datafile_name);

    // We may need XDR support compiled in to read binary .xdr files
    std::string meshfile = data("mesh/input_mesh_name", "Pippo.e");
    mesh.read (&meshfile[0]);

    double scale = 0.1;
    MeshTools::Modification::scale(mesh, scale, scale, scale);
    perf_log.pop("mesh");

    int n_refinements = data("mesh/n_ref", 0);
    for (int k = 0; k < n_refinements; ++k)
    {
        std::cout << "refinement: " << k << std::endl;
        MeshRefinement refinement(mesh);
        refinement.uniformly_refine();
    }

    auto nel = mesh.n_active_elem();
    std::cout << "Mesh has " << nel << " active elements" << std::endl;


    libMesh::EquationSystems es(mesh);
    //es.init();



    std::string pois1 = "poisson1";
    BeatIt::Poisson poisson1(es, pois1);

    std::cout << "Calling setup: ..." << std::flush;
    poisson1.setup(data, pois1);
//    std::cout << " Done!" << std::endl;
//    std::cout << "Calling assemble system: ..." << std::flush;
//    poisson1.assemble_system();
//    std::cout << " Done!" << std::endl;
//    std::cout << "Calling solve system: ..." << std::flush;
//    poisson1.solve_system();
//    std::cout << " Done!" << std::endl;
//    std::cout << "Calling gradient: ..." << std::flush;
//    poisson1.compute_elemental_solution_gradient();
//    std::cout << " Done!" << std::endl;

    std::string pois2 = "poisson2";
    BeatIt::Poisson poisson2(es, pois2);

    std::cout << "Calling setup: ..." << std::flush;
    poisson2.setup(data, pois2);
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


    std::string pois3 = "poisson3";
    BeatIt::Poisson poisson3(es, pois3);

    std::cout << "Calling setup: ..." << std::flush;
    poisson3.setup(data, pois3);
    std::cout << " Done!" << std::endl;
    std::cout << "Calling assemble system: ..." << std::flush;
    poisson3.assemble_system();
    std::cout << " Done!" << std::endl;
    std::cout << "Calling solve system: ..." << std::flush;
    poisson3.solve_system();
    std::cout << " Done!" << std::endl;
    std::cout << "Calling gradient: ..." << std::flush;
    poisson3.compute_elemental_solution_gradient();
    std::cout << " Done!" << std::endl;


    int n_blocks = mesh.n_subdomains ();
//    if(n_blocks > 1)
//    {
//         std::string pois4 = "poisson4";
//         BeatIt::Poisson poisson4(es, pois4);
//
//         std::cout << "Calling setup: ..." << std::flush;
//         poisson4.setup(data, pois4);
//         std::cout << " Done!" << std::endl;
//        std::cout << "Calling assemble system: ..." << std::flush;
//        poisson4.assemble_system();
//        std::cout << " Done!" << std::endl;
//        std::cout << "Calling solve system: ..." << std::flush;
//        poisson4.solve_system();
//        std::cout << " Done!" << std::endl;
//        std::cout << "Calling gradient: ..." << std::flush;
//        poisson4.compute_elemental_solution_gradient();
//        std::cout << " Done!" << std::endl;



        auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").solution;
        // through thickness
        auto& grad2 = poisson2. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois2+"_gradient").solution;
        // short direction
        auto& grad3 = poisson3. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois3+"_gradient").solution;
        //auto& grad4 = poisson4. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois4+"_gradient").solution;

        BeatIt::Util::normalize(*grad2, 0.0, 1.0, 0.0);
        BeatIt::Util::normalize(*grad3, 0.0, 0.0, 1.0);


        libMesh::MeshBase::const_element_iterator el =
             mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el =
             mesh.active_local_elements_end();
        const libMesh::DofMap & dof_map = es.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").get_dof_map();


//  const libMesh::DofMap & dof_map = poisson1.M_equationSystems.get_system<libMesh::ExplicitSystem>("gradient").get_dof_map();
       std::vector<libMesh::dof_id_type> dof_indices;


       //int tissue_block_ID =
          for (; el != end_el; ++el)
         {
    	  libMesh::Elem * elem = *el;
         dof_map.dof_indices(elem, dof_indices);
    	 const auto blockID = elem->subdomain_id();
          
//          if(blockID == 2)
//          {

          double sx = (*grad2)(dof_indices[0]);
    	  double sy = (*grad2)(dof_indices[1]);
    	  double sz = (*grad2)(dof_indices[2]);


          double xfx = (*grad3)(dof_indices[0]);
    	  double xfy = (*grad3)(dof_indices[1]);
    	  double xfz = (*grad3)(dof_indices[2]);

          double fx = xfy * sz - xfz * sy;
          double fy = xfz * sx - xfx * sz;
          double fz = xfx * sy - xfy * sx;

          grad1->set(dof_indices[0],fx);
    	  grad1->set(dof_indices[1],fy);
    	  grad1->set(dof_indices[2],fz);

//    	  grad3->set(dof_indices[0],xfx);
//    	  grad3->set(dof_indices[1],xfy);
//    	  grad3->set(dof_indices[2],xfz);

//           elem->subdomain_id() = 2;
//          }

         }

          BeatIt::Util::normalize(*grad1, 1.0, 0.0, 0.0);

//    }


    // Create Fiber fields
    // potential

    //  long direction
//    auto& grad1 = poisson1. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois1+"_gradient").solution;
    // through thickness
//    auto& grad2 = poisson2. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois2+"_gradient").solution;
    // short direction
//    auto& grad3 = poisson3. M_equationSystems.get_system<libMesh::ExplicitSystem>(pois3+"_gradient").solution;


    poisson1.save_exo("poisson1.exo");

    //get solver type
    std::string solver_type = data("solver", "bidomain");

    auto& fiber_system = es.add_system<libMesh::ExplicitSystem>("fibers");
    fiber_system.add_variable("fibersx", CONSTANT, MONOMIAL);
    fiber_system.add_variable("fibersy", CONSTANT, MONOMIAL);
    fiber_system.add_variable("fibersz", CONSTANT, MONOMIAL);
    fiber_system.init();
    auto& sheets_system = es.add_system<libMesh::ExplicitSystem>("sheets");
    sheets_system.add_variable("sheetsx", CONSTANT, MONOMIAL);
    sheets_system.add_variable("sheetsy", CONSTANT, MONOMIAL);
    sheets_system.add_variable("sheetsz", CONSTANT, MONOMIAL);
    sheets_system.init();
    auto& xfiber_system = es.add_system<libMesh::ExplicitSystem>("xfibers");
    xfiber_system.add_variable("xfibersx", CONSTANT, MONOMIAL);
    xfiber_system.add_variable("xfibersy", CONSTANT, MONOMIAL);
    xfiber_system.add_variable("xfibersz", CONSTANT, MONOMIAL);
    xfiber_system.init();

    // Fibers
    auto& fibers = fiber_system.solution;
    // Sheets
    auto& sheets = sheets_system.solution;
    // XFibers
    auto& xfibers = xfiber_system.solution;

    // Set fibers;
    *fibers = *grad1;
    *sheets = *grad2;
    *xfibers = *grad3;



    double  endo_angle = data("endo_angle", 0.0);
    double  epi_angle = data("epi_angle", 0.0);
    auto& phi = poisson2.get_P0_solution();
    rotate_fibers(*phi,*fibers,*sheets, *xfibers, endo_angle, epi_angle );

    poisson1.deleteSystems();
    poisson2.deleteSystems();
    poisson3.deleteSystems();
//    poisson4.deleteSystems();

    std::string fibers_file_name = data("output_file", "fibers.exo");
    libMesh::ExodusII_IO exporter(mesh);
    exporter.write_equation_systems(fibers_file_name, es);
    exporter.write_element_data(es);
    return 0;
}

