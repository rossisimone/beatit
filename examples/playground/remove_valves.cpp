// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1>Systems Example 1 - Stokes Equations</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example shows how a simple, linear system of equations
// can be solved in parallel.  The system of equations are the familiar
// Stokes equations for low-speed incompressible fluid flow.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>

// Basic include file needed for the mesh functionality.
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/boundary_info.h"
#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/transient_system.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include <libmesh/point_locator_tree.h>
#include "libmesh/getpot.h"
#include "Util/IO/io.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;




int main(int argc, char ** argv)
{
    // Initialize libMesh.
    LibMeshInit init(argc, argv);
    // Read input file
    GetPot commandLine ( argc, argv );
    std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
    GetPot data(datafile_name);

    std::string mesh_file = data("mesh", "NONE");
    ParallelMesh mesh(init.comm());
    mesh.read(mesh_file);
    mesh.print_info();
    ExodusII_IO importer(mesh);
    mesh.prepare_for_use(true);
    mesh.print_info();

    ParallelMesh mesh_no_valves(mesh);

    /*
    std::vector<std::string> vars = importer.get_elem_var_names();
    for(auto & v : vars) std::cout << v << std::endl;
    
    EquationSystems es(mesh);
    auto & sys = es.add_system<ExplicitSystem>("sys");
    for(auto & vn : vars) sys.add_variable(vn, CONSTANT, MONOMIAL);
    es.init();
    es.print_info();

    std::cout << "copying fibers:" << std::endl;
    for(auto & vn : vars) 
    {
        std::cout << vn << std::endl;
        importer.copy_elemental_solution(sys, vn , vn, 1);
    }
    std::cout << "exporting " << std::endl;
     */
    for(auto & elem : mesh_no_valves.active_local_element_ptr_range() )
    {
       bool remove_current_element = false;
       auto blockID = elem->subdomain_id();
       switch(blockID)
       {
          case 13:
          case 14:
          case 24:
          case 50:
          case 60:
          case 661:
          case 662:
          case 663:
          case 664:
          case 665:
          case 667:
            remove_current_element = true;
            break;
          default:
            break;
       }

       if(remove_current_element)
       {
	   mesh_no_valves.delete_elem(elem);
       }

    }
    mesh_no_valves.prepare_for_use(true);
    ExodusII_IO(mesh_no_valves).write("no_valves_test.e");

    // Refine mesh
    //std::string output_name = data("output", "final.e");
    //ExodusII_IO exporter(mesh);
    //ExodusII_IO exporter(mesh);
    //exporter.write_equation_systems(output_name, es);
    //exporter.write_element_data(es);
    return 0;
}
