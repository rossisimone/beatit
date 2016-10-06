#include <iostream>

#include "Util/IO/io.hpp"


#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>

#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/node.h"
#include "libmesh/point.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/unv_io.h"

int main (int argc, char ** argv)
{
    BeatIt::printBanner(std::cout);

    using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);
      libMesh::Mesh mesh(init.comm());
      mesh.set_mesh_dimension(3);
      mesh.set_spatial_dimension(3);
      mesh.reserve_elem(29184);
      mesh.reserve_nodes( 36650 );

      libMesh::Mesh mesh_RA(init.comm());
      mesh_RA.set_mesh_dimension(3);
      mesh_RA.set_spatial_dimension(3);
      mesh_RA.reserve_elem( 4736 );
      mesh_RA.reserve_nodes ( 6155 );

      std::ifstream nodes_list("nodes.txt");
      std::ifstream elements_list("elements.txt");

      std::string line;

      double x, y, z;
      int c = 0;
      int n = 0;

      if(nodes_list.is_open())
      {
			while ( getline (nodes_list,line) )
			{
			    std::istringstream ss(line);
				//line will have
			    ss >> x >> y >> z;
			    libMesh::Node node(x,y,z);
			    mesh.add_point(libMesh::Point(x,y,z),n);
                if(n < 6155 ) mesh.add_point(libMesh::Point(x,y,z),n);
                n++;

//			    libMesh::Elem * elem = mesh.add_elem (new libMesh::NodeElem);
			}
      }

      int n1, n2, n3, n4, n5, n6, n7, n0;
      if(elements_list.is_open())
      {
            while ( getline (elements_list,line) )
            {
                std::istringstream ss(line);


                libMesh::Elem * elem = mesh.add_elem(new libMesh::Hex8);
                // In libMesh:                          In file:
                //
                /*  HEX8: 7        6                    HEX8: 4        6
                *        o--------o                          o--------o
                *       /:       /|                         /:       /|
                *      / :      / |                        / :      / |
                *   4 /  :   5 /  |                     5 /  :   7 /  |
                *    o--------o   |                      o--------o   |
                *    |   o....|...o 2                    |   o....|...o 2
                *    |  .3    |  /                       |  .0    |  /
                *    | .      | /                        | .      | /
                *    |.       |/                         |.       |/
                *    o--------o                          o--------o
                *    0        1                          1        3
                *
                *    In file      0   1   2   3   4   5   6   7
                *    In libmesh  n3  n0  n2  n1  n7  n4  n6  n5
                */
                ss >> n3 >> n0 >> n2 >> n1 >> n7 >> n4 >> n6 >> n5;
                elem->set_node(0) = mesh.node_ptr( n0 );
                elem->set_node(1) = mesh.node_ptr( n1 );
                elem->set_node(2) = mesh.node_ptr( n2 );
                elem->set_node(3) = mesh.node_ptr( n3 );
                elem->set_node(4) = mesh.node_ptr( n4 );
                elem->set_node(5) = mesh.node_ptr( n5 );
                elem->set_node(6) = mesh.node_ptr( n6 );
                elem->set_node(7) = mesh.node_ptr( n7 );
                if(c < 4736 )
                {
                    elem->subdomain_id() = 1;
                    mesh_RA.add_elem(new libMesh::Hex8(elem));
                }
                else if ( (c >= 4736 && c < 12672) || c >= 28864 || ( c >= 26752 && c <= 27071 )) elem->subdomain_id() = 2;
                else elem->subdomain_id() = 3;
                c++;

                // TODO:: Check sides
                //mesh.boundary_info->add_side()
            }
      }

      mesh.prepare_for_use (/*skip_renumber =*/ false);
      mesh.subdomain_name(3) = "H";
      mesh.subdomain_name(1) = "RA";
      mesh.subdomain_name(2) = "LA";

      mesh_RA.prepare_for_use (/*skip_renumber =*/ false);
      std::cout << "num elements: " << c << ", num nodes: " << n << std::endl;
      mesh.print_info();

      nodes_list.close();
      elements_list.close();

      libMesh::ExodusII_IO(mesh).write("heart.e");
      libMesh::ExodusII_IO(mesh_RA).write("RA.e");
      libMesh::GmshIO(mesh).write("heart.msh");
      libMesh::GmshIO(mesh_RA).write("RA.msh");
      libMesh::UNVIO(mesh).write("heart.unv");
      libMesh::UNVIO(mesh_RA).write("RA.unv");
    return 0;
}
