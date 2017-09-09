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
#include "libmesh/gmsh_io.h"

int main (int argc, char ** argv)
{

    using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);
      libMesh::Mesh mesh(init.comm());
      libMesh::GmshIO importer(mesh);
      std::string meshfile = "sphere.msh";
      importer.read(meshfile);
      //libMesh::ExodusII_IO(mesh).write("heart.e");
    return 0;
}
