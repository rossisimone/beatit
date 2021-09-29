// C++ include files that we need
#include <iostream>
// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/boundary_info.h"
#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/getpot.h"

using namespace libMesh;
// The main program.
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
    mesh.read(mesh_file + ".msh");
    mesh.print_info();
    // Refine mesh
    std::string output_name = data("output", mesh_file);
    ExodusII_IO (mesh).write(output_name+"_m0.e");
    int nref = data("nref", 0);
    if(nref > 0)
    {
      for(int nr = 1; nr <= nref; nr++)
      {
          ParallelMesh refined_mesh(mesh);
          MeshRefinement refine(refined_mesh);
          std::cout << "Refine: " << nr << std::endl;
          refine.uniformly_refine(nr);
          std::cout << "Exporting: " << nr << std::endl;
          ExodusII_IO (refined_mesh).write(output_name+"_m"+std::to_string(nr)+".e");
      }
    }
    return 0;
}
