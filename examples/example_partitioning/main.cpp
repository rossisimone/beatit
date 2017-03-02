#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/parallel_mesh.h"


int main(int argc, char** argv)
{
    using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);

      // Create a mesh, with dimension to be overridden later, distributed
      // across the default MPI communicator.
      Mesh mesh(init.comm());

    MeshTools::Generation::build_square (mesh,
                                         10, 10,
                                         0., 10.,
                                         0., 2.,
                                         TRI3);

    mesh.prepare_for_use();
    mesh.update_post_partitioning();

    std::cout << "Is mesh serial: " << mesh.is_serial() << std::endl;
    std::cout << "Num partitions: " << mesh.n_partitions () << std::endl;


    ParallelMesh output(mesh);
    output.write("mesh.exo");
    return 0;
}
