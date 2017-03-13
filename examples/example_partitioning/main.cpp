#include "libmesh/parallel_mesh.h"
#include "libmesh/checkpoint_io.h"

int main(int argc, char** argv)
{
    using namespace libMesh;
    // Initialize libraries, like in example 2.
    LibMeshInit init (argc, argv, MPI_COMM_WORLD);
    ParallelMesh mesh(init.comm());
    CheckpointIO importer(mesh, true);
    importer.read("LA_T4_new_m1.cpr");
    if(mesh.is_prepared() ) std::cout << "Mesh already prepared" << std::endl;
    else mesh.prepare_for_use();
//    mesh.read("LA_T4_new_m1.cpr");
    mesh.print_info(std::cout);
    return 0;
}
