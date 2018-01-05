#include "Util/IO/io.hpp"
#include <sys/stat.h>
#include "libmesh/getpot.h"
#include <iomanip>
#include "libmesh/ghosting_functor.h"
namespace BeatIt
{

void printBanner(std::ostream& cout)
{
    cout << "\n\t .______    _______     ___   .___________.    __  .___________.";
    cout << "\n\t |   _  \\  |   ____|   /   \\  |           |   |  | |           |";
    cout << "\n\t |  |_)  | |  |__     /  ^  \\ `---|  |----`   |  | `---|  |---- ";
    cout << "\n\t |   _  <  |   __|   /  /_\\  \\    |  |        |  |     |  |     ";
    cout << "\n\t |  |_)  | |  |____ /  _____  \\   |  |        |  |     |  |     ";
    cout << "\n\t |______/  |_______/__/     \\__\\  |__|        |__|     |__| 	 ";
    cout << "\n";
}

void saveData(double time, std::vector<double>& var, std::ostream& output)
{
    std::streamsize ss = std::cout.precision();
    output << std::fixed << std::setprecision(15);
    output << time;
    for (const auto& v : var)
    {
        output << " " << v;
    }
    output << std::setprecision(ss);
    output << "\n";
}

bool readList(std::string& list, std::vector<std::string>& container)
{
    auto beg = list.begin();
    auto end = list.end();
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    bool ok = qi::phrase_parse(beg, end, (*~qi::char_(",")) % ',', ascii::blank, container);
    return ok;
}

void createOutputFolder(const libMesh::Parallel::Communicator & comm, std::string& output_folder)
{
    struct stat out_dir;
    if (stat(&output_folder[0], &out_dir) != 0)
    {
        if (comm.rank() == 0)
        {
            mkdir(output_folder.c_str(), 0777);
        }
    }
}

GetPot readInputFile(int argc, char ** argv)
{
    GetPot cl(argc, argv);
    std::string datafile_name = cl.follow("data.beat", 2, "-i", "--input");
    return GetPot(datafile_name);
}

// From: http://stackoverflow.com/a/6417908/2042320
std::string remove_extension(const std::string & filename)
{
    size_t lastdot = filename.find_last_of(".");

    if (lastdot == std::string::npos)
        return filename;

    return filename.substr(0, lastdot);
}

void serial_mesh_partition(const libMesh::Parallel::Communicator & comm, std::string file, libMesh::ParallelMesh * mesh, int n_ref)
{
    //auto& comm = mesh.comm();
    using namespace libMesh;
    int np = comm.size();
    std::cout << "reading mesh to run with " << np << " processors" << std::endl;
    std::ostringstream outputname;
    const bool binary = true;
    std::string filename_no_extension = remove_extension(file);
    if (n_ref > 0)
        filename_no_extension += "_ref_" + std::to_string(n_ref);
    outputname << filename_no_extension << '.' << np << (binary ? ".cpr" : ".cpa");

    int color = comm.rank();
    int key = comm.rank();

    struct stat buffer;
    bool file_already_exists = stat(outputname.str().c_str(), &buffer) == 0;
    if (!file_already_exists)
    {

        libMesh::Parallel::Communicator rank_0_comm;
        comm.split(color, key, rank_0_comm);

        if (comm.rank() == 0)
        {
            ReplicatedMesh mesh_to_split(rank_0_comm);

            mesh_to_split.read(file);

            for (int k = 0; k < n_ref; ++k)
            {
                std::cout << "refinement: " << k << std::endl;
                MeshRefinement refinement(mesh_to_split);
                refinement.uniformly_refine();
            }

            mesh_to_split.print_info();

            DefaultCoupling default_coupling;
            MetisPartitioner partitioner;

            processor_id_type n_procs = np;
            libMesh::out << "\nWriting out files for " << n_procs << " processors...\n\n" << std::endl;
            // Reset the partitioning each time after the first one
            libMesh::out << "Partitioning" << std::endl;
            partitioner.partition(mesh_to_split, n_procs);
            // When running in parallel each processor will write out a portion of the mesh files

            processor_id_type num_chunks = n_procs;
            processor_id_type remaining_chunks = 0;

            processor_id_type my_num_chunks = num_chunks;

            processor_id_type my_first_chunk = 0;

            processor_id_type rank = 0;
            processor_id_type comm_size = 1;

            libMesh::out << "Writing " << my_num_chunks << " Files" << std::endl;

            //const bool binary = !libMesh::on_command_line("--ascii");

            CheckpointIO cpr(mesh_to_split);
            cpr.current_processor_ids().clear();
            for (unsigned int i = my_first_chunk; i < my_first_chunk + my_num_chunks; i++)
                cpr.current_processor_ids().push_back(i);
            cpr.current_n_processors() = n_procs;
            cpr.parallel() = true;
            cpr.binary() = binary;
            std::cout << outputname.str() << std::endl;
            cpr.write(outputname.str());

        }
    }
    else
    {
        std::cout << "File: " << outputname.str() << " alredy exists!!!" << std::endl;
    }
    comm.barrier();
    std::cout << key << ": partitioning done" << std::endl;

    if (mesh)
    {
        CheckpointIO importer(*mesh, true);
        importer.read(outputname.str());
        mesh->skip_partitioning(true);
        mesh->prepare_for_use(true);
    }
    //rank_0_comm.clear();
}

} // namespace BeatIt

