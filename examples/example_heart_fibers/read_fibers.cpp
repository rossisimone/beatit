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
#include "libmesh/mesh_refinement.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/unv_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"

int main(int argc, char ** argv)
{

    using namespace libMesh;
    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    libMesh::Mesh mesh(init.comm());
    libMesh::ExodusII_IO importer(mesh);
    std::string meshfile = "whole_heart_with_fibers_v4.e";
    importer.read(&meshfile[0]);
    mesh.print_info();
    std::cout << "prepare for use ... " << std::endl;
    mesh.prepare_for_use(true);

    std::cout << "Equation Systems ... " << std::endl;
    libMesh::EquationSystems es(mesh);

    typedef libMesh::ExplicitSystem FiberSystem;
    FiberSystem& f_sys = es.add_system < FiberSystem > ("fibers");
    f_sys.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    std::cout << "Initializing fibers system ... " << std::endl;
    f_sys.init();
    auto& f_v = f_sys.solution;

    FiberSystem& s_sys = es.add_system<FiberSystem>("sheets");
    s_sys.add_variable( "sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.add_variable( "sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.add_variable( "sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    std::cout << "Initializing sheets system ... " << std::endl;
    s_sys.init();
    auto& s_v = s_sys.solution;

    FiberSystem& n_sys = es.add_system<FiberSystem>("xfibers");
    n_sys.add_variable( "xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.add_variable( "xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.add_variable( "xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    std::cout << "Initializing cross fibers system ... " << std::endl;
    n_sys.init();
    auto& n_v = n_sys.solution;

    auto vars = importer.get_elem_var_names();
    std::cout << "Outputting variables' names... " << std::endl;
    for (auto && name : vars)
        std::cout << name << std::endl;
    std::cout << "Copying fibersx ... " << std::endl;
    importer.copy_elemental_solution(f_sys, "fibersx", "fibersx");
    std::cout << "Copying fibersy ... " << std::endl;
    importer.copy_elemental_solution(f_sys, "fibersy", "fibersy");
    std::cout << "Copying fibersz ... " << std::endl;
    importer.copy_elemental_solution(f_sys, "fibersz", "fibersz");
    std::cout << "Copying sheetsx ... " << std::endl;
    importer.copy_elemental_solution(s_sys, "sheetsx", "sheetsx");
    std::cout << "Copying sheetsy ... " << std::endl;
    importer.copy_elemental_solution(s_sys, "sheetsy", "sheetsy");
    std::cout << "Copying sheetsz ... " << std::endl;
    importer.copy_elemental_solution(s_sys, "sheetsz", "sheetsz");
    std::cout << "Copying xfibersx ... " << std::endl;
    importer.copy_elemental_solution(n_sys, "xfibersx", "xfibersx");
    std::cout << "Copying xfibersy ... " << std::endl;
    importer.copy_elemental_solution(n_sys, "xfibersy", "xfibersy");
    std::cout << "Copying xfibersz ... " << std::endl;
    importer.copy_elemental_solution(n_sys, "xfibersz", "xfibersz");


    es.print_info();

    std::cout << "Exporting test ... " << std::endl;
    libMesh::ExodusII_IO exporter2(mesh);
    exporter2.write("test.e");
    exporter2.write_element_data(es);

    MeshRefinement mrefine(mesh);
    std::cout << "Uniformly refine ... " << std::endl;
    mrefine.uniformly_refine();
    mesh.print_info();
    std::cout << "Equation systems reinit ... " << std::endl;
    es.reinit();
    es.print_info();

    std::cout << "Exporting ... " << std::endl;
    libMesh::ExodusII_IO exporter(mesh);
    exporter.write("refined_whole_heart_with_fibers_v4.e");
    exporter.write_element_data(es);
    return 0;
}