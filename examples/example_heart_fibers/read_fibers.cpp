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
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"

int main(int argc, char ** argv)
{

    using namespace libMesh;
    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);
    libMesh::Mesh mesh(init.comm());
    libMesh::ExodusII_IO importer(mesh);
    std::string meshfile = "heart.e";
    importer.read(&meshfile[0]);
    mesh.print_info();
    mesh.prepare_for_use(true);

    libMesh::EquationSystems es(mesh);
    typedef libMesh::ExplicitSystem FiberSystem;
    FiberSystem& f_sys = es.add_system < FiberSystem > ("fibers");
    f_sys.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.init();
    auto& f_v = f_sys.solution;
    FiberSystem& s_sys = es.add_system<FiberSystem>("sheets");
    s_sys.add_variable( "sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.add_variable( "sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.add_variable( "sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys.init();
    auto& s_v = s_sys.solution;

    FiberSystem& n_sys = es.add_system<FiberSystem>("xfibers");
    n_sys.add_variable( "xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.add_variable( "xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.add_variable( "xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys.init();
    auto& n_v = n_sys.solution;
    auto vars = importer.get_elem_var_names();
    for (auto && name : vars)
        std::cout << name << std::endl;
    importer.copy_elemental_solution(f_sys, "fibersx", "fibersX");
    importer.copy_elemental_solution(f_sys, "fibersy", "fibersY");
    importer.copy_elemental_solution(f_sys, "fibersz", "fibersZ");
    importer.copy_elemental_solution(s_sys, "sheetsx", "sheetsX");
    importer.copy_elemental_solution(s_sys, "sheetsy", "sheetsY");
    importer.copy_elemental_solution(s_sys, "sheetsz", "sheetsZ");
    importer.copy_elemental_solution(n_sys, "xfibersx", "xfibersX");
    importer.copy_elemental_solution(n_sys, "xfibersy", "xfibersY");
    importer.copy_elemental_solution(n_sys, "xfibersz", "xfibersZ");



    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    libMesh::Mesh mesh2(init.comm());
    libMesh::ExodusII_IO importer2(mesh2);
    importer2.read(&meshfile[0]);
    mesh2.print_info();
    mesh2.prepare_for_use();

    libMesh::EquationSystems es2(mesh2);
    typedef libMesh::ExplicitSystem FiberSystem;
    FiberSystem& f_sys2 = es2.add_system < FiberSystem > ("fibers");
    f_sys2.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys2.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys2.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys2.init();
    auto& f_v2 = f_sys2.solution;
     std::cout << "Copying f:" << std::endl;
    *f_v2 = *f_v;
    FiberSystem& s_sys2 = es2.add_system<FiberSystem>("sheets");
    s_sys2.add_variable( "sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys2.add_variable( "sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys2.add_variable( "sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys2.init();
    auto& s_v2 = s_sys2.solution;
     std::cout << "Copying s:" << std::endl;
    *s_v2 = *s_v;
    FiberSystem& n_sys2 = es2.add_system<FiberSystem>("xfibers");
    n_sys2.add_variable( "xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys2.add_variable( "xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys2.add_variable( "xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys2.init();
     std::cout << "Copying n:" << std::endl;
    auto& n_v2 = n_sys2.solution;
    *n_v2 = *n_v;
     std::cout << "export:" << std::endl;
    libMesh::ExodusII_IO exporter(mesh2);
    std::string output1 = "imported_heart.e";
    exporter.write_equation_systems(output1, es2);
    exporter.write_element_data(es2);

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

/*    libMesh::Mesh mesh3(init.comm());
    libMesh::ExodusII_IO importer3(mesh3);
    importer3.read(output1);
    mesh3.print_info();
    mesh3.prepare_for_use();

    libMesh::EquationSystems es3(mesh3);
    typedef libMesh::ExplicitSystem FiberSystem;
    FiberSystem& f_sys3 = es3.add_system < FiberSystem > ("fibers");
    f_sys3.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys3.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys3.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys3.init();
    auto& f_v3 = f_sys3.solution;
    FiberSystem& s_sys3 = es3.add_system<FiberSystem>("sheets");
    s_sys3.add_variable( "sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys3.add_variable( "sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys3.add_variable( "sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_sys3.init();
    auto& s_v3 = s_sys3.solution;

    FiberSystem& n_sys3 = es3.add_system<FiberSystem>("xfibers");
    n_sys3.add_variable( "xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys3.add_variable( "xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys3.add_variable( "xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_sys3.init();
    auto& n_v3 = n_sys3.solution;
    importer3.copy_elemental_solution(f_sys3, "fibersx", "fibersx");
    importer3.copy_elemental_solution(f_sys3, "fibersy", "fibersy");
    importer3.copy_elemental_solution(f_sys3, "fibersz", "fibersz");
    importer3.copy_elemental_solution(s_sys3, "sheetsx", "sheetsx");
    importer3.copy_elemental_solution(s_sys3, "sheetsy", "sheetsy");
    importer3.copy_elemental_solution(s_sys3, "sheetsz", "sheetsz");
    importer3.copy_elemental_solution(n_sys3, "xfibersx", "xfibersx");
    importer3.copy_elemental_solution(n_sys3, "xfibersy", "xfibersy");
    importer3.copy_elemental_solution(n_sys3, "xfibersz", "xfibersz");

    libMesh::ExodusII_IO exporter3(mesh3);
    std::string output3 = "imported_heart3.e";
    exporter3.write_equation_systems(output3, es3);
    exporter3.write_element_data(es3);

*/

    return 0;
}
