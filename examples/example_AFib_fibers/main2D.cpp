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

/*
 * main.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */
// Basic include files needed for the mesh functionality.
#include "PoissonSolver/Poisson.hpp"
#include "Util/IO/io.hpp"
#include "Util/GenerateFibers.hpp"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/quadrature_gauss.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
#include <iomanip>
#include "libmesh/exodusII_io.h"

#include "Util/Noise.hpp"
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"
#include "libmesh/petsc_matrix.h"

#include "petscmat.h"

using namespace libMesh;




enum class FibrosisModel
{
    NONE,
    // From
    // Patient-specific modeling of atrial fibrosis increases the accuracy of sinus rhythm simulations and may explain maintenance of atrial fibrillation
    // Martin W.Krueger et al.
    KRUEGER, //
    // FROM
    // Epicardial Fibrosis Explains Increased Endo–Epicardial Dissociation and Epicardial Breakthroughs in Human Atrial Fibrillation
    // Ali Gharaviri et al.
    PEZZUTO,
    EPI_ENDO_DISSOCIATION,
    OURS
};

std::map<std::string, FibrosisModel> fibrosis_model_map = { { "none", FibrosisModel::NONE},
                                                            { "krueger", FibrosisModel::KRUEGER},
                                                            { "pezzuto", FibrosisModel::PEZZUTO},
                                                            { "dissociation", FibrosisModel::PEZZUTO},
                                                            { "ours", FibrosisModel::OURS}};

void evaluate(double f[], double s[], double n[], std::vector<double>& phi, std::vector<libMesh::Gradient>& dphi, int blockID, Point& x, GetPot& data, int dim = 2);
void evaluate_atlas(double f[], double s[], double n[], std::vector<double>& phi, std::vector<libMesh::Gradient>& dphi, Point& x, EquationSystems& es, PointLocatorTree& locator);

void assemble_smoother(EquationSystems& es, double D = 1.0, double DD = 1.0, int sector = 1);

void cross(double v1[], double v2[], double cross_product[])
{
    cross_product[0] = v1[1] * v2[2] - v1[2] * v2[1];
    cross_product[1] = v1[2] * v2[0] - v1[0] * v2[2];
    cross_product[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return;
}

void cross(double v1[], int n, double v2[], int m, std::vector<libMesh::Gradient>& dphi, double cross_product[])
{
    auto normalize = [](double& x, double& y, double& z,
            double X, double Y, double Z)
    {
        double norm = std::sqrt( x * x + y * y + z * z);
        if(norm >= 1e-12 )
        {
            x /= norm;
            y /= norm;
            z /= norm;
        }
        else
        {
            x = X;
            y = Y;
            z = Z;
        }
    };

    v1[0] = dphi[n](0);
    v1[1] = dphi[n](1);
    v1[2] = dphi[n](2);
    normalize(v1[0], v1[1], v1[2], 1.0, 0.0, 0.0);

    v2[0] = dphi[m](0);
    v2[1] = dphi[m](1);
    v2[2] = dphi[m](2);
    normalize(v2[0], v2[1], v2[2], 0.0, 1.0, 0.0);
    // compute f as the cross product
    cross(v1, v2, cross_product);
    normalize(cross_product[0], cross_product[1], cross_product[2], 0.0, 0.0, 1.0);
    return;
}

static double laa_threshold = 0.0;
static double fo_threshold = 0.0;
static double lpv_antra_threshold = 0.0;
static double rpv_antra_threshold = 0.0;
static double lpw_threshold = 0.0;
static double rpw_threshold = 0.0;

int main(int argc, char ** argv)
{
    typedef libMesh::ExplicitSystem FiberSystem;

    // Bring in everything from the libMesh namespace
    BeatIt::printBanner(std::cout);
    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    GetPot data = BeatIt::readInputFile(argc, argv);

    std::string ep_data_filename = data("ep_data","ep.beat");
    GetPot data_ep(ep_data_filename);
    // Create the TimeData object
    BeatIt::TimeData datatime;
    // Set it up using the input file
    datatime.setup(data_ep, "");
    // Output on screen the stored variables
    datatime.print();

    // allow us to use higher-order approximation.
    // Create a mesh, with dimension to be overridden later, on the
    // default MPI communicator.
    libMesh::MeshBase * mesh_ptr;
    libMesh::MeshBase * atlas_mesh_ptr;

    // We may need XDR support compiled in to read binary .xdr files
    std::string meshfile = data("mesh", "NONE");
    std::string atlas_meshfile = data("atlas_mesh", "NONE");
    bool read_fibers = data("read_fibers", false);

    libMesh::ExodusII_IO* importer_ptr;
    libMesh::ExodusII_IO* atlas_importer_ptr;
    libMesh::PointLocatorTree* locator_ptr;
    // Read the input mesh.
    if (atlas_meshfile != "NONE")
    {
        std::cout << "reading atlas mesh: " << atlas_meshfile << std::endl;
        atlas_mesh_ptr = new libMesh::ReplicatedMesh(init.comm());
        atlas_importer_ptr = new libMesh::ExodusII_IO(*atlas_mesh_ptr);
        atlas_importer_ptr->read(atlas_meshfile);
        atlas_mesh_ptr->prepare_for_use(true);
        locator_ptr = new PointLocatorTree(*atlas_mesh_ptr);
    }
    if (meshfile != "NONE")
    {
        // We may need XDR support compiled in to read binary .xdr files
        int n_refinements = 0;
        if(read_fibers)
        {
            mesh_ptr = new libMesh::ReplicatedMesh(init.comm());
            importer_ptr = new libMesh::ExodusII_IO(*mesh_ptr);

            importer_ptr->read(meshfile);
            mesh_ptr->prepare_for_use(true);
            std::set<libMesh::subdomain_id_type> bids;
            std::cout << "subdomain IDs:" << std::endl;
            mesh_ptr->subdomain_ids(bids);
            for(auto && bid : bids) std::cout << bid << std::endl;
            std::cout << std::endl;
//            mesh.read(meshfile);
//             //
//
//            libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
//            const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
//            for (; el != end_el; ++el)
//            {
//                libMesh::Elem * elem = *el;
//                elem->subdomain_id() = 1;
//                for (unsigned int side=0; side<elem->n_sides(); side++)
//                {
//                    if (elem->neighbor_ptr(side) != libmesh_nullptr)
//                    {
//                        elem->neighbor_ptr(side)->subdomain_id() = 1;
//
//                    }
//                }
//            }
//            mesh.prepare_for_use(true);
//            std::cout << "N SUBDOMAINS: " << mesh.n_subdomains ()  << std::endl;

        }
        else
        {
            mesh_ptr = new libMesh::ParallelMesh(init.comm());
            BeatIt::serial_mesh_partition(init.comm(), meshfile, dynamic_cast<libMesh::ParallelMesh*>(mesh_ptr), n_refinements);
        }

    }
    else
    {
        std::cout << "ABORTING: MeshFile " << meshfile << " not found!" << std::endl;
        std::cout << "Specify a valid mesh file in the input file under: mesh = " << std::endl;
        throw std::runtime_error("MeshFile not found!");
    }
    libMesh::MeshBase & mesh = *mesh_ptr;

    // Add nodesets for BCs
    bool add_nodesets = data("add_nodesets",false);
    if(add_nodesets)
    {
        std::string ouput_with_nodeset = data("ouput_with_nodeset","mesh_with_nodesets.e");
        std::string nodeset_x_list = data("nodeset_x","");
        std::string nodeset_y_list = data("nodeset_y","");
        std::string nodeset_z_list = data("nodeset_z","");
        std::string nodeset_r_list = data("nodeset_r","");
        std::string nodes_id_list = data("nodes_id","");
        std::string nodeset_id_list = data("nodeset_id","");
        double nodeset_tol = data("nodeset_tolerance", 1e-10);
        std::vector<double> nodeset_x;
        std::vector<double> nodeset_y;
        std::vector<double> nodeset_z;
        std::vector<double> nodeset_r;
        std::vector<int> nodeset_IDs;
        std::vector<int> nodes_IDs;
        BeatIt::readList(nodeset_x_list, nodeset_x);
        BeatIt::readList(nodeset_y_list, nodeset_y);
        BeatIt::readList(nodeset_z_list, nodeset_z);
        BeatIt::readList(nodeset_r_list, nodeset_r);
        BeatIt::readList(nodeset_id_list, nodeset_IDs);
        BeatIt::readList(nodes_id_list, nodes_IDs);

        std::cout << "setup " << nodeset_x.size() << " nodesets" << std::endl;
        for (int k = 0; k < nodeset_x.size() ; ++k)
        {
            std::cout << "(" << nodeset_x[k] << ", " << nodeset_y[k] << ", " << nodeset_z[k] << "), r = " << nodeset_r[k] << ", ID: " << nodeset_IDs[k] << std::endl;
        }
        for(int k = 0; k < nodes_IDs.size(); ++k)
        {
            std::cout << "\nAdding nodeset: " << nodeset_IDs[k] << " at node ID: " << nodes_IDs[k] << std::endl;
            mesh.node_ptr(nodes_IDs[k])->print();
            mesh.get_boundary_info().add_node(mesh.node_ptr(nodes_IDs[k]), nodeset_IDs[k]);
        }
        libMesh::ExodusII_IO(mesh).write(ouput_with_nodeset);
        return 0;
    }

//    std::set<int> added_nodesets;
//    std::cout << "Adding nodesets" << std::endl;
//
//    for(auto& node : mesh.local_node_ptr_range())
//    {
//        libMesh::Point p(*node);
//        for(int k = 0; k < nodeset_x.size(); ++k)
//        {
//            libMesh::Point q(nodeset_x[k], nodeset_y[k], nodeset_z[k]);
//            double r = nodeset_r[k];
//            if(r > 0)
//            {
//                libMesh::Point qmp(p-q);
//                if( qmp.norm() <= r)
//                {
//                    //std::cout << "Adding nodeset: " << nodeset_IDs[k] << " to node:  " << node->id() << std::endl;
//                    mesh.get_boundary_info().add_node(node, nodeset_IDs[k]);
//                    added_nodesets.insert(nodeset_IDs[k]);
//                }
//            }
//        }
//    }
//    if(added_nodesets.size() !=  nodeset_x.size() )
//    {
//        std::cout << "SOMETHING WENT WRONG: may be you are running in parallel!" << std::endl;
//        //return 1;
//    }

    std::cout << "Mesh prepare for use" << std::endl;
    mesh.prepare_for_use();
    //std::cout << "setup " << nodeset_IDs.size() << " nodesets done!" << std::endl;
    //mesh.get_boundary_info().print_info();
    // Create LibMesh equation system
    libMesh::EquationSystems es(mesh);
    std::shared_ptr<libMesh::EquationSystems> atlas_es_ptr;
    if (atlas_meshfile != "NONE")
    {
        std::cout << "reading fibers from: " << atlas_meshfile << std::endl;

        atlas_es_ptr.reset(new libMesh::EquationSystems(*atlas_mesh_ptr) );
        // DEFINE FIBER SYSTEMS:
        FiberSystem& f_sys = atlas_es_ptr->add_system<FiberSystem>("fibers");
        f_sys.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_vector("backup");
        f_sys.init();

        FiberSystem& s_sys = atlas_es_ptr->add_system<FiberSystem>("sheets");
        s_sys.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.init();

        FiberSystem& n_sys = atlas_es_ptr->add_system<FiberSystem>("xfibers");
        n_sys.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.init();

        const int num_fiber_systems = 3;
        std::vector < std::string > fibers(num_fiber_systems);
        fibers[0] = "fibers";
        fibers[1] = "sheets";
        fibers[2] = "xfibers";

        auto& elem_var_names = atlas_importer_ptr->get_elem_var_names();
        auto elem_first = elem_var_names.begin();
        auto elem_end = elem_var_names.end();

        for (int k = 0; k < num_fiber_systems; ++k)
        {
            libMesh::System& system = atlas_es_ptr->get_system(fibers[k]);
            std::string name = system.name();
            std::cout << "Importing System: " << name << std::endl;
            int n_vars = system.n_vars();
            for (int l = 0; l < n_vars; ++l)
            {
                std::string var_name = system.variable_name(l);
                auto elem_it = std::find(elem_first, elem_end, var_name);
                if (elem_it != elem_end)
                {
                    std::cout << "\t elemental variable: " << *elem_it << std::endl;
                    atlas_importer_ptr->copy_elemental_solution(system, *elem_it, *elem_it, 1);
                }
            }
        }
        f_sys.get_vector("backup").close();
        f_sys.solution->close();
        f_sys.get_vector("backup") = *f_sys.solution;

    }
    // Name each system with a number

    unsigned int dim = mesh.mesh_dimension();



    if(!read_fibers)
    {
        lpv_antra_threshold = data("lpv_antra", 0.225);
        rpv_antra_threshold = data("rpv_antra", 0.775);
        lpw_threshold = data("lpw", 0.2);
        rpw_threshold = data("rpw", 0.8);
        //std::cout << "Solving " << N << " Poisson problems" << std::endl;
        std::string input_files_list = data("input_files","data.beat");
        std::vector<std::string> input_files;
        BeatIt::readList(input_files_list, input_files);
        // Solve each Poisson problem
        int N = input_files.size();
        using namespace libMesh;
    //    for(int i = 0; i < input_files.size(); ++i)
    //    {
    //        ExplicitSystem& u_sys = es.add_system<ExplicitSystem>("u"+std::to_string(i));
    //        u_sys.add_variable("u"+std::to_string(i));
    //        u_sys.init();
    //    }

        bool add_LAA_and_FO = data("add_laa_fo", true);
        typedef BeatIt::Poisson Poisson;
        typedef std::shared_ptr<Poisson> PoissonPtr;
        std::vector<PoissonPtr> pv;
        //BeatIt::Poisson poisson(es, "Poisson");

        for (int i = 0; i < input_files.size(); i++)
        {
            std::string pname = "u" + std::to_string(i);
            std::cout << "Solving problem: " << pname << std::endl;
            pv.emplace_back( new  Poisson(es, pname));

            std::cout << "Calling setup: ..." << std::flush;
            GetPot data_i(input_files[i]);
            pv[i]->setup(data_i, "poisson");
            std::cout << " Done!" << std::endl;
            std::cout << "Calling assemble system: ..." << std::flush;
            pv[i]->assemble_system();
            std::cout << " Done!" << std::endl;
            std::cout << "Calling solve system: ..." << std::flush;
            pv[i]->solve_system();
            std::cout << " Done!" << std::endl;
            std::cout << " Done!" << std::endl;
            es.get_system<LinearImplicitSystem>("u"+std::to_string(i)).update();

            if(i == 0 && add_LAA_and_FO)
            {
                laa_threshold = data_i("laa_threshold", 0.22);
                fo_threshold = data_i("fo_threshold", 0.78);
                for(auto & node : mesh.local_node_ptr_range() )
                {
                    unsigned int dn = node->dof_number(es.get_system<ExplicitSystem>("u"+std::to_string(i)).number(), 0, 0);
                    double u1 = es.get_system<ExplicitSystem>("u"+std::to_string(i)).solution->el(dn);
                    if(u1 < laa_threshold)
                    {
                        mesh.get_boundary_info().add_node(node, 123);
                    }
                    else if(u1 > fo_threshold)
                    {
                        mesh.get_boundary_info().add_node(node, 321);
                    }
                }
            }
        }


        // DEFINE FIBER SYSTEMS:
        FiberSystem& f_sys = es.add_system<FiberSystem>("fibers");
        f_sys.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_vector("backup");
        f_sys.init();
        auto& f_v = f_sys.solution;

        FiberSystem& s_sys = es.add_system<FiberSystem>("sheets");
        s_sys.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.init();
        auto& s_v = s_sys.solution;

        FiberSystem& n_sys = es.add_system<FiberSystem>("xfibers");
        n_sys.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.init();
        auto& n_v = n_sys.solution;

        // Define auxiliary vectors used in the computations of the fibers
        double f[3];
        double s[3];
        double n[3];

        std::vector<double> u(N);
        std::vector<libMesh::Gradient> du(N);


        libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        std::cout << "\nGetting dof map:  ... " << std::flush;
        const libMesh::DofMap & dof_map = f_sys.get_dof_map();
        std::cout << " done \n " << std::flush;

        std::vector<libMesh::dof_id_type> dof_indices;
        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;
            dof_map.dof_indices(elem, dof_indices);
            auto elID = elem->id();
            Point centroid = elem->centroid();
            auto blockID = mesh.elem_ptr(elID)->subdomain_id();


            for (int k = 0; k < N; ++k)
            {
                u[k] = es.get_system<LinearImplicitSystem>("u"+std::to_string(k)).point_value(0, centroid, elem);
                du[k] = es.get_system<LinearImplicitSystem>("u"+std::to_string(k)).point_gradient(0, centroid, elem);
            }

            if (atlas_meshfile == "NONE")
                evaluate(f, s, n, u, du, blockID, centroid, data, dim);
            else
                evaluate_atlas(f, s, n, u, du, centroid, *atlas_es_ptr, *locator_ptr);


            // change subdomain ID Pulmonary Veins and Vena Cava:

            f_v->set(dof_indices[0], f[0]);
            f_v->set(dof_indices[1], f[1]);
            f_v->set(dof_indices[2], f[2]);

            s_v->set(dof_indices[0], s[0]);
            s_v->set(dof_indices[1], s[1]);
            s_v->set(dof_indices[2], s[2]);

            n_v->set(dof_indices[0], n[0]);
            n_v->set(dof_indices[1], n[1]);
            n_v->set(dof_indices[2], n[2]);


        }
        f_v->close();
        s_v->close();
        n_v->close();

        // smooth fibers
        bool do_smoothing = data("smoothing", false);
        if(do_smoothing)
        {
    //        f_sys.update();
            std::cout << "Closing 1" << std::endl;
            f_sys.get_vector("backup").close();
            std::cout << "Closing 2" << std::endl;
            f_v->close();
            std::cout << "Copying fibers " << std::endl;
            f_sys.get_vector("backup") = *f_v;
            std::cout << "Smoothing fibers " << std::endl;
    //        double D = data("D", 1.0);
    //        double DD = data("DD", 1.0);
            double f_threshold = data("f_threshold", 0.8);

    //        LinearImplicitSystem& fs2_sys = es.add_system<LinearImplicitSystem>("Smoothing_fibers");
    //        fs2_sys.add_variable("fx", libMesh::FIRST, libMesh::LAGRANGE);
    //        fs2_sys.add_variable("fy", libMesh::FIRST, libMesh::LAGRANGE);
    //        fs2_sys.add_variable("fz", libMesh::FIRST, libMesh::LAGRANGE);
    ////        fs_sys.add_variable("angle", libMesh::FIRST, libMesh::LAGRANGE);
    ////        fs_sys.add_variable("fx", libMesh::CONSTANT, libMesh::MONOMIAL);
    ////        fs_sys.add_variable("fy", libMesh::CONSTANT, libMesh::MONOMIAL);
    ////        fs_sys.add_variable("fz", libMesh::CONSTANT, libMesh::MONOMIAL);
    //        fs2_sys.add_matrix("mass");
    //        fs2_sys.add_matrix("stiffness");
    //        fs2_sys.add_matrix("Dstiffness");
    //        fs2_sys.init();
    //        std::cout << "Assemble smoother fibers " << std::endl;
    //        fs2_sys.assemble_before_solve = false;
    //        fs2_sys.zero_out_matrix_and_rhs = false;
    //        std::cout << "solve 1" << std::endl;
    //        assemble_smoother(es, D, DD, 1);
    //        fs2_sys.solve();
    //        fs2_sys.update();
    //
    //        ExplicitSystem& fs1_sys = es.add_system<ExplicitSystem>("Smoothing_fibers1");
    //        fs1_sys.add_variable("ffx", libMesh::FIRST, libMesh::LAGRANGE);
    //        fs1_sys.add_variable("ffy", libMesh::FIRST, libMesh::LAGRANGE);
    //        fs1_sys.add_variable("ffz", libMesh::FIRST, libMesh::LAGRANGE);
    //        fs1_sys.init();
    //        fs1_sys.solution->close();
    //        *fs1_sys.solution = *fs2_sys.solution;
    //        fs1_sys.update();
    //
    //        std::cout << "solve 2" << std::endl;
    //        assemble_smoother(es, D, DD, 2);
    //        fs2_sys.solve();
    //        fs2_sys.update();

            std::cout << "Reassign fibers" << std::endl;

            //BeatIt::Util::normalize(*fs_sys.solution);
            const libMesh::DofMap & dof_map_f = f_sys.get_dof_map();
            std::vector<libMesh::dof_id_type> dof_indices_f;
    //        auto& fs1_v = fs1_sys.solution;
    //        auto& fs2_v = fs2_sys.solution;

            //reassign the fibers
    //        el = mesh.active_local_elements_begin();
    //        for (; el != end_el; ++el)
    //        {
    //            const libMesh::Elem * elem = *el;
    //            dof_map.dof_indices(elem, dof_indices);
    //            dof_map_f.dof_indices(elem, dof_indices_f);
    //            auto elID = elem->id();
    //            auto blockID = elem->subdomain_id();
    //            Point centroid = elem->centroid();
    //
    //            libMesh::Gradient my_f;
    //            for(int j = 0; j < f_sys.n_vars(); ++j) my_f(j) = f_sys.solution->el(dof_indices_f[j]);
    //            std::vector<libMesh::Gradient> neighbor_f(elem->n_neighbors());
    //            libMesh::Gradient new_f;
    //
    //            for (unsigned int side=0; side<elem->n_sides(); side++)
    //            {
    //                int counter = 0;
    //                if (elem->neighbor_ptr(side) != libmesh_nullptr)
    //                {
    //                    dof_map_f.dof_indices(elem->neighbor_ptr(side), dof_indices_f);
    //                    for(int j = 0; j < f_sys.n_vars(); ++j) neighbor_f[counter](j) = f_sys.solution->el(dof_indices_f[j]);
    //                    counter ++ ;
    //                }
    //            }
    //
    //
    //            new_f = neighbor_f[0];
    //            new_f = new_f.unit();
    //            if( elem->n_neighbors() > 1)
    //            {
    //                double cos_angle = neighbor_f[1] * new_f;
    //                if( cos_angle > 0 ) new_f = new_f + neighbor_f[1];
    //                else new_f = new_f - neighbor_f[1];
    //                new_f = new_f.unit();
    //            }
    //            if( elem->n_neighbors() > 2)
    //            {
    //                double cos_angle = neighbor_f[2] * new_f;
    //                if( cos_angle > 0 ) new_f = new_f + neighbor_f[2];
    //                else new_f = new_f - neighbor_f[2];
    //                new_f = new_f.unit();
    //            }

    //                    auto neighbor_blockID = elem->neighbor_ptr(side)->subdomain_id();
    //                    if(neighbor_blockID != blockID)
    //                    {
    //                        dof_map_f.dof_indices(elem->neighbor_ptr(side), dof_indices_f);
    //                        for(int j = 0; j < f_sys.n_vars(); ++j) neighbor_f(j) = f_sys.solution->el(dof_indices_f[j]);
    //
    //                        double cos_angle = neighbor_f * my_f;
    //                        if( cos_angle > 0 ) new_f = my_f + neighbor_f;
    //                        else new_f = my_f - neighbor_f;
    //
    //                        new_f = new_f.unit();
    //                    }
    //                    break;
    //                }
    //                else
    //                {
    //                    new_f = my_f.unit();
    //                }
    //            fs_v->set(dof_indices[0], new_f(0));
    //            fs_v->set(dof_indices[1], new_f(1));
    //            fs_v->set(dof_indices[2], new_f(2));
    //
    //
    //        }
    //        fs_v->close();
    //        fs_sys.update();

            //reassign the fibers
            el = mesh.active_local_elements_begin();
            for (; el != end_el; ++el)
            {
                const libMesh::Elem * elem = *el;
                dof_map.dof_indices(elem, dof_indices);
                dof_map_f.dof_indices(elem, dof_indices_f);
                auto elID = elem->id();
                Point centroid = elem->centroid();
                auto blockID = elem->subdomain_id();

                libMesh::Gradient old_f;
                libMesh::Gradient new_f;
                for(int j = 0; j < f_sys.n_vars(); ++j) old_f(j) = f_sys.get_vector("backup").el(dof_indices_f[j]);

                double old_angle = std::asin( old_f(1) );

    //            libMesh::Gradient f1; //ff
    //            f1(0) = fs1_sys.point_value(0, centroid, elem);
    //            f1(1) = fs1_sys.point_value(1, centroid, elem);
    //            f1(2) = fs1_sys.point_value(2, centroid, elem);
    //
    //            libMesh::Gradient f2; //f
    //            f2(0) = fs2_sys.point_value(0, centroid, elem);
    //            f2(1) = fs2_sys.point_value(1, centroid, elem);
    //            f2(2) = fs2_sys.point_value(2, centroid, elem);
    //
    //            double magf = f2.norm();
    //            double w = 1.0;
    //            if(magf > f_threshold) w = 0.0;
    //
    //            new_f = f1 * w + (1 - w) * f2;



    //            f[0] = fs_sys.point_value(0, centroid, elem);
    //            f[1] = fs_sys.point_value(1, centroid, elem);
    //            f[2] = fs_sys.point_value(2, centroid, elem);

    //            double mapped_angle = fs_sys.point_value(3, centroid, elem);
    //            double new_angle = mapped_angle - 0.5 * M_PI;
    //            if(old_angle > 0)
    //               new_angle *= -1;
    //            f[0] = std::cos(new_angle);
    //            f[1] = std::sin(new_angle);

                new_f = old_f;
                int  counter = 0;
                for (unsigned int side=0; side<elem->n_sides(); side++)
                {
                    if (elem->neighbor_ptr(side) != libmesh_nullptr)
                    {
                        auto nblockID = elem->neighbor_ptr(side)->subdomain_id();
                        if(nblockID != blockID)
                        {
                            counter ++ ;
                            libMesh::Gradient nf;
                            dof_map_f.dof_indices(elem->neighbor_ptr(side), dof_indices_f);
                            for(int j = 0; j < f_sys.n_vars(); ++j) nf(j) = f_sys.solution->el(dof_indices_f[j]);
                            double cdotn = nf * new_f;
                            if(cdotn >= f_threshold)
                            {
                                new_f = nf + new_f;
                                new_f = new_f.unit();
                            }
                            else
                            {
                                new_f = nf - new_f;
                                new_f = new_f.unit();
                            }
                            //std::cout << counter << std::endl;
                        }
                    }
                }

                new_f = new_f.unit();

                f[0] = new_f(0);
                f[1] = new_f(1);
                f[2] = new_f(2);

                s[0] = 0;
                s[1] = 0;
                s[2] = 1;

                BeatIt::Util::normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);
                cross(s, f, n);
                BeatIt::Util::normalize(n[0], n[1], n[2], 0.0, 0.0, 1.0);


                f_v->set(dof_indices[0], f[0]);
                f_v->set(dof_indices[1], f[1]);
                f_v->set(dof_indices[2], f[2]);
                n_v->set(dof_indices[0], n[0]);
                n_v->set(dof_indices[1], n[1]);
                n_v->set(dof_indices[2], n[2]);
            }

        }
    }

//    {
//        std::cout << "Calling exporter: ..." << ". " << std::flush;
//        ExodusII_IO exporter(mesh);
//        // the solution of the poisson problems
//        bool export_only_fibers = data("only_fibers", false);
//        if (export_only_fibers)
//        {
//            std::vector < std::string > varname(9);
//            varname[0] = "fibersx";
//            varname[1] = "fibersy";
//            varname[2] = "fibersz";
//            varname[3] = "sheetsx";
//            varname[4] = "sheetsy";
//            varname[5] = "sheetsz";
//            varname[6] = "xfibersx";
//            varname[7] = "xfibersy";
//            varname[8] = "xfibersz";
//            exporter.set_output_variables(varname);
//        }
//
//        std::string output = data("output", "test2D.e");
//        exporter.write_equation_systems(output, es);
//        exporter.write_element_data(es);
//    }
//    return 0;
    std::cout << "ATLAS: " << atlas_meshfile << std::endl;
    std::cout << "Create EP ..." << std::endl;

    // Read the model from the input file
    std::string model = data_ep("model", "monowave");
    std::cout << "Create Electrophysiology model ..." << std::endl;
    // Create the model
    BeatIt::ElectroSolver* solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, es);
    // Setup the EP model using the input file
    std::cout << "Calling setup..." << std::endl;
    solver->setup(data_ep, model);
    // Initialize systems
    std::cout << "Calling init ..." << std::endl;
    // Set up initial conditions at time
    solver->init(datatime.M_startTime);

    if(read_fibers)
    {
        // First show the elemental vcariables that can be imported
        auto elemental_variables = importer_ptr->get_elem_var_names();
        for (auto && var : elemental_variables) std::cout << var << std::endl;
        solver->read_fibers(*importer_ptr,1);
    }

    ////////////////////////////////////////////
    // Modify conductivities
    ////////////////////////////////////////////
    bool create_fibrosis = data("create_fibrosis", false);
    std::string fibrosis_model_type = data("fibrosis_model", "none");
    auto it = fibrosis_model_map.find(fibrosis_model_type);
    FibrosisModel fibrosis_model = (it != fibrosis_model_map.end() ) ? it->second : FibrosisModel::NONE;

    if(FibrosisModel::NONE != fibrosis_model)
    {
        std::cout << "Making NOISE! " << std::endl;
        Noise noise;
        std::cout << "setup NOISE! " << std::endl;
        noise.setup(data, "");
        std::cout << "generate NOISE! " << std::endl;
        noise.generate_noise_field(es);

        libMesh::ExplicitSystem& sigma_sys =  es.get_system<ExplicitSystem>("conductivity");
        libMesh::ExplicitSystem& noise_sys = es.get_system<libMesh::ExplicitSystem> ("Noise");
        double fibrosis_threshold1 = data("fib_t1", 0.4);
        double fibrosis_threshold2 = data("fib_t2", 0.55);
        int field = data("fibrosis_field", 0);
        //libMesh::ExplicitSystem& field_system = es.get_system<libMesh::ExplicitSystem> ("u"+std::to_string(field));
        //libMesh::ExplicitSystem& u0_system = es.get_system<libMesh::ExplicitSystem> ("u"+std::to_string(0));

        const libMesh::DofMap & noise_dof_map = noise_sys.get_dof_map();
        const libMesh::DofMap & sigma_dof_map = sigma_sys.get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices_noise;
        std::vector<libMesh::dof_id_type> dof_indices_sigma;

        for(auto & elem : mesh.active_local_element_ptr_range() )
        {
            noise_dof_map.dof_indices(elem,dof_indices_noise);
            sigma_dof_map.dof_indices(elem,dof_indices_sigma);
            double u = 1.0; //field_system.point_value(0, elem->centroid(), elem);
            double u0 = 1.0; //u0_system.point_value(0, elem->centroid(), elem);
            auto blockID = elem->subdomain_id();


            //
            double w = 0;
            if(u > fibrosis_threshold2) w = 1;
            else if(u > fibrosis_threshold1)
            {
                w = (u-fibrosis_threshold1)/std::abs(fibrosis_threshold2-fibrosis_threshold1);
            }
            w = 1;
            libMesh::Gradient sigma;
            for(int j = 0; j < 3; ++j) sigma(j) = sigma_sys.solution->el(dof_indices_sigma[j]);

            // Fossa ovalis
            if(2 == blockID)
            {
                sigma(0) /= 16;
                sigma(1) = sigma(0);
                sigma(2) = sigma(0);
            }
            // everything else
            else
            {
                double noise = noise_sys.solution->el(dof_indices_noise[0]);
                double tn = ( 1.0 - w * noise);

                switch(fibrosis_model)
                {
                    case FibrosisModel::KRUEGER:
                    {
                        int percentage_severe = data("severe", 0);
                        int percentage_moderate = data("moderate", 0);
                        int percentage_mild = data("mild", 0);
                        double noise_threshold = 0.;
                        if(10 == percentage_severe ) noise_threshold = 0.31;
                        if(20 == percentage_severe ) noise_threshold = 0.37;
                        if(30 == percentage_severe ) noise_threshold = 0.42;
                        if(40 == percentage_severe ) noise_threshold = 0.47;
                        if(50 == percentage_severe ) noise_threshold = 0.5;
                        if(60 == percentage_severe ) noise_threshold = 0.55;
                        if(70 == percentage_severe ) noise_threshold = 0.59;
                        if(80 == percentage_severe ) noise_threshold = 0.64;
                        // we are in the severe region
                        if( tn < noise_threshold)
                        {
                            // Anisotropy ratio = 8:1
                            sigma(0) *= 0.16;
                            sigma(1) = sigma(0)/64.0;
                            sigma(2) = sigma(1);
                        }
                        else
                        {
                            noise_threshold = 0.0;
                            if(10 == percentage_moderate ) noise_threshold = 0.31;
                            if(20 == percentage_moderate ) noise_threshold = 0.37;
                            if(30 == percentage_moderate ) noise_threshold = 0.42;
                            if(40 == percentage_moderate ) noise_threshold = 0.47;
                            if(50 == percentage_moderate ) noise_threshold = 0.5;
                            if(60 == percentage_moderate ) noise_threshold = 0.55;
                            if(70 == percentage_moderate ) noise_threshold = 0.59;
                            if(80 == percentage_moderate ) noise_threshold = 0.64;
                            // moderate
                            if( tn < noise_threshold)
                            {
                                // Anisotropy ratio = 6:1
                                sigma(0) *= 0.36;
                                sigma(1) = sigma(0)/36.0;
                                sigma(2) = sigma(1);
                            }
                            else
                            {
                                noise_threshold = 0.0;
                                if(10 == percentage_mild ) noise_threshold = 0.31;
                                if(20 == percentage_mild ) noise_threshold = 0.37;
                                if(30 == percentage_mild ) noise_threshold = 0.42;
                                if(40 == percentage_mild ) noise_threshold = 0.47;
                                if(50 == percentage_mild ) noise_threshold = 0.5;
                                if(60 == percentage_mild ) noise_threshold = 0.55;
                                if(70 == percentage_mild ) noise_threshold = 0.59;
                                if(80 == percentage_mild ) noise_threshold = 0.64;
                                // mild
                                if( tn < noise_threshold)
                                {
                                    // Anisotropy ratio = 4:1
                                    sigma(0) *= 0.64;
                                    sigma(1) = sigma(0)/16.0;
                                    sigma(2) = sigma(1);
                                }
                                else
                                {
                                    // Anisotropy ratio = 1.9:1
                                    // approximate to 2:1 in conduction velocity
                                    sigma(1) = sigma(0)/4.0;
                                    sigma(2) = sigma(1);
                                }
                            }

                        }
                        break;
                    }
                    case FibrosisModel::PEZZUTO:
                    {
                        int percentage = data("noise_percent", 0);
                        double noise_threshold = 0.;
                        if(20 == percentage ) noise_threshold = 0.37;
                        if(30 == percentage ) noise_threshold = 0.42;
                        if(40 == percentage ) noise_threshold = 0.47;
                        if(50 == percentage ) noise_threshold = 0.5;
                        if(60 == percentage ) noise_threshold = 0.55;
                        if(70 == percentage ) noise_threshold = 0.59;
                        if(80 == percentage ) noise_threshold = 0.64;
                        if( tn < noise_threshold)
                        {
                            sigma(1) *= 0.0;
                            sigma(2) *= 0.0;
                        }
                        break;
                    }
                    case FibrosisModel::EPI_ENDO_DISSOCIATION:
                    {
                        int percentage = data("noise_percent", 0);
                        double noise_threshold = 0.;
                        if(20 == percentage ) noise_threshold = 0.37;
                        if(30 == percentage ) noise_threshold = 0.42;
                        if(40 == percentage ) noise_threshold = 0.47;
                        if(50 == percentage ) noise_threshold = 0.5;
                        if(60 == percentage ) noise_threshold = 0.55;
                        if(70 == percentage ) noise_threshold = 0.59;
                        if(80 == percentage ) noise_threshold = 0.64;
                        if( noise < noise_threshold)
                        {
                            sigma(1) *= 0.0;
                        }
                        break;
                    }

                    case FibrosisModel::OURS:
                    default:
                    {
                        sigma(1) *= tn * tn;
                        sigma(2) *= tn * tn;
                        break;
                    }
                }
            }
            for(int j = 0; j < 3; ++j) sigma_sys.solution->set(dof_indices_sigma[j], sigma(j) );

        }

    }

    std::cout << "Calling exporter: ..." << ". " << std::flush;
    ExodusII_IO exporter(mesh);
    // the solution of the poisson problems
    bool export_only_fibers = data("only_fibers", false);
    if (export_only_fibers)
    {
        std::vector < std::string > varname(9);
        varname[0] = "fibersx";
        varname[1] = "fibersy";
        varname[2] = "fibersz";
        varname[3] = "sheetsx";
        varname[4] = "sheetsy";
        varname[5] = "sheetsz";
        varname[6] = "xfibersx";
        varname[7] = "xfibersy";
        varname[8] = "xfibersz";
        exporter.set_output_variables(varname);
    }

    std::string output = data("output", "test2D.e");
    exporter.write_equation_systems(output, es);
    exporter.write_element_data(es);
    ////////////////////////////////////////////
    // keep going with the EP
    ////////////////////////////////////////////


    std::cout << "Assembling matrices" << std::endl;
    solver->assemble_matrices(datatime.M_dt);
    // output file counter
    int save_iter = 0;
    // Export initial condition at time
  //    solver->save_exo_timestep(save_iter, datatime.M_time);
    solver->save_potential(save_iter, datatime.M_startTime);
    // Parameters to save the activation times
    // A node is activated if the transmembrane potential > threshold
    double threshold = data_ep("threshold", -10.0);
    // We export the activation times every at_save_iter iterations
    int at_save_iter = data_ep("at_save_iter", 25);

    // Old parameters that define the method, not to be changed
    // TODO: clean this part
    std::string system_mass = data_ep(model + "/diffusion_mass", "mass");
    std::string iion_mass = data_ep(model + "/reaction_mass", "lumped_mass");
    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;

    // Start loop in time
    std::cout << "Time loop starts:" << std::endl;
    // Control the time loop using the TimeData object
    for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
    {
        // We are doing a new iteration
        // let's update first the information in the TimeData object
        datatime.advance();
        // Ouput to screen the current time and the current iteration
        std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter << std::endl;
        // Advance the solution in time: u_n <- u_n+1
        solver->advance();
        // Solve ionic model and evaluate ionic currents
        solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        // Solve monodomain model
        solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        // Update the activation times
        solver->update_activation_time(datatime.M_time, threshold);
        //Export the solution if at the right timestep
        if (0 == datatime.M_iter % datatime.M_saveIter)
        {
            // update output file counter
            save_iter++;
            // export current solution
            solver->save_potential(save_iter, datatime.M_time);
  //            solver->save_exo_timestep(save_iter, datatime.M_time);
        }
        // export the activation times if at the corresponding timestep
        if (0 == datatime.M_iter % (at_save_iter * datatime.M_saveIter))
        {
            // export activation times
            solver->save_activation_times_nemesis(save_iter);
        }

    }


    std::cout << "laa_threshold: " << laa_threshold << std::endl;
    std::cout << "fo_threshold: " << fo_threshold << std::endl;
    std::cout << "lpv_antra_threshold: " << lpv_antra_threshold << std::endl;
    std::cout << "rpv_antra_threshold: " << rpv_antra_threshold << std::endl;
    std::cout << "lpw_threshold: " << lpw_threshold << std::endl;
    std::cout << "rpw_threshold: " << rpw_threshold << std::endl;
    // // export activation times again
    solver->save_activation_times_nemesis(save_iter);
    //solver->save_activation_times(save_iter);
    // delete solver before ending the simulation
    // avoiding memory leaks
    delete solver;



    return 0;
}

void evaluate(double f[], double s[], double n[], std::vector<double>& phi, std::vector<libMesh::Gradient>& dphi, int blockID, Point& x, GetPot& data, int dim)
{
    auto normalize = [](double& x, double& y, double& z,
            double X, double Y, double Z)
    {
        double norm = std::sqrt( x * x + y * y + z * z);
        if(norm >= 1e-12 )
        {
            x /= norm;
            y /= norm;
            z /= norm;
        }
        else
        {
            x = X;
            y = Y;
            z = Z;
        }
    };


    std::string endo_epi = data("layer", "NONE");
    enum class Layer { NONE, Endo, Epi, Milano };
    double u0 = phi[0];
    Layer layer = Layer::NONE;
    if( endo_epi == "epi" ) layer = Layer::Epi;
    if( endo_epi == "endo" ) layer = Layer::Endo;
    if( endo_epi == "milano" ) layer = Layer::Milano;

    enum class NPVs { THREE, FOUR, FIVE, SIX };
    NPVs npvs = NPVs::FOUR;
    int n_pvs = data("n_pvs", 4);
    if( 3 == n_pvs  ) npvs = NPVs::THREE;
    if( 5 == n_pvs  ) npvs = NPVs::FIVE;
    if( 6 == n_pvs  ) npvs = NPVs::SIX;

    s[0] = 0;  s[1] = 0;  s[2] = 1;
    if(dim == 3)
    {
        auto m = dphi.size();
        s[0] = dphi[m-1](0);  s[1] = dphi[m-1](1);  s[m-1] = dphi[m-1](2);
        normalize(s[0],s[1],s[2],  1.0,  0.0 , 0.0);
        layer = Layer::NONE;
    }

//    for(auto & ui : phi) std::cout << ui << " " << std::flush;

    switch(layer)
    {
        case Layer::Milano:
        {
            double tau_mv = data("tau_mv", 0.65);
            double tau_lpv = data("tau_lpv", 0.65);
            double tau_rpv = data("tau_rpv", 0.1);
            if( phi[0] > tau_mv )
            {
                n[0] = dphi[0](0);  n[1] = dphi[0](1);  n[2] = 0.0;
                normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                cross(s, n, f);
            }
            else
            {
                if( phi[1] > tau_lpv || phi[1] < tau_rpv )
                {
                    n[0] = dphi[0](0);  n[1] = dphi[0](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);

                }
                else
                {
                    n[0] = dphi[2](0);  n[1] = dphi[2](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);
                }
            }
            break;
        }
        case Layer::Endo:
        {
            // Left PVs + // Right PVs
            if(phi[4] < lpv_antra_threshold || phi[4] > rpv_antra_threshold)
            {
//                std::cout << "PVs: "  << std::endl;
//                std::cout << "lpv_antra_threshold: " << lpv_antra_threshold << std::endl;
//                std::cout << "rpv_antra_threshold: " << rpv_antra_threshold << std::endl;

                n[0] = dphi[2](0);  n[1] = dphi[2](1);  n[2] = 0.0;
                normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                cross(s, n, f);
                if( NPVs::THREE == npvs && phi[4] < lpv_antra_threshold)
                {
                    n[0] = dphi[4](0);  n[1] = dphi[4](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);
                }
            }
            else
            {
                // LAA
                if(phi[0] < laa_threshold)
                {
                    std::cout << "LAA"  << std::endl;

                    n[0] = dphi[0](0);  n[1] = dphi[0](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);
                }
                // FO e

                else if(phi[0] > fo_threshold)
                {
                    std::cout << "FO"  << std::endl;
                    f[0] = dphi[4](0);  f[1] = dphi[4](1);  n[2] = 0.0;
                    normalize(f[0],f[1],f[2],  1.0,  0.0 , 0.0);
                    cross(s, f, n);
                }
                else
                {
                        // Posterior Wall
                    //(u2_phi>0.5)*(u1_phi>0.2)*(u1_phi<0.8)*(u3_phi>0.1)
                    //+ (u2_phi<0.5)*(u1_phi>0.2)*(u1_phi<0.8)*(u3_phi>0.7);
                        bool pw = (phi[2]>0.5)*(phi[1]>lpw_threshold)*(phi[1]<rpw_threshold)*(phi[3]>0.1)
                                + (phi[2]<0.5)*(phi[1]>lpw_threshold)*(phi[1]<rpw_threshold)*(phi[3]>0.7);
                        // (u2_phi<0.5)*(u1_phi>0.2)*(u1_phi<0.85)*(u3_phi>0.1)*(u3_phi<0.7)
                        bool aw = (phi[2]<0.5)*(phi[1]>0.2)*(phi[1]<0.85)*(phi[3]>0.1)*(phi[3]<0.7);
                        if( pw == true)
                        {
                            std::cout << "PW"  << std::endl;

                            n[0] = dphi[1](0);  n[1] = dphi[1](1);  n[2] = 0.0;
                            normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                            cross(s, n, f);
                        }
                        else if (true == aw)
                        {
                            std::cout << "AW"  << std::endl;

                            n[0] = dphi[5](0);  n[1] = dphi[5](1);  n[2] = 0.0;
                            normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                            cross(s, n, f);
                        }
                        else
                        {
                            std::cout << "ELSE"  << std::endl;

                            n[0] = dphi[3](0);  n[1] = dphi[3](1);  n[2] = 0.0;
                            normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                            cross(s, n, f);
                        }
                }
            }
            break;
        }
        case Layer::Epi:
        {
            // Left PVs
            if(phi[4] < lpv_antra_threshold )
            {
                n[0] = dphi[2](0);  n[1] = dphi[2](1);  n[2] = 0.0;
                normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                cross(s, n, f);
                if( NPVs::THREE == npvs && phi[4] < lpv_antra_threshold)
                {
                    n[0] = dphi[4](0);  n[1] = dphi[4](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);
                }
            }
            // Right PVs
            else if(phi[4] > rpv_antra_threshold )
            {
                // Right Superior PVs
                if(phi[2] < 0.5 )
                {
                    n[0] = dphi[2](0);  n[1] = dphi[2](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);
                }
                else
                {
                    // Right Inferior PVs
                    f[0] = dphi[4](0);  f[1] = dphi[4](1);  f[2] = 0.0;
                    normalize(f[0],f[1],f[2],  1.0,  0.0 , 0.0);
                    cross(s, f, n);
                }
            }
            else
            {
                // LAA
                if(phi[0] < laa_threshold)
                {
                    n[0] = dphi[0](0);  n[1] = dphi[0](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);
                }
                else
                {
                        // Posterior Wall
                        bool pw = (phi[2]>0.5)*(phi[1]>lpw_threshold)*(phi[1]<2)*(phi[3]>0.1)
                                 + (phi[2]<0.5)*(phi[1]>lpw_threshold)*(phi[1]<rpw_threshold)*(phi[3]>0.35);
                        if( pw == true)
                        {
                            // septum
                            if(phi[1] > rpw_threshold )
                            {
                                    f[0] = dphi[4](0);  f[1] = dphi[4](1);  f[2] = 0.0;
                                    normalize(f[0],f[1],f[2],  1.0, -1.0 , 0.0);
                                    normalize(f[0],f[1],f[2],  1.0,  0.0 , 0.0);
                                    cross(s, f, n);
                            }
                            else
                            {
                                //  (u3_phi<0.6)*(u2_phi<0.2)*(u1_phi<0.5)
                                //bool aw = (phi[3]<0.6)*(phi[2]<0.2)*(phi[1]<0.5);
                                bool aw = (phi[3]<0.6)*(phi[1]>lpw_threshold)*(phi[1]<0.5)*(phi[2]<0.5);
                                // anterior wall left part
                                if(aw == true)
                                {
//                                    double w = 90/(0.5-0.2) * (phi[1] - 0.5);
//                                    f[0] = dphi[1](1);  f[1] = dphi[1](0);  f[2] = 0.0;
//                                    normalize(f[0],f[1],f[2],  1.0,  0.0 , 0.0);
//                                    cross(s, f, n);
                                    n[0] = dphi[5](0);  n[1] = dphi[5](1);  n[2] = 0.0;
                                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                                    cross(s, n, f);


//                                    double theta = M_PI / 180.0*w + M_PI / 2;
//                                    f[0] = std::cos(theta)*n[0] - std::sin(theta)*n[1];
//                                    f[1] = std::sin(theta)*n[0] + std::cos(theta)*n[1];
//                                    f[2] = 0.0;
//                                    normalize(f[0],f[1],f[2],  1.0,  0.0 , 0.0);
//                                    cross(s, f, n);
                                }
                                // posterior wall front
                                else if (phi[2] < 0.5)
                                {
                                    n[0] = dphi[4](0);  n[1] = dphi[4](1);  n[2] = 0.0;
                                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                                    cross(s, n, f);

                                }
                                // posterior wall back
                                else
                                {
                                    n[0] = dphi[1](0);  n[1] = dphi[1](1);  n[2] = 0.0;
                                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                                    cross(s, n, f);
                                }
                            }
                        }
                        else
                        {
                            n[0] = dphi[3](0);  n[1] = dphi[3](1);  n[2] = 0.0;
                            normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                            cross(s, n, f);
                        }
                }
            }
            break;
        }
        case Layer::NONE:
        default:
        {
            int el = blockID-1;
            auto df = dphi[el];
            df(2) = 0;
            df = df.unit();
            f[0] = df(0);  f[1] = df(1);  f[2] = 0;

            if( f[0] < 0 )
            {
                f[0] *= -1;
                f[1] *= -1;
            }
            cross(s, f, n);

            break;
        }
    }
}


void assemble_smoother(EquationSystems& es, double D, double DD, int sector)
{
    const libMesh::MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

     LinearImplicitSystem& system  =  es.get_system<LinearImplicitSystem>("Smoothing_fibers");
     int n_vars   =  system.n_vars();
     ExplicitSystem& f_sys  =  es.get_system<ExplicitSystem>("fibers");
     ExplicitSystem& s_sys  =  es.get_system<ExplicitSystem>("sheets");
    system.matrix->zero();
    system.rhs->zero();
         // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const libMesh::DofMap & dof_map = system.get_dof_map();
    const libMesh::DofMap & dof_map_f = f_sys.get_dof_map();
         // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    libMesh::FEType fe_type = dof_map.variable_type(0);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr<libMesh::FEBase> fe_qp(libMesh::FEBase::build(dim, fe_type));
     // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule_stiffness(dim, libMesh::FOURTH);
     // Tell the finite element object to use our quadrature rule.
    fe_qp->attach_quadrature_rule(&qrule_stiffness);

     // The element Jacobian * quadrature weight at each integration point.
    const std::vector<libMesh::Real> & JxW_qp = fe_qp->get_JxW();
     // The physical XY locations of the quadrature points on the element.
    // These might be useful for evaluating spatially varying material
    // properties at the quadrature points.
    const std::vector<libMesh::Point> & q_point_qp = fe_qp->get_xyz();

    const std::vector<std::vector<libMesh::Real> > & phi_qp = fe_qp->get_phi();
    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<libMesh::RealGradient> > & dphi_qp = fe_qp->get_dphi();
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    libMesh::DenseMatrix<libMesh::Number> Se;
    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseMatrix<libMesh::Number> Me;
     libMesh::DenseVector<libMesh::Number> Fe;

          // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_f;
         // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    //
    // Element iterators are a nice way to iterate through all the
    // elements, or all the elements that have some property.  The
    // iterator el will iterate from the first to// the last element on
    // the local processor.  The iterator end_el tells us when to stop.
    // It is smart to make this one const so that we don't accidentally
    // mess it up!  In case users later modify this program to include
    // refinement, we will be safe and will only consider the active
    // elements; hence we use a variant of the active_elem_iterator.
    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();

    std::cout << "assemble smoother" << std::endl;
    double ftxyzi = 0.0;

    std::vector<double> myf(n_vars);

    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        dof_map_f.dof_indices(elem,dof_indices_f);
        libMesh::Gradient frhs;
        for(int j = 0; j < f_sys.n_vars(); ++j) frhs(j) = f_sys.solution->el(dof_indices_f[j]);
        libMesh::Gradient s;
        for(int j = 0; j < s_sys.n_vars(); ++j) s(j) = s_sys.solution->el(dof_indices_f[j]);

        myf[0] = frhs(0);
        myf[1] = frhs(1);
        if( sector != 1 && myf[1] < - 0 )
        {
            myf[0] = -frhs(0);
            myf[1] = -frhs(1);
        }
        myf[2] = frhs(2);

        double angle = std::asin( myf[1] );
        double map = 0.5 * M_PI - std::abs( angle );
        map = -2.0 / M_PI * angle * angle + M_PI * 0.5;
        map = 0.5*(1+std::cos(2*angle));
        map = std::sin(2*angle);
        if(myf.size() == 4) myf[3] = map;

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_qp->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        Se.resize(dof_indices.size(), dof_indices.size());
        Ke.resize(dof_indices.size(), dof_indices.size());
        Me.resize(dof_indices.size(), dof_indices.size());
         Fe.resize(dof_indices.size());


         int n_b = phi_qp.size();
       //  std::cout << "Assembling Poisson ... " << std::endl;
        for (unsigned int qp = 0; qp < qrule_stiffness.n_points(); qp++)
        {
            for (int idim = 0; idim < n_vars; ++idim)
            {
                for (unsigned int i = 0; i < n_b; i++)
                {
                    for (unsigned int j = 0; j <n_b; j++)
                    {
                        // stiffness term
                        Se(i + idim * n_b, j + idim * n_b) += JxW_qp[qp] *  dphi_qp[i][qp] * ( D * dphi_qp[j][qp]);
                        Se(i + idim * n_b, i + idim * n_b) += DD * JxW_qp[qp] * phi_qp[i][qp] * phi_qp[j][qp];
                        Ke(i + idim * n_b, j + idim * n_b) += JxW_qp[qp] *  dphi_qp[i][qp] * ( dphi_qp[j][qp]);
                        Me(i + idim * n_b, i + idim * n_b) +=  JxW_qp[qp] * phi_qp[i][qp] * phi_qp[j][qp];
                    }
                    Fe(i + idim * n_b) += JxW_qp[qp] * myf[idim] * phi_qp[i][qp];
                }

            }
        }
        system.matrix->add_matrix(Se, dof_indices);
        system.get_matrix("stiffness").add_matrix(Ke, dof_indices);
        system.get_matrix("mass").add_matrix(Me, dof_indices);
        system.rhs   ->add_vector(Fe, dof_indices);

    }
    system.rhs->close();
    system.matrix->close();
    system.get_matrix("mass").close();
    system.get_matrix("stiffness").close();
    system.get_matrix("Dstiffness").close();

//    system.matrix->zero();
//    system.get_matrix("Dstiffness").zero();
//    system.get_matrix("Dstiffness").add(D, system.get_matrix("stiffness"));
//    typedef libMesh::PetscMatrix<libMesh::Number> PetscMat;
//    PetscMat * mat1 = dynamic_cast<PetscMat*>(system.request_matrix("Dstiffness"));
//    PetscMat * mat2 = dynamic_cast<PetscMat*>(system.request_matrix("stiffness"));
//    PetscMat * mat3 = dynamic_cast<PetscMat*>(system.matrix);
//    Mat mat4 =  mat3->mat();
//    MatSetOption(mat3->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//    Mat mat5;
//    MatMatMult(mat1->mat(), mat2->mat(), MAT_INITIAL_MATRIX,1.0, &mat5);
//    PetscMat mat6(mat5, mesh.comm());
//    //double norm = mat6.linfty_norm();
//    mat3->add(1.0, mat6);
////    std::cout << " NORM: " << norm << std::endl;
//    system.matrix->add(1.0, system.get_matrix("mass"));
}



void evaluate_atlas(double f[], double s[], double n[], std::vector<double>& phi, std::vector<libMesh::Gradient>& dphi, Point& x, EquationSystems& es, PointLocatorTree& locator)
{
    typedef libMesh::ExplicitSystem FiberSystem;

    auto normalize = [](double& x, double& y, double& z,
            double X, double Y, double Z)
    {
        double norm = std::sqrt( x * x + y * y + z * z);
        if(norm >= 1e-12 )
        {
            x /= norm;
            y /= norm;
            z /= norm;
        }
        else
        {
            x = X;
            y = Y;
            z = Z;
        }
    };

    auto m = dphi.size();
    s[0] = dphi[m-1](0);  s[1] = dphi[m-1](1);  s[m-1] = dphi[m-1](2);
    normalize(s[0],s[1],s[2],  1.0,  0.0 , 0.0);
    // rescale phi 6
    double phi_min = 0.10508553192315474;
    double phi_max = 0.9999999997209817;
    double R = 5.0* ( phi[m-2] - phi_min ) / (phi_max-phi_min);
    double X = x(0);
    double Y = x(1);
    // (atan((coordsY)/(coordsX))*180/acos(-1)+90)*(coordsX>0)+
    // (atan((coordsY)/(coordsX))*180/acos(-1)+270)*(coordsX<0)
    libMesh::Point x2D;
    if(X < 1e-10 && X > -1e-10 )
    {
        x2D(0) = 0;
        x2D(1) = R;
        x2D(2) = 0;
    }
    else
    {
        double angle = std::atan(Y/X) * ( X >= 0 ) * ( Y >= 0 )
                     + ( 2.0 * M_PI + std::atan(Y/X) ) * ( X >= 0 ) * ( Y < 0 )
                     + ( std::atan(Y/X) + M_PI ) * (X < 0);

        //angle = std::atan(Y/X);
        double X2D =  R * std::cos(angle);
        double Y2D =  R * std::sin(angle);

        x2D(0) = X2D;
        x2D(1) = Y2D;
        x2D(2) = 0;
    }
    const Elem *atlas_element = locator(x2D);
//    std::cout << " Point: " << std::flush;
//    x.print();
//    std::cout << " maps to: " << std::flush;
//    x2D.print();
//    std::cout << ": fibers = " << std::flush;
    if(atlas_element)
    {
        FiberSystem& f_sys = es.get_system<FiberSystem>("fibers");
        const libMesh::DofMap & dof_map_f = f_sys.get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices_f;
        dof_map_f.dof_indices(atlas_element, dof_indices_f);
        libMesh::Gradient f2D;

        for(int j = 0; j < f_sys.n_vars(); ++j) f2D(j) = f_sys.get_vector("backup").el(dof_indices_f[j]);

        // ROTATE
        auto S = dphi[m-1];
        S = S.unit();
        auto S2D = libMesh::Point(0,0,1);
        S2D.unit();

        libMesh::Point u = S.cross(S2D);
        u = u.unit();
        // rotation angle
        double cos_tetha = S * S2D;
        double tetha = -std::acos(cos_tetha);
        double sin_tetha = std::sin(tetha);

        TensorValue<Number> R;
        R(0,0) =          cos_tetha + u(0) * u(0) * ( 1 - cos_tetha);
        R(0,1) = - u(2) * sin_tetha + u(0) * u(1) * ( 1 - cos_tetha);
        R(0,2) =   u(1) * sin_tetha + u(0) * u(2) * ( 1 - cos_tetha);
        R(1,0) =   u(2) * sin_tetha + u(1) * u(0) * ( 1 - cos_tetha);
        R(1,1) =          cos_tetha + u(1) * u(1) * ( 1 - cos_tetha);
        R(1,2) = - u(0) * sin_tetha + u(1) * u(2) * ( 1 - cos_tetha);
        R(2,0) = - u(1) * sin_tetha + u(2) * u(0) * ( 1 - cos_tetha);
        R(2,1) =   u(0) * sin_tetha + u(2) * u(1) * ( 1 - cos_tetha);
        R(2,2) =          cos_tetha + u(2) * u(2) * ( 1 - cos_tetha);
        libMesh::Point Rf =  R * f2D;


        f[0] = Rf(0);
        f[1] = Rf(1);
        f[2] = Rf(2);
//        f2D.print();
//        std::cout << std::endl;
        cross(s, f, n);
    }
    else
    {
        std::cout << "point not found!!!!" << std::endl;
        f[0] = -1.0;
        f[1] = -1.0;
        f[2] = -1.0;
        n[0] = -1.0;
        n[1] = -1.0;
        n[2] = -1.0;
    }

}
