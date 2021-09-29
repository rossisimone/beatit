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


using namespace libMesh;

enum class Layer { Endo, Mid, Epi };


enum class FibrosisModel
{
    NONE,
    // From
    // Patient-specific modeling of atrial fibrosis increases the accuracy of sinus rhythm simulations and may explain maintenance of atrial fibrillation
    // Martin W.Krueger et al.
    KRUEGER, //
    // FROM
    // Epicardial Fibrosis Explains Increased Endo�Epicardial Dissociation and Epicardial Breakthroughs in Human Atrial Fibrillation
    // Ali Gharaviri et al.
    PEZZUTO,
    EPI_ENDO_DISSOCIATION,
    OURS
};

struct Thresholds
{
    Thresholds (GetPot& data)
    {
        t_floor = data("t_floor", 0.65);
        t_laa = data("t_laa", 0.35);
        t_fo = data("t_fo", 0.14);
        t_antra_rspv = data("t_antra_rspv", 0.75);
        t_antra_ripv = data("t_antra_ripv", 0.75);
        t_lateral = data("t_lateral", 0.38);
        t_lateral_u5 = data("t_lateral_u5", 0.1);
        t_antra_left = data("t_antra_left", 0.2);
        t_lower_antra_left = data("t_lower_antra_left", 0.2);

        t_septum_rspv_1 = data("t_septum_rspv_1", 0.17);
        t_septum_rspv_2 = data("t_septum_rspv_2", 0.5);
        t_septum_ripv_1 = data("t_septum_ripv_1", 0.17);
        t_septum_ripv_2 = data("t_septum_ripv_2", 0.37);

        t_anterior_bottom_1 = data("t_anterior_bottom_1", 0.4);
        t_anterior_bottom_2 = data("t_anterior_bottom_2", 0.5);
        t_anterior_top_1 = data("t_anterior_top_1", 0.499);
        t_anterior_top_2 = data("t_anterior_top_2", 0.3);

        t_posterior_bottom = data("t_posterior_bottom", 0.4);
        t_posterior_top = data("t_posterior_top", 0.0);

        t_left_right = data("t_left_right", 0.5);
        t_anterior_posterior = data("t_anterior_posterior", 0.5);
        t_anterior_floor = data("t_anterior_floor", 0.85);
        t_bachmann = data("t_bachmann", 0.45);
        t_anterior_lateral = data("t_anterior_lateral", 0.35);
        t_left_lateral_ridge = data("t_left_lateral_ridge", 0.1);

    }

    double t_floor;
    double t_laa;
    double t_fo;
    double t_antra_rspv;
    double t_antra_ripv;
    double t_lateral;
    double t_lateral_u5;
    double t_antra_left;

    double t_septum_rspv_1;
    double t_septum_rspv_2;
    double t_septum_ripv_1;
    double t_septum_ripv_2;

    double t_anterior_bottom_1;
    double t_anterior_bottom_2;
    double t_anterior_top_1;
    double t_anterior_top_2;


    double t_posterior_bottom;
    double t_posterior_top;

    double t_left_right;
    double t_anterior_posterior;
    double t_anterior_floor;
    double t_anterior_lateral;
    double t_bachmann;
    double t_left_lateral_ridge;
    double t_lower_antra_left;


};

Thresholds * thresholds;

std::map<std::string, FibrosisModel> fibrosis_model_map = { { "none", FibrosisModel::NONE},
                                                            { "krueger", FibrosisModel::KRUEGER},
                                                            { "pezzuto", FibrosisModel::PEZZUTO},
                                                            { "dissociation", FibrosisModel::PEZZUTO},
                                                            { "ours", FibrosisModel::OURS}};

void evaluate(double f[], double s[], double n[], std::vector<double>& phi, std::vector<libMesh::Gradient>& dphi, int blockID, Point& x, GetPot& data);

void assemble_smoother(EquationSystems& es, double D = 1.0);

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

int main(int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace
    BeatIt::printBanner(std::cout);
    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    GetPot data = BeatIt::readInputFile(argc, argv);

    thresholds = new Thresholds(data);
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
    libMesh::ParallelMesh mesh(init.comm());

    // We may need XDR support compiled in to read binary .xdr files
    std::string meshfile = data("mesh", "NONE");

    // Read the input mesh.
    if (meshfile != "NONE")
    {
        // We may need XDR support compiled in to read binary .xdr files
        int n_refinements = data("n_refinements", 0);
        BeatIt::serial_mesh_partition(init.comm(), meshfile, &mesh, n_refinements);

    }
    else
    {
        std::cout << "ABORTING: MeshFile " << meshfile << " not found!" << std::endl;
        std::cout << "Specify a valid mesh file in the input file under: mesh = " << std::endl;
        throw std::runtime_error("MeshFile not found!");
    }
    // Add nodesets for BCs
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
//    mesh.node_ptr(55164)->print();
//    for(int k = 0; k < nodes_IDs.size(); ++k)
//    {
//        mesh.get_boundary_info().add_node(mesh.node_ptr(nodes_IDs[k]), nodeset_IDs[k]);
//    }
    std::set<int> added_nodesets;
    std::cout << "Adding nodesets" << std::endl;

    for(auto& node : mesh.local_node_ptr_range())
    {
        libMesh::Point p(*node);
        for(int k = 0; k < nodeset_x.size(); ++k)
        {
            libMesh::Point q(nodeset_x[k], nodeset_y[k], nodeset_z[k]);
            double r = nodeset_r[k];
            if(r > 0)
            {
                libMesh::Point qmp(p-q);
                if( qmp.norm() <= r)
                {
                    //std::cout << "Adding nodeset: " << nodeset_IDs[k] << " to node:  " << node->id() << std::endl;
                    mesh.get_boundary_info().add_node(node, nodeset_IDs[k]);
                    added_nodesets.insert(nodeset_IDs[k]);
                }
            }
        }
    }
    if(added_nodesets.size() !=  nodeset_x.size() )
    {
        std::cout << "SOMETHING WENT WRONG: may be you are running in parallel!" << std::endl;
        //return 1;
    }

    std::cout << "Mesh prepare for use" << std::endl;
    mesh.prepare_for_use();
    std::cout << "setup " << nodeset_IDs.size() << " nodesets done!" << std::endl;
    //mesh.get_boundary_info().print_info();
    // Create LibMesh equation system
    libMesh::EquationSystems es(mesh);
    // Name each system with a number

    //std::cout << "Solving " << N << " Poisson problems" << std::endl;
    std::string input_files_list = data("input_files","data.beat");
    std::vector<std::string> input_files;
    BeatIt::readList(input_files_list, input_files);
    // Solve each Poisson problem
    int N = input_files.size();
    using namespace libMesh;
    for(int i = 0; i < input_files.size(); ++i)
    {
        ExplicitSystem& u_sys = es.add_system<ExplicitSystem>("u"+std::to_string(i));
        u_sys.add_variable("u"+std::to_string(i));
        u_sys.init();
    }


    BeatIt::Poisson poisson(es, "Poisson");
    for (int i = 0; i < input_files.size(); i++)
    {
        std::cout << "Solving problem: " << i << std::endl;
        std::cout << "Calling setup: ..." << std::flush;
        GetPot data_i(input_files[i]);
        poisson.setup(data_i, "poisson");
        std::cout << " Done!" << std::endl;
        std::cout << "Calling assemble system: ..." << std::flush;
        poisson.assemble_system();
        std::cout << " Done!" << std::endl;
        std::cout << "Calling solve system: ..." << std::flush;
        poisson.solve_system();
        std::cout << " Done!" << std::endl;
        std::cout << "Copy solution: ..." << std::flush;
        (*es.get_system<ExplicitSystem>("u"+std::to_string(i)).solution) = (*es.get_system<LinearImplicitSystem>("Poisson").solution);
        std::cout << " Done!" << std::endl;

        // let's set up the LAA nodeset
        if(i == 1)
        {
            double laa_threshold = thresholds->t_laa; //data_i("laa_threshold", 0.22);
            double fo_threshold = thresholds->t_fo; //data_i("fo_threshold", 0.1);
            for(auto & node : mesh.local_node_ptr_range() )
            {
                unsigned int dn = node->dof_number(es.get_system<ExplicitSystem>("u"+std::to_string(i)).number(), 0, 0);
                double u1 = es.get_system<ExplicitSystem>("u"+std::to_string(i)).solution->el(dn);
                if(u1 > laa_threshold)
                {
                    mesh.get_boundary_info().add_node(node, 123);
                }
                else if(u1 < fo_threshold)
                {
                    mesh.get_boundary_info().add_node(node, 321);
                }
            }
        }

        // let's set up the LAA nodeset
        if(i == 3)
        {
            //double threshold = data_i("lpv_threshold", 0.3);
            double threshold = thresholds->t_antra_left;
            for(auto & node : mesh.local_node_ptr_range() )
            {
                unsigned int dn = node->dof_number(es.get_system<ExplicitSystem>("u"+std::to_string(i)).number(), 0, 0);
                double u3 = es.get_system<ExplicitSystem>("u"+std::to_string(i)).solution->el(dn);
                if(u3 < threshold)
                {
                    mesh.get_boundary_info().add_node(node, 333);
                }
            }
        }

        es.get_system<ExplicitSystem>("u"+std::to_string(i)).update();
    }

    // DEFINE FIBER SYSTEMS:
    typedef libMesh::ExplicitSystem FiberSystem;
    FiberSystem& f_sys = es.add_system<FiberSystem>("fibers");
    f_sys.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_sys.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
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
            u[k] = es.get_system<ExplicitSystem>("u"+std::to_string(k)).point_value(0, centroid, elem);
            du[k] = es.get_system<ExplicitSystem>("u"+std::to_string(k)).point_gradient(0, centroid, elem);
        }

        evaluate(f, s, n, u, du, blockID, centroid, data);
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

        // Fossa Ovalis
        if(blockID == 2)
            mesh.elem_ptr(elID)->subdomain_id() = blockID;

    }

    std::cout << "Calling exporter: ..." << ". " << std::flush;
    ExodusII_IO exporter(mesh);
    // the solution of the poisson problems
    bool export_only_fibers = data("only_fibers", false);
    std::vector < std::string > varname;
    if (export_only_fibers)
    {
        varname.push_back("fibersx");
        varname.push_back("fibersy");
        varname.push_back("fibersz");
        varname.push_back("sheetsx");
        varname.push_back("sheetsy");
        varname.push_back("sheetsz");
        varname.push_back("xfibersx");
        varname.push_back("xfibersy");
        varname.push_back("xfibersz");
        exporter.set_output_variables(varname);
    }
    bool export_fields = data("fields", false);
    if (export_fields)
    {
        varname.push_back("u0");
        varname.push_back("u1");
        varname.push_back("u2");
        varname.push_back("u3");
        varname.push_back("u4");
        varname.push_back("u5");
        varname.push_back("u6");
    }
    exporter.set_output_variables(varname);

    std::string output = data("output", "test.e");
    exporter.write_equation_systems(output, es);
    exporter.write_element_data(es);

    // smooth fibers
    bool do_smoothing = data("smoothing", false);
    if(do_smoothing)
    {
        std::cout << "Smoothign fibers " << std::endl;
        double D = data("D", 1.0);
        LinearImplicitSystem& fs_sys = es.add_system<LinearImplicitSystem>("Smoothing_fibers");
        fs_sys.add_variable("fx", libMesh::FIRST, libMesh::LAGRANGE);
        fs_sys.add_variable("fy", libMesh::FIRST, libMesh::LAGRANGE);
        fs_sys.add_variable("fz", libMesh::FIRST, libMesh::LAGRANGE);
        fs_sys.init();
        std::cout << "Assemble smoother fibers " << std::endl;
        fs_sys.assemble_before_solve = false;
        fs_sys.zero_out_matrix_and_rhs = false;
        assemble_smoother(es, D);
        std::cout << "solve" << std::endl;
        fs_sys.solve();
        fs_sys.update();

        std::cout << "Reassign fibers" << std::endl;

        BeatIt::Util::normalize(*fs_sys.solution);

        //reassign the fibers
        el = mesh.active_local_elements_begin();
        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;
            dof_map.dof_indices(elem, dof_indices);
            auto elID = elem->id();
            Point centroid = elem->centroid();

            f[0] = fs_sys.point_value(0, centroid, elem);
            f[1] = fs_sys.point_value(1, centroid, elem);
            f[2] = fs_sys.point_value(2, centroid, elem);


            s[0] = s_v->el(dof_indices[0]);
            s[1] = s_v->el(dof_indices[1]);
            s[2] = s_v->el(dof_indices[2]);

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
        libMesh::ExplicitSystem& field_system = es.get_system<libMesh::ExplicitSystem> ("u"+std::to_string(field));
        libMesh::ExplicitSystem& u0_system = es.get_system<libMesh::ExplicitSystem> ("u"+std::to_string(0));

        const libMesh::DofMap & noise_dof_map = noise_sys.get_dof_map();
        const libMesh::DofMap & sigma_dof_map = sigma_sys.get_dof_map();
        std::vector<libMesh::dof_id_type> dof_indices_noise;
        std::vector<libMesh::dof_id_type> dof_indices_sigma;

        for(auto & elem : mesh.active_local_element_ptr_range() )
        {
            noise_dof_map.dof_indices(elem,dof_indices_noise);
            sigma_dof_map.dof_indices(elem,dof_indices_sigma);
            double u = field_system.point_value(0, elem->centroid(), elem);
            double u0 = u0_system.point_value(0, elem->centroid(), elem);
            auto blockID = elem->subdomain_id();


            //
            double w = 0;
            if(u > fibrosis_threshold2) w = 1;
            else if(u > fibrosis_threshold1)
            {
                w = (u-fibrosis_threshold1)/std::abs(fibrosis_threshold2-fibrosis_threshold1);
            }
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

    // // export activation times again
    solver->save_activation_times_nemesis(save_iter);
    //solver->save_activation_times(save_iter);
    // delete solver before ending the simulation
    // avoiding memory leaks
    delete thresholds;
    delete solver;



    return 0;
}

void evaluate(double f[], double s[], double n[], std::vector<double>& phi, std::vector<libMesh::Gradient>& dphi, int blockID, Point& x, GetPot& data)
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

    double sx = dphi[0](0);
    double sy = dphi[0](1);
    double sz = dphi[0](2);
    normalize(sx, sy, sz, 1.0, 0.0, 0.0);
    s[0] = sx;  s[1] = sy;  s[2] = sz;

//
//    double fx, fy, fz;
//    double nx, ny, nz;

    double u0 = phi[0];
    Layer layer = Layer::Endo;
    std::string which_layer = data("layer", "all");
    if(which_layer == "all")
    {
        if( u0 > 0.5 ) layer = Layer::Epi;
    }
    if(which_layer == "epi")
    {
        layer = Layer::Epi;
    }

    bool common_left_trunk = data("common_left_trunk", false);

    //else if ( u0 > 0.25 ) layer = Layer::Mid;

    // FLOOR
    // ANTERIOR FLOOR
    // thresholds->t_anterior_floor
    // M75 = 0.85
    // M73 = 0.85
    if( phi[4] > thresholds->t_anterior_floor && phi[3] > thresholds->t_anterior_posterior)
    {
        cross(s, 0, n, 4, dphi, f);
    }
    // POSTERIOR FLOOR
    else if( phi[4] > thresholds->t_floor && phi[3] < thresholds->t_anterior_posterior )
    {
        cross(s, 0, n, 4, dphi, f);

        /*
        if( Layer::Epi == layer || phi[4] > thresholds->t_anterior_floor)
        {
            cross(s, 0, n, 4, dphi, f);
        }
        else
        {
            cross(s, 0, n, 2, dphi, f);
        }
        */
    }
    else
    {
        // LAA
        // thresholds->t_laa
        // M73 = 0.32
        if( phi[1] > thresholds->t_laa)
        {
            cross(s, 0, n, 1, dphi, f);
        }
        // FO
        // thresholds->t_fo
        // M73 = 0.2
        else if( phi[1] < thresholds->t_fo)
        {
            cross(s, 0, f, 5, dphi, n);
            if( Layer::Epi == layer)
            {
                cross(s, 0, n, 5, dphi, f);
            }
        }
        else
        {
            // RIGHT ANTRA
            // thresholds->t_antra_rspv
            // M75 = 0.78
            // M73 = 0.78
            // thresholds->t_antra_ripv
            // M75 = 0.78
            // M75 = 0.78
            if( phi[6] > thresholds->t_antra_rspv ||
                phi[6] > thresholds->t_antra_ripv)
            {
                cross(s, 0, n, 3, dphi, f);
                if( Layer::Epi == layer)
                {
                    if( phi[3] < thresholds->t_anterior_posterior )
                    {
                        cross(s, 0, f, 4, dphi, n);
                    }
                }
            }
            // thresholds->t_septum_ripv_2
            // M75 = 0.37
            // M73 = 0.32
            // thresholds->t_left_right
            // M75 = 0.5
            // M73 = 0.65
            // thresholds->t_anterior_posterior
            // M75 = 0.5
            // M73 = 0.42
            else if ( phi[3] < thresholds->t_septum_ripv_2 &&
                      phi[2] > thresholds->t_left_right &&
                      Layer::Epi == layer )
            {
                cross(s, 0, f, 4, dphi, n);
            }
            // LEFT ANTRA
            // thresholds->t_antra_left
            // M75 = 0.2
            // M73 = 0.27
            else if( phi[6] < thresholds->t_antra_left)
            {
                if(common_left_trunk) cross(s, 0, n, 6, dphi, f);
                else cross(s, 0, n, 3, dphi, f);
            }
            else
            {
                // ANTERIOR
                if( phi[3] > thresholds->t_anterior_posterior)
                {
                    // thresholds->t_anterior_lateral
                    // M75 = 0.35
                    // M73 = 0.5
                    // thresholds->t_left_lateral_ridge
                    // M75 = 0.1
                    // M73 = 0.2
                    if( phi[5] < thresholds->t_left_lateral_ridge &&
                        phi[6] < thresholds->t_anterior_lateral &&
                        phi[4] < thresholds->t_anterior_top_1 )
                    //if( phi[2] < thresholds->t_lower_antra_left )
                    {
                        // ENDO use U4
                        cross(s, 0, n, 4, dphi, f);
                    }
                    else if
                    ( phi[4] > thresholds->t_anterior_top_2 &&
                      phi[5] < thresholds->t_lateral_u5 &&
                      phi[6] < thresholds->t_anterior_bottom_1 )
                    {
                        // ENDO use U4
                        cross(s, 0, n, 4, dphi, f);
                    }
                    else
                    {
                        // ENDO use U5
                        cross(s, 0, n, 5, dphi, f);
                    }
                    if( Layer::Epi == layer)
                    {
                        // thresholds->t_bachmann
                        // M75 = 0.45
                        // M73 = 0.450.
                        if( phi[4] > thresholds->t_bachmann )
                        {
                            cross(s, 0, n, 4, dphi, f);
                        }
                        else
                        {
                            if( phi[2] > thresholds->t_left_right )
                            {
                                cross(s, 0, n, 2, dphi, f);
                            }
                            else
                            {
                                //if( phi[6] <  0.3 && phi[4] > 0.3)
                                //if(phi[5] < 0.1 && phi[6] <  0.35)
                                if( phi[5] < thresholds->t_left_lateral_ridge &&
                                    phi[6] < thresholds->t_anterior_lateral &&
                                    phi[4] < thresholds->t_anterior_top_1 )
//                                if( phi[2] < thresholds->t_lower_antra_left )
                                {
                                    cross(s, 0, n, 4, dphi, f);
                                }
                                else
                                {
                                    cross(s, 0, n, 5, dphi, f);
                                }
                            }
                        }
                    }

                }
                // POSTERIOR
                else
                {
                    // ENDO use U2
                    cross(s, 0, n, 2, dphi, f);
                    // thresholds->t_lateral
                    // M73 = 0.295
                    if( phi[1] > thresholds->t_lateral )
                    {
                        cross(s, 0, n, 4, dphi, f);
                    }

                    if( Layer::Epi == layer)
                    {

                        // M75 M1 = 0.21
                        // M75 M3 = 0.245
                        if( phi[1] < thresholds->t_septum_ripv_1 )
                        {
                            cross(s, 0, f, 4, dphi, n);
                        }
                    }
                }
            }
        }
    }

    return;


    switch(layer)
    {
    case Layer::Epi:
    {
        double u1 = phi[1];
        // LAA
        double t1 = 0.15;
        double t2 = 0.005;
        double t3 = 0.1;
        double t4 = 0.52;
        double t5 = 0.5;
        double t6 = 0.7;
        if(u1 > t1)
        {
            // s = du0
            // n = du1
            // f = du0 x du1
            cross(s, 0, n, 1, dphi, f);
            blockID = 1;
        }
        // Fossa Ovalis
        else if (u1 < t2)
        {
            // s = du0
            // f = du1
            // n = du0 x du1
            cross(s, 0, n, 1, dphi, f);
            blockID = 2;
        }
        else
        {
            double u2 = phi[2];
            // Right pulmonary veins
            if(u2 < t4) // u8 > 0.7
            {
                // s = du0
                // n = du4
                // f = du0 x du4
                cross(s, 0, n, 2, dphi, f)  ;
                blockID = 3;
            }
            else
            {
                double u3 = phi[3];
                if(u3 < t5)
                {
                    // s = du0
                    // n = du7
                    // f = du0 x du7
                    cross(s, 0, n, 3, dphi, f);
                    blockID = 4;
                }
                else
                {
                    double u8 = phi[8];
                    if(u8 > t6)
                    {
                        // s = du0
                        // n = du7
                        // f = du0 x du7
                        cross(s, 0, n, 4, dphi, f);
                        blockID = 4;
                    }
                    else
                    {
                        cross(s, 0, n, 8, dphi, f);
                    }

                }
            }
        }

        break;
    }
        case Layer::Endo:
        {
            double u1 = phi[1];
            // LAA
            double t1 = 0.15;
            double t2 = 0.005;
            double t3 = 0.1;
            double t4 = 0.77;
            double t5 = 0.78;
            double t6 = 1.3;
            if(u1 > t1)
            {
                // s = du0
                // n = du1
                // f = du0 x du1
                cross(s, 0, f, 1, dphi, n);
                blockID = 1;
            }
            // Fossa Ovalis
            else if (u1 < t2)
            {
                // s = du0
                // f = du1
                // n = du0 x du1
                cross(s, 0, n, 1, dphi, f);
                blockID = 2;
            }
            else
            {
                double u8 = phi[8];
                // Right pulmonary veins
                if(u8 > t4) // u8 > 0.7
                {
                    // s = du0
                    // n = du4
                    // f = du0 x du4
                    cross(s, 0, n, 4, dphi, f)  ;
                    blockID = 3;
                }
                else
                {
                    double u7 = phi[7];
                    if(u7 > t5)
                    {
                        // s = du0
                        // n = du7
                        // f = du0 x du7
                        cross(s, 0, n, 8, dphi, f);
                        blockID = 4;
                    }
                    else
                    {
                        double u6 = phi[6];
                        if(u6 > t6)
                        {
                            // s = du0
                            // n = du7
                            // f = du0 x du7
                            cross(s, 0, n, 6, dphi, f);
                            blockID = 4;
                        }
                        else
                        {
                            cross(s, 0, n, 5, dphi, f);
                        }

                    }
                }
            }

            break;
        }
        case Layer::Mid:
        {
            double u1 = phi[1];
            // LAA
            double t1 = 0.2;
            double t2 = 0.005;
            double t3 = 0.1;
            double t4 = 0.7;
            double t5 = 0.6;
            if(u1 > t1)
            {
                // s = du0
                // n = du1
                // f = du0 x du1
                cross(s, 0, n, 1, dphi, f);
                blockID = 1;
            }
            // Fossa Ovalis
            else if (u1 < t2)
            {
                // s = du0
                // f = du1
                // n = du0 x du1
                cross(s, 0, n, 1, dphi, f);
                blockID = 2;
            }
            else if(u1 > t3)
            {
                // s = du0
                // n = du7
                // f = du0 x du7
                double w = (u1 - t3 ) / (t1-t3);
                cross(s, 0, n, 1, dphi, f);
                double nt[3];
                double ft[3];
                cross(s, 0, nt, 7, dphi, ft);
                double dot = f[0] * ft[0] + f[1] * ft[1] + f[2] * ft[2];
                double angle = std::acos(dot);
                if( std::abs(angle) > M_PI / 2  )
                {
                    ft[0] = -ft[0];
                    ft[1] = -ft[1];
                    ft[2] = -ft[2];
                }
                for(int j = 0; j < 3; ++j)
                {
                    f[j] = w * f[j] + (1-w) * ft[j];
                    n[j] = w * n[j] + (1-w) * nt[j];
                }
                normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);
                normalize(n[0], n[1], n[2], 0.0, 0.0, 1.0);
                blockID = 1;

            }
            else
            {
                double u8 = phi[8];
                // Right pulmonary veins
                if(u8 > t4) // u8 > 0.7
                {
                    // s = du0
                    // n = du4
                    // f = du0 x du4
                    cross(s, 0, n, 4, dphi, f)  ;
                    blockID = 3;
                }
                else if(u8 > t5) // u8 > 0.6
                {
                    // s = du0
                    // n = du7
                    // f = du0 x du7
                    double w = (u8 - t5 ) / (t4-t5);
                    cross(s, 0, n, 4, dphi, f);
                    double nt[3];
                    double ft[3];
                    cross(s, 0, nt, 7, dphi, ft);
                    for(int j = 0; j < 3; ++j)
                    {
                        f[j] = w * f[j] + (1-w) * ft[j];
                        n[j] = w * n[j] + (1-w) * nt[j];
                    }
                    normalize(f[0], f[1], f[2], 1.0, 0.0, 0.0);
                    normalize(n[0], n[1], n[2], 0.0, 0.0, 1.0);
                    blockID = 4;

                }
                else
                {
                    // s = du0
                    // n = du7
                    // f = du0 x du7
                    cross(s, 0, n, 7, dphi, f);
                    blockID = 4;
                }
            }
            break;
        }
        default:
        {
            std::cout << "We should not be here" << std::endl;
            throw std::runtime_error("Error");
            break;
        }
    }



}


void assemble_smoother(EquationSystems& es, double D)
{
    const libMesh::MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

     LinearImplicitSystem& system  =  es.get_system<LinearImplicitSystem>("Smoothing_fibers");
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
    libMesh::DenseMatrix<libMesh::Number> Ke;
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

    double ftxyzi = 0.0;
    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        dof_map_f.dof_indices(elem,dof_indices_f);
        libMesh::Gradient frhs;
        for(int j = 0; j < dim; ++j) frhs(j) = f_sys.solution->el(dof_indices_f[j]);
        libMesh::Gradient s;
        for(int j = 0; j < dim; ++j) s(j) = s_sys.solution->el(dof_indices_f[j]);

        libMesh::TensorValue<double> H;
        for(int i = 0; i < dim; ++i)
        {
            H(i,i) += 1.0;
            for(int j = 0; j < dim; ++j)
            {
                H(i,j) -= s(i) * s(j);
            }
        }

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
        Ke.resize(dof_indices.size(), dof_indices.size());
         Fe.resize(dof_indices.size());


         int n_b = phi_qp.size();
       //  std::cout << "Assembling Poisson ... " << std::endl;
        for (unsigned int qp = 0; qp < qrule_stiffness.n_points(); qp++)
        {
            for (int idim = 0; idim < dim; ++idim)
            {
                for (unsigned int i = 0; i < n_b; i++)
                {
                    for (unsigned int j = 0; j <n_b; j++)
                    {
                        // stiffness term
                        Ke(i + idim * n_b, j + idim * n_b) += JxW_qp[qp] *  dphi_qp[i][qp] * ( H * dphi_qp[j][qp]);
                        Ke(i + idim * n_b, j + idim * n_b) += JxW_qp[qp] * phi_qp[i][qp] * phi_qp[j][qp];
                    }
                    Fe(i + idim * n_b) += JxW_qp[qp] * frhs(idim) * phi_qp[i][qp];
                }

            }
        }
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs   ->add_vector(Fe, dof_indices);

    }

    system.matrix->close();
    system.rhs->close();
}


