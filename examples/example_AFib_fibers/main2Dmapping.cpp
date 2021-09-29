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

void cross(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& cross_product)
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
    libMesh::ReplicatedMesh flat_mesh(init.comm());
    libMesh::ReplicatedMesh atlas_mesh(init.comm());
    libMesh::ReplicatedMesh la_mesh(init.comm());
    const unsigned int dim = flat_mesh.mesh_dimension();

    // We may need XDR support compiled in to read binary .xdr files
    std::string flat_meshfile = data("flat_mesh", "NONE");
    std::string atlas_meshfile = data("atlas_mesh", "NONE");
    std::string la_meshfile = data("la_mesh", "NONE");

    atlas_mesh.read(atlas_meshfile);
    flat_mesh.read(flat_meshfile);
    la_mesh.read(la_meshfile);

    libMesh::PointLocatorTree* locator_ptr;
    libMesh::PointLocatorTree atlas_locator(atlas_mesh);

    //std::cout << "setup " << nodeset_IDs.size() << " nodesets done!" << std::endl;
    //mesh.get_boundary_info().print_info();
    // Create LibMesh equation system
    libMesh::EquationSystems atlas_es(atlas_mesh);
    libMesh::EquationSystems flat_es(flat_mesh);
    libMesh::EquationSystems la_es(la_mesh);

    int N = 0;
    if(true)
    {
        lpv_antra_threshold = data("lpv_antra", 0.225);
        rpv_antra_threshold = data("rpv_antra", 0.775);
        lpw_threshold = data("lpw", 0.2);
        rpw_threshold = data("rpw", 0.8);

        std::string input_files_list = data("input_files","data.beat");
        std::vector<std::string> input_files;
        BeatIt::readList(input_files_list, input_files);
        // Solve each Poisson problem
        N = input_files.size();
        using namespace libMesh;

        bool add_LAA_and_FO = data("add_laa_fo", true);
        typedef BeatIt::Poisson Poisson;
        typedef std::shared_ptr<Poisson> PoissonPtr;
        std::vector<PoissonPtr> pv;
        std::vector<PoissonPtr> flat_pv;
        std::vector<PoissonPtr> la_pv;
        //BeatIt::Poisson poisson(es, "Poisson");

        for (int i = 0; i < input_files.size(); i++)
        {
            std::string pname = "u" + std::to_string(i);
            std::cout << "Solving problem: " << pname << std::endl;
            pv.emplace_back( new  Poisson(atlas_es, pname));
            flat_pv.emplace_back( new  Poisson(flat_es, pname));
            la_pv.emplace_back( new  Poisson(la_es, pname));

            std::cout << "Calling setup: ..." << std::flush;
            GetPot data_i(input_files[i]);
            pv[i]->setup(data_i, "poisson");
            flat_pv[i]->setup(data_i, "poisson");
            la_pv[i]->setup(data_i, "poisson");
            std::cout << " Done!" << std::endl;
            std::cout << "Calling assemble system: ..." << std::flush;
            pv[i]->assemble_system();
            std::cout << " Done!" << std::endl;
            std::cout << "Calling solve system: ..." << std::flush;
            pv[i]->solve_system();
            std::cout << " Done!" << std::endl;
            std::cout << " Done!" << std::endl;
            atlas_es.get_system<LinearImplicitSystem>("u"+std::to_string(i)).update();

            if(i == 0 && add_LAA_and_FO)
            {
                laa_threshold = data_i("laa_threshold", 0.22);
                fo_threshold = data_i("fo_threshold", 0.78);
                for(auto & node : atlas_mesh.local_node_ptr_range() )
                {
                    unsigned int dn = node->dof_number(atlas_es.get_system<ExplicitSystem>("u"+std::to_string(i)).number(), 0, 0);
                    double u1 = atlas_es.get_system<ExplicitSystem>("u"+std::to_string(i)).solution->el(dn);
                    if(u1 < laa_threshold)
                    {
                        atlas_mesh.get_boundary_info().add_node(node, 123);
                    }
                    else if(u1 > fo_threshold)
                    {
                        atlas_mesh.get_boundary_info().add_node(node, 321);
                    }
                }
            }
            // Copy solution to flat and la systems
        }


        // DEFINE FIBER SYSTEMS:
        std::vector<libMesh::EquationSystems * > es_ptr(3);
        es_ptr[0] = &atlas_es;
        es_ptr[1] = &flat_es;
        es_ptr[2] = &la_es;

        for(int k = 1; k < 3; ++k)
        {
            FiberSystem& f_sys = es_ptr[k]->add_system<FiberSystem>("fibers");
            f_sys.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
            f_sys.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
            f_sys.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
            f_sys.add_vector("backup");
            f_sys.init();

            FiberSystem& s_sys = es_ptr[k]->add_system<FiberSystem>("sheets");
            s_sys.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
            s_sys.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
            s_sys.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
            s_sys.init();

            FiberSystem& n_sys = es_ptr[k]->add_system<FiberSystem>("xfibers");
            n_sys.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
            n_sys.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
            n_sys.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
            n_sys.init();

        }
        FiberSystem& f_sys = atlas_es.add_system<FiberSystem>("fibers");
        f_sys.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
        f_sys.add_vector("backup");
        f_sys.init();
        auto& f_v = f_sys.solution;

        FiberSystem& s_sys = atlas_es.add_system<FiberSystem>("sheets");
        s_sys.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
        s_sys.init();
        auto& s_v = s_sys.solution;

        FiberSystem& n_sys = atlas_es.add_system<FiberSystem>("xfibers");
        n_sys.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
        n_sys.init();
        auto& n_v = n_sys.solution;


        // Define auxiliary vectors used in the computations of the fibers
        double f[3];
        double s[3];
        double n[3];

        std::vector<libMesh::Gradient> du(N);


        libMesh::MeshBase::const_element_iterator el = atlas_mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el = atlas_mesh.active_local_elements_end();
        std::cout << "\nGetting dof map:  ... " << std::flush;
        const libMesh::DofMap & dof_map = f_sys.get_dof_map();
        std::cout << " done \n " << std::flush;

        std::vector<libMesh::dof_id_type> dof_indices;
        std::vector<double> u(N);

        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;
            dof_map.dof_indices(elem, dof_indices);
            auto elID = elem->id();
            Point centroid = elem->centroid();
            auto blockID = atlas_mesh.elem_ptr(elID)->subdomain_id();


            for (int k = 0; k < N; ++k)
            {
                u[k] = atlas_es.get_system<LinearImplicitSystem>("u"+std::to_string(k)).point_value(0, centroid, elem);
                du[k] = atlas_es.get_system<LinearImplicitSystem>("u"+std::to_string(k)).point_gradient(0, centroid, elem);
            }

            //if (atlas_meshfile == "NONE")
                evaluate(f, s, n, u, du, blockID, centroid, data, dim);
            //else
            //    evaluate_atlas(f, s, n, u, du, centroid, *atlas_es_ptr, *locator_ptr);


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
            // copy fibers to flat and la
        }
        f_v->close();
        s_v->close();
        n_v->close();
    }


    std::cout << "ATLAS: " << atlas_meshfile << std::endl;
    //flat_es.print_info();
    // REMAP
    int missed_nodes = 0;
    for(auto & node : flat_mesh.node_ptr_range() )
    {
        libMesh::Point p(*node);
        const Elem *atlas_element = atlas_locator(p);
        std::vector<double> u(N, -666.0);


        if(atlas_element)
        {
            for (int k = 0; k < N; ++k)
            {
                u[k] = atlas_es.get_system<LinearImplicitSystem>("u"+std::to_string(k)).point_value(0, p, atlas_element);
            }
        }
        else
        {
            const Elem *atlas_element2 = atlas_locator.perform_linear_search(p, nullptr, true, 1e-5);
            if(atlas_element2)
            {
                for (int k = 0; k < N; ++k)
                {
                    u[k] = atlas_es.get_system<LinearImplicitSystem>("u"+std::to_string(k)).point_value(0, p, atlas_element2);
                }
            }
            else
            {
                missed_nodes ++;
            }
        }

        for (int k = 0; k < N; ++k)
        {
            LinearImplicitSystem & system = flat_es.get_system<LinearImplicitSystem>("u"+std::to_string(k));
            const libMesh::DofMap & dof_map = system.get_dof_map();
            std::vector<libMesh::dof_id_type> dof_indices;
            dof_map.dof_indices(node, dof_indices);
            system.solution->set(dof_indices[0], u[k]);
        }
    }

    int missed_elem = 0;

    for(auto & elem : flat_mesh.element_ptr_range() )
    {
        std::vector<double> ff(3, 0.0);
        ff[0] = 1.0; ff[1] = 0.0; ff[2] = 0.0;
        std::vector<double> ss(3, 0.0);
        std::vector<double> xf(3, 0.0);
        xf[0] = 1.0; xf[1] = 0.0; xf[2] = 0.0;
        libMesh::Point p( elem->centroid() );
        const Elem *atlas_element = atlas_locator(p);
        FiberSystem & atlas_fsystem = atlas_es.get_system<FiberSystem>("fibers");
        const libMesh::DofMap & atlas_f_dof_map = atlas_fsystem.get_dof_map();
        std::vector<libMesh::dof_id_type> atlas_dof_indices;


        if(atlas_element)
        {
            atlas_f_dof_map.dof_indices(atlas_element, atlas_dof_indices);
            atlas_fsystem.solution->localize(ff, atlas_dof_indices);
        }
        else
        {
            const Elem *atlas_element2 = atlas_locator.perform_linear_search(p, nullptr, true, 1e-5);
            if(atlas_element2)
            {
                atlas_f_dof_map.dof_indices(atlas_element2, atlas_dof_indices);
                atlas_fsystem.solution->localize(ff, atlas_dof_indices);
            }
            else
            {
                missed_elem ++;
            }
        }

        FiberSystem & fsystem = flat_es.get_system<FiberSystem>("fibers");
        const libMesh::DofMap & f_dof_map = fsystem.get_dof_map();
        std::vector<libMesh::dof_id_type> f_dof_indices;
        f_dof_map.dof_indices(elem, f_dof_indices);
        fsystem.solution->set(f_dof_indices[0], ff[0]);
        fsystem.solution->set(f_dof_indices[1], ff[1]);
        fsystem.solution->set(f_dof_indices[2], ff[2]);
    }

    FiberSystem & f_system_2D = flat_es.get_system<FiberSystem>("fibers");
    FiberSystem & f_system_3D = la_es.get_system<FiberSystem>("fibers");
    f_system_2D.solution->close();
    f_system_3D.solution->close();
    *f_system_3D.solution = *f_system_2D.solution;

    for(auto & elem : la_mesh.element_ptr_range() )
    {
        std::vector<double> ff(3, 0.0);
        ff[0] = 1.0; ff[1] = 0.0; ff[2] = 0.0;
        std::vector<double> ss(3, 0.0);
        std::vector<double> xf(3, 0.0);
        xf[0] = 1.0; xf[1] = 0.0; xf[2] = 0.0;
        FiberSystem & la_fsystem = la_es.get_system<FiberSystem>("fibers");
        const libMesh::DofMap & la_f_dof_map = la_fsystem.get_dof_map();
        std::vector<libMesh::dof_id_type> la_dof_indices;
        la_f_dof_map.dof_indices(elem, la_dof_indices);
        la_fsystem.solution->localize(ff, la_dof_indices);

        libMesh::Gradient f2D(ff[0], ff[1], ff[2]);

        libMesh::Point p0(*elem->node_ptr(0));
        libMesh::Point p1(*elem->node_ptr(1));
        libMesh::Point p2(*elem->node_ptr(2));

        libMesh::Point v1(p1-p0);
        libMesh::Point v2(p2-p0);

        libMesh::Point S = v1.cross(v2);
        S = -S.unit();
        auto S2D = libMesh::Point(0,0,1);
        S2D.unit();
        libMesh::Point u = S.cross(S2D);
        u = u.unit();
        // rotation angle
        double cos_tetha = S * S2D;
        double tetha = -std::acos(cos_tetha);
        double sin_tetha = std::sin(tetha);

        libMesh::TensorValue<libMesh::Number> RR;
        RR(0,0) =          cos_tetha + u(0) * u(0) * ( 1 - cos_tetha);
        RR(0,1) = - u(2) * sin_tetha + u(0) * u(1) * ( 1 - cos_tetha);
        RR(0,2) =   u(1) * sin_tetha + u(0) * u(2) * ( 1 - cos_tetha);
        RR(1,0) =   u(2) * sin_tetha + u(1) * u(0) * ( 1 - cos_tetha);
        RR(1,1) =          cos_tetha + u(1) * u(1) * ( 1 - cos_tetha);
        RR(1,2) = - u(0) * sin_tetha + u(1) * u(2) * ( 1 - cos_tetha);
        RR(2,0) = - u(1) * sin_tetha + u(2) * u(0) * ( 1 - cos_tetha);
        RR(2,1) =   u(0) * sin_tetha + u(2) * u(1) * ( 1 - cos_tetha);
        RR(2,2) =          cos_tetha + u(2) * u(2) * ( 1 - cos_tetha);
        libMesh::Point Rf =  RR * f2D;
        Rf = Rf.unit();
        ff[0] = Rf(0);
        ff[1] = Rf(1);
        ff[2] = Rf(2);
        ss[0] = S(0);
        ss[1] = S(1);
        ss[2] = S(2);
        cross(ss, ff, xf);
        la_es.get_system<FiberSystem>("fibers").solution->set(la_dof_indices[0], ff[0]);
        la_es.get_system<FiberSystem>("fibers").solution->set(la_dof_indices[1], ff[1]);
        la_es.get_system<FiberSystem>("fibers").solution->set(la_dof_indices[2], ff[2]);
    }

    for (int k = 0; k < N; ++k)
    {
        LinearImplicitSystem & system_2D = flat_es.get_system<LinearImplicitSystem>("u"+std::to_string(k));
        LinearImplicitSystem & system_3D = la_es.get_system<LinearImplicitSystem>("u"+std::to_string(k));
        system_2D.solution->close();
        system_3D.solution->close();
        *system_3D.solution = *system_2D.solution;
    }


    std::cout << "Missed nodes: " << missed_nodes << std::endl;
    std::cout << "Missed elements: " << missed_elem << std::endl;
    std::cout << "Calling exporter: ..." << ". " << std::flush;
    ExodusII_IO flat_exporter(flat_mesh);
    ExodusII_IO la_exporter(la_mesh);

    ExodusII_IO atlas_exporter(atlas_mesh);
    atlas_exporter.write_equation_systems("atlas.e", atlas_es);
    atlas_exporter.write_element_data(atlas_es);



//    // the solution of the poisson problems
//    bool export_only_fibers = data("only_fibers", false);
//    if (export_only_fibers)
//    {
//        std::vector < std::string > varname(9);
//        varname[0] = "fibersx";
//        varname[1] = "fibersy";
//        varname[2] = "fibersz";
//        varname[3] = "sheetsx";
//        varname[4] = "sheetsy";
//        varname[5] = "sheetsz";
//        varname[6] = "xfibersx";
//        varname[7] = "xfibersy";
//        varname[8] = "xfibersz";
//        exporter.set_output_variables(varname);
//    }

    flat_exporter.write_equation_systems("flat.e", flat_es);
    flat_exporter.write_element_data(flat_es);
    la_exporter.write_equation_systems("la.e", la_es);
    la_exporter.write_element_data(la_es);

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

    //std::cout << "N: " << phi.size() << std::endl;
    //for ( auto && u : phi) std::cout << u << " " << std::flush;
    //std::cout << std::endl;
    //std::cout << "laa_threshold: " << laa_threshold << std::endl;
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
                std::cout << "lpv_antra_threshold: " << lpv_antra_threshold << std::endl;
                std::cout << " LPV " << std::endl;
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
                std::cout << "phi[4]: " << phi[4] << std::endl;
                std::cout << "rpv_antra_threshold: " << rpv_antra_threshold << std::endl;
                // Right Superior PVs
                if(phi[2] < 0.5 )
                {
                    std::cout << " RSPV " << std::endl;
                    n[0] = dphi[2](0);  n[1] = dphi[2](1);  n[2] = 0.0;
                    normalize(n[0],n[1],n[2],  1.0,  0.0 , 0.0);
                    cross(s, n, f);
                }
                else
                {
                    std::cout << " RIPV " << std::endl;
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
                    std::cout << " LAA " << std::endl;
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
                            //std::cout << " PW " << std::endl;
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
                            std::cout << " ELSE " << std::endl;

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
