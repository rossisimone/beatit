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

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
#include <iomanip>
#include "libmesh/exodusII_io.h"

using namespace libMesh;

const int N = 4;
void evaluate(double f[], double s[], double n[], double phi[], double dphi[][3], int blockID, Point& x, GetPot& data);
void cross(double v1[], double v2[], double cross_product[])
{
    cross_product[0] = v1[1] * v2[2] - v1[2] * v2[1];
    cross_product[1] = v1[2] * v2[0] - v1[0] * v2[2];
    cross_product[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return;
}
void cross(double v1[], int n, double v2[], int m, double dphi[][3], double cross_product[])
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

    v1[0] = dphi[n][0];
    v1[1] = dphi[n][1];
    v1[2] = dphi[n][2];
    normalize(v1[0], v1[1], v1[2], 1.0, 0.0, 0.0);

    v2[0] = dphi[m][0];
    v2[1] = dphi[m][1];
    v2[2] = dphi[m][2];
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
    // allow us to use higher-order approximation.
    // Create a mesh, with dimension to be overridden later, on the
    // default MPI communicator.
    libMesh::Mesh mesh(init.comm());

    // We may need XDR support compiled in to read binary .xdr files
    std::string meshfile = data("mesh", "NONE");

    // Read the input mesh.
    if (meshfile != "NONE")
    {
        mesh.read(&meshfile[0]);
    }
    else
    {
        std::cout << "ABORTING: MeshFile " << meshfile << " not found!" << std::endl;
        std::cout << "Specify a valid mesh file in the input file under: mesh = " << std::endl;
        throw std::runtime_error("MeshFile not found!");
    }

    // Create LibMesh equation system
    libMesh::EquationSystems es(mesh);
    // Name each system with a number
    std::string pois[N];
    for (int i = 0; i < N; i++)
        pois[i] = "poisson" + std::to_string(i);
    // create solver holders
    typedef BeatIt::Poisson Poisson;
    typedef std::shared_ptr<Poisson> PoissonPtr;
    Poisson* pv[N];
    // Solve each Poisson problem

    for (int i = 0; i < N; i++)
    {
        std::cout << "Solving Poisson " << i << std::endl;
        const int N = 1;
        std::string pois1 = pois[i];
        BeatIt::Poisson poisson(es, pois1);
        std::cout << "Calling setup: ..." << std::flush;
        poisson.setup(data, pois1);
        std::cout << " Done!" << std::endl;
        std::cout << "Calling assemble system: ..." << std::flush;
        poisson.assemble_system();
        std::cout << " Done!" << std::endl;
        std::cout << "Calling solve system: ..." << std::flush;
        poisson.solve_system();
        std::cout << " Done!" << std::endl;
        std::cout << "Calling gradient: ..." << std::flush;
        poisson.compute_elemental_solution_gradient();
        std::cout << " Done!" << std::endl;
        pv[i] = &poisson;
        //normalize the gradients
        std::cout << "Normalize: ..." << std::flush;
        auto& grad = es.get_system<libMesh::ExplicitSystem>(pois[i] + "_gradient").solution;
        BeatIt::Util::normalize(*grad, 1.0, 0.0, 0.0);
        std::cout << " Done!" << std::endl;
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
    double fx = 0, fy = 0, fz = 0;
    double f[3];
    double sx = 0, sy = 0, sz = 0;
    double s[3];
    double nx = 0, ny = 0, nz = 0;
    double n[3];

    double phi[N];
    double dphi[N][3];

    // Loop over each element and based on the blockID
    // assign the fibers in some way
    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    std::cout << "\nGetting dof map:  ... " << std::flush;
    const libMesh::DofMap & dof_map = f_sys.get_dof_map();
    const libMesh::DofMap & dof_map_p = es.get_system<libMesh::ExplicitSystem>(pois[0] + "_P0").get_dof_map();
    std::cout << " done \n " << std::flush;
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_p;
    std::cout << "\nLooping over elements: ... \n" << std::flush;
    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        dof_map_p.dof_indices(elem, dof_indices_p);
        auto elID = elem->id();
        Point centroid = elem->centroid();
        const auto blockID = mesh.elem_ptr(elID)->subdomain_id();

        for (int k = 0; k < N; ++k)
        {
            phi[k] = (*es.get_system<libMesh::ExplicitSystem>(pois[k] + "_P0").solution)(dof_indices_p[0]);
            dphi[k][0] = (*es.get_system<libMesh::ExplicitSystem>(pois[k] + "_gradient").solution)(dof_indices[0]);
            dphi[k][1] = (*es.get_system<libMesh::ExplicitSystem>(pois[k] + "_gradient").solution)(dof_indices[1]);
            dphi[k][2] = (*es.get_system<libMesh::ExplicitSystem>(pois[k] + "_gradient").solution)(dof_indices[2]);
        }

        evaluate(f, s, n, phi, dphi, blockID, centroid, data);
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

        mesh.elem_ptr(elID)->subdomain_id() = 1;

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
    std::string output = data("output", "test.e");
    exporter.write_equation_systems(output, es);
    exporter.write_element_data(es);

    return 0;
}

void evaluate(double f[], double s[], double n[], double phi[], double dphi[][3], int blockID, Point& x, GetPot& data)
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

    double potential = phi[0];
    double sx = dphi[0][0];
    double sy = dphi[0][1];
    double sz = dphi[0][2];

    normalize(sx, sy, sz, 1.0, 0.0, 0.0);

    double fx, fy, fz;
    double nx, ny, nz;


//    if (blockID == 117)
//    {
//        potential = phi[1];
//    }

    //read centerline
    double cx = data("blockID"+std::to_string(blockID)+"/cx",1.0);
    double cy = data("blockID"+std::to_string(blockID)+"/cy",0.0);
    double cz = data("blockID"+std::to_string(blockID)+"/cz",0.0);

    switch (blockID)
    {
        // floor
        case 1:
        case 2:
        case 3:
        case 8:
        case 9:
        {
            double cdot = cx * sx + cy * sy + cz * sz;
            nx = cx - cdot * sx;
            ny = cy - cdot * sy;
            nz = cz - cdot * sz;

            normalize(nx, ny, nz, 0.0, 0.0, 1.0);

            fx = sy * nz - sz * ny;
            fy = sz * nx - sx * nz;
            fz = sx * ny - sy * nx;
            normalize(fx, fy, fz, 1.0, 0.0, 0.0);

            break;
        }
        case 4:
        case 6:
        case 7:
        case 10:
        {
            nx = dphi[1][0];
            ny = dphi[1][1];
            nz = dphi[1][2];
            normalize(nx, ny, nz, 0.0, 0.0, 1.0);

            fx = sy * nz - sz * ny;
            fy = sz * nx - sx * nz;
            fz = sx * ny - sy * nx;
            normalize(fx, fy, fz, 1.0, 0.0, 0.0);
            break;
        }
        case 5:
        case 11:
        case 12:
        {
            nx = dphi[2][0];
            ny = dphi[2][1];
            nz = dphi[2][2];
            normalize(nx, ny, nz, 0.0, 0.0, 1.0);

            fx = sy * nz - sz * ny;
            fy = sz * nx - sx * nz;
            fz = sx * ny - sy * nx;
            normalize(fx, fy, fz, 1.0, 0.0, 0.0);
            break;
        }
        default:
        {
            fx = -dphi[0][1];
            fy = dphi[0][0];
            fz = 0.0;
            normalize(fx, fy, fz, 1.0, 0.0, 0.0);

            nx = fy * sz - fz * sy;
            ny = fz * sx - fx * sz;
            nz = fx * sy - fy * sx;
            normalize(nx, ny, nz, 0.0, 0.0, 1.0);
            break;
        }
    }

    f[0] = fx;  f[1] = fy;  f[2] = fz;
    s[0] = sx;  s[1] = sy;  s[2] = sz;
    n[0] = nx;  n[1] = ny;  n[2] = nz;


}

