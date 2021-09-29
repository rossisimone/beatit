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
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "Electromechanics/Electromechanics.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"


#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "Util/IO/io.hpp"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/quadrature_gauss.h"

#include <random>
#include "libmesh/vtk_io.h"
#include <iomanip>      // std::setw

enum class Subdomain : unsigned int{ TISSUE = 1,
                                     BATH  = 2};

struct IonicModel
{
    void setup(GetPot& data, std::string section = "")
    {
        v0 = data("v0", 0.0);
        v1 = data("v1", 0.1);
        v2 = data("v2", 1.0);
        k = data("k", 8.0);
    }

    double iion(double vn)
    {
        return k * (vn-v0) * (vn-v1) * (vn-v2);
    }

    void print()
    {
        std::cout << "Ionic model parameters:" << std::endl;
        std::cout << "k: " << k << std::endl;
        std::cout << "v0: " << v0 << std::endl;
        std::cout << "v1: " << v1 << std::endl;
        std::cout << "v2: " << v2 << std::endl;

    }

    double k;
    double v0, v1, v2;
};


struct Pacing
{
    void setup(GetPot& data, std::string section = "")
    {
        duration = data(section+"/duration", 2.0);
        amplitude = data(section+"/amplitude", 50.0);
        start_time = data(section+"/start_time", 0.0);
        double maxx = data(section+"/maxx", 1.0);
        double maxy = data(section+"/maxy", 1.0);
        double maxz = data(section+"/maxz", 1.0);
        double minx = data(section+"/minx", 0.0);
        double miny = data(section+"/miny", 0.0);
        double minz = data(section+"/minz", 0.0);
        box = libMesh::BoundingBox(libMesh::Point(minx, miny, minz), libMesh::Point(maxx, maxy, maxz));
    }

    double istim(const libMesh::Point& p, double time)
    {
        double I = 0.0;
        if(time > start_time && time <= start_time+ duration)
        {
            if(box.contains_point(p))
            {
                I = amplitude;
            }
        }
        return I;
    }

    void print()
    {
        std::cout << "Pacing parameters:" << std::endl;
        std::cout << "start_time: " << start_time << std::endl;
        std::cout << "duration: " << duration << std::endl;
        std::cout << "amplitude: " << amplitude << std::endl;
        std::cout << "min: " << std::endl;
        box.min().print();
        std::cout << "\nmax: " << std::endl;
        box.max().print();
        std::cout << std::endl;
    }

    double duration;
    double amplitude;
    double start_time;
    libMesh::BoundingBox box;
};

int main(int argc, char ** argv)
{
    // Read input file
    GetPot data = BeatIt::readInputFile(argc, argv);
    // Use libMesh stuff without writing libMesh:: everytime
    using namespace libMesh;
    // Initialize libMesh
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // save output in folder
    std::string output_folder = data("output_folder", "Output");
    BeatIt::createOutputFolder(init.comm(), output_folder);

    // Create empty mesh
    Mesh mesh(init.comm());
    // Create Mesh
    // number of elements in the x,y and z direction
    int nelx = data("nelx", 10);
    int nely = data("nely", 10);
    // if nelz = 0 we run in 2D
    int nelz = data("nelz", 10);

    // the cube dimensions are defined as [minx, maxx] x [miny, maxy] x [minz, maxz]
    double maxx = data("maxx", 1.0);
    double maxy = data("maxy", 1.0);
    double maxz = data("maxz", 1.0);

    double minx = data("minx", 0.0);
    double miny = data("miny", 0.0);
    double minz = data("minz", 0.0);

    // Get mesh parameters
    std::string eltype = data("eltype", "simplex");
    std::map<std::string, ElemType> elem_type_map = { {"TRI3", TRI3},
                                                 {"QUAD4", QUAD4},
                                                 {"TET4", TET4},
                                                 {"HEX8", HEX8} };
    auto elType = TRI3;
    auto elem_type_it = elem_type_map.find(eltype);
    if(elem_type_it != elem_type_map.end()) elType = elem_type_it->second;

    std::cout << "Creating the cube [" << minx << ", " << maxx << "] x ["
                                       << miny << ", " << maxy << "] x ["
                                       << minx << ", " << maxx << "] " << std::endl;
    std::cout << "Using " << nelx << " x " << nely << " x " << nelz << " elements." << std::endl;
    std::cout << "Element type " << elem_type_it->first << std::endl;

    // Create mesh
    MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, minx, maxx, miny, maxy, minz, maxz, elType);

    // setup subdomains
    // tissue domain box
    double tissue_maxx = data("tissue_maxx", 1.0);
    double tissue_maxy = data("tissue_maxy", 1.0);
    double tissue_maxz = data("tissue_maxz", 1.0);
    double tissue_minx = data("tissue_minx", 0.0);
    double tissue_miny = data("tissue_miny", 0.0);
    double tissue_minz = data("tissue_minz", 0.0);
    libMesh::BoundingBox box(libMesh::Point(tissue_minx, tissue_miny, tissue_minz), libMesh::Point(tissue_maxx, tissue_maxy, tissue_maxz));

    // right_interface
    for( auto el : mesh.element_ptr_range() )
    {
        auto c = el->centroid();
        // Are we in the tissue region?
        if( box.contains_point(c) )
        {
            el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
        }
        // we are not
        else
        {
            el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::BATH);
        }
    }
    // output the details about the mesh
    mesh.print_info();
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();


    // Information about the time, such as
    // current iteration, timestep, final time etc
    // can be conveniently read from the input file
    // and stored in a TimeData abject
    // In the input file specify
    //
    //    [time]
    //        # Timestep
    //        dt = 0.125      # Default: 1.0
    //        # Simulation initial time
    //        init_time = 0.0 # Default: 0.0
    //        # Simulation end time
    //        final_time = 10  # Default: 1.0
    //        # Maximum number of timesteps
    //        max_iter = 200000000   # Default: 99999999
    //        # Export the solution every save_iter iterations
    //        save_iter = 8   # Default: 1
    //    [../]
    //
    // Create the TimeData object
    BeatIt::TimeData datatime;
    // Set it up using the input file
    datatime.setup(data, "");
    // Output on screen the stored variables
    datatime.print();

    // Create libMesh Equations systems
    // This will hold the mesh and create the corresponding
    // finite element spaces
    libMesh::EquationSystems es(mesh);

    libMesh::LinearImplicitSystem& elliptic = es.add_system<libMesh::LinearImplicitSystem>("elliptic");
    elliptic.add_variable("Ve", FIRST, LAGRANGE);

    libMesh::LinearImplicitSystem& parabolic = es.add_system<libMesh::LinearImplicitSystem>("parabolic");
    std::set< libMesh::subdomain_id_type > tissue_subdomains;
    tissue_subdomains.insert(static_cast<libMesh::subdomain_id_type>( Subdomain::TISSUE) );
    parabolic.add_variable("V", FIRST, LAGRANGE, &tissue_subdomains);
    parabolic.add_vector("Ve");
    parabolic.add_vector("iion");
    parabolic.add_vector("ML");
    parabolic.add_matrix("M");
    parabolic.add_matrix("Ke");
    parabolic.add_matrix("KeFace");
    parabolic.add_matrix("Ki");
    parabolic.add_vector("KiV");

    es.init();

    elliptic.zero_out_matrix_and_rhs = false;
    elliptic.assemble_before_solve = false;

    const DofMap & parabolic_dof_map = parabolic.get_dof_map();
    const DofMap & elliptic_dof_map = elliptic.get_dof_map();

    // ASSEMBLE MATRICES
    // right_interface
    {
        // zero out stuff
        elliptic.matrix->zero();
        parabolic.get_vector("ML").zero();
        parabolic.get_matrix("M").zero();
        parabolic.get_matrix("Ke").zero();
        parabolic.get_matrix("KeFace").zero();
        parabolic.get_matrix("Ki").zero();
        // get parameters
        double sigma_il = data("sigma_il", 1.0);
        double sigma_it = data("sigma_it", 1.0);
        double sigma_el = data("sigma_el", 1.0);
        double sigma_et = data("sigma_et", 1.0);
        double sigma_bb = data("sigma_b", 1.0);
        double chi = data("chi", 1000.0);
        double Cm = data("Cm", 1.0);

        double fx = data("fx", 1.0);
        double fy = data("fy", 0.0);
        double fz = data("fz", 0.0);



        FEType fe_type = parabolic_dof_map.variable_type(0);
        std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
        libMesh::QGauss qrule (dim, FIFTH);
        fe->attach_quadrature_rule (&qrule);

        std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
        QGauss qface(dim-1, FIFTH);
        fe_face->attach_quadrature_rule (&qface);


        const std::vector<Real> & JxW = fe->get_JxW();
        const std::vector<Point> & q_point = fe->get_xyz();
        const std::vector<std::vector<Real>> & phi = fe->get_phi();
        const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

        const std::vector<Real> & JxW_face = fe_face->get_JxW();
        const std::vector<Point> & q_point_face = fe_face->get_xyz();
        const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
        const std::vector<std::vector<RealGradient>> & dphi_face = fe_face->get_dphi();
        const std::vector<Point> & normal = fe_face->get_normals ();
        DenseMatrix<Number> Kie;
        DenseMatrix<Number> Ke;
        DenseMatrix<Number> Ke_face;
        DenseMatrix<Number> Ki;
        DenseMatrix<Number> Kip;
        DenseMatrix<Number> M;
        DenseVector<Number> M_lumped;

        std::vector<dof_id_type> parabolic_dof_indices;
        std::vector<dof_id_type> elliptic_dof_indices;

        for( auto elem : mesh.element_ptr_range() )
        {
            parabolic_dof_map.dof_indices (elem, parabolic_dof_indices);
            elliptic_dof_map.dof_indices (elem, elliptic_dof_indices);

            fe->reinit(elem);

            int elliptic_ndofs = elliptic_dof_indices.size();
            Kie.resize (elliptic_ndofs, elliptic_ndofs);
            int parabolic_ndofs = parabolic_dof_indices.size();
            Ke.resize (parabolic_ndofs, parabolic_ndofs);
            Ke_face.resize (parabolic_ndofs, parabolic_ndofs);
            Ki.resize (parabolic_ndofs, parabolic_ndofs);
            M.resize (parabolic_ndofs, parabolic_ndofs);
            M_lumped.resize (parabolic_ndofs);

            libMesh::TensorValue<double> eye(1.0, 0.0, 0.0,
                                             0.0, 1.0, 0.0,
                                             0.0, 0.0, 1.0);
            libMesh::TensorValue<double> fof(fx*fx, fx*fy, fx*fz,
                                             fy*fx, fy*fy, fy*fz,
                                             fz*fx, fz*fy, fz*fz);

            libMesh::TensorValue<double> sigma_i = sigma_il * fof + sigma_it * ( eye - fof );
            libMesh::TensorValue<double> sigma_e = sigma_el * fof + sigma_et * ( eye - fof );
            libMesh::TensorValue<double> sigma_b = sigma_bb * eye;

            libMesh::TensorValue<double> sigma_elliptic = sigma_bb * eye;
            if(elliptic_ndofs == parabolic_ndofs)
            {
                sigma_elliptic = sigma_i + sigma_e;
            }

            for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
            {
                for (unsigned int i = 0; i != elliptic_ndofs; i++)
                {
                    for (unsigned int j = 0; j != elliptic_ndofs; j++)
                    {
                        Kie(i, j) += JxW[qp] * dphi[i][qp] * (sigma_elliptic * dphi[j][qp]);
                    }
                }
                // we are in the tissue
                if (elliptic_ndofs == parabolic_ndofs)
                {
                    for (unsigned int i = 0; i != parabolic_ndofs; i++)
                    {
                        for (unsigned int j = 0; j != parabolic_ndofs; j++)
                        {
                            Ki(i, j) += JxW[qp] * dphi[i][qp] * (sigma_i * dphi[j][qp]);
                            Ke(i, j) += JxW[qp] * dphi[i][qp] * (sigma_e * dphi[j][qp]);
                            M(i, j) += JxW[qp] * chi * phi[i][qp] * phi[j][qp];
                            M_lumped(i) += JxW[qp] * Cm * chi / datatime.M_dt * phi[i][qp] * phi[j][qp];
                        }
                    }
                }
            }

            auto subdomain_id = elem->subdomain_id();
            if ( subdomain_id == static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE) )
            {
                for (auto side : elem->side_index_range())
                {
                    if (elem->neighbor_ptr(side) != nullptr)
                    {
                        auto neighbor_subdomain_id = elem->neighbor_ptr(side)->subdomain_id();
                        if (neighbor_subdomain_id != subdomain_id)
                        {
                            fe_face->reinit(elem, side);
                            for (unsigned int qp=0; qp<qface.n_points(); qp++)
                            {
                                for (unsigned int i=0; i != parabolic_ndofs; i++)
                                  for (unsigned int j=0; j != parabolic_ndofs; j++)
                                      Ke_face(i,j) -= JxW_face[qp] * ( sigma_e * dphi_face[j][qp] * normal[qp]) * phi_face[i][qp];
                            }
                        }
                    }
                }
            }
            elliptic.matrix->add_matrix (Kie, elliptic_dof_indices, elliptic_dof_indices);
            parabolic.get_vector("ML").add_vector(M_lumped, parabolic_dof_indices);
            parabolic.get_matrix("M").add_matrix(M, parabolic_dof_indices, parabolic_dof_indices);
            parabolic.get_matrix("Ke").add_matrix(Ke, parabolic_dof_indices, parabolic_dof_indices);
            parabolic.get_matrix("KeFace").add_matrix(Ke_face, parabolic_dof_indices, parabolic_dof_indices);
            parabolic.get_matrix("Ki").add_matrix(Ki, parabolic_dof_indices, parabolic_dof_indices);

        }
        elliptic.matrix->close();
        parabolic.get_vector("ML").close();
        parabolic.get_matrix("M").close();
        parabolic.get_matrix("Ke").close();
        parabolic.get_matrix("KeFace").close();
        parabolic.get_matrix("Ki").close();

//        elliptic.matrix->print();
//        parabolic.get_matrix("M").print();
//        parabolic.get_matrix("Ke").print();
//        parabolic.get_matrix("Ki").print();
        //elliptic.matrix->set(0, 0, 1e10);
    }

    // setup pacing protocol
    Pacing pacing;
    pacing.setup(data, "pacing");
    pacing.print();
    IonicModel ionic_model;
    ionic_model.setup(data, "ionic_model");
    ionic_model.print();
    // setup initial condition
    parabolic.solution->close();
    parabolic.solution->zero();
    parabolic.solution->add(ionic_model.v0);

    // save initial time:
    int save_iter = 0;
    {
        std::cout << "Exporting  ... " << std::endl;
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << save_iter;
        std::string step_str = ss.str();
        libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/bath_" + step_str + ".pvtu", es);
        std::cout << "done " << std::endl;
    }

    // method
    // 0: solve for Cm dV/dt = - nabla ( sigma_e / chi \nabla Ve) - Iion
    // 1: solve for Cm dV/dt = nabla ( sigma_i / chi \nabla V) + nabla ( sigma_i / chi \nabla Ve) - Iion
    int method = data("method", 0);
    // Start loop in time
    std::cout << "Time loop starts:" << std::endl;
    std::vector<dof_id_type> parabolic_dof_indices;
    std::vector<dof_id_type> elliptic_dof_indices;
    // Control the time loop using the TimeData object
    for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
    {
        // We are doing a new iteration
        // let's update first the information in the TimeData object
        datatime.advance();
        // Ouput to screen the current time and the current iteration
        std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter << std::endl;

        // assemble I ion and get Ve
        {
            parabolic.get_vector("iion").zero();

            for(auto node : mesh.node_ptr_range() )
            {
                parabolic_dof_map.dof_indices (node, parabolic_dof_indices);
                elliptic_dof_map.dof_indices (node, elliptic_dof_indices);
                int elliptic_ndofs = elliptic_dof_indices.size();
                int parabolic_ndofs = parabolic_dof_indices.size();

                if(elliptic_ndofs == parabolic_ndofs)
                {
                    double stimulus = pacing.istim(*node, datatime.M_time);
                    double vn = (*parabolic.solution)(parabolic_dof_indices[0]);
                    double i_ion = ionic_model.iion(vn);
                    double i_tot = - i_ion - stimulus;
                    parabolic.get_vector("iion").set(parabolic_dof_indices[0], i_tot);

                    double ven = (*elliptic.current_local_solution)(elliptic_dof_indices[0]);
                    parabolic.get_vector("Ve").set(parabolic_dof_indices[0], ven);
                }
            }
            parabolic.get_vector("iion").close();
            //std::cout << "iion" << std::endl;
            //parabolic.get_vector("iion").print();
            parabolic.get_vector("Ve").close();
            //std::cout << "Ve" << std::endl;
            //parabolic.get_vector("Ve").print();
        }


        parabolic.rhs->zero();
        // V^n+1 = V^n + ML^-1 * ( Ke * Ve - M * Iion)
        if(method == 0)
        {
            parabolic.get_matrix("KeFace").vector_mult_add(*parabolic.rhs, parabolic.get_vector("Ve"));
           // *parabolic.rhs *= -1.0;
            parabolic.get_matrix("Ke").vector_mult_add(*parabolic.rhs, parabolic.get_vector("Ve"));
        }
        // V^n+1 = V^n + ML^-1 * ( - Ki * Ve - Ki * V - M * Iion)
        else if(method == 1)
        {
            parabolic.get_matrix("Ki").vector_mult_add(*parabolic.rhs, *parabolic.solution);
            parabolic.get_matrix("Ki").vector_mult_add(*parabolic.rhs, parabolic.get_vector("Ve"));
            *parabolic.rhs *= -1.0;
        }

        parabolic.get_matrix("M").vector_mult_add(*parabolic.rhs, parabolic.get_vector("iion"));

        *parabolic.rhs /= parabolic.get_vector("ML");
        *parabolic.solution += *parabolic.rhs;

        //std::cout << "V n+1" << std::endl;
        //parabolic.solution->print();

        //std::cout << "KiV" << std::endl;
        //parabolic.get_vector("KiV").print();

        parabolic.add_vector("KiV").zero();
        parabolic.get_matrix("Ki").vector_mult_add(parabolic.add_vector("KiV"), *parabolic.solution);
        // setup rhs for elliptic problem
        {
            elliptic.rhs->zero();
            for(auto node : mesh.node_ptr_range() )
            {
                parabolic_dof_map.dof_indices (node, parabolic_dof_indices);
                elliptic_dof_map.dof_indices (node, elliptic_dof_indices);
                int elliptic_ndofs = elliptic_dof_indices.size();
                int parabolic_ndofs = parabolic_dof_indices.size();

                if(elliptic_ndofs == parabolic_ndofs)
                {
                    double rhs = parabolic.get_vector("KiV")(parabolic_dof_indices[0]);
                    elliptic.rhs->set(elliptic_dof_indices[0], -rhs);
                }
            }
            elliptic.rhs->close();
        }

        //std::cout << "- Ki V" << std::endl;
        //elliptic.rhs->print();
        elliptic.solve();
        //std::cout << "Ve n+1" << std::endl;
        //elliptic.solution->print();

        parabolic.update();
        elliptic.update();

        //std::cout << "Done ... " << std::endl;

        //Export the solution if at the right timestep
        if (0 == datatime.M_iter % datatime.M_saveIter)
        {
            std::cout << "Exporting  ... " << std::endl;
            // update output file counter
            save_iter++;

//            std::cout << "Exporting EXODUS  ... " << std::endl;
//            std::set<std::string> exported_variables;
//            exported_variables.insert("monowave");
//            exported_variables.insert("wave");
//            exported_variables.insert("istim");
//            exported_variables.insert(em.M_elasticity->M_myName);
//            exported_variables.insert("I4f");
//            exported_variables.insert("activation");
//
            std::ostringstream ss;
            ss << std::setw(4) << std::setfill('0') << save_iter;
            std::string step_str = ss.str();
            libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/bath_" + step_str + ".pvtu", es);
            std::cout << "done " << std::endl;
        }
    }
    // The end
    std::cout << "Good luck with your simulation :P" << std::endl;
    return 0;
}
