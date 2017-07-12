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
 * DynamicElasticity.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: srossi
 */

#include "Elasticity/DynamicElasticity.hpp"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/type_tensor.h"

// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

//
#include "libmesh/vector_value.h"
#include "libmesh/linear_solver.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"

#include "libmesh/perf_log.h"

#include <sys/stat.h>
#include "BoundaryConditions/BCData.hpp"
#include "Util/IO/io.hpp"
#include "Util/MapsToLibMeshTypes.hpp"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"

#include "Util/GenerateFibers.hpp"

//Materials
#include "Elasticity/Materials/LinearMaterial.hpp"
#include "Elasticity/Materials/IsotropicMaterial.hpp"
#include "Elasticity/Materials/BenNeohookean.hpp"

namespace BeatIt
{

typedef libMesh::TransientLinearImplicitSystem LinearSystem;

DynamicElasticity::DynamicElasticity(libMesh::EquationSystems& es, std::string system_name)
                : super(es, system_name)
{
    // TODO Auto-generated constructor stub
}

DynamicElasticity::~DynamicElasticity()
{
    // TODO Auto-generated destructor stub
}

void DynamicElasticity::setupSystem(std::string section)
{
    int dimension = M_equationSystems.get_mesh().mesh_dimension();
    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // DISPLACEMENT PART
    std::cout
                    << "* DYNAMIC ELASTICITY: Creating new System for the DYNAMIC ELASTICITY Solver"
                    << std::endl;
    LinearSystem& system = M_equationSystems.add_system<LinearSystem>(M_myName);

    std::string disp_name = "displacement";
    std::string vel_name = "velocity";

    int ord = M_datafile(section + "/order", 1);
    std::cout << "* DYNAMIC ELASTICITY: Reading displacement field - order: " << ord << std::flush;
    auto order_it = BeatIt::libmesh_order_map.find(ord);
    std::cout << " ... " << std::flush;
    libMesh::Order order =
                    (order_it != BeatIt::libmesh_order_map.end()) ?
                                    order_it->second :
                                    throw std::runtime_error("Order Explode!!!");
    std::cout << "order found!" << std::endl;

    std::string fam = M_datafile(section + "/fefamily", "lagrange");
    std::cout << "* DYNAMIC ELASTICITY: Reading displacement field - family: " << fam << std::flush;
    auto fefamily_it = BeatIt::libmesh_fefamily_map.find(fam);
    std::cout << " ... " << std::flush;
    libMesh::FEFamily fefamily =
                    (fefamily_it != BeatIt::libmesh_fefamily_map.end()) ?
                                    fefamily_it->second :
                                    throw std::runtime_error("FeFamily Explode!!!");
    std::cout << "fefamily found!" << std::endl;
    std::cout
                    << "* DYNAMIC ELASTICITY: Setting up displacement field - order: "
                    << ord
                    << ", fefamily = "
                    << fam
                    << std::endl;
    system.add_variable(vel_name + "x", order, fefamily);
    if (dimension > 1)
        system.add_variable(vel_name + "y", order, fefamily);
    if (dimension > 2)
        system.add_variable(vel_name + "z", order, fefamily);

    LinearSystem& disp_system = M_equationSystems.add_system<LinearSystem>(disp_name);
    disp_system.add_variable(disp_name + "x", order, fefamily);
    if (dimension > 1)
        disp_system.add_variable(disp_name + "y", order, fefamily);
    if (dimension > 2)
        disp_system.add_variable(disp_name + "z", order, fefamily);

    LinearSystem& pos_system = M_equationSystems.add_system<LinearSystem>("pos");
    pos_system.add_variable("posx", order, fefamily);
    if (dimension > 1)
        pos_system.add_variable("posy", order, fefamily);
    if (dimension > 2)
        pos_system.add_variable("posz", order, fefamily);


    // PRESSURE SYSTEM
    std::string formulation = M_datafile(section + "/formulation", "primal");
    std::cout << "* DYNAMIC ELASTICITY: Using a " << formulation << " formulation" << std::endl;
    std::string pressure_name = "pressure";

    if (formulation == "mixed")
    {
        throw std::runtime_error("MIXED FORMULATION NOT CODED IN DYNAMIC CASE");
        M_solverType = ElasticSolverType::Mixed;
        ord = M_datafile(section + "/p_order", 1);
        fam = M_datafile(section + "/p_fefamily", "lagrange");
        std::cout
                        << "* DYNAMIC ELASTICITY: Setting up pressure  field - order: "
                        << ord
                        << ", fefamily = "
                        << fam
                        << std::endl;
        libMesh::Order p_order = BeatIt::libmesh_order_map.find(ord)->second;
        libMesh::FEFamily p_fefamily = BeatIt::libmesh_fefamily_map.find(fam)->second;
        system.add_variable(pressure_name, p_order, p_fefamily);
    }
    else
    {
        M_solverType = ElasticSolverType::Primal;
        LinearSystem& p_system = M_equationSystems.add_system<LinearSystem>("Pressure_Projection");
        p_system.add_variable(pressure_name, libMesh::FIRST, libMesh::LAGRANGE);
        p_system.init();
    }

    system.add_vector("residual");
    system.add_vector("step");
    system.add_matrix("mass");

    std::cout << "* ELASTICITY: Reading BC ... " << std::endl;
    M_bch.readBC(M_datafile, section);
    M_bch.showMe();
    std::cout << "* ELASTICITY: Setup Homogenuous Dirichlet BC ... " << std::endl;

    for (auto&& bc_ptr : M_bch.M_bcs)
    {
        auto bc_type = bc_ptr->get_type();
        if (bc_type == BCType::Dirichlet)
        {
            std::set<libMesh::boundary_id_type> dirichlet_boundary_ids;
            auto bc_mode = bc_ptr->get_mode();
            std::vector<unsigned int> variables;

            switch (bc_mode)
            {
                default:
                case BCMode::Full:
                {
                    for (unsigned int k = 0; k < dimension; ++k)
                        variables.push_back(k);
                    break;
                }
                case BCMode::Component:
                {
                    auto component = bc_ptr->get_component();
                    switch (component)
                    {
                        case BCComponent::X:
                        {
                            variables.push_back(0);
                            break;
                        }
                        case BCComponent::Y:
                        {
                            variables.push_back(1);
                            break;
                        }
                        case BCComponent::Z:
                        {
                            variables.push_back(2);
                            break;
                        }
                        default:
                        {
                            break;
                        }
                    }
                    break;
                }
            }

			auto num_flags = bc_ptr->size();
			for(int nflag = 0; nflag < num_flags; nflag++)
			{
				dirichlet_boundary_ids.insert( bc_ptr->get_flag(nflag) );
			}

//            for (auto&& flag : bc_ptr->M_flag)
//            {
//                // Create a ZeroFunction to initialize dirichlet_bc
//                dirichlet_boundary_ids.insert(flag);
//            }

            libMesh::ZeroFunction<> zf;
            libMesh::DirichletBoundary dirichlet_bc(dirichlet_boundary_ids, variables, &zf);

            system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

        } // end if Dirichlet
    } // end for loop on BC
    system.init();
    system.get_dof_map().print_info(std::cout);
    auto ptr = system.get_dof_map().get_dirichlet_boundaries ();
    std::cout << "NUMBER DIRICHLET BC: " << ptr->size() << std::endl;
    std::cout << "* DYNAMIC ELASTICITY: init " << std::endl;
    disp_system.init();
    pos_system.init();

    // PRESSURE SYSTEM
    std::string time_integrator = M_datafile(section + "/time_integrator", "explicit");
    if ("implicit" == time_integrator)
        M_timeIntegratorType = DynamicTimeIntegratorType::Implicit;
    else
        M_timeIntegratorType = DynamicTimeIntegratorType::Explicit;
    if (fefamily != libMesh::FEFamily::LAGRANGE)
    {
        M_usingDG = true;
        DynamicTimeIntegratorType::Explicit;
    }
    else
        M_usingDG = false;
}

void DynamicElasticity::update_displacements(double dt)
{
    LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem& disp_system = M_equationSystems.get_system<LinearSystem>("displacement");
    *disp_system.solution = *system.solution;
    disp_system.solution->scale(dt);
    *disp_system.solution += *disp_system.old_local_solution;
    system.update();
    disp_system.update();
}

void DynamicElasticity::advance()
{
    LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem& disp_system = M_equationSystems.get_system<LinearSystem>("displacement");

    *system.old_local_solution = *system.solution;
    *system.older_local_solution = *system.solution;

    *disp_system.old_local_solution = *disp_system.solution;
    *disp_system.older_local_solution = *disp_system.solution;
}

void DynamicElasticity::assemble_residual(
                double dt,
                libMesh::NumericVector<libMesh::Number>* activation_ptr)
{
//  std::cout << "* DYNAMIC ELASTICITY: assembling ... " << std::endl;

    using libMesh::UniquePtr;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int max_dim = 3;
    // Get a reference to the LinearImplicitSystem we are solving
    LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
    double time = system.time;
    LinearSystem& disp_system = M_equationSystems.get_system<LinearSystem>("displacement");
    LinearSystem& pos_system = M_equationSystems.get_system<LinearSystem>("pos");
    ParameterSystem& fiber_system = M_equationSystems.get_system<ParameterSystem>("fibers");
    ParameterSystem& sheets_system = M_equationSystems.get_system<ParameterSystem>("sheets");
    ParameterSystem& xfiber_system = M_equationSystems.get_system<ParameterSystem>("xfibers");

//    std::cout << "Disp! " <<std::endl;
//    disp_system.solution->print(std::cout);
    system.get_vector("residual").zero();
    system.get_vector("step").zero();
    system.get_matrix("mass").zero();
    system.rhs->zero();
    system.matrix->zero();
//    system.solution->print(std::cout);
//    disp_system.solution->print(std::cout);

    unsigned int ux_var = system.variable_number("velocityx");
    unsigned int uy_var, uz_var;

    if (dim > 1)
        uy_var = system.variable_number("velocityy");
    if (dim > 2)
        uz_var = system.variable_number("velocityz");

    const libMesh::DofMap & dof_map = system.get_dof_map();
//    dof_map.print_info(std::cout);
//    libMesh::FEType fe_disp = dof_map.variable_type(ux_var);
    libMesh::FEType fe_disp(1, libMesh::LAGRANGE);

    UniquePtr<libMesh::FEBase> fe_u(libMesh::FEBase::build(dim, fe_disp));
    UniquePtr<libMesh::FEBase> fe_elem_face(libMesh::FEBase::build(dim, fe_disp));
    UniquePtr<libMesh::FEBase> fe_neighbor_face(libMesh::FEBase::build(dim, fe_disp));

    auto order = fe_u->get_order();

    libMesh::QGauss qrule_1(dim, libMesh::SECOND);

    fe_u->attach_quadrature_rule(&qrule_1);

    const std::vector<libMesh::Real> & JxW_u = fe_u->get_JxW();
    const std::vector<libMesh::Point> & q_point_u = fe_u->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<libMesh::Real> > & phi_u = fe_u->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & dphi_u = fe_u->get_dphi();

    libMesh::DenseVector<libMesh::Number> Fe;
    libMesh::DenseVector<libMesh::Number> Fee;
    libMesh::DenseVector<libMesh::Number> Fen;
    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseMatrix<libMesh::Number> Kee;
    libMesh::DenseMatrix<libMesh::Number> Ken;
    libMesh::DenseMatrix<libMesh::Number> Kne;
    libMesh::DenseMatrix<libMesh::Number> Knn;
//    libMesh::TensorValue<libMesh::Number> Mtest;
//    libMesh::TensorValue<libMesh::Number> Mtrial;

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_ux;
    std::vector<libMesh::dof_id_type> dof_indices_uy;
    std::vector<libMesh::dof_id_type> dof_indices_uz;

    // Grad U
    std::vector<double> disp_k;
    std::vector<double> vel_k;
    std::vector<double> vel_n;
    libMesh::VectorValue<libMesh::Number> Uk;
    libMesh::VectorValue<libMesh::Number> Xk;
    libMesh::VectorValue<libMesh::Number> Xk_neighbor;
    libMesh::VectorValue<libMesh::Number> Vk;
    libMesh::VectorValue<libMesh::Number> Vn;
    libMesh::TensorValue<libMesh::Number> dU;
    // Grad U
    libMesh::TensorValue<libMesh::Number> dUk;
    // Strain
    libMesh::TensorValue<libMesh::Number> E;
    // Strain
    libMesh::TensorValue<libMesh::Number> Ek;
    // Identity
    libMesh::TensorValue<libMesh::Number> id;
    id(0, 0) = 1.0;
    id(1, 1) = 1.0;
    id(2, 2) = 1.0;
    // Stress
    libMesh::TensorValue<libMesh::Number> Sk;
    libMesh::TensorValue<libMesh::Number> Jk;
    libMesh::TensorValue<libMesh::Number> Jkn;
    // Stress
    libMesh::TensorValue<libMesh::Number> S;
    // Grad W (test function)
    libMesh::TensorValue<libMesh::Number> dW;

    libMesh::VectorValue<libMesh::Number> g;
    libMesh::VectorValue<libMesh::Number> trial;
    libMesh::VectorValue<libMesh::Number> test;
    libMesh::TensorValue<libMesh::Number> dtrial;

    libMesh::TensorValue<libMesh::Number> dUk_neighbor;
    libMesh::VectorValue<libMesh::Number> Vk_neighbor;
    libMesh::TensorValue<libMesh::Number> Sk_neighbor;
    libMesh::TensorValue<libMesh::Number> Skaverage;

    libMesh::VectorValue<libMesh::Number> Uk_neighbor;
    libMesh::TensorValue<libMesh::Number> dS;
    libMesh::TensorValue<libMesh::Number> trialN;
    libMesh::TensorValue<libMesh::Number> testN;
    libMesh::TensorValue<libMesh::Number> UkN;
//    libMesh::VectorValue<libMesh::Number> JVk;



    double body_force[3];
    const std::vector<libMesh::Point> & q_point = fe_u->get_xyz();

    double gamma_f;
    double gamma_s;
    double gamma_n;
    std::vector<double> gamma_f_k;
    libMesh::TensorValue<libMesh::Number> FA;

    // On the boundary
    libMesh::UniquePtr<libMesh::FEBase> fe_face(libMesh::FEBase::build(dim, fe_disp));
    libMesh::QGauss qface(dim - 1, libMesh::SECOND);
    libMesh::QGauss qface_neighbor(dim - 1, libMesh::SECOND);
    fe_face->attach_quadrature_rule(&qface);
    fe_neighbor_face->attach_quadrature_rule(&qface);

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    double rho;
    libMesh::RealGradient f0;
    libMesh::RealGradient s0;
    libMesh::RealGradient n0;
    std::vector<libMesh::dof_id_type> neighbor_dof_indices;

    std::vector<double> disp_k_neighbor;
    std::vector<double> vel_k_neighbor;
    std::vector<double> x_k;
    std::vector<double> x_k_neighbor;

    // Data for surface integrals on the neighbor boundary
    const std::vector<std::vector<libMesh::Real> > & phi_neighbor_face =
                    fe_neighbor_face->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > & dphi_neighbor_face = fe_neighbor_face
                    ->get_dphi();
//	    std::cout << "* ELASTICITY:loop starts ... " << std::endl;

    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        rho = M_materialMap[0]->M_density;
        auto elID = elem->id();

//        for(unsigned int nn = 0; nn < 3; nn++)
//        {
//            const libMesh::Node& node = elem->node_ref(nn);
//            node.print_info(std::cout);
//        }
        dof_map.dof_indices(elem, dof_indices);

        dof_map.dof_indices(elem, dof_indices_ux, ux_var);
        if (dim > 1)
            dof_map.dof_indices(elem, dof_indices_uy, uy_var);
        if (dim > 2)
            dof_map.dof_indices(elem, dof_indices_uz, uz_var);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_ux_dofs = dof_indices_ux.size();
        fe_u->reinit(elem);
        Fe.resize(n_dofs);
        Ke.resize(n_dofs, n_dofs);
        Me.resize(n_dofs, n_dofs);

        // get uk
        disp_k.resize(n_dofs);
        vel_k.resize(n_dofs);
        vel_n.resize(n_dofs);

        disp_system.current_local_solution->get(dof_indices, disp_k);
        system.current_local_solution->get(dof_indices, vel_k);
        system.old_local_solution->get(dof_indices, vel_n);

        if (activation_ptr)
        {
            gamma_f_k.resize(n_ux_dofs);
            for (int nd = 0; nd < n_ux_dofs; nd++)
            {
                int index = dof_indices_ux[nd];
                gamma_f_k[nd] = (*activation_ptr)(index);
            }
        }
        // local solution for a triangle contains
        /*
         *                        solution_k = [ ux1, ux2, ux3, uy1, uy2, uy3, uz1, uz2, uz3, p1, p2, p3];
         */
        /*
         *            Jacobian =  | K_ux_wx        K_uy_wx        K_uz_wx        B_p_wx   |
         *                               | K_ux_wy        K_uy_wy        K_uz_wy        B_p_wy   |
         *                               | K_ux_wz        K_uy_wz        K_uz_wz        B_p_wz  |
         *                               | D_ux_q          D_uy_q          D_uz_q          C_p_q     |
         *
         *
         *				residual=  |  R_wx  |
         *								 |  R_wy  |
         *								 |  R_wx  |
         *								 |  R_q    |
         *
         */

        // Block K
        int index = 0;
        for (unsigned int qp = 0; qp < qrule_1.n_points(); qp++)
        {

            dUk *= 0.0;
            Ek *= 0.0;
            Sk *= 0.0;
            gamma_f *= 0.0;
            FA *= 0.0;
            Uk *= 0.0;
            Vk *= 0.0;
            Vn *= 0.0;
            const unsigned int n_phi = phi_u.size();
            for (unsigned int l = 0; l < n_phi; l++)
            {
                for (int idim = 0; idim < dim; idim++)
                {
                    Uk(idim) += phi_u[l][qp] * disp_k[l + idim * n_phi];
                    Vk(idim) += phi_u[l][qp] * vel_k[l + idim * n_phi];
                    Vn(idim) += phi_u[l][qp] * vel_n[l + idim * n_phi];
                    for (int jdim = 0; jdim < dim; jdim++)
                    {
                        dUk(idim, jdim) += dphi_u[l][qp](jdim) * disp_k[l + idim * n_phi];
                    }
                }
            }

            if (activation_ptr)
            {
                for (unsigned int l = 0; l < n_phi; ++l)
                {
                    gamma_f += phi_u[l][qp] * gamma_f_k[l];
                }
                gamma_s = 1.0 / std::sqrt(1.0 + gamma_f) - 1.0;
                gamma_n = gamma_s;
                for (int idim = 0; idim < dim; idim++)
                {
                    FA(idim, idim) += 1.0;
                    for (int jdim = 0; jdim < dim; jdim++)
                    {
                        FA(idim, jdim) += gamma_f * f0(idim) * f0(jdim) + gamma_s * s0(idim) * s0(
                                        jdim) + gamma_n * n0(idim) * n0(jdim);
                    }
                }
            }

            M_materialMap[0]->M_gradU = dUk;
            M_materialMap[0]->M_FA = FA;
            M_materialMap[0]->updateVariables();
            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
            Sk = M_materialMap[0]->M_PK1;
            // Residual
            const double x = q_point[qp](0);
            const double y = q_point[qp](1);
            const double z = q_point[qp](2);
            const double time = 0.0;
            for (int idim = 0; idim < dim; idim++)
                body_force[idim] = M_rhsFunction(time, x, y, z, idim);

            for (unsigned int n = 0; n < phi_u.size(); ++n)
            {
                for (int jdim = 0; jdim < dim; jdim++)
                {
                    dW *= 0.0;
                    dW(jdim, 0) = JxW_u[qp] * dphi_u[n][qp](0);
                    dW(jdim, 1) = JxW_u[qp] * dphi_u[n][qp](1);
                    dW(jdim, 2) = JxW_u[qp] * dphi_u[n][qp](2);
                    // Compute  - \nabla \cdot \sigma + f
                    Fe(n + jdim * n_ux_dofs) -= Sk.contract(dW);
                    Fe(n + jdim * n_ux_dofs) += JxW_u[qp] * rho * body_force[jdim] * phi_u[n][qp];
                    Fe(n + jdim * n_ux_dofs) +=
                                    JxW_u[qp] * rho * (Vn(jdim) - Vk(jdim)) / dt * phi_u[n][qp];
                }
            }

            // Matrix
            // for each dimension of the test function
            for (int jdim = 0; jdim < dim; jdim++)
            {
                for (unsigned int n = 0; n < phi_u.size(); ++n)
                {
                    for (unsigned int m = 0; m < phi_u.size(); ++m)
                    {
                        int index1 = n + jdim * n_ux_dofs;
                        int index2 = m + jdim * n_ux_dofs;
                        if (M_timeIntegratorType == DynamicTimeIntegratorType::Explicit)
                        {
//                              index2 = index1;
                        }
                        Me(index1, index2) += rho / dt * JxW_u[qp] * phi_u[n][qp] * phi_u[m][qp];
                    }
                }
            }

            if (M_timeIntegratorType == DynamicTimeIntegratorType::Implicit)
            {
                // For each test function
                for (unsigned int n = 0; n < phi_u.size(); ++n)
                {
                    // for each dimension of the test function
                    for (int jdim = 0; jdim < dim; jdim++)
                    {
                        dW *= 0.0;
                        dW(jdim, 0) = JxW_u[qp] * dphi_u[n][qp](0);
                        dW(jdim, 1) = JxW_u[qp] * dphi_u[n][qp](1);
                        dW(jdim, 2) = JxW_u[qp] * dphi_u[n][qp](2);

                        // for each trial function
                        for (unsigned int m = 0; m < phi_u.size(); ++m)
                        {
                            // for each dimension of the trial function
                            for (int idim = 0; idim < dim; idim++)
                            {
                                dU *= 0.0;
                                E *= 0.0;
                                S *= 0.0;

                                dU(idim, 0) = dphi_u[m][qp](0);
                                dU(idim, 1) = dphi_u[m][qp](1);
                                dU(idim, 2) = dphi_u[m][qp](2);
                                // dV = dt * dU
                                dU *= dt;
                                M_materialMap[0]->evaluateJacobian(dU, 0.0);
                                S = M_materialMap[0]->M_total_jacobian;
                                auto int_1 = n + jdim * n_ux_dofs;
                                auto int_2 = m + idim * n_ux_dofs;
                                //                              std::cout << int_1 << ", " << int_2 << ", SdW: " << S.contract(dW) << std::endl;
                                Ke(n + jdim * n_ux_dofs, m + idim * n_ux_dofs) += S.contract(dW);
                            }
                        }
                    }
                }
            }
        }

        Ke += Me;

        apply_BC(elem, Ke, Fe, fe_face, qface, mesh, n_ux_dofs);
//        Fe.print(std::cout);
//        for(auto && id : dof_indices) std::cout << id << ", ";
//        std::cout << ""<< std::endl;
//
//        for(auto && id : dof_indices_ux) std::cout << id << ", ";
//        std::cout << ""<< std::endl;
        //
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            double mu = 96.0;
            double gamma_mu = M_datafile("gm",6.0);
            double kappa = 208.0;
            double gamma_kappa = M_datafile("gk",5.0);
            // On the boundary we use strong boundary conditions for the time being
            if (elem->neighbor(side) == libmesh_nullptr)
            {

                const unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, side);

                auto bc = M_bch.get_bc(boundary_id);


                if (bc)
                {
                    const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();
                    const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();
                    const std::vector<libMesh::Point> & qface_point = fe_face->get_xyz();
                    const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_face->get_dphi();
                    const std::vector<libMesh::Point>& normals = fe_face->get_normals();

                    auto bc_type = bc->get_type();
                    fe_face->reinit(elem, side);
                    int n_phi = phi_face.size();

                    if (BCType::NitscheSymmetric == bc_type)
                    {
                        double length = elem->side(side)->volume();
                        double area = elem->volume();
                        double hE = area / length;
                        //Fee.resize(n_dofs);

                        // loop over quadrature nodes
                        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                        {
//                            std::cout << "qp: " << qp << std::endl;
//                            for(int n = 0; n < n_phi; ++n) std::cout << "n: " << n  << ", phi: " << phi_face[n][qp] << std::endl;
                            Vk *= 0.0;
                            Uk *= 0.0;
                            dUk *= 0.0;
//                            std::cout << "Uk" << std::endl;
//                            Uk.print(std::cout);
//                            std::cout << "" << std::endl;

                            // Evaluate displacements velocities and gradients at the quadrature node
                            for (unsigned int l = 0; l < phi_face.size(); l++)
                            {
                                for (int idim = 0; idim < dim; idim++)
                                {
//                                    std::cout << "l: " << l << ", phi: " << phi_face[l][qp] << ", val: " << phi_face[l][qp] * disp_k[l + idim * n_phi] << ", index: " << l + idim * n_phi << std::endl;
                                    Vk(idim) += phi_face[l][qp] * vel_k[l + idim * n_phi];
                                    Uk(idim) += phi_face[l][qp] * disp_k[l + idim * n_phi];
                                    for (int jdim = 0; jdim < dim; jdim++)
                                    {
                                        dUk(idim, jdim) +=
                                                        dphi_face[l][qp](jdim) * disp_k[l + idim * n_phi];
                                    }
                                }
                            }
//                            std::cout << "disp_k: " << std::flush;
//                            for(auto && dk : disp_k) std::cout << dk << ", ";
//                            std::cout << "Uk after" << std::endl;
//                            Uk.print(std::cout);
//                            std::cout << "" << std::endl;
                            // Position of the quadrature node
                            const double xq = qface_point[qp](0);
                            const double yq = qface_point[qp](1);
                            const double zq = qface_point[qp](2);

                            //Fill with the value of bc
                            g *= 0.0;
                            for (int pdim = 0; pdim < dim; ++pdim)
                                g(pdim) = bc->get_function()(0.0, xq, yq, zq, pdim);

//                            std::cout << "Middle: " << hE << std::endl;
//                            Fe.print(std::cout);

                            for (int idim = 0; idim < dim; ++idim)
                            {
                                for (unsigned int i = 0; i < phi_face.size(); i++)
                                {
                                    test *= 0.0;
                                    test(idim) = JxW_face[qp] * phi_face[i][qp];
                                    // \int g
                                    auto val = mu * gamma_mu / hE * (g - Uk) * test;
                                    auto val2 = Uk(idim) ;
//                                    std::cout << "Fee mu: " << val << ", val2: " << val2 << ", gm: " << gamma_mu << std::endl;
                                    Fe(i + idim * n_ux_dofs) +=
                                                    mu * gamma_mu / hE * (g - Uk) * test;

                                    double gdotn = g.contract(normals[qp]);
                                    double ukdotn = Uk.contract(normals[qp]);
                                    val = kappa * gamma_kappa / hE * (gdotn - ukdotn) * (test * normals[qp]);
//                                    std::cout << "Fee kappa: " << val << std::endl;
                                    Fe(i + idim * n_ux_dofs) +=
                                                    kappa * gamma_kappa / hE * (gdotn - ukdotn) * (test * normals[qp]);


                                    dW *= 0.0;
                                    for (int pdim = 0; pdim < dim; ++pdim)
                                        dW(idim, pdim) = JxW_face[qp] * dphi_face[i][qp](pdim);
//                                    std::cout << "dW: " << std::endl;
//                                    M_materialMap[0]->M_gradU = dW;
//                                    M_materialMap[0]->updateVariables();
//                                    M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
//                                    dS = M_materialMap[0]->M_PK1;
//                                    Fe(i + idim * n_ux_dofs) -= (dS * normals[qp]) * (g - Uk);
//                                    std::cout << "dUk: " << std::endl;
                                    M_materialMap[0]->M_gradU = dUk;
                                    M_materialMap[0]->updateVariables();
                                    M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
                                    Sk = M_materialMap[0]->M_PK1;
                                    Fe(i + idim * n_ux_dofs) += 0.0*(Sk * normals[qp]) * test;

                                    for (int jdim = 0; jdim < dim; ++jdim)
                                    {
                                        for (unsigned int j = 0; j < phi_face.size(); j++)
                                        {
                                            trial *= 0.0;
                                            trial(jdim) = dt * phi_face[j][qp];
                                            // \int g
                                            Ke(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) +=
                                                            mu * gamma_mu / hE * trial * test;

                                            double wdotn = trial.contract(normals[qp]);
                                            Ke(i + idim * n_ux_dofs, j + jdim * n_ux_dofs)  +=
                                                            kappa * gamma_kappa / hE *  wdotn * (test * normals[qp]);

//                                            Ke(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) -= (dS * normals[qp]) * trial;

                                            dU *= 0.0;
                                            S *= 0.0;

                                            dU(jdim, 0) = dphi_face[j][qp](0);
                                            dU(jdim, 1) = dphi_face[j][qp](1);
                                            dU(jdim, 2) = dphi_face[j][qp](2);
                                            // dU = dt * dV
                                            dU *= dt;
                                            M_materialMap[0]->evaluateJacobian(dU, 0.0);
                                            S = M_materialMap[0]->M_total_jacobian;

                                            Ke(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) -= 0.0*(S * normals[qp]) * test;
                                        } // end loop over trial functions
                                    } // end loop over dimension trial functions
                                }// end loop over test
//                                std::cout << "Fe: " << std::endl;
//                                Fe.print(std::cout);
                            } // end loop over dimension test function
                        }// end loop over qp
//                        std::cout << "After: " << std::endl;
//                        Fee.print(std::cout);
                        //Fe += Fee;

                    }
                }
            }
            else // interior egde
            {
                if (M_usingDG)
                {

                    double eta = M_datafile("eta", 1.0);
                    double beta = M_datafile("beta", 1.0);

                    // Store a pointer to the neighbor we are currently
                    // working on.
                    const libMesh::Elem * neighbor = elem->neighbor(side);

                    // Get the global id of the element and the neighbor
                    const unsigned int elem_id = elem->id();
                    const unsigned int neighbor_id = neighbor->id();


                    // If the neighbor has the same h level and is active
                    // perform integration only if our global id is bigger than our neighbor id.
                    // We don't want to compute twice the same contributions.
                    // If the neighbor has a different h level perform integration
                    // only if the neighbor is at a lower level.
                    if ((neighbor->active() && (neighbor->level() == elem->level()) && (elem_id < neighbor_id)) || (neighbor
                                    ->level() < elem->level()))
                    {
                        // Pointer to the element side
                        libMesh::UniquePtr<libMesh::Elem> elem_side(elem->build_side(side));

                        double length = elem->side(side)->volume();
                        double area = elem->volume();
                        double hE = area / length;

                        // The quadrature point locations on the neighbor side
                        std::vector<libMesh::Point> qface_neighbor_point;

                        // The quadrature point locations on the element side


                        const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();
                        const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();
                        const std::vector<libMesh::Point> & qface_point = fe_face->get_xyz();
                        const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_face->get_dphi();
                        const std::vector<libMesh::Point>& normals = fe_face->get_normals();
                        // Reinitialize shape functions on the element side
                        fe_face->reinit(elem, side);

                        // Find their locations on the neighbor
                        unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);

                        libMesh::FEInterface::inverse_map(elem->dim(), fe_u->get_fe_type(),
                                        neighbor,
                                        qface_point, qface_neighbor_point);


                        // Calculate the neighbor element shape functions at those locations
                        fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);

                        // Get the degree of freedom indices for the
                        // neighbor.  These define where in the global
                        // matrix this neighbor will contribute to.
                        dof_map.dof_indices(neighbor, neighbor_dof_indices);
                        const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
                        Fee.resize(n_dofs);
                        Fen.resize(n_neighbor_dofs);
                        Kee.resize(n_dofs,n_dofs);
                        Kne.resize(n_neighbor_dofs,n_dofs);
                        Ken.resize(n_dofs,n_neighbor_dofs);
                        Knn.resize(n_neighbor_dofs,n_neighbor_dofs);
//                        std::cout << "USING DG" << std::endl;

                        disp_k_neighbor.resize(n_neighbor_dofs);
                        vel_k_neighbor.resize(n_neighbor_dofs);
                        x_k.resize(n_dofs);
                        x_k_neighbor.resize(n_neighbor_dofs);

                        disp_system.current_local_solution->get(neighbor_dof_indices,
                                        disp_k_neighbor);
                        system.current_local_solution->get(neighbor_dof_indices, vel_k_neighbor);
                        //pos_system.solution->get(neighbor_dof_indices,
//                        x_k_neighbor);
//                        pos_system.solution->get(dof_indices, x_k);


                        double Vk_normal_jump = 0;

                        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                        {
                            // Position of the quadrature node
                            const double xq = qface_point[qp](0);
                            const double yq = qface_point[qp](1);
                            if(xq > 4-1e-12 && xq < 4+1e-12 && yq > 10 - time && time < 4)
                            {
//                                std::cout << "skipping side: " << side << " on element " << elID << std::endl;
                                break;
                            }
                            if(xq > 4-1e-12 && xq < 4+1e-12 && yq > 6 && time >= 4)
                            {
//                                std::cout << "skipping side: " << side << " on element " << elID << std::endl;
                                break;
                            }
                            if( std::abs(yq - xq-2.0) < 1e-12 && yq < 6 && yq > 4 && time >= 5)
                            {
//                                std::cout << "skipping side: " << side << " on element " << elID << std::endl;
                                break;
                            }
                            if(xq > 2-1e-12 && xq < 2+1e-12 && yq > 2 && yq < 4 && time >= 6)
                            {
//                                std::cout << "skipping side: " << side << " on element " << elID << std::endl;
                                break;
                            }
//                            if(yq > 2-1e-12 && yq < 2+1e-12 && xq > 2.0  && xq < 2.0+2.0*(time-7.5) && time >= 7.5 && time < 9.5)
//                            {
////                                std::cout << "skipping side: " << side << " on element " << elID << std::endl;
//                                break;
//                            }
                            if(yq > 2-1e-12 && yq < 2+1e-12 && xq > 2.0  && xq < 6.0 && time >= 8.5)
                            {
//                                std::cout << "skipping side: " << side << " on element " << elID << std::endl;
                                break;
                            }
                            if( std::abs(yq - xq+4.0) < 1e-12 && yq < 2.0+2.0*(time-11) && yq > 2 && time >= 11)
                            {
//                                std::cout << "skipping side: " << side << " on element " << elID << std::endl;
                                break;
                            }

                            const double zq = qface_point[qp](2);
                            Uk *= 0.0;
                            Vk *= 0.0;
                            dUk *= 0.0;
                            Uk_neighbor *= 0.0;
                            dUk_neighbor *= 0.0;
                            Vk_neighbor *= 0.0;
                            dW *= 0.0;
                            Sk *= 0.0;
                            Sk_neighbor *= 0.0;
                            dS *= 0.0;
                            Xk *= 0.0;
                            Xk_neighbor *= 0.0;
                            Jk *= 0.0;
                            Jkn *= 0.0;
                            Skaverage *= 0.0;
                            const unsigned int n_phi_face = phi_face.size();
                            const unsigned int n_phi = phi_u.size();

                            for (unsigned int l = 0; l < phi_face.size(); l++)
                            {
                                for (int idim = 0; idim < dim; idim++)
                                {
                                    Vk(idim) += phi_face[l][qp] * vel_k[l + idim * n_phi];
                                    Uk(idim) += phi_face[l][qp] * disp_k[l + idim * n_phi];
                                    //Xk(idim) += phi_face[l][qp] * x_k[l + idim * n_phi];
                                    for (int jdim = 0; jdim < dim; jdim++)
                                    {
                                        dUk(idim, jdim) +=
                                                        dphi_face[l][qp](jdim) * disp_k[l + idim * n_phi];
                                    }
                                }
                            }
//                            JVk = Vk - Vk_neighbor;
//                            JVk.print(std::cout);
                            M_materialMap[0]->M_gradU = dUk;
                            M_materialMap[0]->updateVariables();
                            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
                            Sk = M_materialMap[0]->M_PK1;
                            Skaverage += 0.5 * Sk;
                            for (unsigned int l = 0; l < phi_neighbor_face.size(); l++)
                            {
                                for (int idim = 0; idim < dim; idim++)
                                {
                                    Vk_neighbor(idim) +=
                                                    phi_neighbor_face[l][qp] * vel_k_neighbor[l + idim * n_phi];
                                    Uk_neighbor(idim) +=
                                                    phi_neighbor_face[l][qp] * disp_k_neighbor[l + idim * n_phi];
                                   // Xk_neighbor(idim) += phi_neighbor_face[l][qp] * x_k_neighbor[l + idim * n_phi];
                                    for (int jdim = 0; jdim < dim; jdim++)
                                    {
                                        dUk_neighbor(idim, jdim) +=
                                                        dphi_neighbor_face[l][qp](jdim) * disp_k_neighbor[l + idim * n_phi];
                                    }
                                }
                            }


                            M_materialMap[0]->M_gradU = dUk_neighbor;
                            M_materialMap[0]->updateVariables();
                            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
                            Sk_neighbor = M_materialMap[0]->M_PK1;
                            Skaverage += 0.5 * Sk_neighbor;

                            double ukdotn = Uk.contract(normals[qp]);
                            double ukndotn = Uk_neighbor.contract(normals[qp]);

                            // HANSBO & CO. // NOT WORKING
//                            std::cout << "Assembling IP on side: " << side << " on element " << elID << std::endl;
                            for (int idim = 0; idim < dim; ++idim)
                            {
                                for (unsigned int i = 0; i < phi_face.size(); i++)
                                {
                                    test *= 0.0;
                                    test(idim) = JxW_face[qp] * phi_face[i][qp];
                                    // \int g
                                    Fee(i + idim * n_ux_dofs) += mu * gamma_mu / hE * (Uk_neighbor - Uk) * test;

                                    Fee(i + idim * n_ux_dofs) += kappa * gamma_kappa / hE * (ukndotn - ukdotn) * (test * normals[qp]);

                                    Fee(i + idim * n_ux_dofs) += (0.5 * (Sk + Sk_neighbor) * normals[qp]) * test;

//                                    dW *= 0.0;
//                                    dW(idim, 0) = JxW_face[qp] * dphi_face[i][qp](0);
//                                    dW(idim, 1) = JxW_face[qp] * dphi_face[i][qp](1);
//                                    dW(idim, 2) = JxW_face[qp] * dphi_face[i][qp](2);

//                                    M_materialMap[0]->M_gradU = dW;
//                                    M_materialMap[0]->updateVariables();
//                                    M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
//                                    dS = M_materialMap[0]->M_PK1;

//                                    Fe(i + idim * n_ux_dofs) -= 0.5 * (dS * normals[qp]) * (Uk_neighbor - Uk);

                                    for (int jdim = 0; jdim < dim; ++jdim)
                                    {
                                        for (unsigned int j = 0; j < phi_face.size(); j++)
                                        {
                                            // ++
                                            trial *= 0.0;
                                            trial(jdim) = dt * phi_face[j][qp];
                                            Kee(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) +=
                                                            mu * gamma_mu / hE * trial * test;
                                            double wdotn = trial.contract(normals[qp]);
                                            Kee(i + idim * n_ux_dofs, j + jdim * n_ux_dofs)  +=
                                                            kappa * gamma_kappa / hE *  wdotn * (test * normals[qp]);

                                            dU *= 0.0;
                                            S *= 0.0;
                                            dU(jdim, 0) = dphi_face[j][qp](0);
                                            dU(jdim, 1) = dphi_face[j][qp](1);
                                            dU(jdim, 2) = dphi_face[j][qp](2);
                                            // dU = dt * dV
                                            dU *= dt;
                                            M_materialMap[0]->M_gradU = dUk;
                                            M_materialMap[0]->updateVariables();
                                            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
                                            M_materialMap[0]->evaluateJacobian(dU, 0.0);
                                            S = M_materialMap[0]->M_total_jacobian;
                                            Kee(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) -= 0.5 * (S * normals[qp]) * test;

                                            // +-
                                            trial *= 0.0;
                                            trial(jdim) = dt * phi_neighbor_face[j][qp];
                                            Ken(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) -=
                                                            mu * gamma_mu / hE * trial * test;
                                            double wndotn = trial.contract(normals[qp]);
                                            Ken(i + idim * n_ux_dofs, j + jdim * n_ux_dofs)  -=
                                                            kappa * gamma_kappa / hE *  wndotn * (test * normals[qp]);

                                            dU *= 0.0;
                                            S *= 0.0;
                                            dU(jdim, 0) = dphi_neighbor_face[j][qp](0);
                                            dU(jdim, 1) = dphi_neighbor_face[j][qp](1);
                                            dU(jdim, 2) = dphi_neighbor_face[j][qp](2);
                                            // dU = dt * dV
                                            dU *= dt;
                                            M_materialMap[0]->M_gradU = dUk_neighbor;
                                            M_materialMap[0]->updateVariables();
                                            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
                                            M_materialMap[0]->evaluateJacobian(dU, 0.0);
                                            S = M_materialMap[0]->M_total_jacobian;
                                            Ken(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) -= 0.5 * (S * normals[qp]) * test;

                                        }
                                    }

                                    test *= 0.0;
                                    test(idim) = JxW_face[qp] * phi_neighbor_face[i][qp];
                                    // \int g
                                    Fen(i + idim * n_ux_dofs) -= mu * gamma_mu / hE * (Uk_neighbor - Uk) * test;

                                    Fen(i + idim * n_ux_dofs) -= kappa * gamma_kappa / hE * (ukndotn - ukdotn) * (test * normals[qp]);
//
                                    Fen(i + idim * n_ux_dofs) -= (0.5 * (Sk + Sk_neighbor) * normals[qp]) * test;

                                    for (int jdim = 0; jdim < dim; ++jdim)
                                    {
                                        for (unsigned int j = 0; j < phi_face.size(); j++)
                                        {
                                            // ++
                                            trial *= 0.0;
                                            trial(jdim) = dt * phi_face[j][qp];
                                            Kne(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) -=  mu * gamma_mu / hE * trial * test;
                                            double wdotn = trial.contract(normals[qp]);
                                            Kne(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) -=  kappa * gamma_kappa / hE * wdotn * (test * normals[qp]);

                                            dU *= 0.0;
                                            S *= 0.0;
                                            dU(jdim, 0) = dphi_face[j][qp](0);
                                            dU(jdim, 1) = dphi_face[j][qp](1);
                                            dU(jdim, 2) = dphi_face[j][qp](2);
                                            // dU = dt * dV
                                            dU *= dt;
                                            M_materialMap[0]->M_gradU = dUk;
                                            M_materialMap[0]->updateVariables();
                                            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
                                            M_materialMap[0]->evaluateJacobian(dU, 0.0);
                                            S = M_materialMap[0]->M_total_jacobian;
                                            Kne(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) += 0.5 * (S * normals[qp]) * test;

                                            // +-
                                            trial *= 0.0;
                                            trial(jdim) = dt * phi_neighbor_face[j][qp];
                                            Knn(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) += mu * gamma_mu / hE * trial * test;
                                            double wndotn = trial.contract(normals[qp]);
                                            Knn(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) +=  kappa * gamma_kappa / hE * wndotn * (test * normals[qp]);

                                            dU *= 0.0;
                                            S *= 0.0;
                                            dU(jdim, 0) = dphi_neighbor_face[j][qp](0);
                                            dU(jdim, 1) = dphi_neighbor_face[j][qp](1);
                                            dU(jdim, 2) = dphi_neighbor_face[j][qp](2);
                                            // dU = dt * dV
                                            dU *= dt;
                                            M_materialMap[0]->M_gradU = dUk_neighbor;
                                            M_materialMap[0]->updateVariables();
                                            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
                                            M_materialMap[0]->evaluateJacobian(dU, 0.0);
                                            S = M_materialMap[0]->M_total_jacobian;
                                            Knn(i + idim * n_ux_dofs, j + jdim * n_ux_dofs) += 0.5 * (S * normals[qp]) * test;

                                        }
                                    }
                                }
                            }

//                            UkN *= 0.0;
//                            for(int p = 0; p < dim; ++p)
//                            {
//                                for(int q = 0; q < dim; ++q)
//                                {
//                                    UkN(p,q) = (Uk(p) - Uk_neighbor(p)) * normals[qp](q);
//                                }
//                            }
//
//                            M_materialMap[0]->M_gradU = 0.5*(dUk + dUk_neighbor);
//                            M_materialMap[0]->updateVariables();
//                            M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
//                            M_materialMap[0]->evaluateJacobian(UkN, 0.0);
//                            S = M_materialMap[0]->M_total_jacobian;

//                            // RADOVITSKY & NOELS
//                            for (int idim = 0; idim < dim; ++idim)
//                            {
//                                for (unsigned int i = 0; i < phi_face.size(); i++)
//                                {
//
//                                    test *= 0.0;
//                                    test(idim) = JxW_face[qp] * phi_face[i][qp];
//                                    Fee(i + idim * n_ux_dofs) += eta *(Skaverage * normals[qp]) * test;
//                                    testN *= 0.0;
//                                    for(int p = 0; p < dim; ++p)
//                                    {
//                                        for(int q = 0; q < dim; ++q)
//                                        {
//                                            testN(p,q) = test(p) * normals[qp](q);
//                                        }
//                                    }
//                                    Fee(i + idim * n_ux_dofs) += beta / hE * testN.contract(S);
//
//                                   test *= 0.0;
//                                   test(idim) = JxW_face[qp] * phi_neighbor_face[i][qp];
//                                   Fen(i + idim * n_ux_dofs) -= eta * (Skaverage * normals[qp]) * test;                                   testN *= 0.0;
//                                   for(int p = 0; p < dim; ++p)
//                                   {
//                                       for(int q = 0; q < dim; ++q)
//                                       {
//                                           testN(p,q) = test(p) * normals[qp](q);
//                                       }
//                                   }
//                                   Fen(i + idim * n_ux_dofs) -= beta / hE * testN.contract(S);
//
//
//                                }
//                            }

//                            // IP
//                            for (int idim = 0; idim < dim; ++idim)
//                            {
//                                for (unsigned int i = 0; i < phi_face.size(); i++)
//                                {
//                                    double eta = M_datafile("eta", 33.33);
//
//                                    test *= 0.0;
//                                    test(idim) = JxW_face[qp] * phi_face[i][qp];
//                                    Fee(i + idim * n_ux_dofs) -= eta / hE * (Uk_neighbor - Uk) * test;
//                                    testN *= 0.0;
//                                    for(int p = 0; p < dim; ++p)
//                                    {
//                                        for(int q = 0; q < dim; ++q)
//                                        {
//                                            testN(p,q) = test(p) * normals[qp](q);
//                                        }
//                                    }
//                                    Fee(i + idim * n_ux_dofs) += testN.contract(0.5 * (Sk + Sk_neighbor));
//
//                                    dW(idim, 0) = JxW_face[qp] * dphi_face[i][qp](0);
//                                    dW(idim, 1) = JxW_face[qp] * dphi_face[i][qp](1);
//                                    dW(idim, 2) = JxW_face[qp] * dphi_face[i][qp](2);
//                                    M_materialMap[0]->M_gradU = dW;
//                                    M_materialMap[0]->updateVariables();
//                                    M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
//                                    dS = M_materialMap[0]->M_PK1;
//                                    Fee(i + idim * n_ux_dofs) += 0.5 * (dS * normals[qp])  * (Uk_neighbor - Uk);
//
//
//
//                                    test *= 0.0;
//                                    test(idim) = JxW_face[qp] * phi_neighbor_face[i][qp];
//                                    Fen(i + idim * n_ux_dofs) += eta / hE * (Uk_neighbor - Uk) * test;
//                                    testN *= 0.0;
//                                    for(int p = 0; p < dim; ++p)
//                                    {
//                                        for(int q = 0; q < dim; ++q)
//                                        {
//                                            testN(p,q) = test(p) * normals[qp](q);
//                                        }
//                                    }
//                                    Fen(i + idim * n_ux_dofs) -= testN.contract(0.5 * (Sk + Sk_neighbor));
//
//                                    dW(idim, 0) = JxW_face[qp] * dphi_neighbor_face[i][qp](0);
//                                    dW(idim, 1) = JxW_face[qp] * dphi_neighbor_face[i][qp](1);
//                                    dW(idim, 2) = JxW_face[qp] * dphi_neighbor_face[i][qp](2);
//                                    M_materialMap[0]->M_gradU = dW;
//                                    M_materialMap[0]->updateVariables();
//                                    M_materialMap[0]->evaluateStress(ElasticSolverType::Primal);
//                                    dS = M_materialMap[0]->M_PK1;
//                                    Fen(i + idim * n_ux_dofs) -= 0.5 * (dS * normals[qp])  * (Uk_neighbor - Uk);
//
//                                }
//                            } // end IP

                        } // end dg qp loop

                        dof_map.constrain_element_vector(Fee, dof_indices);
                        dof_map.constrain_element_vector(Fen, neighbor_dof_indices);

                        system.rhs->add_vector(Fee, dof_indices);
                        system.rhs->add_vector(Fen, neighbor_dof_indices);
//                        std::cout << "BEFORE ADDING KEE" << std::endl;
                        system.matrix->add_matrix(Kee, dof_indices);
//                        std::cout << "BEFORE ADDING KEN" << std::endl;
                        system.matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
                        system.matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
                        system.matrix->add_matrix(Knn, neighbor_dof_indices);
//                        std::cout << "AFTER ADDING KEN" << std::endl;
                        //system.rhs->add_vector(Fen, dof_indices);
                    } // end if active element

                } // end M_usingDG
            } // end boundary if

        } // end for loop on edges
//        std::cout << "Constrainig" << std::endl;
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
      dof_map.constrain_element_matrix_and_vector (Me, Fe, dof_indices);
//      Me.print(std::cout);
        system.matrix->add_matrix(Ke, dof_indices);
        system.get_matrix("mass").add_matrix(Me, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);


    } // end loop over elements
    system.matrix->close();
    system.get_matrix("mass").close();
    system.rhs->close();

}

void DynamicElasticity::init_exo_output(const std::string& output_filename)
{
    if (M_usingDG)
    {
        M_saveStep = 0;
        M_exporter->write_discontinuous_exodusII(M_outputFolder + output_filename +"-s.0000",
                        M_equationSystems);
        M_saveStep++;

        //M_exporter->append(true);
    }
    else
        super::init_exo_output(output_filename);
}

void
DynamicElasticity::save_exo(
        const std::string& output_filename,
        int step,
        double time)
{

    if(M_usingDG)
    {
        std::cout << "* DYNAMIC ELASTICITY: EXODUSII::Exporting  Discontinuous data time " << time << ", step " << M_saveStep << " in: " << M_outputFolder << " ... " << std::flush;
        M_exporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );
        std::string save_step = std::to_string(M_saveStep);
        std::string pre = "-s.";
        int size = save_step.size();
        for(int j = 0; j < 4-size; ++j) pre += "0";
        M_exporter->write_discontinuous_exodusII(M_outputFolder + output_filename+pre+save_step,
                        M_equationSystems);
        M_saveStep++;
    }
    else super::save_exo(output_filename, step, time);
}
void
DynamicElasticity::save(const std::string& output_filename, int step)
{
    if (M_usingDG)
    {
//        M_GMVexporter->write_discontinuous_gmv  ( M_outputFolder+output_filename+"."+std::to_string(step),
//                                                       M_equationSystems, true);
        M_VTKexporter->write_equation_systems( M_outputFolder+output_filename+"."+std::to_string(step), M_equationSystems);
    }
    else
        super::save(output_filename, step);

}


void DynamicElasticity::project_pressure()
{
    std::cout << "* ELASTICITY: projecting pressure ... " << std::endl;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearSystem & system = M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem & disp_system = M_equationSystems.get_system<LinearSystem>("displacement");
    LinearSystem & system_p = M_equationSystems.get_system<LinearSystem>("Pressure_Projection");
    system.update();
    unsigned int ux_var = disp_system.variable_number("displacementx");
    unsigned int uy_var, uz_var, p_var;

    if (dim > 1)
        uy_var = disp_system.variable_number("displacementy");
    if (dim > 2)
        uz_var = disp_system.variable_number("displacementz");
    p_var = system_p.variable_number("pressure");

    const libMesh::DofMap & dof_map = system.get_dof_map();
    libMesh::FEType fe_type = dof_map.variable_type(0);
    libMesh::UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule(dim, fe->get_order());
    fe->attach_quadrature_rule(&qrule);

    const libMesh::DofMap & dof_map_p = system_p.get_dof_map();
    libMesh::FEType fe_type_p = dof_map_p.variable_type(0);
    libMesh::UniquePtr<libMesh::FEBase> fe_p(libMesh::FEBase::build(dim, fe_type_p));
    libMesh::QGauss qrule_p(dim, libMesh::SECOND);
    fe_p->attach_quadrature_rule(&qrule_p);

    const std::vector<libMesh::Real> & JxW = fe_p->get_JxW();
    const std::vector<std::vector<libMesh::RealGradient> > & dphi = fe_p->get_dphi();
    const std::vector<std::vector<libMesh::Real> > & phi = fe_p->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > & dphi_u = fe->get_dphi();
    const std::vector<std::vector<libMesh::Real> > & phi_u = fe->get_phi();

    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseVector<libMesh::Number> Fe;
    // Grad U
    // Grad U
    std::vector<double> solution_k;

    libMesh::TensorValue<libMesh::Number> dUk;

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_ux;
    std::vector<libMesh::dof_id_type> dof_indices_uy;
    std::vector<libMesh::dof_id_type> dof_indices_uz;
    std::vector<libMesh::dof_id_type> dof_indices_p;

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;

        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_ux, ux_var);
        if (dim > 1)
        {
            dof_map.dof_indices(elem, dof_indices_uy, uy_var);
        }
        if (dim > 2)
        {
            dof_map.dof_indices(elem, dof_indices_uz, uz_var);
        }
        dof_map_p.dof_indices(elem, dof_indices_p, p_var);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_ux_dofs = dof_indices_ux.size();
        const unsigned int n_uy_dofs = dof_indices_uy.size();
        const unsigned int n_uz_dofs = dof_indices_uz.size();
        const unsigned int n_p_dofs = dof_indices_p.size();

        fe->reinit(elem);
        fe_p->reinit(elem);

        Ke.resize(dof_indices_p.size(), dof_indices_p.size());
        Fe.resize(dof_indices_p.size());

        // get uk
        solution_k.resize(n_dofs);
        //        for(auto && di : dof_indices)  std::cout << "dof id: " << di << std::endl;

        system.current_local_solution->get(dof_indices, solution_k);

        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
//              std::cout << "* ELASTICITY: evaluate dUk ... " << std::endl;

            dUk *= 0.0;
            const unsigned int n_phi = phi_u.size();
            for (unsigned int l = 0; l < n_phi; ++l)
            {
                for (int idim = 0; idim < dim; idim++)
                {
                    for (int jdim = 0; jdim < dim; jdim++)
                    {
                        dUk(idim, jdim) += dphi_u[l][qp](jdim) * solution_k[l + idim * n_phi];
                    }
                }
            }
//               std::cout << "* ELASTICITY: evaluate p ... " << std::endl;

            M_materialMap[0]->M_gradU = dUk;
            double p = M_materialMap[0]->evaluatePressure();

            for (unsigned int i = 0; i < phi.size(); i++)
            {
                Fe(i) += JxW[qp] * p * phi[i][qp];

                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    // mass term
                    //
                    //Ke(i, j) += JxW[qp] * phi[i][qp] * phi[j][qp];
                    //lumping
                    Ke(i, i) += JxW[qp] * phi[i][qp] * phi[j][qp];
                }
            }
        }
//           std::cout << "* ELASTICITY: add matrices ... " << std::endl;

        system_p.matrix->add_matrix(Ke, dof_indices_p);
        system_p.rhs->add_vector(Fe, dof_indices_p);

    }
//    std::cout << "* ELASTICITY: solving p ... " << std::endl;

    system_p.matrix->close();
    system_p.rhs->close();

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);

    rval = M_projectionsLinearSolver->solve(*system_p.matrix, nullptr, *system_p.solution,
                    *system_p.rhs, tol, max_iter);
}




} /* namespace BeatIt */
