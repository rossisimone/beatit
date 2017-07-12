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
 * MixedDynamicElasticity.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: srossi
 */

#include "Elasticity/MixedDynamicElasticity.hpp"
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
typedef libMesh::TransientExplicitSystem ExplicitSystem;

MixedDynamicElasticity::MixedDynamicElasticity(libMesh::EquationSystems& es, std::string system_name)
                : super(es, system_name)
{
    // TODO Auto-generated constructor stub
}

MixedDynamicElasticity::~MixedDynamicElasticity()
{
    // TODO Auto-generated destructor stub
}

void MixedDynamicElasticity::setupSystem(std::string section)
{
    int dimension = M_equationSystems.get_mesh().mesh_dimension();
    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // DISPLACEMENT PART
    std::cout
                    << "* MIXED DYNAMIC ELASTICITY: Creating new System for the DYNAMIC ELASTICITY Solver"
                    << std::endl;
    LinearSystem& system = M_equationSystems.add_system<LinearSystem>(M_myName);

    std::string disp_name = "displacement";
    std::string vel_name = "velocity";

    int ord = M_datafile(section + "/order", 1);
    std::cout << "* MIXED DYNAMIC ELASTICITY: Reading displacement field - order: " << ord << std::flush;
    auto order_it = BeatIt::libmesh_order_map.find(ord);
    std::cout << " ... " << std::flush;
    libMesh::Order order =
                    (order_it != BeatIt::libmesh_order_map.end()) ?
                                    order_it->second :
                                    throw std::runtime_error("Order Explode!!!");
    std::cout << "order found!" << std::endl;

    std::string fam = M_datafile(section + "/fefamily", "lagrange");
    std::cout << "* MIXED DYNAMIC ELASTICITY: Reading displacement field - family: " << fam << std::flush;
    auto fefamily_it = BeatIt::libmesh_fefamily_map.find(fam);
    std::cout << " ... " << std::flush;
    libMesh::FEFamily fefamily =
                    (fefamily_it != BeatIt::libmesh_fefamily_map.end()) ?
                                    fefamily_it->second :
                                    throw std::runtime_error("FeFamily Explode!!!");
    std::cout << "fefamily found!" << std::endl;
    std::cout
                    << "* MIXED DYNAMIC ELASTICITY: Setting up displacement field - order: "
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

    // PRESSURE SYSTEM
    std::string formulation = M_datafile(section + "/formulation", "primal");
    std::cout << "* MIXED DYNAMIC ELASTICITY: Using a " << formulation << " formulation" << std::endl;
    std::string pressure_name = "pressure";

    M_solverType = ElasticSolverType::Mixed;
    ord = M_datafile(section + "/p_order", 1);
    fam = M_datafile(section + "/p_fefamily", "lagrange");
    std::cout
                    << "* MIXED DYNAMIC ELASTICITY: Setting up pressure  field - order: "
                    << ord
                    << ", fefamily = "
                    << fam
                    << std::endl;
    libMesh::Order p_order = BeatIt::libmesh_order_map.find(ord)->second;
    libMesh::FEFamily p_fefamily = BeatIt::libmesh_fefamily_map.find(fam)->second;
    M_solverType = ElasticSolverType::Primal;
    LinearSystem& p_system = M_equationSystems.add_system<LinearSystem>("pressure");
    p_system.add_variable(pressure_name, p_order, p_fefamily);
    p_system.add_vector("midpoint_pr", true, libMesh::GHOSTED);
    p_system.add_vector("step");
    p_system.init();


    system.add_vector("midpoint_vel", true, libMesh::GHOSTED);
    system.add_vector("rhs_vel");
    system.add_vector("step");

    disp_system.add_vector("midpoint_disp", true, libMesh::GHOSTED);
    disp_system.add_vector("rhs_disp");

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

            libMesh::ZeroFunction<> zf;
            libMesh::DirichletBoundary dirichlet_bc(dirichlet_boundary_ids, variables, &zf);

            system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

        } // end if Dirichlet
    } // end for loop on BC
    system.init();
    system.get_dof_map().print_info(std::cout);
    auto ptr = system.get_dof_map().get_dirichlet_boundaries ();
    std::cout << "NUMBER DIRICHLET BC: " << ptr->size() << std::endl;
    std::cout << "* MIXED DYNAMIC ELASTICITY: init " << std::endl;
    disp_system.init();

    // PRESSURE SYSTEM
    std::string time_integrator = M_datafile(section + "/time_integrator", "explicit");
    if ("implicit" == time_integrator)
    {
        std::cout << "implicit method not coded for mixed dynamic elasticity" << std::endl;
        throw std::runtime_error("implicit method not coded for mixed dynamic elasticity");
    }
    else
    {
        M_timeIntegratorType = DynamicTimeIntegratorType::Explicit;
    }

}

void
MixedDynamicElasticity::update_displacements(double dt)
{
    LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem& disp_system = M_equationSystems.get_system<LinearSystem>("displacement");

    //system.get_vector("midpoint_velocity") = 0.5 * (*system.solution + *system.old_local_solution);

    // u^n+1 = u^n + dt * (v^n+1 + v^n ) / 2
    system.get_vector("midpoint_vel").zero();
    system.get_vector("midpoint_vel").add(0.5, *system.current_local_solution);
    system.get_vector("midpoint_vel").add(0.5, *system.old_local_solution);
    disp_system.solution->add(dt, system.get_vector("midpoint_vel") );
}

void MixedDynamicElasticity::advance()
{
    LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem& disp_system = M_equationSystems.get_system<LinearSystem>("displacement");
    LinearSystem& p_system = M_equationSystems.add_system<LinearSystem>("pressure");

    std::cout << "update velocity" << std::endl;
    *system.older_local_solution  = *system.old_local_solution;
    *system.old_local_solution = *system.solution;
    std::cout << "update pressure" << std::endl;
    *p_system.older_local_solution = *p_system.old_local_solution;
    *p_system.old_local_solution = *p_system.solution;
    std::cout << "update disp 1" << std::endl;
    *disp_system.older_local_solution = *disp_system.old_local_solution;
    std::cout << "update disp 2" << std::endl;
    *disp_system.old_local_solution = *disp_system.solution;
    std::cout << "done" << std::endl;


}

void MixedDynamicElasticity::assemble_residual(
                double dt,
                libMesh::NumericVector<libMesh::Number>* activation_ptr)
{
  std::cout << "* MIXED DYNAMIC ELASTICITY: assembling ... " << std::endl;

    using libMesh::UniquePtr;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int max_dim = 3;
    // Get a reference to the LinearImplicitSystem we are solving
    LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
    double time = system.time;

    // Get Systems
    LinearSystem& disp_system = M_equationSystems.get_system<LinearSystem>("displacement");
    ParameterSystem& fiber_system = M_equationSystems.get_system<ParameterSystem>("fibers");
    ParameterSystem& sheets_system = M_equationSystems.get_system<ParameterSystem>("sheets");
    ParameterSystem& xfiber_system = M_equationSystems.get_system<ParameterSystem>("xfibers");
    LinearSystem& p_system = M_equationSystems.add_system<LinearSystem>("pressure");
    ParameterSystem& dummy_system       = M_equationSystems.get_system<ParameterSystem>("dumb");

    system.update();
    disp_system.update();
    p_system.update();


    // Zero out the vectors
//    system.get_vector("residual").zero();
    system.get_vector("step").zero();
    p_system.get_vector("step").zero();
    //system.get_matrix("mass").zero();
    system.rhs->zero();
    system.matrix->zero();


    // Define the midpoint vectors
    // at which we evaluate the residual
    disp_system.get_vector("midpoint_disp").zero();
    disp_system.get_vector("midpoint_disp").close();
    disp_system.get_vector("midpoint_disp").add(0.5, *disp_system.current_local_solution);
    disp_system.get_vector("midpoint_disp").add(0.5, *disp_system.old_local_solution);
    p_system.get_vector("midpoint_pr").zero();
    disp_system.get_vector("midpoint_disp").close();
    p_system.get_vector("midpoint_pr").add(0.5, *p_system.current_local_solution);
    p_system.get_vector("midpoint_pr").add(0.5, *p_system.old_local_solution);
    system.get_vector("midpoint_vel").zero();
    disp_system.get_vector("midpoint_disp").close();
    system.get_vector("midpoint_vel").add(0.5, *system.current_local_solution);
    system.get_vector("midpoint_vel").add(0.5, *system.old_local_solution);

    // Set variable numbers
    unsigned int ux_var = system.variable_number("velocityx");
    unsigned int uy_var, uz_var;

    if (dim > 1)
        uy_var = system.variable_number("velocityy");
    if (dim > 2)
        uz_var = system.variable_number("velocityz");

    // Velocity/Displacement DOF map
    const libMesh::DofMap & dof_map = system.get_dof_map();
    // Pressure DOF map
    const libMesh::DofMap & p_dof_map = p_system.get_dof_map();
    // Activation
    const libMesh::DofMap & dof_map_activation = dummy_system.get_dof_map();
    // Fibers
    const libMesh::DofMap & dof_map_fibers= fiber_system.get_dof_map();

    // Velocity/Displacement FE Type
    libMesh::FEType fe_disp(1, libMesh::LAGRANGE);
    // Pressure FE Type
    libMesh::FEType fe_pressure(1, libMesh::LAGRANGE);

    UniquePtr<libMesh::FEBase> fe_u(libMesh::FEBase::build(dim, fe_disp));
    UniquePtr<libMesh::FEBase> fe_p(libMesh::FEBase::build(dim, fe_pressure));

//    UniquePtr<libMesh::FEBase> fe_elem_face(libMesh::FEBase::build(dim, fe_disp));
//    UniquePtr<libMesh::FEBase> fe_neighbor_face(libMesh::FEBase::build(dim, fe_disp));


    auto order = fe_u->get_order();
    // Quadrature rule
    libMesh::QGauss qrule_1(dim, libMesh::SECOND);
    fe_u->attach_quadrature_rule(&qrule_1);
    fe_p->attach_quadrature_rule(&qrule_1);

    const std::vector<libMesh::Real> & JxW_u = fe_u->get_JxW();
    const std::vector<libMesh::Point> & q_point_u = fe_u->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<libMesh::Real> > & phi_u = fe_u->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > & dphi_u = fe_u->get_dphi();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<libMesh::Real> & JxW_p = fe_p->get_JxW();
    const std::vector<std::vector<libMesh::Real> > & phi_p = fe_p->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > & dphi_p = fe_p->get_dphi();


    libMesh::DenseVector<libMesh::Number> Fe;
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseVector<libMesh::Number> Fep;
    libMesh::DenseMatrix<libMesh::Number> Mep;

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_ux;
    std::vector<libMesh::dof_id_type> dof_indices_uy;
    std::vector<libMesh::dof_id_type> dof_indices_uz;
    std::vector<libMesh::dof_id_type> dof_indices_p;
    std::vector<libMesh::dof_id_type> dof_indices_activation;

    // Grad U
    std::vector<double> disp_k;
    std::vector<double> p_k_midpoint;
    std::vector<double> p_n;
    std::vector<double> p_k;
    std::vector<double> vel_k_midpoint;
    std::vector<double> vel_k;
    std::vector<double> vel_n;
    libMesh::VectorValue<libMesh::Number> Uk;
    libMesh::VectorValue<libMesh::Number> Xk;
    libMesh::VectorValue<libMesh::Number> Vk;
    libMesh::VectorValue<libMesh::Number> Vn;
    libMesh::TensorValue<libMesh::Number> dU;
    // Grad U
    libMesh::TensorValue<libMesh::Number> dUk;
    // Grad V
    libMesh::TensorValue<libMesh::Number> dVk;
    // Strain
    libMesh::TensorValue<libMesh::Number> E;
    // Strain
    libMesh::TensorValue<libMesh::Number> Ek;
    // Strain
    libMesh::TensorValue<libMesh::Number> Hk;
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

    // Pressure
    double pn, p_midpoint, pk, p, q, dp;
    // Pressure Gradient
    libMesh::VectorValue <libMesh::Number> grad_pk;


    // Residual
    libMesh::VectorValue <libMesh::Number> res_k;
    libMesh::VectorValue <libMesh::Number> d2U_H_gradq;
    libMesh::VectorValue <libMesh::Number>  body_force;

    // On the boundary
    libMesh::UniquePtr<libMesh::FEBase> fe_face(libMesh::FEBase::build(dim, fe_disp));
    libMesh::QGauss qface(dim - 1, libMesh::SECOND);
    fe_face->attach_quadrature_rule(&qface);

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // density
    double rho;
    // tau / h
    double tau;
    //
    bool stabilize = M_stabilize;

    // anisotropy
    libMesh::RealGradient f0;
    libMesh::RealGradient s0;
    libMesh::RealGradient n0;
    std::vector<libMesh::dof_id_type> dof_indices_fibers;

    std::vector<double> x_k;

    // Activation
    double gamma_f;
    double gamma_s;
    double gamma_n;
    std::vector<double> gamma_f_k;
    libMesh::TensorValue<libMesh::Number> FA;


    // Loop over elements
    for (; el != end_el; ++el)
    {
        // Get element
        const libMesh::Elem * elem = *el;
        // Element ID
        auto elID = elem->id();
        // Get density of the material in this element
        rho = M_materialMap[0]->M_density;
        // element size
        double h = elem->hmin();
        // stabilization parameter
        tau = M_materialMap[0]->M_tau;


        // Velocity/Displacement DOFS
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_ux, ux_var);
        p_dof_map.dof_indices(elem, dof_indices_p);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_ux_dofs = dof_indices_ux.size();
        const unsigned int n_dofs_p = dof_indices_p.size();

        fe_u->reinit(elem);
        fe_p->reinit(elem);
        Fe.resize(n_dofs);
        Me.resize(n_dofs, n_dofs);
        Fep.resize(n_dofs_p);
        Mep.resize(n_dofs_p, n_dofs_p);

        // get uk
        disp_k.resize(n_dofs);
        vel_k_midpoint.resize(n_dofs);
        vel_k.resize(n_dofs);
        vel_n.resize(n_dofs);
        p_k.resize(n_dofs_p);
        p_n.resize(n_dofs_p);
        p_k_midpoint.resize(n_dofs_p);

        disp_system.get_vector("midpoint_disp").get(dof_indices, disp_k);
        system.get_vector("midpoint_vel").get(dof_indices, vel_k_midpoint);
        system.current_local_solution->get(dof_indices, vel_k);
        system.old_local_solution->get(dof_indices, vel_n);
        p_system.get_vector("midpoint_pr").get(dof_indices_p, p_k_midpoint);
        p_system.current_local_solution->get(dof_indices_p, p_k);
        p_system.old_local_solution->get(dof_indices_p, p_n);

        // If the material is active, get activation
        if (activation_ptr)
        {
            dof_map_activation.dof_indices (elem, dof_indices_activation);
            const unsigned int n_dofs_activation   = dof_indices_activation.size();
            gamma_f_k.resize(n_dofs_activation);
            activation_ptr->get(dof_indices_activation, gamma_f_k);
        }


        // Get Fibers
        dof_map_fibers.dof_indices(elem, dof_indices_fibers);
        // fiber direction
        f0(0) = (*fiber_system.solution)(dof_indices_fibers[0]);
        f0(1) = (*fiber_system.solution)(dof_indices_fibers[1]);
        f0(2) = (*fiber_system.solution)(dof_indices_fibers[2]);
        // sheet direction
        s0(0) = (*sheets_system.solution)(dof_indices_fibers[0]);
        s0(1) = (*sheets_system.solution)(dof_indices_fibers[1]);
        s0(2) = (*sheets_system.solution)(dof_indices_fibers[2]);
        // crossfiber direction
        n0(0) = (*xfiber_system.solution)(dof_indices_fibers[0]);
        n0(1) = (*xfiber_system.solution)(dof_indices_fibers[1]);
        n0(2) = (*xfiber_system.solution)(dof_indices_fibers[2]);

        // Block K
        int index = 0;
        for (unsigned int qp = 0; qp < qrule_1.n_points(); qp++)
        {

            dUk *= 0.0;
            Ek *= 0.0;
            gamma_f *= 0.0;
            FA *= 0.0;
            Uk *= 0.0;
            Vk *= 0.0;
            Vn *= 0.0;
            pk = 0.0;
            grad_pk *= 0.0;

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
                        dVk(idim, jdim) += dphi_u[l][qp](jdim) * vel_k_midpoint[l + idim * n_phi];
                    }
                }
            } // end grad Uk on qp

            const unsigned int n_phi_p = phi_p.size();
            for(unsigned int l = 0; l < n_phi_p; ++l )
            {
                pk += phi_p[l][qp] * p_k[l];
                pn += phi_p[l][qp] * p_n[l];
                p_midpoint += phi_p[l][qp] * p_k_midpoint[l];
                pk += phi_p[l][qp] * p_k[l];
                grad_pk += dphi_p[l][qp] * p_k_midpoint[l];
            }// end grad pk on qp


            if (activation_ptr)
            {
                FA = Material::M_identity;

                gamma_f *= 0.0;
                for(unsigned int l = 0; l < n_phi; ++l )
                {
                    gamma_f += phi_p[l][qp] * gamma_f_k[l];
                }
                // We need to separate the dim = 2 and dim = 3 cases
                // in order to maintain the plane strain condition
                if(dim == 3)
                {
                    gamma_s = 1.0 / std::sqrt(1.0+gamma_f) - 1.0;
                    gamma_n = gamma_s;
                }
                else
                {
                    gamma_s = 1.0 / (1.0+gamma_f) - 1.0;
                    gamma_n = 0.0;
                }
                for(int idim = 0; idim < dim; idim++)
                {
                    for(int jdim = 0; jdim < dim; jdim++)
                    {
                        FA(idim, jdim) += gamma_f * f0(idim) * f0(jdim)
                                        + gamma_s * s0(idim) * s0(jdim)
                                        + gamma_n * n0(idim) * n0(jdim);
                    }
                }
                M_materialMap[0]->M_FA = FA;
                M_materialMap[0]->M_FAinv = FA.inverse();
                M_materialMap[0]->M_CAinv = M_materialMap[0]->M_FAinv  * M_materialMap[0]->M_FAinv;
            }

            M_materialMap[0]->M_gradU = dUk;
            M_materialMap[0]->M_pressure = p_midpoint;
            M_materialMap[0]->M_f0 = f0;
            M_materialMap[0]->updateVariables();
            M_materialMap[0]->evaluateStress(ElasticSolverType::Mixed);
            Sk = M_materialMap[0]->M_PK1;

            // Residual
            const double x = q_point_u[qp](0);
            const double y = q_point_u[qp](1);
            const double z = q_point_u[qp](2);
            const double time = 0.0;
            for (int idim = 0; idim < dim; idim++)
                body_force(idim) = M_rhsFunction(time, x, y, z, idim);


            // the strog residual is
            // - div \sigma - \rho f
            // for linear elements and nonlinear material div \sigma = H grad p
            // for linear elements and linear material div \sigma = grad p
            // for linear materials we set H = I
            Hk = M_materialMap[0]->H();
            res_k = rho / dt * (Vk - Vn) - Hk *  grad_pk - rho * body_force;


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
                    Fe(n + jdim * n_ux_dofs) += JxW_u[qp] * rho * body_force(jdim) * phi_u[n][qp];
                    Fe(n + jdim * n_ux_dofs) +=
                                    JxW_u[qp] * rho * (Vn(jdim) - Vk(jdim)) / dt * phi_u[n][qp];
                }
            }

            p =  M_materialMap[0]->evaluatePressureResidual();
            for(unsigned int n = 0; n < phi_p.size(); ++n )
            {
                // In this way we can do the compressible and the incopressible case together
                // p-dot = k div v
                // p^n+1 - p^n = k div v^n
                // p^k + delta_p - p^n = dt * k div v^n
                // delta_p/dt  = (p^n-p^k)/dt + k div v^n

                // if assemble_pressure_mass == 0 then we are solving the incompressible case
                Fep(n) += JxW_p[qp] * (pn - pk) / dt * phi_p[n][qp];
                double K = M_materialMap[0]->d2U(M_materialMap[0]->M_Jk);
                Fep(n) += JxW_p[qp] * K * Hk.contract(dVk) * phi_p[n][qp];
                if(stabilize)
                {
                    // the stabilized version is
                    // q * pdot - q * k div v - q * k div v' = 0
                    // q * (pk-pn)/dt - q * k div v + grad q * v'
                    // v'  = - tau Res_k
                    // q * (pk-pn)/dt - q * k div v - grad q * tau * Res_k
                    // q * (pk-pn)/dt - q * k div v - tau * ( rho (vk-vn) / dt - grad p - rho f)  * grad q

                    // we use += since we are using the same convetion as in the line above for
                    // q * pk
                    // These 2 terms should have the opposite sign
                    // q * pdot - q * k div v - tau Res_k * grad q
                    // Fe(n+dim*n_ux_dofs) -= tau * JxW_p[qp] * M_materialMap[0]->d2U() * grad_pk * dphi_p[n][qp];



                    // the stabilized version is
                    // q * pdot - q * U'' * ( H : grad v ) - q * U'' * ( H : grad v' )
                    // q * pdot - q * U'' * ( H : grad v ) - q * U'' * ( H : grad v' )
                    // q * pdot - q * U'' * ( H : grad v )  + U'' * ( H * grad q ) * v'
                    // v'  = - tau * Res_k
                    // q * pdot - q * U'' * ( H : grad v ) + U'' * tau * Res_k * ( H * grad q )

                    // we use += since we are using the same convetion as in the line above for
                    // q * pk
                    // These 2 terms should have the opposite sign
                    // q * pk - q * U' - U'' * ( tau Res_k ) * ( H * grad q )
                    //Fe(n+dim*n_ux_dofs) += JxW_p[qp] * tau * M_materialMap[0]->d2U() * res_k * ( Hk * dphi_p[n][qp] );
                    Fep(n) += JxW_p[qp] * tau * K * h * res_k * ( Hk * dphi_p[n][qp] );
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
                        //lumped mass matrix
                        Me(index1, index2) += rho / dt * JxW_u[qp] * phi_u[n][qp] * phi_u[m][qp];
                    }
                }
            }
            for (unsigned int n = 0; n < phi_p.size(); ++n)
            {
                for (unsigned int m = 0; m < phi_p.size(); ++m)
                {
                    //lumped mass matrix
                    Mep(n, m) += JxW_p[qp]  / dt * phi_p[n][qp] * phi_p[m][qp];
                }
            }

        }

        apply_BC(elem, Me, Fe, fe_face, qface, mesh, n_ux_dofs);

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
                            Vk *= 0.0;
                            Uk *= 0.0;
                            dUk *= 0.0;

                            // Evaluate displacements velocities and gradients at the quadrature node
                            for (unsigned int l = 0; l < phi_face.size(); l++)
                            {
                                for (int idim = 0; idim < dim; idim++)
                                {
                                    Vk(idim) += phi_face[l][qp] * vel_k[l + idim * n_phi];
                                    Uk(idim) += phi_face[l][qp] * disp_k[l + idim * n_phi];
                                    for (int jdim = 0; jdim < dim; jdim++)
                                    {
                                        dUk(idim, jdim) +=
                                                        dphi_face[l][qp](jdim) * disp_k[l + idim * n_phi];
                                    }
                                }
                            }
                            // Position of the quadrature node
                            const double xq = qface_point[qp](0);
                            const double yq = qface_point[qp](1);
                            const double zq = qface_point[qp](2);

                            //Fill with the value of bc
                            g *= 0.0;
                            for (int pdim = 0; pdim < dim; ++pdim)
                                g(pdim) = bc->get_function()(0.0, xq, yq, zq, pdim);

                            for (int idim = 0; idim < dim; ++idim)
                            {
                                for (unsigned int i = 0; i < phi_face.size(); i++)
                                {
                                    test *= 0.0;
                                    test(idim) = JxW_face[qp] * phi_face[i][qp];
                                    // \int g
                                    auto val = mu * gamma_mu / hE * (g - Uk) * test;
                                    auto val2 = Uk(idim) ;
                                    Fe(i + idim * n_ux_dofs) +=
                                                    mu * gamma_mu / hE * (g - Uk) * test;

                                    double gdotn = g.contract(normals[qp]);
                                    double ukdotn = Uk.contract(normals[qp]);
                                    val = kappa * gamma_kappa / hE * (gdotn - ukdotn) * (test * normals[qp]);
                                    Fe(i + idim * n_ux_dofs) +=
                                                    kappa * gamma_kappa / hE * (gdotn - ukdotn) * (test * normals[qp]);
                                }// end loop over test
                            } // end loop over dimension test function
                        }// end loop over qp
                    }
                }
            }
        } // end for loop on edges

        dof_map.constrain_element_matrix_and_vector (Me, Fe, dof_indices);
        system.matrix->add_matrix(Me, dof_indices);
        p_system.matrix->add_matrix(Mep, dof_indices_p);
        system.rhs->add_vector(Fe, dof_indices);
        p_system.rhs->add_vector(Fep, dof_indices_p);


    } // end loop over elements
    system.matrix->close();
    //system.get_matrix("mass").close();
    system.rhs->close();
    p_system.matrix->close();
    p_system.rhs->close();
    std::cout << "Assembly completed" << std::endl;

}



void
MixedDynamicElasticity::solve_system()
{
        LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);
        LinearSystem& p_system  =  M_equationSystems.get_system<LinearSystem>("pressure");

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0,0.0);
    rval = M_linearSolver->solve (*system.matrix, nullptr,
                                                        system.get_vector("step"),
                                                        *system.rhs, tol, max_iter);

    std::cout << "Solving for pressure" << std::endl;
    rval = M_projectionsLinearSolver->solve (*p_system.matrix, nullptr,
                                                        p_system.get_vector("step"),
                                                        *p_system.rhs, tol, max_iter);
    std::cout << "Adding vectors" << std::endl;

    *p_system.solution += p_system.get_vector("step");
}


} /* namespace BeatIt */
