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
 * MixedElasticity.cpp
 *
 *  Created on: Oct 24, 2016
 *      Author: srossi
 */

#include "Elasticity/MixedElasticity.hpp"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/type_tensor.h"

// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
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

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

#include "libmesh/exodusII_io.h"

#include "libmesh/perf_log.h"

#include <sys/stat.h>
#include "BoundaryConditions/BCData.hpp"
#include "Util/IO/io.hpp"
#include "Util/MapsToLibMeshTypes.hpp"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"


//Materials
#include "Elasticity/Materials/LinearMaterial.hpp"
#include "Elasticity/Materials/Neohookean.hpp"
#include "Elasticity/Materials/IsotropicMaterial.hpp"
#include "Elasticity/Materials/HolzapfelOgden.hpp"

#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"

namespace BeatIt {


MixedElasticity::MixedElasticity( libMesh::EquationSystems& es, std::string system_name ) :
		Elasticity(es, system_name)
{
	// TODO Auto-generated constructor stub

}

MixedElasticity::~MixedElasticity() {
	// TODO Auto-generated destructor stub
}



void
MixedElasticity::assemble_residual(double dt , libMesh::NumericVector<libMesh::Number>* activation_ptr)
{
	typedef libMesh::LinearImplicitSystem    LinearSystem;

	std::cout << "* MIXED ELASTICITY: assembling ... " << std::flush;

    using std::unique_ptr;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int max_dim = 3;
    // Get a reference to the LinearImplicitSystem we are solving
    LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);
    system.update();
    double time = system.time;

    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");
    ParameterSystem& dummy_system       = M_equationSystems.get_system<ParameterSystem>("dumb");

    system.get_vector("residual").zero();
    system.rhs->zero();
    system.matrix->zero();
    system.update();
    unsigned int ux_var = system.variable_number ("displacementx");
    unsigned int uy_var, uz_var, p_var;

   if(dim>1)uy_var =  system.variable_number ("displacementy");
   if(dim>2)uz_var =  system.variable_number ("displacementz");
	p_var =  system.variable_number ("pressure");

    const libMesh::DofMap & dof_map = system.get_dof_map();
    //dof_map.create_dof_constraints(mesh, time);
    const libMesh::DofMap & dof_map_fibers= fiber_system.get_dof_map();
    const libMesh::DofMap & dof_map_activation = dummy_system.get_dof_map();

    libMesh::FEType fe_disp= dof_map.variable_type(ux_var);
    libMesh::FEType fe_pr = dof_map.variable_type(p_var);

    std::unique_ptr<libMesh::FEBase> fe_u(libMesh::FEBase::build(dim, fe_disp) );
    std::unique_ptr<libMesh::FEBase> fe_p(libMesh::FEBase::build(dim, fe_pr) );

    auto order_u = fe_u->get_order();
    auto order_p = fe_u->get_order();

//    libMesh::QGauss qrule_1(dim,  libMesh::SECOND );
//    libMesh::QGauss qrule_2(dim,  libMesh::SECOND );
    libMesh::QGauss qrule_1(dim,  libMesh::FIRST );
    libMesh::QGauss qrule_2(dim,  libMesh::FIRST );

    fe_u->attach_quadrature_rule(&qrule_1);
    fe_p->attach_quadrature_rule(&qrule_2);

    const libMesh::FEMap& fe_u_map = fe_u->get_fe_map ();

    const std::vector<libMesh::Real> & JxW_u = fe_u->get_JxW();
    const std::vector<libMesh::Real> & JxW_p = fe_p->get_JxW();

    const std::vector<libMesh::Point> & q_point_u = fe_u->get_xyz();
    const std::vector<libMesh::Point> & q_point_p = fe_p->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<libMesh::Real> > & phi_u = fe_u->get_phi();
    const std::vector<std::vector<libMesh::Real> > & phi_p = fe_p->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & dphi_u = fe_u->get_dphi();
    const std::vector<std::vector<libMesh::RealGradient> > & dphi_p = fe_p->get_dphi();

    libMesh::DenseVector<libMesh::Number> Fe;
    libMesh::DenseMatrix<libMesh::Number>  Ke;

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_ux;
	std::vector<libMesh::dof_id_type> dof_indices_uy;
	std::vector<libMesh::dof_id_type> dof_indices_uz;
	std::vector<libMesh::dof_id_type> dof_indices_p;
    std::vector<libMesh::dof_id_type> dof_indices_activation;

	// Grad U
	std::vector<double> solution_k;
    libMesh::TensorValue <libMesh::Number> dU;
	// Grad U
    libMesh::TensorValue <libMesh::Number> dUk;
    // Strain
    libMesh::TensorValue <libMesh::Number> E;
    // Strain
    libMesh::TensorValue <libMesh::Number> Ek;
    // Identity
    libMesh::TensorValue <libMesh::Number> id;
    id(0,0) = 1.0; id(1,1) = 1.0; id(2,2) = 1.0;
    // Stress
    libMesh::TensorValue <libMesh::Number> Sk;
    // Stress
    libMesh::TensorValue <libMesh::Number> S;
    // Grad W (test function)
    libMesh::TensorValue <libMesh::Number> dW;
    // Cofactor
    libMesh::TensorValue <libMesh::Number> Hk;
    // Linearization Cofactor
    libMesh::TensorValue <libMesh::Number> dH;
    // Jacobian
    libMesh::TensorValue <libMesh::Number> h_jac;
    // Pressure
    double pk, p, q, dp;
    // Pressure Gradient
    libMesh::VectorValue <libMesh::Number> grad_pk;
    libMesh::VectorValue <libMesh::Number> grad_q;
    libMesh::VectorValue <libMesh::Number> Uk;
    // Dirichlet function
    libMesh::VectorValue<libMesh::Number> g;
    libMesh::VectorValue<libMesh::Number> trial;
    libMesh::VectorValue<libMesh::Number> test;
    libMesh::TensorValue<libMesh::Number> dtrial;

    // Residual
	libMesh::VectorValue <libMesh::Number> res_k;
	libMesh::VectorValue <libMesh::Number> d2U_H_gradq;
	libMesh::VectorValue <libMesh::Number>  body_force;
    const std::vector<libMesh::Point> & q_point = fe_u->get_xyz();
	// On the boundary
    std::unique_ptr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_disp));
    libMesh::QGauss qface(dim-1,  libMesh::SECOND);
    fe_face->attach_quadrature_rule (&qface);

	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    int assemble_pressure_mass = 1;
    bool stabilize = M_stabilize;
    // density
    double rho;
    // tau / h
    double tau;
    // anisotropy
    libMesh::RealGradient f0;
    libMesh::RealGradient s0;
    libMesh::RealGradient n0;
    std::vector<libMesh::dof_id_type> dof_indices_fibers;

    // Active Material
    double gamma_f;
    double gamma_s;
    double gamma_n;
    std::vector<double> gamma_f_k;
    libMesh::TensorValue <libMesh::Number> FA;


    for ( ; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        auto elID = elem->id();
        auto blockID = elem->subdomain_id();
        double h = elem->hmin();
        rho = M_materialMap[blockID]->M_density;
        tau = M_materialMap[blockID]->M_tau;
        auto h2 = h*h;

        dof_map.dof_indices (elem, dof_indices);


//        for(auto && v : dof_indices) std::cout << v << std::endl;
        dof_map.dof_indices (elem, dof_indices_ux, ux_var);
        if(dim>1)dof_map.dof_indices (elem, dof_indices_uy, uy_var);
        if(dim>2)dof_map.dof_indices (elem, dof_indices_uz, uz_var);

        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_ux_dofs = dof_indices_ux.size();

        fe_u->reinit (elem);

//
////        std::cout << "Slice 3 Done ... Printing" << std::endl;
//        h2.print(std::cout);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
        fe_p->reinit (elem);
        Fe.resize (n_dofs);
        Ke.resize(n_dofs, n_dofs);

        // get uk
        solution_k.resize(n_dofs);

        system.current_local_solution->get(dof_indices, solution_k);

        if(activation_ptr)
        {
            dof_map_activation.dof_indices (elem, dof_indices_activation);
            const unsigned int n_dofs_activation   = dof_indices_activation.size();
            gamma_f_k.resize(n_dofs_activation);
            //for(auto && nd : dof_indices_activation) std::cout << nd << std::endl;
            //std::cout << "getting " << n_dofs_activation << " dofs for  gf ... " << std::flush;
            activation_ptr->get(dof_indices_activation, gamma_f_k);
            //std::cout << " done" << std::endl;
        }

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
        assemble_pressure_mass = M_materialMap[blockID]->isIncompressible() ? 0 : 1;

        for (unsigned int qp=0; qp<qrule_1.n_points(); qp++)
        {

            dUk *= 0.0;
            Ek *= 0.0;
            Sk *= 0.0;
            pk = 0.0;
            grad_pk *= 0.0;
            const unsigned int n_phi = phi_u.size();
            for(unsigned int l = 0; l < n_phi; ++l )
            {
                for(int idim = 0; idim < dim; idim++)
                {
                    for(int jdim = 0; jdim < dim; jdim++)
                    {
                        dUk(idim, jdim) += dphi_u[l][qp](jdim) *solution_k[l+idim*n_phi];
                    }
                }
            } // end grad Uk on qp

            M_materialMap[blockID]->M_gradU = dUk;
            M_materialMap[blockID]->M_f0 = f0;
            M_materialMap[blockID]->M_s0 = s0;

            const unsigned int n_phi_p = phi_p.size();
            for(unsigned int l = 0; l < n_phi_p; ++l )
            {
                pk += phi_p[l][qp] * solution_k[l+dim*n_phi];
                grad_pk += dphi_p[l][qp] * solution_k[l+dim*n_phi];
            }// end grad pk on qp

            M_materialMap[blockID]->M_pressure = pk;

            // Residual
            const double x = q_point[qp](0);
            const double y = q_point[qp](1);
            const double z = q_point[qp](2);
            const double time = 0.0;
            for(int idim = 0; idim < dim; idim++)
            {
                body_force(idim) =M_rhsFunction(time,x,y,z,idim);
            }

            if(activation_ptr)
            {
                FA = Material::M_identity;

                gamma_f *= 0.0;
                for(unsigned int l = 0; l < n_phi; ++l )
                {
                    gamma_f += phi_p[l][qp] * gamma_f_k[l];
                }

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
                M_materialMap[blockID]->M_FA = FA;
                M_materialMap[blockID]->M_FAinv = FA.inverse();
                M_materialMap[blockID]->M_CAinv = M_materialMap[blockID]->M_FAinv  * M_materialMap[blockID]->M_FAinv;
            }

            M_materialMap[blockID]->updateVariables();
            M_materialMap[blockID]->evaluateStress( ElasticSolverType::Mixed);
            Sk = M_materialMap[blockID]->M_PK1;



            // the strog residual is
            // - div \sigma - \rho f
            // for linear elements and nonlinear material div \sigma = H grad p
            // for linear elements and linear material div \sigma = grad p
            // for linear materials we set H = I
            Hk = M_materialMap[blockID]->H();
            res_k = - Hk *  grad_pk - rho * body_force;

            for(unsigned int n = 0; n < phi_u.size(); ++n )
            {
                for(int jdim= 0; jdim < dim; jdim++)
                {
                    dW *= 0.0;
                    dW(jdim, 0) = JxW_u[qp]* dphi_u[n][qp](0);
                    dW(jdim, 1) = JxW_u[qp]* dphi_u[n][qp](1);
                    dW(jdim, 2) = JxW_u[qp]* dphi_u[n][qp](2);
                    // Compute  - \nabla \cdot \sigma + f
                    Fe(n+jdim*n_ux_dofs) -= Sk.contract(dW);
                    Fe(n+jdim*n_ux_dofs) += JxW_u[qp] * rho * body_force(jdim) * phi_u[n][qp];

                }
            }

            p =  M_materialMap[blockID]->evaluatePressureResidual();
            for(unsigned int n = 0; n < phi_p.size(); ++n )
            {
                // In this way we can do the compressible and the incopressible case together
                // pk + p = pk - k div u
                // if assemble_pressure_mass == 0 then we are solving the incompressible case
                // Fe(n+dim*n_ux_dofs) -=   JxW_p[qp] * ( assemble_pressure_mass * pk + p ) * phi_p[n][qp];
                // Changing sign to make it symmetric!!!
                Fe(n+dim*n_ux_dofs) +=   JxW_p[qp] * ( assemble_pressure_mass * pk + p ) * phi_p[n][qp];
                if(stabilize)
                {
                    // the stabilized version is
                    // q * pk - q * k div u - q * k div u'
                    // q * pk - q * k div u + grad q * u'
                    // u'  = - tau Res_k
                    // q * pk - q * k div u - tau Res_k * grad q
                    // q * pk - q * k div u - tau ( - grad p - rho f) * grad q
                    // q * pk - q * k div u - tau ( - grad p - rho f) * grad q
                    // q * pk - q * k div u +  tau * grad q * grad p + tau * rho * f * grad q

                    // we use += since we are using the same convetion as in the line above for
                    // q * pk
                    // These 2 terms should have the opposite sign
                    // q * pk - q * k div u - tau Res_k * grad q
                    // Fe(n+dim*n_ux_dofs) -= tau * JxW_p[qp] * M_materialMap[0]->d2U() * grad_pk * dphi_p[n][qp];



                    // the stabilized version is
                    // q * pk - q * U' - q * U'' * ( H : grad u' )
                    // q * pk - q * U' + U'' * ( H * grad q ) * u'
                    // u'  = - tau * Res_k
                    // q * pk - q * U' - U'' * ( H * grad q ) * ( tau Res_k )

                    // we use += since we are using the same convetion as in the line above for
                    // q * pk
                    // These 2 terms should have the opposite sign
                    // q * pk - q * U' - U'' * ( tau Res_k ) * ( H * grad q )
                    //Fe(n+dim*n_ux_dofs) += JxW_p[qp] * tau * M_materialMap[0]->d2U() * res_k * ( Hk * dphi_p[n][qp] );

                    //Fe(n+dim*n_ux_dofs) += JxW_p[qp] * tau * M_materialMap[0]->d2U( M_materialMap[0]->M_Jk) * h2 * res_k * ( Hk * dphi_p[n][qp] );
                    Fe(n+dim*n_ux_dofs) -= JxW_p[qp] * tau * M_materialMap[blockID]->d2U( M_materialMap[blockID]->M_Jk) * h2 * res_k * ( Hk * dphi_p[n][qp] );
                }
            }

            // Matrix
            // For each test function
            for(unsigned int n = 0; n < phi_u.size(); ++n )
            {

                // for each dimension of the test function
                for(int jdim= 0; jdim < dim; jdim++)
                {
                    dW *= 0.0;
                    dW(jdim, 0) = JxW_u[qp]* dphi_u[n][qp](0);
                    dW(jdim, 1) = JxW_u[qp]* dphi_u[n][qp](1);
                    dW(jdim, 2) = JxW_u[qp]* dphi_u[n][qp](2);

                    // for each trial function
                    for(unsigned int m = 0; m < phi_u.size(); ++m )
                    {
                        // for each dimension of the trial function
                        for(int idim = 0; idim < dim; idim++)
                        {
                            dU*= 0.0;
                            E*= 0.0;
                            S*= 0.0;

                            dU(idim, 0) =  dphi_u[m][qp](0);
                            dU(idim, 1) =  dphi_u[m][qp](1);
                            dU(idim, 2) =  dphi_u[m][qp](2);

                            M_materialMap[blockID]->evaluateDeviatoricJacobian(dU, 0.0);
                            S = M_materialMap[blockID]->M_deviatoric_jacobian;
                            //dU.print(std::cout);
                            Ke(n+jdim*n_ux_dofs,m+idim*n_ux_dofs) += S.contract(dW);
                        }
                    }
                    // for each pressure  trial function
                    for(unsigned int m = 0; m < phi_p.size(); ++m )
                    {
                        q = phi_p[m][qp];
                        M_materialMap[blockID]->evaluateVolumetricJacobian(dU, q);
                        S = M_materialMap[blockID]->M_volumetric_jacobian;
                        Ke(n+jdim*n_ux_dofs,m+dim*n_ux_dofs) += S.contract(dW);
                    }
                }
            }

            // For each pressure test function
            for(unsigned int n = 0; n < phi_p.size(); ++n )
            {
                q =  JxW_p[qp]*phi_p[n][qp];
            	grad_q = JxW_p[qp]*dphi_p[n][qp];

                // For each displacement trial function
                for(unsigned int m = 0; m < phi_u.size(); ++m )
                {
                    // for each dimension of the trial function
                    for(int idim = 0; idim < dim; idim++)
                    {
                        dU*= 0.0;
                        E*= 0.0;
                        S*= 0.0;

                        dU(idim, 0) =  dphi_u[m][qp](0);
                        dU(idim, 1) =  dphi_u[m][qp](1);
                        dU(idim, 2) =  dphi_u[m][qp](2);

                        dp = M_materialMap[blockID]->dpdF(dU);
                        Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) -= dp * q;

                        if(stabilize)
                        {
                            M_materialMap[blockID]->dH(dU, dH);
                            // the stabilized version is
                            // q * pk - q * U' - q * U'' * ( H : grad u' )
                            // q * pk - q * U' + U'' * ( H * grad q ) * u'
                            // u'  = - tau * Res_k
                            // q * pk - q * U' - U'' * ( H * grad q ) * ( tau Res_k )

                            // we use += since we are using the same convetion as in the line above for
                            // q * pk
                            // These 2 terms should have the opposite sign
                            // q * pk - q * U' - U'' * ( tau Res_k ) * ( H * grad q )

                            // Res_k =  - H  grad_p - \rho f
                            // - U'' * ( tau Res_k ) * ( H * grad q ) = + U'' * tau * ( H  grad_p + \rho f ) * ( H * grad q )
                            // The source term part does not enter the matrix
                            //  U'' * tau * ( H  grad_p ) * ( H * grad q )
                            // The linearization is composed of 4 terms
                            //  1) U''' * (H: dF ) * tau * ( H  grad pk ) * ( H * grad q )

                            Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) -= tau * M_materialMap[blockID]->d3U( M_materialMap[blockID]->M_Jk ) * Hk.contract(dU) * ( Hk *  grad_pk ) * ( Hk * grad_q );
                            //  2) U'' * tau * ( dH  grad pk ) * ( H * grad q )
                            Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) -= tau * M_materialMap[blockID]->d2U( M_materialMap[blockID]->M_Jk ) * ( dH *  grad_pk ) * ( Hk * grad_q );
                            //  3) U'' * tau * ( H  grad pk ) * ( dH * grad q )
                            Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) -= tau * M_materialMap[blockID]->d2U( M_materialMap[blockID]->M_Jk ) * ( Hk *  grad_pk ) * ( dH * grad_q );
                            // The last term is assembled afterwards

                        }
                    }
                }
                // For each pressure trial function
                for(unsigned int m = 0; m < phi_p.size(); ++m )
                {
                    Ke(n+dim*n_ux_dofs,m+dim*n_ux_dofs) -= assemble_pressure_mass * q * phi_p[m][qp];
                    if(stabilize)
                    {
                        // the stabilized version is
                        // q * pk - q * k div u - q * k div u'
                        // q * pk - q * k div u + grad q * u'
                        // u'  = - tau Res_k
                        // q * pk - q * k div u - tau ( - grad p - rho f) * grad q
                        // q * pk - q * k div u - tau ( - grad p - rho f) * grad q
                        // q * pk - q * k div u +  tau * grad q * grad p + tau * rho * f * grad q

                        // Same sign as the mass term

                        //Ke(n+dim*n_ux_dofs,m+dim*n_ux_dofs) += JxW_p[qp] * tau * M_materialMap[0]->d2U()
                        //                                         * (dphi_p[n][qp])
                        //                                        * (dphi_p[m][qp]);
                        // the stabilized version is
                        // q * pk - q * U' - q * U'' * ( H : grad u' )
                        // q * pk - q * U' + U'' * ( H * grad q ) * u'
                        // u'  = - tau * Res_k
                        // q * pk - q * U' - U'' * ( H * grad q ) * ( tau Res_k )

                        // we use += since we are using the same convetion as in the line above for
                        // q * pk
                        // These 2 terms should have the opposite sign
                        // q * pk - q * U' - U'' * ( tau Res_k ) * ( H * grad q )

                        // Res_k =  - H  grad_p - \rho f
                        // - U'' * ( tau Res_k ) * ( H * grad q ) = + U'' * tau * ( H  grad_p + \rho f ) * ( H * grad q )
                        // The source term part does not enter the matrix
                        //  U'' * tau * ( H  grad_p ) * ( H * grad q )
                        // The linearization is composed of 4 terms
                        //  4) U'' * tau * ( H  grad dp ) * ( H * grad q )
//                        Ke(n+dim*n_ux_dofs,m+dim*n_ux_dofs) += JxW_p[qp] * tau * M_materialMap[0]->d2U()
//                                                                         * (Hk * dphi_p[n][qp])
//                                                                         * (Hk * dphi_p[m][qp]);
                        Ke(n+dim*n_ux_dofs,m+dim*n_ux_dofs) -= h2 * tau * M_materialMap[blockID]->d2U( M_materialMap[blockID]->M_Jk)
                                                                         * (Hk * dphi_p[m][qp])
                                                                         * (Hk * grad_q);
                    }
                }
            }
        }



        Elasticity::apply_BC(elem, Ke, Fe, fe_face, qface, mesh, n_ux_dofs, nullptr, 0.0, time, &solution_k);

        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);

    }

    std::cout << " closing. " << std::endl;
    system.matrix->close();

    system.rhs->close();
    std::cout << " done. " << std::endl;
//    {
//            std::cout << "* MIXED ELASTICITY: Assigning field split information ... " << std::flush;
//
//            IS is_u_local;
//            IS is_p_global;
//            IS is_p_local;
//            std::vector<libMesh::dof_id_type> p_indices;
//            int var_num = dim;
//            dof_map.local_variable_indices(p_indices, mesh, var_num);
//
//            ISCreateGeneral(PETSC_COMM_SELF, p_indices.size(), reinterpret_cast<int*>(&p_indices[0]),PETSC_COPY_VALUES,&is_p_local);
//            ISAllGather(is_p_local, &is_p_global);
//            int nmin =  dof_map.first_dof();
//            int nmax = dof_map.end_dof();
//            ISComplement(is_p_local, nmin, nmax, &is_u_local);
//            typedef libMesh::PetscMatrix<libMesh::Number> PetscMatrix;
//            M_linearSolver->init(dynamic_cast<PetscMatrix *>(system.matrix));
//
//           PCFieldSplitSetIS(M_linearSolver->pc(),"p",is_p_local);
//           PCFieldSplitSetIS(M_linearSolver->pc(),"u",is_u_local);
//           std::cout << "  done" << std::endl;
//
//    //       PetscBool isMatSymm;
//    //       MatIsSymmetric(dynamic_cast<PetscMatrix *>(bidomain_system.matrix)->mat(), 1e-8,&isMatSymm);
//    //       std::cout << "PETSC, is the system matrix symmetric? " << isMatSymm << std::endl;
//           int size;
//           ISGetSize(is_u_local, &size);
//           std::cout << "is_u_local size: " << size << std::endl;
//           ISGetSize(is_p_local, &size);
//           std::cout << "is_p_local size: " << size << std::endl;
//           ISGetSize(is_p_global, &size);
//           std::cout << "is_p_global size: " << size << std::endl;
//
//        }

}




} /* namespace BeatIt */
