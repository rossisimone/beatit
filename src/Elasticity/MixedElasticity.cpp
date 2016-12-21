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
#include "Elasticity/Materials/IsotropicMaterial.hpp"

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
MixedElasticity::assemble_residual()
{
	typedef libMesh::LinearImplicitSystem    LinearSystem;

	std::cout << "* MIXED ELASTICITY: assembling ... " << std::endl;

    using libMesh::UniquePtr;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int max_dim = 3;
    // Get a reference to the LinearImplicitSystem we are solving
    LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);
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
    libMesh::FEType fe_disp= dof_map.variable_type(ux_var);
    libMesh::FEType fe_pr = dof_map.variable_type(p_var);

    UniquePtr<libMesh::FEBase> fe_u(libMesh::FEBase::build(dim, fe_disp) );
    UniquePtr<libMesh::FEBase> fe_p(libMesh::FEBase::build(dim, fe_pr) );

    //libMesh::QGauss qrule_1(dim, libMesh::FOURTH  );
    //libMesh::QGauss qrule_2(dim,  libMesh::FOURTH   );
    libMesh::QGauss qrule_1(dim, libMesh::SECOND  );
    libMesh::QGauss qrule_2(dim,  libMesh::SECOND   );

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
    // Residual
	libMesh::VectorValue <libMesh::Number> res_k;
	libMesh::VectorValue <libMesh::Number> d2U_H_gradq;
	libMesh::VectorValue <libMesh::Number>  body_force;
    const std::vector<libMesh::Point> & q_point = fe_u->get_xyz();
	// On the boundary
    libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_disp));
    libMesh::QGauss qface(dim-1,  libMesh::FOURTH);
    fe_face->attach_quadrature_rule (&qface);

	libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    int assemble_pressure_mass = 1;
    bool stabilize = M_stabilize;
    // density
    double rho;
    // tau / h
    double tau;

    for ( ; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        auto elID = elem->id();
        double h = elem->hmin();
        rho = M_materialMap[0]->M_density;
        tau = h * h * M_materialMap[0]->M_tau;
        tau = M_materialMap[0]->M_tau;
        dof_map.dof_indices (elem, dof_indices);
//        for(auto && v : dof_indices) std::cout << v << std::endl;
        dof_map.dof_indices (elem, dof_indices_ux, ux_var);
        if(dim>1)dof_map.dof_indices (elem, dof_indices_uy, uy_var);
        if(dim>2)dof_map.dof_indices (elem, dof_indices_uz, uz_var);

        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_ux_dofs = dof_indices_ux.size();

        fe_u->reinit (elem);
        auto& dxdxi   = fe_u_map.get_dxyzdxi();
        auto& dxdeta  = fe_u_map.get_dxyzdeta();
        auto& dxdzeta  = fe_u_map.get_dxyzdzeta();

        std::cout << "max element length  sqaured is: " << h*h << ", while the jacobian is " <<  dxdzeta.size() << std::endl;
        h_jac.slice(0) = dxdxi[0];
//        std::cout << "Slice 1 Done ... " << std::flush;
        h_jac.slice(1) = dxdeta[0];
//        std::cout << "Slice 2 Done ... " << std::flush;
        //h_jac.slice(2) = dxdzeta[0];
        h_jac *= 0.0;
        h_jac(0,0) = 0.5;
        h_jac(0,1) = 0.5;
        h_jac(1,1) = 0.5;
        h_jac(2,2) = 0.;
        auto h2 = h_jac*h_jac;

//        std::cout << "Slice 3 Done ... Printing" << std::endl;
        h2.print(std::cout);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
        fe_p->reinit (elem);
        Fe.resize (n_dofs);
        Ke.resize(n_dofs, n_dofs);

        // get uk
        solution_k.resize(n_dofs);
        //  	  for(auto && di : dof_indices)  std::cout << "dof id: " << di << std::endl;

        system.current_local_solution->get(dof_indices, solution_k);
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
        assemble_pressure_mass = M_materialMap[0]->isIncompressible() ? 0 : 1;

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

            M_materialMap[0]->M_gradU = dUk;

            const unsigned int n_phi_p = phi_p.size();
            for(unsigned int l = 0; l < n_phi_p; ++l )
            {
                pk += phi_p[l][qp] * solution_k[l+dim*n_phi];
                grad_pk += dphi_p[l][qp] * solution_k[l+dim*n_phi];
            }// end grad pk on qp

            M_materialMap[0]->M_pressure = pk;
            M_materialMap[0]->evaluateStress( ElasticSolverType::Mixed);
            Sk = M_materialMap[0]->M_PK1;

            // Residual
            const double x = q_point[qp](0);
            const double y = q_point[qp](1);
            const double z = q_point[qp](2);
            const double time = 0.0;
            for(int idim = 0; idim < dim; idim++)
            {
                body_force(idim) =M_rhsFunction(time,x,y,z,idim);
            }

            // the strog residual is
            // - div \sigma - \rho f
            // for linear elements and nonlinear material div \sigma = H grad p
            // for linear elements and linear material div \sigma = grad p
            // for linear materials we set H = I
            Hk = M_materialMap[0]->H();
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
                    Fe(n+jdim*n_ux_dofs) += rho * body_force(jdim) * phi_u[n][qp];

                }
            }

            p =  M_materialMap[0]->evaluatePressureResidual();
            for(unsigned int n = 0; n < phi_p.size(); ++n )
            {
                // In this way we can do the compressible and the incopressible case together
                // pk + p = pk - k div u
                // if assemble_pressure_mass == 0 then we are solving the incompressible case
                Fe(n+dim*n_ux_dofs) -=   JxW_p[qp] * ( assemble_pressure_mass * pk + p ) * phi_p[n][qp];
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
                    Fe(n+dim*n_ux_dofs) += JxW_p[qp] * tau * M_materialMap[0]->d2U() * h2 * res_k * ( Hk * dphi_p[n][qp] );
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

                            M_materialMap[0]->evaluateDeviatoricJacobian(dU, 0.0);
                            S = M_materialMap[0]->M_deviatoric_jacobian;
                            Ke(n+jdim*n_ux_dofs,m+idim*n_ux_dofs) += S.contract(dW);
                        }
                    }
                    // for each pressure  trial function
                    for(unsigned int m = 0; m < phi_p.size(); ++m )
                    {
                        q = phi_p[m][qp];
                        M_materialMap[0]->evaluateVolumetricJacobian(dU, q);
                        S = M_materialMap[0]->M_volumetric_jacobian;
                        Ke(n+jdim*n_ux_dofs,m+dim*n_ux_dofs) += S.contract(dW);
                    }
                }
            }

            // For each pressure test function
            for(unsigned int n = 0; n < phi_p.size(); ++n )
            {
                q =  JxW_p[qp]*phi_p[n][qp];
                // For each pressure trial function
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

                        dp = M_materialMap[0]->dpdF(dU);
                        Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) += dp * q;

                        if(stabilize)
                        {
                            M_materialMap[0]->dH(dU, dH);
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
                            Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) += tau * M_materialMap[0]->d3U() * Hk.contract(dU) * ( Hk *  grad_pk ) * ( Hk * dphi_p[n][qp] );
                            //  2) U'' * tau * ( dH  grad pk ) * ( H * grad q )
                            Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) += tau * M_materialMap[0]->d2U() * ( dH *  grad_pk ) * ( Hk * dphi_p[n][qp] );
                            //  3) U'' * tau * ( H  grad pk ) * ( dH * grad q )
                            Ke(n+dim*n_ux_dofs,m+idim*n_ux_dofs) += tau * M_materialMap[0]->d2U() * ( Hk *  grad_pk ) * ( dH * dphi_p[n][qp] );
                            // The last term is assembled afterwards

                        }
                    }
                }
                // For each pressure trial function
                for(unsigned int m = 0; m < phi_p.size(); ++m )
                {
                    Ke(n+dim*n_ux_dofs,m+dim*n_ux_dofs) += assemble_pressure_mass * q * phi_p[m][qp];
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
                        Ke(n+dim*n_ux_dofs,m+dim*n_ux_dofs) += JxW_p[qp] * tau * M_materialMap[0]->d2U()
                                                                         * (h2 * Hk * dphi_p[n][qp])
                                                                         * (Hk * dphi_p[m][qp]);
                    }
                }
            }
        }

        Elasticity::apply_BC(elem, Ke, Fe, fe_face, qface, mesh, n_ux_dofs);
//        Ke.print();
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
//        Ke.print();

        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);

    }

    system.matrix->close();
//    system.matrix->print(std::cout);
    system.rhs->close();
//    system.rhs->print(std::cout);
}


} /* namespace BeatIt */
