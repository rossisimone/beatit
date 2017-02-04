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

#include "Util/GenerateFibers.hpp"

//Materials
#include "Elasticity/Materials/LinearMaterial.hpp"
#include "Elasticity/Materials/IsotropicMaterial.hpp"

namespace BeatIt {

typedef libMesh::TransientLinearImplicitSystem    LinearSystem;

DynamicElasticity::DynamicElasticity( libMesh::EquationSystems& es, std::string system_name )
	: super(es, system_name)
{
	// TODO Auto-generated constructor stub
}

DynamicElasticity::~DynamicElasticity() {
	// TODO Auto-generated destructor stub
}

void
DynamicElasticity::setupSystem(std::string section )
{
    int dimension = M_equationSystems.get_mesh().mesh_dimension();
	// ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // DISPLACEMENT PART
    std::cout << "* DYNAMIC ELASTICITY: Creating new System for the DYNAMIC ELASTICITY Solver" << std::endl;
    LinearSystem& system  =  M_equationSystems.add_system<LinearSystem>(M_myName);

    std::string disp_name = "displacement";
    std::string vel_name = "velocity";

    int ord = M_datafile(section+"/order",  1);
    std::cout << "* DYNAMIC ELASTICITY: Reading displacement field - order: " << ord  << std::flush;
    auto order_it =  BeatIt::libmesh_order_map.find(ord);
    std::cout << " ... " << std::flush;
    libMesh::Order  order  = ( order_it !=  BeatIt::libmesh_order_map.end() ) ? order_it->second : throw std::runtime_error("Order Explode!!!");
    std::cout << "order found!" << std::endl;

    std::string fam = M_datafile(section+"/fefamily",  "lagrange");
    std::cout << "* DYNAMIC ELASTICITY: Reading displacement field - family: " << fam  << std::flush;
   auto  fefamily_it  =BeatIt::libmesh_fefamily_map.find(fam);
    std::cout << " ... " << std::flush;
    libMesh::FEFamily  fefamily = ( fefamily_it !=  BeatIt::libmesh_fefamily_map.end() ) ? fefamily_it->second : throw std::runtime_error("FeFamily Explode!!!");
    std::cout << "fefamily found!" << std::endl;
    std::cout << "* DYNAMIC ELASTICITY: Setting up displacement field - order: " << ord << ", fefamily = " << fam  << std::endl;
    system.add_variable(vel_name+"x", order, fefamily);
    if(dimension > 1) system.add_variable(vel_name+"y", order, fefamily);
    if(dimension > 2) system.add_variable(vel_name+"z", order, fefamily);

	LinearSystem& disp_system  =  M_equationSystems.add_system<LinearSystem>(disp_name);
    disp_system.add_variable(disp_name+"x", order, fefamily);
    if(dimension > 1) disp_system.add_variable(disp_name+"y", order, fefamily);
    if(dimension > 2) disp_system.add_variable(disp_name+"z", order, fefamily);


    // PRESSURE SYSTEM
    std::string formulation = M_datafile(section+"/formulation",  "primal");
    std::cout << "* DYNAMIC ELASTICITY: Using a " << formulation << " formulation" << std::endl;
    std::string pressure_name = "pressure";

    if(formulation == "mixed")
    {
		throw std::runtime_error("MIXED FORMULATION NOT CODED IN DYNAMIC CASE");
        M_solverType = ElasticSolverType::Mixed;
		ord = M_datafile(section+"/p_order",  1);
		fam = M_datafile(section+"/p_fefamily",  "lagrange");
	    std::cout << "* DYNAMIC ELASTICITY: Setting up pressure  field - order: " << ord << ", fefamily = " << fam  << std::endl;
	    libMesh::Order  p_order  = BeatIt::libmesh_order_map.find(ord)->second;
	    libMesh::FEFamily  p_fefamily   = BeatIt::libmesh_fefamily_map.find(fam)->second;
	    system.add_variable(pressure_name, p_order, p_fefamily);
    }
    else
    {
        M_solverType = ElasticSolverType::Primal;
        LinearSystem& p_system  =  M_equationSystems.add_system<LinearSystem>("Pressure_Projection");
        p_system.add_variable(pressure_name, libMesh::FIRST, libMesh::LAGRANGE);
        p_system.init();
    }

    system.add_vector("residual");
    system.add_vector("step");
    std::cout << "* DYNAMIC ELASTICITY: init "  << std::endl;
    system.init();
    disp_system.init();

}


void
DynamicElasticity::update_displacements(double dt)
{
    std::cout << "* DYNAMIC ELASTICITY: update_displacements: " << std::endl;
    LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem& disp_system  =  M_equationSystems.get_system<LinearSystem>("displacement");

    *disp_system.solution = *system.solution;
    disp_system.solution->scale(dt);
    *disp_system.solution += *disp_system.old_local_solution;

}


void
DynamicElasticity::advance()
{
    LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem& disp_system  =  M_equationSystems.get_system<LinearSystem>("displacement");

    *system.old_local_solution    = *system.solution;
    *system.older_local_solution  = *system.solution;

    *disp_system.old_local_solution   = *disp_system.solution;
    *disp_system.older_local_solution = *disp_system.solution;
}

void
DynamicElasticity::assemble_residual(double dt, libMesh::NumericVector<libMesh::Number>* activation_ptr)
{
  std::cout << "* DYNAMIC ELASTICITY: assembling ... " << std::endl;

    using libMesh::UniquePtr;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int max_dim = 3;
    // Get a reference to the LinearImplicitSystem we are solving
    LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem& disp_system  =  M_equationSystems.get_system<LinearSystem>("displacement");
    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");

    system.get_vector("residual").zero();
    system.rhs->zero();
    system.matrix->zero();
    system.update();
    disp_system.update();
    unsigned int ux_var = disp_system.variable_number ("displacementx");
    unsigned int uy_var, uz_var;

   if(dim>1)uy_var =  disp_system.variable_number ("displacementy");
   if(dim>2)uz_var =  disp_system.variable_number ("displacementz");

   const libMesh::DofMap & dof_map = system.get_dof_map();
    libMesh::FEType fe_disp= dof_map.variable_type(ux_var);

    UniquePtr<libMesh::FEBase> fe_u(libMesh::FEBase::build(dim, fe_disp) );
    auto order = fe_u->get_order();

    libMesh::QGauss qrule_1(dim, order);

    fe_u->attach_quadrature_rule(&qrule_1);

    const std::vector<libMesh::Real> & JxW_u = fe_u->get_JxW();
    const std::vector<libMesh::Point> & q_point_u = fe_u->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<libMesh::Real> > & phi_u = fe_u->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> > & dphi_u = fe_u->get_dphi();

    libMesh::DenseVector<libMesh::Number> Fe;
    libMesh::DenseMatrix<libMesh::Number>  Ke;

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_ux;
	std::vector<libMesh::dof_id_type> dof_indices_uy;
	std::vector<libMesh::dof_id_type> dof_indices_uz;

	// Grad U
	std::vector<double> disp_k;
    std::vector<double> vel_k;
    std::vector<double> vel_n;
    libMesh::VectorValue <libMesh::Number> Uk;
    libMesh::VectorValue <libMesh::Number> Vk;
    libMesh::VectorValue <libMesh::Number> Vn;
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

	  double body_force[3];
      const std::vector<libMesh::Point> & q_point = fe_u->get_xyz();

      double gamma_f;
      double gamma_s;
      double gamma_n;
      std::vector<double> gamma_f_k;
      libMesh::TensorValue <libMesh::Number> FA;

	    // On the boundary
     libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_disp));
     libMesh::QGauss qface(dim-1,  libMesh::FIRST);
     fe_face->attach_quadrature_rule (&qface);

     libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	 const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	    double rho;
	     libMesh::RealGradient f0;
         libMesh::RealGradient s0;
         libMesh::RealGradient n0;


//	    std::cout << "* ELASTICITY:loop starts ... " << std::endl;


  for ( ; el != end_el; ++el)
  {
      const libMesh::Elem * elem = *el;
      rho = M_materialMap[0]->M_density;
      auto elID = elem->id();

      dof_map.dof_indices (elem, dof_indices);

      dof_map.dof_indices (elem, dof_indices_ux, ux_var);
      if(dim>1)dof_map.dof_indices (elem, dof_indices_uy, uy_var);
      if(dim>2)dof_map.dof_indices (elem, dof_indices_uz, uz_var);

	const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_ux_dofs = dof_indices_ux.size();

      fe_u->reinit (elem);
      Fe.resize (n_dofs);
      Ke.resize(n_dofs, n_dofs);

      // get uk
      disp_k.resize(n_dofs);
//  	  for(auto && di : dof_indices)  std::cout << "dof id: " << di << std::endl;

      disp_system.current_local_solution->get(dof_indices, disp_k);
      system.current_local_solution->get(dof_indices, vel_k);
      system.old_local_solution->get(dof_indices, vel_n);

      if(activation_ptr)
      {
          gamma_f_k.resize(n_ux_dofs);
          for(int nd = 0; nd < n_ux_dofs; nd++ )
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

      for (unsigned int qp=0; qp<qrule_1.n_points(); qp++)
      {

          dUk *= 0.0;
          Ek *= 0.0;
          Sk *= 0.0;
          gamma_f *= 0.0;
          FA *= 0.0;
          Uk *= 0.0;
          Vk *= 0.0;
          const unsigned int n_phi = phi_u.size();
          for(unsigned int l = 0; l < n_phi; ++l )
          {
              for(int idim = 0; idim < dim; idim++)
              {
                  Uk(idim) += phi_u[l][qp]*disp_k[l+idim*n_phi];
                  Vk(idim) += phi_u[l][qp]*vel_k[l+idim*n_phi];
                  Vn(idim) += phi_u[l][qp]*vel_n[l+idim*n_phi];
                  for(int jdim = 0; jdim < dim; jdim++)
                  {
                      dUk(idim, jdim) += dphi_u[l][qp](jdim) *disp_k[l+idim*n_phi];
                  }
              }
          }

          if(activation_ptr)
          {
              for(unsigned int l = 0; l < n_phi; ++l )
              {
                  gamma_f += phi_u[l][qp] * gamma_f_k[l];
              }
              gamma_s = 1.0 / std::sqrt(1.0+gamma_f) - 1.0;
              gamma_n = gamma_s;
              for(int idim = 0; idim < dim; idim++)
              {
                  FA(idim, idim) += 1.0;
                  for(int jdim = 0; jdim < dim; jdim++)
                  {
                      FA(idim, jdim) += gamma_f * f0(idim) * f0(jdim)
                                      + gamma_s * s0(idim) * s0(jdim)
                                      + gamma_n * n0(idim) * n0(jdim);
                  }
              }
          }



    		M_materialMap[0]->M_gradU = dUk;
            M_materialMap[0]->M_FA = FA;
            M_materialMap[0]->updateVariables();
            M_materialMap[0]->evaluateStress( ElasticSolverType::Primal);
            Sk = M_materialMap[0]->M_PK1;
    	  	  // Residual
    		const double x = q_point[qp](0);
            const double y = q_point[qp](1);
            const double z = q_point[qp](2);
            const double time = 0.0;
            for(int idim = 0; idim < dim; idim++) body_force[idim] =M_rhsFunction(time,x,y,z,idim);

    	  	  for(unsigned int n = 0; n < phi_u.size(); ++n )
    	  	  {
					  for(int jdim= 0; jdim < dim; jdim++)
					  {
						  dW *= 0.0;
						  dW(jdim, 0) = JxW_u[qp]* dphi_u[n][qp](0);
						  dW(jdim, 1) = JxW_u[qp]* dphi_u[n][qp](1);
						  dW(jdim, 2) = JxW_u[qp]* dphi_u[n][qp](2);
						  // Compute  - \nabla \cdot \sigma + f
						  Fe(n+jdim*n_ux_dofs) -= dt * Sk.contract(dW);
						  Fe(n+jdim*n_ux_dofs) += dt * rho * body_force[jdim] * phi_u[n][qp];
                          Fe(n+jdim*n_ux_dofs) += rho * ( Vn(jdim) - Vk(jdim) ) * phi_u[n][qp];
					  }
    	  	  }

//                // Matrix
//                // For each test function
//                for(unsigned int n = 0; n < phi_u.size(); ++n )
//                {
//                    // for each trial function
//                    for(unsigned int m = 0; m < phi_u.size(); ++m )
//                    {
//                        // for each dimension of the test function
//                        for(int jdim= 0; jdim < dim; jdim++)
//                        {
//                        // Lumped mass matrix
//                        Ke(n+jdim*n_ux_dofs,m+jdim*n_ux_dofs) += JxW_u[qp] * phi_u[n][qp] * phi_u[m][qp];
//                        }
//                    }
//                }
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

							  M_materialMap[0]->evaluateJacobian(dU, 0.0);
							  S = M_materialMap[0]->M_total_jacobian;
							  auto int_1 = n+jdim*n_ux_dofs;
                              auto int_2 = m+idim*n_ux_dofs;
//	                            std::cout << int_1 << ", " << int_2 << ", SdW: " << S.contract(dW) << std::endl;

						  }
					   }
				  }
    	  	  }
      }

      apply_BC(elem, Ke, Fe, fe_face, qface, mesh, n_ux_dofs);
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
  }
    system.matrix->close();
    system.rhs->close();

}



void
DynamicElasticity::project_pressure()
{
    std::cout << "* ELASTICITY: projecting pressure ... " << std::endl;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearSystem & system = M_equationSystems.get_system<LinearSystem>(M_myName);
    LinearSystem & disp_system = M_equationSystems.get_system<LinearSystem>("displacement");
    LinearSystem & system_p = M_equationSystems.get_system<LinearSystem>("Pressure_Projection");
    system.update();
    unsigned int ux_var = disp_system.variable_number ("displacementx");
    unsigned int uy_var, uz_var, p_var;

   if(dim>1)uy_var =  disp_system.variable_number ("displacementy");
   if(dim>2)uz_var =  disp_system.variable_number ("displacementz");
    p_var = system_p.variable_number ("pressure");

    const libMesh::DofMap & dof_map = system.get_dof_map();
    libMesh::FEType fe_type = dof_map.variable_type(0);
    libMesh::UniquePtr<libMesh::FEBase> fe (libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule (dim, fe->get_order() );
    fe->attach_quadrature_rule (&qrule);

    const libMesh::DofMap & dof_map_p = system_p.get_dof_map();
    libMesh::FEType fe_type_p = dof_map_p.variable_type(0);
    libMesh::UniquePtr<libMesh::FEBase> fe_p (libMesh::FEBase::build(dim, fe_type_p));
    libMesh::QGauss qrule_p (dim,libMesh::SECOND);
    fe_p->attach_quadrature_rule (&qrule_p);


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

      libMesh::TensorValue <libMesh::Number> dUk;


    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> dof_indices_ux;
    std::vector<libMesh::dof_id_type> dof_indices_uy;
    std::vector<libMesh::dof_id_type> dof_indices_uz;
    std::vector<libMesh::dof_id_type> dof_indices_p;

    libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for ( ; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;

        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_ux, ux_var);
        if(dim > 1)
        {
            dof_map.dof_indices (elem, dof_indices_uy, uy_var);
        }
        if(dim > 2)
        {
            dof_map.dof_indices (elem, dof_indices_uz, uz_var);
        }
        dof_map_p.dof_indices (elem, dof_indices_p, p_var);

        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_ux_dofs = dof_indices_ux.size();
        const unsigned int n_uy_dofs = dof_indices_uy.size();
        const unsigned int n_uz_dofs = dof_indices_uz.size();
        const unsigned int n_p_dofs = dof_indices_p.size();

        fe->reinit (elem);
        fe_p->reinit (elem);

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
               for(unsigned int l = 0; l < n_phi; ++l )
               {
                   for(int idim = 0; idim < dim; idim++)
                   {
                       for(int jdim = 0; jdim < dim; jdim++)
                       {
                           dUk(idim, jdim) += dphi_u[l][qp](jdim) *solution_k[l+idim*n_phi];
                       }
                   }
               }
//               std::cout << "* ELASTICITY: evaluate p ... " << std::endl;

               M_materialMap[0]->M_gradU = dUk;
               double p = M_materialMap[0]->evaluatePressure();

               for (unsigned int i = 0; i < phi.size(); i++)
               {
                   Fe(i) +=  JxW[qp] * p* phi[i][qp];

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

                 system_p.matrix->add_matrix (Ke, dof_indices_p);
                 system_p.rhs->add_vector    (Fe, dof_indices_p);

      }
//    std::cout << "* ELASTICITY: solving p ... " << std::endl;

    system_p.matrix->close();
    system_p.rhs->close();

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

    rval = M_projectionsLinearSolver->solve (*system_p.matrix, nullptr,
                                    *system_p.solution,
                                                        *system_p.rhs, tol, max_iter);
}



} /* namespace BeatIt */
