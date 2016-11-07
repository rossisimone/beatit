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

    libMesh::QGauss qrule_1(dim, libMesh::FOURTH  );
    libMesh::QGauss qrule_2(dim,  libMesh::FOURTH   );

    fe_u->attach_quadrature_rule(&qrule_1);
    fe_p->attach_quadrature_rule(&qrule_2);

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
	  // Pressure
	  double pk, p, q, dp;

      double body_force[3];
      const std::vector<libMesh::Point> & q_point = fe_u->get_xyz();
	    // On the boundary
    libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_disp));
    libMesh::QGauss qface(dim-1,  libMesh::FOURTH);
    fe_face->attach_quadrature_rule (&qface);

	  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();


//	    std::cout << "* ELASTICITY:loop starts ... " << std::endl;
	  int assemble_pressure_mass = 1;

  for ( ; el != end_el; ++el)
  {
      const libMesh::Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);

      dof_map.dof_indices (elem, dof_indices_ux, ux_var);
      if(dim>1)dof_map.dof_indices (elem, dof_indices_uy, uy_var);
      if(dim>2)dof_map.dof_indices (elem, dof_indices_uz, uz_var);

	const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_ux_dofs = dof_indices_ux.size();

      fe_u->reinit (elem);
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

//    	  	  	    std::cout << "* ELASTICITY: evaluate dUk ... " << std::endl;
//
//    	  	  for(auto && di : dof_indices)  std::cout << "dof id: " << di << std::endl;
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

		  M_materialMap[0]->M_gradU = dUk;
			const unsigned int n_phi_p = phi_p.size();
			for(unsigned int l = 0; l < n_phi_p; ++l )
			{
				pk += phi_p[l][qp] * solution_k[l+dim*n_phi];
			}

			M_materialMap[0]->M_pressure = pk;

            M_materialMap[0]->evaluateStress( ElasticSolverType::Mixed);
            Sk = M_materialMap[0]->M_total_stress;
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
						  Fe(n+jdim*n_ux_dofs) -= Sk.contract(dW);
						  Fe(n+jdim*n_ux_dofs) += body_force[jdim] * phi_u[n][qp];
					  }
    	  	  }

			p =  M_materialMap[0]->evaluatePressureResidual();
				  for(unsigned int n = 0; n < phi_p.size(); ++n )
				  {
						  Fe(n+dim*n_ux_dofs) -=   JxW_p[qp] * p * phi_p[n][qp];
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
						  }
				  }
         	  	  // For each pressure trial function
				  for(unsigned int m = 0; m < phi_p.size(); ++m )
				  {
							  Ke(n+dim*n_ux_dofs,m+dim*n_ux_dofs) += assemble_pressure_mass * q * phi_p[m][qp];
				  }
    	  	  }

      }

//      Ke.print(std::cout);
//      std::cout << std::endl;
//      std::cout << std::endl;

      Elasticity::apply_BC(elem, Ke, Fe, fe_face, qface, mesh, n_ux_dofs);

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);


  }
    system.matrix->close();
    system.rhs->close();
}


} /* namespace BeatIt */
