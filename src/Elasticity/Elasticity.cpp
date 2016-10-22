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
 * Elasticity.cpp
 *
 *  Created on: Oct 19, 2016
 *      Author: srossi
 */

#include "Elasticity/Elasticity.hpp"
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

namespace BeatIt {

typedef libMesh::LinearImplicitSystem    LinearSystem;


Elasticity::Elasticity( libMesh::EquationSystems& es, std::string system_name )
    : M_equationSystems(es)
    , M_exporter()
    , M_outputFolder()
    , M_datafile()
    , M_linearSolver()
    , M_bch()
    , M_rhsFunction()
    , M_myName(system_name)
{
	// TODO Auto-generated constructor stub

}

Elasticity::~Elasticity()
{
	// TODO Auto-generated destructor stub
}

void
Elasticity::write_equation_system(const std::string& es)
{
    M_equationSystems.write(es, libMesh::WRITE);
}

void
Elasticity::read_equation_system(const std::string& es)
{
    M_equationSystems.clear ();
    M_equationSystems.read(es, libMesh::READ);
}

void
Elasticity::init_exo_output(const std::string& output_filename)
{
    M_exporter->write_equation_systems (M_outputFolder+output_filename, M_equationSystems);
    M_exporter->append(true);
}


void
Elasticity::save_exo(const std::string& output_filename, int step, double time)
{
    std::cout << "* ELASTICITY: EXODUSII::Exporting in: "  << M_outputFolder << " ... " << std::flush;
    M_exporter->write_timestep (M_outputFolder+output_filename, M_equationSystems, step, time);
//	M_exporter->write_element_data(M_equationSystems);
    std::cout << "done " << std::endl;
}

void
Elasticity::deleteSystems()
{
	M_equationSystems.delete_system(M_myName);
}


void
Elasticity::setup(const GetPot& data, std::string section )
{
    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Read Input File
    M_datafile = data;
    //Read output folder from datafile
    std::string output_folder = M_datafile(section+"/output_folder",  "Output");
    M_outputFolder = "./" + output_folder + "/";
    // Create folder to save output
    BeatIt::createOutputFolder(M_equationSystems.get_mesh().comm(), M_outputFolder);

    int dimension = M_equationSystems.get_mesh().mesh_dimension();
	// ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // 1) ADR
    std::cout << "* ELASTICITY: Creating new System for the ELASTICITY Solver" << std::endl;
    LinearSystem& system  =  M_equationSystems.add_system<LinearSystem>(M_myName);

    std::string disp_name = "displacement";

    int ord = M_datafile(section+"/order",  1);
    std::string fam = M_datafile(section+"/fefamily",  "lagrange");
    std::cout << "* ELASTICITY: Setting up displacement field - order: " << ord << ", fefamily = " << fam  << std::endl;
    libMesh::Order  order  = BeatIt::libmesh_order_map.find(ord)->second;
    libMesh::FEFamily  fefamily  =BeatIt::libmesh_fefamily_map.find(fam)->second;
    system.add_variable(disp_name+"x", order, fefamily);
    if(dimension > 1) system.add_variable(disp_name+"y", order, fefamily);
    if(dimension > 2) system.add_variable(disp_name+"z", order, fefamily);

    // PRESSURE SYSTEM
    std::string formulation = M_datafile(section+"/formulation",  "primal");

    if(formulation == "mixed")
    {
		std::string pressure_name = "pressure";
		ord = M_datafile(section+"/p_order",  1);
		fam = M_datafile(section+"/p_fefamily",  "lagrange");
		std::cout << "* ELASTICITY: Setting up pressure  field - order: " << ord << ", fefamily = " << fam  << std::endl;
		libMesh::Order  p_order  = BeatIt::libmesh_order_map.find(ord)->second;
		libMesh::FEFamily  p_fefamily   = BeatIt::libmesh_fefamily_map.find(fam)->second;
		system.add_variable(pressure_name, p_order, p_fefamily);
    }
    system.add_vector("residual");
    system.add_vector("step");



  std::cout << "* ELASTICITY: Reading BC ... " << std::endl;
    M_bch.readBC(data, section);
    M_bch.showMe();
  std::cout << "* ELASTICITY: Setup Homogenuous Dirichlet BC ... " << std::endl;

  for(auto&& bc_ptr : M_bch.M_bcs)
  {
	  auto bc_type = bc_ptr->get_type();
	  if(bc_type == BCType::Dirichlet)
	  {
		  std::set< libMesh::boundary_id_type> dirichlet_boundary_ids;
		  auto bc_mode = bc_ptr->get_mode();
		  std::vector<unsigned int> variables;

		  switch(bc_mode)
		  {
		  	  default:
			  case BCMode::Full:
			  {
				  for(unsigned int k = 0; k < dimension; ++k) variables.push_back(k);
				  break;
			  }
			  case BCMode::Component:
			  {
				  auto component = bc_ptr->get_component();
				  switch(component)
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

		  for(auto&& flag : bc_ptr->M_flag)
		  {
			  // Create a ZeroFunction to initialize dirichlet_bc
			  dirichlet_boundary_ids.insert(flag);
		  }

		  libMesh::ZeroFunction<> zf;
		  libMesh::DirichletBoundary dirichlet_bc(dirichlet_boundary_ids,
																				 variables,
																				 &zf);

		  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

	  } // end if Dirichlet
  } // end for loop on BC

    /// Init system
  std::cout << "* ELASTICITY: Init System ... " << std::endl;
	system.init();
    M_equationSystems.print_info();
    /// Setting up BC

  std::cout << "* ELASTICITY: Setup body forces ... " << std::endl;

	std::string rhs = M_datafile(section+"/rhs", "no_rhs");
    M_rhsFunction.read(rhs);
    M_rhsFunction.showMe();

  std::cout << "* ELASTICITY: Setup linear solvers ... " << std::endl;

	M_linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( M_equationSystems.comm() );
    M_linearSolver->set_solver_type(libMesh::GMRES);
    M_linearSolver->set_preconditioner_type( libMesh::AMG_PRECOND);
    M_linearSolver->init();
	M_projectionsLinearSolver =  libMesh::LinearSolver<libMesh::Number>::build( M_equationSystems.comm() );
    M_projectionsLinearSolver->set_solver_type(libMesh::GMRES);
    M_projectionsLinearSolver->set_preconditioner_type( libMesh::SOR_PRECOND);
    M_projectionsLinearSolver->init();


    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Setup Exporters
    M_exporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );

    std::string mats = data(section+"/materials", "NO_MATERIAL");
    std::vector<std::string> materials;
    BeatIt::readList(mats, materials);
    std::string base_path = section+"/materials/";
    for(auto&& material : materials)
    {
    	std::string path = base_path + material;
    	unsigned int materialID = data(path+"/matID", 0);
    // Setup material
    	M_materialMap[materialID].reset( Material::MaterialFactory::Create(material));
    	M_materialMap[materialID]->setup(M_datafile, path);
    }


}


void
Elasticity::assemble_residual()
{
  std::cout << "* ELASTICITY: assembling ... " << std::endl;

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
   p_var = 0;
//   p_var =  system.variable_number ("pressure");

    const libMesh::DofMap & dof_map = system.get_dof_map();
    libMesh::FEType fe_disp= dof_map.variable_type(ux_var);
    libMesh::FEType fe_pr = dof_map.variable_type(p_var);

    UniquePtr<libMesh::FEBase> fe_u(libMesh::FEBase::build(dim, fe_disp) );
    UniquePtr<libMesh::FEBase> fe_p(libMesh::FEBase::build(dim, fe_pr) );

    libMesh::QGauss qrule_1(dim,libMesh::FIRST  );
    libMesh::QGauss qrule_2(dim, fe_pr.default_quadrature_order()  );

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
	  double pk;

	  double body_force[3];
  const std::vector<libMesh::Point> & q_point = fe_u->get_xyz();
	    // On the boundary
    libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_disp));
    libMesh::QGauss qface(dim-1, fe_disp.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);


	  libMesh::MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	  const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

//	    std::cout << "* ELASTICITY:loop starts ... " << std::endl;


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
      fe_p->reinit (elem);
      Fe.resize (n_dofs);
      Ke.resize(n_dofs, n_dofs);

      // get uk
      solution_k.resize(n_dofs);
//  	  for(auto && di : dof_indices)  std::cout << "dof id: " << di << std::endl;

      system.current_local_solution->get(dof_indices, solution_k);
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

//    	  	  	    std::cout << "* ELASTICITY: evaluate Sk ... " << std::endl;
    	  	  M_materialMap[0]->M_gradU = dUk;
    	  	  M_materialMap[0]->evaluateStress();
//    	  	  Ek = 0.5 * ( dUk + dUk.transpose() );
//    	  	  double trEk = Ek(0,0)+Ek(1,1)+Ek(2,2);
//    	  	  Sk = 2.0 * mu * ( Ek - 1.0 / 3.0 * trEk * id ) + kappa * trEk * id;
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
						  Fe(n+jdim*n_ux_dofs) -= Sk.contract(dW);
						  Fe(n+jdim*n_ux_dofs) += body_force[jdim] * phi_u[n][qp];
					  }

    	  	  }
    	  	  // Matrix
//    	  	  	    std::cout << "* ELASTICITY: evaluate J ... " << std::endl;
    	  	  for(unsigned int m = 0; m < phi_u.size(); ++m )
    	  	  {
				  for(unsigned int n = 0; n < phi_u.size(); ++n )
				  {
					  for(int idim = 0; idim < dim; idim++)
					  {
						  dU*= 0.0;
						  E*= 0.0;
						  S*= 0.0;

						  dU(idim, 0) =  dphi_u[m][qp](0);
						  dU(idim, 1) =  dphi_u[m][qp](1);
						  dU(idim, 2) =  dphi_u[m][qp](2);

						  M_materialMap[0]->evaluateJacobian(dU, 0.0);
//						  E = 0.5 * (dU + dU.transpose() );
//						  double trE = E(0,0)+E(1,1)+E(2,2);
//						  S = 2.0 * mu * (E - 1.0/3.0 * trE * id) + kappa * trE * id;
						  S = M_materialMap[0]->M_total_jacobian;
						  for(int jdim= 0; jdim < dim; jdim++)
						  {
							  dW *= 0.0;
							  dW(jdim, 0) = JxW_u[qp]* dphi_u[n][qp](0);
							  dW(jdim, 1) = JxW_u[qp]* dphi_u[n][qp](1);
							  dW(jdim, 2) = JxW_u[qp]* dphi_u[n][qp](2);
//    	  	  	    std::cout << "* ELASTICITY: contract ... " << std::endl;
//    	  	  	    double val = ;
//    	  	  	    std::cout << "* ELASTICITY: Add val ... " << std::endl;
							  Ke(n+jdim*n_ux_dofs,m+idim*n_ux_dofs) += S.contract(dW);
						  }
					   }
				  }
    	  	  }

      }

//      std::cout << "\n Fe1: " << mu <<", " << kappa <<  std::endl;
//      Fe.print(std::cout);

      apply_BC(elem, Ke, Fe, fe_face, qface, mesh, n_ux_dofs);
//            std::cout << "\n Fe2: " << mu <<", " << kappa <<  std::endl;
//      Fe.print(std::cout);

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

//      std::cout << "\n Ke: " << mu <<", " << kappa <<  std::endl;
//      Ke.print(std::cout);
////    	  	  	    std::cout << "* ELASTICITY: BC ... " << std::endl;
//
//
//
////    	  	  	    std::cout << "* ELASTICITY: Constrain ... . " << std::endl;
//
//      std::cout << "\n Fe3: " << mu <<", " << kappa <<  std::endl;
//      Fe.print(std::cout);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);





  }
    system.matrix->close();
    system.rhs->close();

}


void
Elasticity::apply_BC( const libMesh::Elem*& elem,
                   libMesh::DenseMatrix<libMesh::Number>& Ke,
                   libMesh::DenseVector<libMesh::Number>& Fe,
                   libMesh::UniquePtr<libMesh::FEBase>& fe_face,
                   libMesh::QGauss& qface,
                   const libMesh::MeshBase& mesh, int n_ux_dofs)
{
	const unsigned int dim = mesh.mesh_dimension();
	for (unsigned int side=0; side<elem->n_sides(); side++)
	{
		  if (elem->neighbor(side) == libmesh_nullptr)
    	  {
    			  const unsigned int boundary_id =
                mesh.boundary_info->boundary_id (elem, side);

                auto bc = M_bch.get_bc(boundary_id);

                if(bc)
                {
					auto bc_type = bc->get_type();
                    switch(bc_type)
                    {
						case BCType::Neumann:
						{
							const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();
							const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();
							const std::vector<libMesh::Point> & qface_point = fe_face->get_xyz();
							const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_face->get_dphi();
							const std::vector<libMesh::Point>&  normals = fe_face->get_normals();
							fe_face->reinit(elem, side);
						    for (unsigned int qp=0; qp<qface.n_points(); qp++)
						    {
								const double xq = qface_point[qp](0);
                                const double yq = qface_point[qp](1);
                                const double zq = qface_point[qp](2);
						    	for(int idim = 0; idim < dim; ++idim)
						    	{
									const double traction = bc->get_function()(0.0, xq, yq, zq, idim);
									for (unsigned int i=0; i< n_ux_dofs; i++)
									{
										Fe(i+idim*n_ux_dofs) += JxW_face[qp] * traction * phi_face[i][qp];
									}
						    	}
						    }
							break;
						}
						default:
						{
							break;
						}
                    }

                }
		 }
	}
}

void
Elasticity::solve_system()
{
	    LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

    rval = M_linearSolver->solve (*system.matrix, nullptr,
														system.get_vector("step"),
														*system.rhs, tol, max_iter);


}


void
Elasticity::newton()
{
	    LinearSystem& system  =  M_equationSystems.get_system<LinearSystem>(M_myName);
	    assemble_residual();
	    auto res_norm = system.rhs->linfty_norm ();
		double tol = 1e-9;
		int max_iter = 20;
		int iter = 0;

		std::cout << "* ELASTICITY: Performing Newton solve: " << std::endl;
		std::cout << "\t\t\t  iter: " << iter << ", residual: " << res_norm << std::endl;

 	    while(res_norm > tol && iter < max_iter)
 	    {
             iter++;
 	    	solve_system();
 	    	(*system.solution) += system.get_vector("step");

 	      	assemble_residual();
             res_norm = system.rhs->linfty_norm ();
			std::cout << "\t\t\t  iter: " << iter << ", residual: " << res_norm << std::endl;
 	    }

		std::cout << "* ELASTICITY: Newton solve completed. " << std::endl;

}

} /* namespace BeatIt */
