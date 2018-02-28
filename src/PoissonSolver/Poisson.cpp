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
 * Poisson.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */


#include "PoissonSolver/Poisson.hpp"

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

namespace BeatIt {

typedef libMesh::LinearImplicitSystem    PoissonSystem;

//Poisson::Poisson(libMesh::MeshBase & mesh)
//    : M_equationSystems(mesh)
//    , M_exporter()
//    , M_outputFolder()
//    , M_datafile()
//    , M_linearSolver()
//    , M_bch()
//	, M_rhsFunction()
//{
//
//}


Poisson::Poisson( libMesh::EquationSystems& es, std::string system_name )
    : M_equationSystems(es)
    , M_exporter()
    , M_outputFolder()
    , M_datafile()
    , M_linearSolver()
    , M_bch()
    , M_rhsFunction()
    , M_myName(system_name)
	, M_myNameGradient(system_name+"_gradient")
	, M_myNameP0(system_name+"_P0")
{

}

Poisson::~Poisson() {
	// TODO Auto-generated destructor stub
}

void
Poisson::deleteSystems()
{
	M_equationSystems.delete_system(M_myName);
	M_equationSystems.delete_system(M_myNameGradient);
	M_equationSystems.delete_system(M_myNameP0);
}

void
Poisson::setup(const GetPot& data, std::string section )
{
    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Read Input File
    M_datafile = data;
    //Read output folder from datafile
    std::string output_folder = M_datafile(section+"/output_folder",  "Output");
    M_outputFolder = "./" + output_folder + "/";
    std::cout << "* POISSON: Output folder: " << M_outputFolder << std::endl;
        // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // 1) ADR
    std::cout << "* POISSON: Creating new System for the Poisson Solver" << std::endl;
    PoissonSystem& system  =  M_equationSystems.add_system<PoissonSystem>(M_myName);
    // TO DO: Generalize to higher order
    system.add_variable( M_myName+"_phi", libMesh::FIRST);
    libMesh::ExplicitSystem& grad_system  =  M_equationSystems.add_system<libMesh::ExplicitSystem>(M_myNameGradient);
    grad_system.add_variable( M_myName+"_dphix", libMesh::CONSTANT, libMesh::MONOMIAL);
    grad_system.add_variable( M_myName+"_dphiy", libMesh::CONSTANT, libMesh::MONOMIAL);
    grad_system.add_variable( M_myName+"_dphiz", libMesh::CONSTANT, libMesh::MONOMIAL);

    libMesh::ExplicitSystem& sol_system  =  M_equationSystems.add_system<libMesh::ExplicitSystem>(M_myNameP0);
    sol_system.add_variable( M_myName+"_phi_p0", libMesh::CONSTANT, libMesh::MONOMIAL);

	system.init();
	grad_system.init();
	sol_system.init();
    M_equationSystems.print_info();
    /// Setting up BC
    M_bch.readBC(data, section);
    M_bch.showMe();

    std::string rhs = data(section+"/rhs", "0.0");
    M_rhsFunction.read(rhs);
    M_rhsFunction.showMe();
//    std::map<std::string, libMesh::SolverType> solver_map;
//    solver_map["cg"] = libMesh::CG;
//    solver_map["cgs"] = libMesh::CGS;
//    solver_map["gmres"] = libMesh::GMRES;
//    std::map<std::string, libMesh::PreconditionerType> prec_map;
//    prec_map["jacobi"] =  libMesh::JACOBI_PRECOND;
//    prec_map["sor"] =  libMesh::SOR_PRECOND;
//    prec_map["ssor"] =  libMesh::SSOR_PRECOND;
//    prec_map["cholesky"] =  libMesh::CHOLESKY_PRECOND;
//    prec_map["lu"] =  libMesh::LU_PRECOND;
//    prec_map["ilu"] =  libMesh::ILU_PRECOND;
//    prec_map["amg"] =  libMesh::AMG_PRECOND;
//    M_linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( M_equationSystems.comm() );
    typedef libMesh::PetscLinearSolver<libMesh::Number> PetscSolver;

    M_linearSolver.reset(new PetscSolver(M_equationSystems.comm()));
    M_linearSolver->set_solver_type(libMesh::CG);
    M_linearSolver->set_preconditioner_type( libMesh::AMG_PRECOND);
    M_linearSolver->init();
    KSPSetOptionsPrefix(M_linearSolver->ksp(),"poisson_");
    PCSetOptionsPrefix(M_linearSolver->pc(),"poisson_");

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Setup Exporters
    M_exporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );

	struct stat out_dir;
    if( stat(&M_outputFolder[0],&out_dir) != 0  )
    {
        if ( system.get_mesh().comm().rank() == 0 )
        {
            mkdir ( M_outputFolder.c_str(), 0777 );
        }
    }

}


double
Poisson::get_solution_norm()
{
    PoissonSystem& system  =  M_equationSystems.get_system<PoissonSystem>(M_myName);
    return system.solution->l1_norm();
}


void
Poisson::solve_system()
{
	    PoissonSystem& system  =  M_equationSystems.get_system<PoissonSystem>(M_myName);

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

    rval = M_linearSolver->solve (*system.matrix, *system.solution,
														*system.rhs, tol, max_iter);
}


void
Poisson::save_exo(const std::string& output_filename)
{
    std::cout << "* POISSON: EXODUSII::Exporting in: "  << M_outputFolder << " ... " << std::flush;
    M_exporter->write_equation_systems (M_outputFolder+output_filename, M_equationSystems);
    M_exporter->write_element_data(M_equationSystems);
    std::cout << "done " << std::endl;
}

void
Poisson::write_equation_system(const std::string& es)
{
    M_equationSystems.write(es, libMesh::WRITE);
}

void
Poisson::read_equation_system(const std::string& es)
{
    M_equationSystems.clear ();
    M_equationSystems.read(es, libMesh::READ);
}


const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Poisson::get_solution()
{
    libMesh::ExplicitSystem& p_system  =  M_equationSystems.get_system<PoissonSystem>(M_myName);
    return p_system.solution;
}

const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Poisson::get_P0_solution()
{
    libMesh::ExplicitSystem& s_system  =  M_equationSystems.get_system<libMesh::ExplicitSystem>(M_myNameP0);
    return s_system.solution;
}

const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Poisson::get_gradient()
{
    libMesh::ExplicitSystem& g_system  =  M_equationSystems.get_system<libMesh::ExplicitSystem>(M_myNameGradient);
    return g_system.solution;
}

void
Poisson::assemble_system()
{
     using libMesh::UniquePtr;

     const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
     const unsigned int dim = mesh.mesh_dimension();

      PoissonSystem& system  =  M_equationSystems.get_system<PoissonSystem>(M_myName);
     system.matrix->zero();
          // A reference to the  DofMap object for this system.  The  DofMap
     // object handles the index translation from node and element numbers
     // to degree of freedom numbers.  We will talk more about the  DofMap
     // in future examples.
     const libMesh::DofMap & dof_map = system.get_dof_map();
          // Get a constant reference to the Finite Element type
     // for the first (and only) variable in the system.
     libMesh::FEType fe_type = dof_map.variable_type(0);

     // Build a Finite Element object of the specified type.  Since the
     // FEBase::build() member dynamically creates memory we will
     // store the object as a UniquePtr<FEBase>.  This can be thought
     // of as a pointer that will clean up after itself.  Introduction Example 4
     // describes some advantages of  UniquePtr's in the context of
     // quadrature rules.
     UniquePtr<libMesh::FEBase> fe_qp(libMesh::FEBase::build(dim, fe_type));
      // A 5th order Gauss quadrature rule for numerical integration.
     libMesh::QGauss qrule_stiffness(dim, libMesh::FOURTH);
      // Tell the finite element object to use our quadrature rule.
     fe_qp->attach_quadrature_rule(&qrule_stiffness);

       // Declare a special finite element object for
  // boundary integration.
  UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));

  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  libMesh::QGauss qface(dim-1, libMesh::FOURTH);

  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);
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
     libMesh::DenseMatrix<libMesh::Number> Ke;
      libMesh::DenseVector<libMesh::Number> Fe;

           // This vector will hold the degree of freedom indices for
     // the element.  These define where in the global system
     // the element degrees of freedom get mapped.
     std::vector<libMesh::dof_id_type> dof_indices;
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

     double ftxyzi = 0.0;

     for (; el != end_el; ++el)
     {
         const libMesh::Elem * elem = *el;
         dof_map.dof_indices(elem, dof_indices);
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
         Ke.resize(dof_indices.size(), dof_indices.size());
          Fe.resize(dof_indices.size());

        //  std::cout << "Assembling Poisson ... " << std::endl;
         for (unsigned int qp = 0; qp < qrule_stiffness.n_points(); qp++)
         {
             for (unsigned int i = 0; i < phi_qp.size(); i++)
             {
                 for (unsigned int j = 0; j < phi_qp.size(); j++)
                 {
                     // stiffness term
                     Ke(i, j) += JxW_qp[qp] * dphi_qp[i][qp] * dphi_qp[j][qp];
                 }
                 ftxyzi = M_rhsFunction(0.0,  q_point_qp[qp](0),  q_point_qp[qp](0),  q_point_qp[qp](0), 0);
                 Fe(i) +=  JxW_qp[qp] * ftxyzi * phi_qp[i][qp];
             }
         }

          apply_BC(elem, Ke, Fe, fe_face, qface, mesh);


         system.matrix->add_matrix(Ke, dof_indices);
         system.rhs->add_vector    (Fe, dof_indices);

     }
}

void
Poisson::apply_BC( const libMesh::Elem*& elem,
                   libMesh::DenseMatrix<libMesh::Number>& Ke,
                   libMesh::DenseVector<libMesh::Number>& Fe,
                   libMesh::UniquePtr<libMesh::FEBase>& fe_face,
                   libMesh::QGauss& qface,
                   const libMesh::MeshBase& mesh)
{
        for (unsigned int side=0; side<elem->n_sides(); side++)
        {

//            if (elem->neighbor(side) == libmesh_nullptr)

            {
                const unsigned int boundary_id =
                mesh.boundary_info->boundary_id (elem, side);
                //std::cout << "BID: " << boundary_id << std::endl;

                auto bc = M_bch.get_bc(boundary_id);

                if(bc)
                {
                    // The Jacobian * Quadrature Weight at the quadrature
                    // points on the face.
                    const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();
                    const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();
                    const std::vector<libMesh::Point> & qface_point = fe_face->get_xyz();
                    const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_face->get_dphi();
                    const std::vector<libMesh::Point>&  normals = fe_face->get_normals();
                    fe_face->reinit(elem, side);

                    auto bc_type = bc->get_type();
                    switch(bc_type)
                    {
                        case BCType::Dirichlet:
                        {
                            double beta = 1e6;


                            for (unsigned int qp=0; qp<qface.n_points(); qp++)
                            {
                                // The location on the boundary of the current
                                // face quadrature point.
                                const double xf = qface_point[qp](0);
                                const double yf = qface_point[qp](1);
                                const double zf = qface_point[qp](2);
                                const double value = bc->get_function()(0.0, xf, yf, zf, 0);
                                for (unsigned int i=0; i<phi_face.size(); i++)
                                {
                                    double nabla_test_normal = normals[qp]* dphi_face[i][qp];
                                    for (unsigned int j=0; j<phi_face.size(); j++)
                                    {
                                        double nabla_trial_normal = normals[qp] * dphi_face[j][qp];
                                        Ke(i,j) += JxW_face[qp] * beta * phi_face[i][qp] * phi_face[j][qp];
                                    }
                                    Fe(i) += JxW_face[qp] * beta * value * phi_face[i][qp];

                                }
                            }
                            break;
                        }
                    case BCType::NitscheSymmetric:
                    {
                        double hE = elem->side(side)->hmax();
                        double beta = 1e4 / hE;

                        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                        {
                            // The location on the boundary of the current
                            // face quadrature point.
                            const double xf = qface_point[qp](0);
                            const double yf = qface_point[qp](1);
                            const double zf = qface_point[qp](2);
                            const double value = bc->get_function()(0.0, xf, yf,
                                    zf, 0);
                            for (unsigned int i = 0; i < phi_face.size(); i++)
                            {
                                double nabla_test_normal = normals[qp]
                                        * dphi_face[i][qp];
                                for (unsigned int j = 0; j < phi_face.size();
                                        j++)
                                {
                                    double nabla_trial_normal = normals[qp]
                                            * dphi_face[j][qp];
                                    Ke(i, j) += JxW_face[qp] * beta
                                            * phi_face[i][qp] * phi_face[j][qp];
                                    Ke(i, j) -= JxW_face[qp] * nabla_test_normal
                                            * phi_face[j][qp];
                                    Ke(i, j) -= JxW_face[qp]
                                            * nabla_trial_normal
                                            * phi_face[i][qp];
                                }
                                Fe(i) += JxW_face[qp] * beta * value
                                        * phi_face[i][qp];
                                Fe(i) -= JxW_face[qp] * value
                                        * nabla_test_normal;

                            }
                        }
                        break;
                    }
                        case BCType::NitscheUnsymmetric:
                        {
                            const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_face->get_dphi();

                            // The XYZ locations (in physical space) of the
                            // quadrature points on the face.  This is where
                            // we will interpolate the boundary value function.
                            const std::vector<libMesh::Point>&  normals = fe_face->get_normals();

                            for (unsigned int qp=0; qp<qface.n_points(); qp++)
                            {
                                // The location on the boundary of the current
                                // face quadrature point.
                                const double xf = qface_point[qp](0);
                                const double yf = qface_point[qp](1);
                                const double zf = qface_point[qp](2);
                                const double value = bc->get_function()(0.0, xf, yf, zf, 0);
                                libMesh::RealVectorValue normal( normals[qp]);

                                // Matrix contribution of the L2 projection.
                                for (unsigned int i=0; i<phi_face.size(); i++)
                                {
                                    double nabla_test_normal = normals[qp] * dphi_face[i][qp];
                                    for (unsigned int j=0; j<phi_face.size(); j++)
                                    {
                                        double nabla_trial_normal = normals[qp] * dphi_face[j][qp];
                                        Ke(i,j) += JxW_face[qp] * nabla_test_normal      * phi_face[j][qp];
                                        Ke(i,j) -= JxW_face[qp] * nabla_trial_normal     * phi_face[i][qp];
                                    }
                                    Fe(i) += JxW_face[qp]        * value * nabla_test_normal;

                                }
                            }
                            break;
                        }
                        case BCType::Neumann:
                        {
                            for (unsigned int qp=0; qp<qface.n_points(); qp++)
                            {
                                // The location on the boundary of the current
                                // face quadrature point.
                                const double xf = qface_point[qp](0);
                                const double yf = qface_point[qp](1);
                                const double zf = qface_point[qp](2);
                                const double value = bc->get_function()(0.0, xf, yf, zf, 0);

                                // Matrix contribution of the L2 projection.
                                for (unsigned int i=0; i<phi_face.size(); i++)
                                {
                                    Fe(i) += JxW_face[qp] * value * phi_face[i][qp];
                                }
                            }
                            break;
                        }
                        default:
                        {
                            throw std::runtime_error("- Error - Wrong BC Type for Poisson");
                            break;
                        }
                    }


                }
            }
        }
}


void
Poisson::compute_elemental_solution_gradient()
{
    std::cout << "* POISSON: Evaluating the gradient ... " << std::flush;
    libMesh::ExplicitSystem& g_system  =  M_equationSystems.get_system<libMesh::ExplicitSystem>(M_myNameGradient);
    libMesh::ExplicitSystem& s_system  =  M_equationSystems.get_system<libMesh::ExplicitSystem>(M_myNameP0);
    PoissonSystem& p_system  =  M_equationSystems.get_system<PoissonSystem>(M_myName);
    p_system.update();
    using libMesh::UniquePtr;

    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap & g_dof_map = g_system.get_dof_map();
    const libMesh::DofMap & p_dof_map = p_system.get_dof_map();
    const libMesh::DofMap & s_dof_map = s_system.get_dof_map();
    libMesh::FEType g_fe_type = g_dof_map.variable_type(0);
    libMesh::FEType p_fe_type = p_dof_map.variable_type(0);
    UniquePtr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, p_fe_type));
    libMesh::QGauss qrule(dim, libMesh::FIRST);
    if( qrule.get_order() != libMesh::FIRST )
    {
        std::cout << "\n- Error - computing gradients assumes single elemental value, a single quad point" << std::endl;;
        throw std::runtime_error("- Error - Poisson: wrong number of quadrature points!");
    }
    fe->attach_quadrature_rule(&qrule);
    const std::vector<std::vector<libMesh::Real > > & phi = fe->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > & dphi = fe->get_dphi();

    std::vector<libMesh::dof_id_type> p_dof_indices;
    std::vector<libMesh::dof_id_type> g_dof_indices;
    std::vector<libMesh::dof_id_type> s_dof_indices;


    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();

    libMesh::RealGradient grad;

    std::vector<double> sol(1,0.0);
    std::vector<double> elemental_solution;

    std::vector<double> grad_vector(3,0.0);
    double elSol = 0.0;
    double nodeSol = 0.0;
    for (; el != end_el; ++el)
    {
        const libMesh::Elem * elem = *el;
        p_dof_map.dof_indices(elem, p_dof_indices);
        g_dof_map.dof_indices(elem, g_dof_indices);
        s_dof_map.dof_indices(elem, s_dof_indices);

        fe->reinit(elem);

        grad *= 0.0;
        sol[0] = 0.0;
        p_system.current_local_solution->get(p_dof_indices, elemental_solution);
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            // Compute the old solution & its gradient.
            for (unsigned int l=0; l<elemental_solution.size(); l++)
            {
            	//nodeSol =  (*p_system.current_local_solution)(p_dof_indices[l]);
            	nodeSol = elemental_solution[l];
				sol[0] +=nodeSol * phi[l][qp];
                grad.add_scaled ( dphi[l][qp], nodeSol );
            }
        }

        //copy from RealGradient to std::vector
        int cnt = 0;
        for(auto&& v : grad_vector)
        {
            v = grad(cnt);
            ++cnt;
        }
        g_system.solution->add_vector(grad_vector, g_dof_indices);
        s_system.solution->add_vector(sol, s_dof_indices);
    }
    g_system.solution->close();
    s_system.solution->close();
    std::cout << " done" << std::endl;

}

} /* namespace BeatIt */
