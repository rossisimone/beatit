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
#include "libmesh/petsc_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"

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
#include "Elasticity/Materials/BenNeohookean.hpp"
#include "Elasticity/Materials/IsotropicMaterial.hpp"
#include "Elasticity/Materials/HolzapfelOgden.hpp"

#include "Util/Timer.hpp"

namespace libMesh
{
class GMVIO;
}

namespace BeatIt
{

typedef libMesh::TransientLinearImplicitSystem LinearSystem;

Elasticity::Elasticity(libMesh::EquationSystems& es, std::string system_name) :
		M_equationSystems(es), M_exporter(), M_outputFolder(), M_datafile(), M_linearSolver(), M_bch(), M_rhsFunction(), M_myName(
				system_name), M_JacIsAssembled(false), M_stabilize(false), M_currentNewtonIter(
				0)
{
	// TODO Auto-generated constructor stub

}

Elasticity::~Elasticity()
{
	// TODO Auto-generated destructor stub
}

void Elasticity::write_equation_system(const std::string& es)
{
	M_equationSystems.write(es, libMesh::WRITE);
}

void Elasticity::read_equation_system(const std::string& es)
{
	M_equationSystems.clear();
	M_equationSystems.read(es, libMesh::READ);
}

void Elasticity::init_exo_output(const std::string& output_filename)
{
	M_exporter->write_equation_systems(M_outputFolder + output_filename,
			M_equationSystems);
	M_exporter->append(true);
}

void Elasticity::save(const std::string& output_filename, int step)
{
	M_GMVexporter->write_equation_systems(
			M_outputFolder + output_filename + "." + std::to_string(step),
			M_equationSystems);
}

void Elasticity::save_exo(const std::string& output_filename, int step,
		double time)
{
	if (M_solverType == ElasticSolverType::Primal)
	{
		project_pressure();
	}
	std::cout << "* ELASTICITY: EXODUSII::Exporting  time " << time << " in: "
			<< M_outputFolder << " ... " << std::flush;
	M_exporter->write_timestep(M_outputFolder + output_filename,
			M_equationSystems, step, time);
	M_exporter->write_element_data(M_equationSystems);

//	M_exporter->write_element_data(M_equationSystems);
	std::cout << "done " << std::endl;
}

void Elasticity::deleteSystems()
{
	M_equationSystems.delete_system(M_myName);
}

void Elasticity::setup(const GetPot& data, std::string section)
{
	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	// Read Input File
	M_datafile = data;
	//Read output folder from datafile
	std::string output_folder = M_datafile(section + "/output_folder",
			"Output");
	M_outputFolder = "./" + output_folder + "/";
	// Create folder to save output
	BeatIt::createOutputFolder(M_equationSystems.get_mesh().comm(),
			M_outputFolder);

	std::cout << "Setting up system " << std::endl;
	setupSystem(section);
	setupPositionVector();
	std::cout << "Setting up parameters " << std::endl;
	setupParameters(section);
}

void Elasticity::setupSystem(std::string section)
{
	int dimension = M_equationSystems.get_mesh().mesh_dimension();
	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	// Starts by creating the equation systems
	// DISPLACEMENT PART
	std::cout << "* ELASTICITY: Creating new System for the ELASTICITY Solver"
			<< std::endl;
	LinearSystem& system = M_equationSystems.add_system<LinearSystem>(M_myName);

	std::string disp_name = "displacement";

	int ord = M_datafile(section + "/order", 1);
	std::cout << "* ELASTICITY: Reading displacement field - order: " << ord
			<< std::flush;
	auto order_it = BeatIt::libmesh_order_map.find(ord);
	std::cout << " ... " << std::flush;
	libMesh::Order order =
			(order_it != BeatIt::libmesh_order_map.end()) ?
					order_it->second :
					throw std::runtime_error("Order Explode!!!");
	std::cout << "order found!" << std::endl;

	std::string fam = M_datafile(section + "/fefamily", "lagrange");
	std::cout << "* ELASTICITY: Reading displacement field - family: " << fam
			<< std::flush;
	auto fefamily_it = BeatIt::libmesh_fefamily_map.find(fam);
	std::cout << " ... " << std::flush;
	libMesh::FEFamily fefamily =
			(fefamily_it != BeatIt::libmesh_fefamily_map.end()) ?
					fefamily_it->second :
					throw std::runtime_error("FeFamily Explode!!!");
	std::cout << "fefamily found!" << std::endl;
	std::cout << "* ELASTICITY: Setting up displacement field - order: " << ord
			<< ", fefamily = " << fam << std::flush;
	std::cout << "Adding variable_x!" << std::endl;
	system.add_variable(disp_name + "x", order, fefamily);
	std::cout << "Adding variable_y!" << std::endl;
	if (dimension > 1)
		system.add_variable(disp_name + "y", order, fefamily);
	if (dimension > 2)
		system.add_variable(disp_name + "z", order, fefamily);
	std::cout << " ... done " << std::endl;

	// PRESSURE SYSTEM
	std::string formulation = M_datafile(section + "/formulation", "primal");
	std::cout << "* ELASTICITY: Using a " << formulation << " formulation"
			<< std::endl;
	std::string pressure_name = "pressure";

	if (formulation == "mixed")
	{
		M_solverType = ElasticSolverType::Mixed;
		ord = M_datafile(section + "/p_order", 1);
		fam = M_datafile(section + "/p_fefamily", "lagrange");
		std::cout << "* ELASTICITY: Setting up pressure  field - order: " << ord
				<< ", fefamily = " << fam << std::flush;
		libMesh::Order p_order = BeatIt::libmesh_order_map.find(ord)->second;
		libMesh::FEFamily p_fefamily =
				BeatIt::libmesh_fefamily_map.find(fam)->second;
		system.add_variable(pressure_name, p_order, p_fefamily);
		std::cout << " ... done " << std::endl;

	}
	else
	{
		M_solverType = ElasticSolverType::Primal;
		std::cout << "* ELASTICITY: Setting up pressure  field ... "
				<< std::flush;
		LinearSystem& p_system = M_equationSystems.add_system<LinearSystem>(
				"Pressure_Projection");
		p_system.add_variable(pressure_name, libMesh::FIRST, libMesh::LAGRANGE);
		p_system.init();
		std::cout << " done " << std::endl;

	}

	std::cout << "* ELASTICITY: Adding residual and step" << std::flush;
	system.add_vector("residual");
	system.add_vector("step");
	std::cout << " done " << std::endl;

	std::cout << "* ELASTICITY: Reading BC ... " << std::flush;
	M_bch.readBC(M_datafile, section);
	M_bch.showMe();
	std::cout << " done " << std::endl;

	std::cout << "* ELASTICITY: Setup Homogenuous Dirichlet BC ... "
			<< std::flush;

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
				std::cout << "Setting mode FULL" << std::endl;
				for (unsigned int k = 0; k < dimension; ++k)
					variables.push_back(k);
				break;
			}
			case BCMode::Component:
			{
				std::cout << "Setting mode COMPONENT" << std::flush;
				auto component = bc_ptr->get_component();
				switch (component)
				{
				case BCComponent::X:
				{
					std::cout << " X" << std::endl;
					variables.push_back(0);
					break;
				}
				case BCComponent::Y:
				{
					std::cout << " Y" << std::endl;
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

			for (auto&& flag : bc_ptr->M_flag)
			{
				// Create a ZeroFunction to initialize dirichlet_bc
				dirichlet_boundary_ids.insert(flag);
			}
//
//            libMesh::ZeroFunction<> zf;
//            libMesh::DirichletBoundary dirichlet_bc(dirichlet_boundary_ids,
//                                                                                   variables,
//                                                                                   &zf );

			libMesh::DirichletBoundary dirichlet_bc(dirichlet_boundary_ids,
					variables, &(bc_ptr->get_function()));

			system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

		} // end if Dirichlet
	} // end for loop on BC
	  //std::cout << " done " << std::endl;
	  //system.get_dof_map().print_info();

	// Add position X
	system.init();

}

void Elasticity::setupPositionVector()
{
	LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);

	system.add_vector("X");
	// FILL VECTOR X
	std::cout << "* ELASTICITY: Filling position vector ... " << std::flush;
	auto it = M_equationSystems.get_mesh().active_nodes_begin();
	auto it_end = M_equationSystems.get_mesh().active_nodes_end();
	auto sys_num = system.number();
	for (; it != it_end; it++)
	{
		auto * node = *it;
		auto n_vars = node->n_vars(sys_num);
		for (auto v = 0; v != n_vars; ++v)
		{
			double value = 0.0;
			if (v < 3)
				value = (*node)(v);

			auto nc = node->n_comp(sys_num, v);
			for (auto c = 0; c != nc; ++c)
			{
				auto dof_id = node->dof_number(sys_num, v, c);
				system.get_vector("X").set(dof_id, value);
			}
		}
	}
	system.get_vector("X").close();
}

void Elasticity::setupParameters(std::string section)
{
	int dimension = M_equationSystems.get_mesh().mesh_dimension();
	LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);

	//Add vectors for anisotropic case
	std::cout << "Fibers? " << M_equationSystems.has_system("fibers")
			<< std::endl;
	std::cout << "Sheets? " << M_equationSystems.has_system("sheets")
			<< std::endl;
	std::cout << "XFibers? " << M_equationSystems.has_system("xfibers")
			<< std::endl;
	if (!M_equationSystems.has_system("fibers"))
	{
		ParameterSystem& fiber_system = M_equationSystems.add_system<
				ParameterSystem>("fibers");
		fiber_system.add_variable("fibersx", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		fiber_system.add_variable("fibersy", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		fiber_system.add_variable("fibersz", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		fiber_system.init();
		std::cout << "* ELASTICITY: Setup fibers: " << std::endl;
		std::string fibers_data = M_datafile(section + "/fibers",
				"1.0, 0.0, 0.0");
		std::cout << "*             f: " << fibers_data << std::endl;
		Util::project_function(fibers_data, fiber_system);
	}
	if (!M_equationSystems.has_system("sheets"))
	{
		ParameterSystem& sheets_system = M_equationSystems.add_system<
				ParameterSystem>("sheets");
		sheets_system.add_variable("sheetsx", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		sheets_system.add_variable("sheetsy", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		sheets_system.add_variable("sheetsz", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		sheets_system.init();
		std::string sheets_data = M_datafile(section + "/sheets",
				"0.0, 1.0, 0.0");
		std::cout << "*             s: " << sheets_data << std::endl;
		Util::project_function(sheets_data, sheets_system);
	}
	if (!M_equationSystems.has_system("xfibers"))
	{
		ParameterSystem& xfiber_system = M_equationSystems.add_system<
				ParameterSystem>("xfibers");
		xfiber_system.add_variable("xfibersx", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		xfiber_system.add_variable("xfibersy", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		xfiber_system.add_variable("xfibersz", libMesh::CONSTANT,
				libMesh::MONOMIAL);
		xfiber_system.init();
		std::string xfibers_data = M_datafile(section + "/xfibers",
				"0.0, 0.0, 1.0");
		std::cout << "*             n: " << xfibers_data << std::endl;
		Util::project_function(xfibers_data, xfiber_system);
	}
	ParameterSystem& dummy_system =
			M_equationSystems.add_system<ParameterSystem>("dumb");
	dummy_system.add_variable("dumb", libMesh::FIRST, libMesh::LAGRANGE);
	ParameterSystem& I4f_system = M_equationSystems.add_system<ParameterSystem>(
			"I4f");
	I4f_system.add_variable("I4f", libMesh::FIRST, libMesh::LAGRANGE);
	I4f_system.add_vector("patch_area");
	I4f_system.add_vector("F11");
	I4f_system.add_vector("F12");
	I4f_system.add_vector("F13");
	I4f_system.add_vector("F21");
	I4f_system.add_vector("F22");
	I4f_system.add_vector("F23");
	I4f_system.add_vector("F31");
	I4f_system.add_vector("F32");
	I4f_system.add_vector("F33");
	I4f_system.add_vector("nodal_fibersx");
	I4f_system.add_vector("nodal_fibersy");
	I4f_system.add_vector("nodal_fibersz");

	/// Init system
	std::cout << "* ELASTICITY: Init Systems ... " << std::endl;
	dummy_system.init();
	std::cout << "* ELASTICITY: I4f Systems ... " << std::endl;
	I4f_system.init();
	std::cout << "* ELASTICITY: print info Systems ... " << std::endl;
	M_equationSystems.print_info();
	/// Setting up BC

	std::cout << "* ELASTICITY: Setup body forces ... " << std::endl;

	std::string rhs = M_datafile(section + "/rhs", "no_rhs");
	M_rhsFunction.read(rhs);
	M_rhsFunction.showMe();

	std::cout << "* ELASTICITY: Setup linear solvers ... " << std::endl;

	M_linearSolver = libMesh::LinearSolver<libMesh::Number>::build(
			M_equationSystems.comm());
//    M_linearSolver->set_solver_type(libMesh::GMRES);
	M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
	M_linearSolver->init();
	M_projectionsLinearSolver = libMesh::LinearSolver<libMesh::Number>::build(
			M_equationSystems.comm());
	M_projectionsLinearSolver->set_solver_type(libMesh::GMRES);
	M_projectionsLinearSolver->set_preconditioner_type(libMesh::SOR_PRECOND);
	M_projectionsLinearSolver->init();

	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	// Setup Exporters
	M_exporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
	M_GMVexporter.reset(new Exporter(M_equationSystems.get_mesh()));
	M_VTKexporter.reset(new VTKExporter(M_equationSystems.get_mesh()));
	std::string mats = M_datafile(section + "/materials", "NO_MATERIAL");
	std::vector<std::string> materials;
	BeatIt::readList(mats, materials);
	std::string base_path = section + "/materials/";
	for (auto&& material : materials)
	{
		std::string path = base_path + material;
		unsigned int materialID = M_datafile(path + "/matID", 0);
		// Setup material
		M_materialMap[materialID].reset(
				Material::MaterialFactory::Create(material));

		M_materialMap[materialID]->setup(M_datafile, path, dimension);
	}

	M_newtonData.tol = M_datafile(section + "/newton/tolerance", 1e-9);
	M_newtonData.max_iter = M_datafile(section + "/newton/max_iter", 20);
//    if(M_solverType == ElasticSolverType::Primal)
//    {
//    	std::cout << "Solver Type = PRIMAL" << std::endl;
//    }
//    else     	std::cout << "Solver Type = MIXED" << std::endl;

	M_stabilize = M_datafile(section + "/stabilize", false);
	std::cout << "* ELASTICITY: Using stabilization: " << M_stabilize
			<< std::endl;
}

void Elasticity::setTime(double time)
{
	LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
	system.time = time;
}

void Elasticity::assemble_residual(double /* dt */,
		libMesh::NumericVector<libMesh::Number>* activation_ptr)
{
	std::cout << "* ELASTICITY: assembling ... " << std::endl;

	using libMesh::UniquePtr;

	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	const unsigned int max_dim = 3;
	// Get a reference to the LinearImplicitSystem we are solving
	LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
	double time = system.time;
	ParameterSystem& fiber_system =
			M_equationSystems.get_system<ParameterSystem>("fibers");
	ParameterSystem& sheets_system = M_equationSystems.get_system<
			ParameterSystem>("sheets");
	ParameterSystem& xfiber_system = M_equationSystems.get_system<
			ParameterSystem>("xfibers");

	system.get_vector("residual").zero();
	system.rhs->zero();
	system.matrix->zero();
	system.update();
	unsigned int ux_var = system.variable_number("displacementx");
	unsigned int uy_var, uz_var;

	if (dim > 1)
		uy_var = system.variable_number("displacementy");
	if (dim > 2)
		uz_var = system.variable_number("displacementz");

	const libMesh::DofMap & dof_map = system.get_dof_map();
	libMesh::FEType fe_disp = dof_map.variable_type(ux_var);
	const libMesh::DofMap & dof_map_fibers = fiber_system.get_dof_map();

	UniquePtr<libMesh::FEBase> fe_u(libMesh::FEBase::build(dim, fe_disp));
	auto order = fe_u->get_order();

	libMesh::QGauss qrule_1(dim, order);

	fe_u->attach_quadrature_rule(&qrule_1);

	const std::vector<libMesh::Real> & JxW_u = fe_u->get_JxW();
	const std::vector<libMesh::Point> & q_point_u = fe_u->get_xyz();

	// The element shape functions evaluated at the quadrature points.
	const std::vector<std::vector<libMesh::Real> > & phi_u = fe_u->get_phi();

	const std::vector<std::vector<libMesh::RealGradient> > & dphi_u =
			fe_u->get_dphi();

	libMesh::DenseVector<libMesh::Number> Fe;
	libMesh::DenseMatrix<libMesh::Number> Ke;

	std::vector<libMesh::dof_id_type> dof_indices;
	std::vector<libMesh::dof_id_type> dof_indices_ux;
	std::vector<libMesh::dof_id_type> dof_indices_uy;
	std::vector<libMesh::dof_id_type> dof_indices_uz;

	// Grad U
	std::vector<double> solution_k;
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
	// Stress
	libMesh::TensorValue<libMesh::Number> S;
	// Grad W (test function)
	libMesh::TensorValue<libMesh::Number> dW;

	double body_force[3];
	const std::vector<libMesh::Point> & q_point = fe_u->get_xyz();

	double gamma_f;
	double gamma_s;
	double gamma_n;
	std::vector<libMesh::dof_id_type> dof_indices_fibers;

	std::vector<double> gamma_f_k;
	libMesh::TensorValue<libMesh::Number> FA;

	// On the boundary
	libMesh::UniquePtr<libMesh::FEBase> fe_face(
			libMesh::FEBase::build(dim, fe_disp));
	libMesh::QGauss qface(dim - 1, libMesh::FIRST);
	fe_face->attach_quadrature_rule(&qface);

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	double rho;
	libMesh::RealGradient f0;
	libMesh::RealGradient s0;
	libMesh::RealGradient n0;

//	    std::cout << "* ELASTICITY:loop starts ... " << std::endl;

	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		rho = M_materialMap[0]->M_density;
		auto elID = elem->id();

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

		// get uk
		solution_k.resize(n_dofs);
//  	  for(auto && di : dof_indices)  std::cout << "dof id: " << di << std::endl;

		system.current_local_solution->get(dof_indices, solution_k);

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
			const unsigned int n_phi = phi_u.size();
			for (unsigned int l = 0; l < n_phi; ++l)
			{
				for (int idim = 0; idim < dim; idim++)
				{
					for (int jdim = 0; jdim < dim; jdim++)
					{
						dUk(idim, jdim) += dphi_u[l][qp](jdim)
								* solution_k[l + idim * n_phi];
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
						FA(idim, jdim) += gamma_f * f0(idim) * f0(jdim)
								+ gamma_s * s0(idim) * s0(jdim)
								+ gamma_n * n0(idim) * n0(jdim);
					}
				}
			}

			M_materialMap[0]->M_f0 = f0;
			M_materialMap[0]->M_s0 = s0;

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
					Fe(n + jdim * n_ux_dofs) += JxW_u[qp] * rho
							* body_force[jdim] * phi_u[n][qp];
				}
			}

			// Matrix
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

							M_materialMap[0]->evaluateJacobian(dU, 0.0);
							S = M_materialMap[0]->M_total_jacobian;
							auto int_1 = n + jdim * n_ux_dofs;
							auto int_2 = m + idim * n_ux_dofs;
//	                            std::cout << int_1 << ", " << int_2 << ", SdW: " << S.contract(dW) << std::endl;

							Ke(n + jdim * n_ux_dofs, m + idim * n_ux_dofs) +=
									S.contract(dW);
						}
					}
				}
			}
		}

		apply_BC(elem, Ke, Fe, fe_face, qface, mesh, n_ux_dofs, nullptr, 0.0,
				time);
		dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
//      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

		system.matrix->add_matrix(Ke, dof_indices);
		system.rhs->add_vector(Fe, dof_indices);
	}
	system.matrix->close();
	system.rhs->close();

}

void Elasticity::apply_BC(const libMesh::Elem*& elem,
		libMesh::DenseMatrix<libMesh::Number>& Ke,
		libMesh::DenseVector<libMesh::Number>& Fe,
		libMesh::UniquePtr<libMesh::FEBase>& fe_face, libMesh::QGauss& qface,
		const libMesh::MeshBase& mesh, int n_ux_dofs, MaterialPtr /* mat */,
		double /* dt */, double time, std::vector<double>* solk)
{
	const unsigned int dim = mesh.mesh_dimension();
	for (unsigned int side = 0; side < elem->n_sides(); side++)
	{
		if (elem->neighbor_ptr(side) == libmesh_nullptr)
		{
			const unsigned int boundary_id = mesh.boundary_info->boundary_id(
					elem, side);

			auto bc = M_bch.get_bc(boundary_id);

			if (bc)
			{
				const std::vector<libMesh::Real> & JxW_face =
						fe_face->get_JxW();
				const std::vector<std::vector<libMesh::Real> > & phi_face =
						fe_face->get_phi();
				int n_phi = phi_face.size();

				const std::vector<libMesh::Point> & qface_point =
						fe_face->get_xyz();
				const std::vector<std::vector<libMesh::RealGradient> > & dphi_face =
						fe_face->get_dphi();
				const std::vector<libMesh::Point>& normals =
						fe_face->get_normals();
				fe_face->reinit(elem, side);

				auto bc_type = bc->get_type();
				switch (bc_type)
				{
				case BCType::Neumann:
				{

					auto mode = bc->get_mode();
					for (unsigned int qp = 0; qp < qface.n_points(); qp++)
					{
						const double xq = qface_point[qp](0);
						const double yq = qface_point[qp](1);
						const double zq = qface_point[qp](2);
						if (BCMode::Full == mode)
						{
							for (int idim = 0; idim < dim; ++idim)
							{
								const double traction = bc->get_function()(time,
										xq, yq, zq, idim);
								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									Fe(i + idim * n_ux_dofs) += JxW_face[qp]
											* traction * phi_face[i][qp];
								}
							}
						}
						else if (BCMode::Component == mode)
						{
							int idim = 0;
							if (BCComponent::Y == bc->get_component())
								idim = 1;
							if (BCComponent::Z == bc->get_component())
								idim = 2;
							const double traction = bc->get_function()(time, xq,
									yq, zq, 0);
							for (unsigned int i = 0; i < n_ux_dofs; i++)
							{
								Fe(i + idim * n_ux_dofs) += JxW_face[qp]
										* traction * phi_face[i][qp];
							}
						}
						else if (BCMode::Normal == mode)
						{
							const double traction = bc->get_function()(time, xq,
									yq, zq, 0);
							for (int idim = 0; idim < dim; ++idim)
							{
								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									Fe(i + idim * n_ux_dofs) += JxW_face[qp]
											* traction * normals[qp](idim)
											* phi_face[i][qp];
								}
							}
						}

					}

					break;
				}
				case BCType::Robin:
				{
					auto& solution_k = *solk;
					auto mode = bc->get_mode();
					double beta = 1e6;

					for (unsigned int qp = 0; qp < qface.n_points(); qp++)
					{
						double uk = 0.0;
						// The location on the boundary of the current
						// face quadrature point.
						const double xq = qface_point[qp](0);
						const double yq = qface_point[qp](1);
						const double zq = qface_point[qp](2);

						const unsigned int n_phi = phi_face.size();

						if (BCMode::Component == mode)
						{
							int idim = 0;
							if (BCComponent::Y == bc->get_component())
								idim = 1;
							if (BCComponent::Z == bc->get_component())
								idim = 2;

							for (unsigned int l = 0; l < n_phi; ++l)
							{
								uk += phi_face[l][qp]
										* solution_k[l + idim * n_phi];
							} // end grad Uk on qp

							const double value = bc->get_function()(time, xq,
									yq, zq, 0);
							// Add the constant
							beta = bc->get_function()(time, xq,
									yq, zq, 1);
							for (unsigned int i = 0; i < n_ux_dofs; i++)
							{
								Fe(i + idim * n_ux_dofs) += JxW_face[qp] *
										* (value - beta * uk) * phi_face[i][qp];
							}

							for (unsigned int i = 0; i < n_ux_dofs; i++)
							{
								for (unsigned int j = 0; j < n_ux_dofs; j++)
								{
									Ke(i + idim * n_ux_dofs,
											j + idim * n_ux_dofs) +=
											JxW_face[qp] * beta
													* phi_face[i][qp]
													* phi_face[j][qp];
								}
							}

						}
						else if (BCMode::Full == mode)
						{
							beta = bc->get_function()(time, xq,
									yq, zq, dim+1);
							for (int idim = 0; idim < dim; ++idim)
							{

								for (unsigned int l = 0; l < n_phi; ++l)
								{
									uk += phi_face[l][qp]
											* solution_k[l + idim * n_phi];
								} // end grad Uk on qp

								const double value = bc->get_function()(time,
										xq, yq, zq, idim);
								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									Fe(i + idim * n_ux_dofs) += JxW_face[qp]
											* (value - beta * uk)
											* phi_face[i][qp];
								}

								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									for (unsigned int j = 0; j < n_ux_dofs; j++)
									{
										Ke(i + idim * n_ux_dofs,
												j + idim * n_ux_dofs) +=
												JxW_face[qp] * beta
														* phi_face[i][qp]
														* phi_face[j][qp];
									}
								}
							}
						}
						else if (BCMode::Normal == mode)
						{
							beta = bc->get_function()(time, xq,
									yq, zq, 1);
							// Defines a region on which to apply the BC
							const double value = bc->get_function()(time,
									xq, yq, zq, 0);
							double un = 0.0;
							double u1 = 0.0;
							double u2 = 0.0;
							double u3 = 0.0;
							for (unsigned int l = 0; l < n_phi; ++l)
							{
								u1 += phi_face[l][qp] * solution_k[l];
								if (dim > 1)
									u2 += phi_face[l][qp]
											* solution_k[l + n_phi];
								if (dim > 2)
									u3 += phi_face[l][qp]
											* solution_k[l + 2 * n_phi];
							} // end grad Uk on qp

							un = u1 * normals[qp](0);
							if (dim > 1)
								un += u2 * normals[qp](1);
							if (dim > 2)
								un += u3 * normals[qp](2);

							for (int idim = 0; idim < dim; ++idim)
							{
								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									//Fe(i+idim*n_ux_dofs) += JxW_face[qp] * beta * (value - ukn ) * phi_face[i][qp] * normals[qp](idim);
									//Fe(i+idim*n_ux_dofs) += JxW_face[qp] * beta * (value - ukn ) * phi_face[i][qp] * normals[qp](idim);
									//Fe(i+idim*n_ux_dofs) -= JxW_face[qp] * u1 * normals[qp](0) * phi_face[i][qp] * normals[qp](idim);
									//Fe(i+idim*n_ux_dofs) -= JxW_face[qp] * u2 * normals[qp](1) * phi_face[i][qp] * normals[qp](idim);
									Fe(i + idim * n_ux_dofs) -= JxW_face[qp]
											* (un - beta * value)
											* phi_face[i][qp]
											* normals[qp](idim);
								}
							}

							for (int idim = 0; idim < dim; ++idim)
							{
								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									for (int jdim = 0; jdim < dim; ++jdim)
									{
										for (unsigned int j = 0;
												j < n_ux_dofs; j++)
										{
											//Ke(i+idim*n_ux_dofs,j+jdim*n_ux_dofs) += JxW_face[qp] * beta * phi_face[i][qp] * normals[qp](idim) * phi_face[j][qp] * normals[qp](jdim);
											Ke(i + idim * n_ux_dofs,
													j + jdim * n_ux_dofs) +=
													JxW_face[qp] * beta
															* phi_face[i][qp]
															* normals[qp](
																	idim)
															* phi_face[j][qp]
															* normals[qp](
																	jdim);
										}
									}
								}
							}
						}

					}

					break;
				}
				case BCType::Penalty:
				{
					auto& solution_k = *solk;
					auto mode = bc->get_mode();
					double beta = 1e6;

					for (unsigned int qp = 0; qp < qface.n_points(); qp++)
					{
						double uk = 0.0;
						// The location on the boundary of the current
						// face quadrature point.
						const double xq = qface_point[qp](0);
						const double yq = qface_point[qp](1);
						const double zq = qface_point[qp](2);

						const unsigned int n_phi = phi_face.size();

						if (BCMode::Component == mode)
						{
							int idim = 0;
							if (BCComponent::Y == bc->get_component())
								idim = 1;
							if (BCComponent::Z == bc->get_component())
								idim = 2;

							for (unsigned int l = 0; l < n_phi; ++l)
							{
								uk += phi_face[l][qp]
										* solution_k[l + idim * n_phi];
							} // end grad Uk on qp

							const double value = bc->get_function()(time, xq,
									yq, zq, 0);
							for (unsigned int i = 0; i < n_ux_dofs; i++)
							{
								Fe(i + idim * n_ux_dofs) += JxW_face[qp] * beta
										* (value - uk) * phi_face[i][qp];
							}

							for (unsigned int i = 0; i < n_ux_dofs; i++)
							{
								for (unsigned int j = 0; j < n_ux_dofs; j++)
								{
									Ke(i + idim * n_ux_dofs,
											j + idim * n_ux_dofs) +=
											JxW_face[qp] * beta
													* phi_face[i][qp]
													* phi_face[j][qp];
								}
							}

						}
						else if (BCMode::Full == mode)
						{
							for (int idim = 0; idim < dim; ++idim)
							{

								for (unsigned int l = 0; l < n_phi; ++l)
								{
									uk += phi_face[l][qp]
											* solution_k[l + idim * n_phi];
								} // end grad Uk on qp

								const double value = bc->get_function()(time,
										xq, yq, zq, idim);
								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									Fe(i + idim * n_ux_dofs) += JxW_face[qp]
											* beta * (value - uk)
											* phi_face[i][qp];
								}

								for (unsigned int i = 0; i < n_ux_dofs; i++)
								{
									for (unsigned int j = 0; j < n_ux_dofs; j++)
									{
										Ke(i + idim * n_ux_dofs,
												j + idim * n_ux_dofs) +=
												JxW_face[qp] * beta
														* phi_face[i][qp]
														* phi_face[j][qp];
									}
								}
							}
						}
						else if (BCMode::Normal == mode)
						{
//							beta = 1e3;
							beta = 1e6;
							// Defines a region on which to apply the BC
							int function_size = bc->get_function().size();
							std::vector<double> center(function_size);
							double radius = 0;
							double distance = 0;
							if (function_size >= 2)
								radius = bc->get_function()(time, xq, yq, zq,
										1);
							for (int fs = 2; fs < center.size(); ++fs)
							{
								center[fs] = bc->get_function()(time, xq, yq,
										zq, fs);
							}
							if (function_size >= 3)
								distance += (center[2] - xq) * (center[2] - xq);
							if (function_size >= 4)
								distance += (center[3] - yq) * (center[3] - yq);
							if (function_size >= 5)
								distance += (center[4] - zq) * (center[4] - zq);

							if (distance <= radius * radius
									|| function_size <= 1)
							{

								const double value = bc->get_function()(time,
										xq, yq, zq, 0);
								double un = 0.0;
								double u1 = 0.0;
								double u2 = 0.0;
								double u3 = 0.0;
								for (unsigned int l = 0; l < n_phi; ++l)
								{
									u1 += phi_face[l][qp] * solution_k[l];
									if (dim > 1)
										u2 += phi_face[l][qp]
												* solution_k[l + n_phi];
									if (dim > 2)
										u3 += phi_face[l][qp]
												* solution_k[l + 2 * n_phi];
								} // end grad Uk on qp

								un = u1 * normals[qp](0);
								if (dim > 1)
									un += u2 * normals[qp](1);
								if (dim > 2)
									un += u3 * normals[qp](2);

								for (int idim = 0; idim < dim; ++idim)
								{
									for (unsigned int i = 0; i < n_ux_dofs; i++)
									{
										//Fe(i+idim*n_ux_dofs) += JxW_face[qp] * beta * (value - ukn ) * phi_face[i][qp] * normals[qp](idim);
										//Fe(i+idim*n_ux_dofs) += JxW_face[qp] * beta * (value - ukn ) * phi_face[i][qp] * normals[qp](idim);
										//Fe(i+idim*n_ux_dofs) -= JxW_face[qp] * u1 * normals[qp](0) * phi_face[i][qp] * normals[qp](idim);
										//Fe(i+idim*n_ux_dofs) -= JxW_face[qp] * u2 * normals[qp](1) * phi_face[i][qp] * normals[qp](idim);
										Fe(i + idim * n_ux_dofs) -= JxW_face[qp]
												* beta * (un - value)
												* phi_face[i][qp]
												* normals[qp](idim);
									}
								}

								for (int idim = 0; idim < dim; ++idim)
								{
									for (unsigned int i = 0; i < n_ux_dofs; i++)
									{
										for (int jdim = 0; jdim < dim; ++jdim)
										{
											for (unsigned int j = 0;
													j < n_ux_dofs; j++)
											{
												//Ke(i+idim*n_ux_dofs,j+jdim*n_ux_dofs) += JxW_face[qp] * beta * phi_face[i][qp] * normals[qp](idim) * phi_face[j][qp] * normals[qp](jdim);
												Ke(i + idim * n_ux_dofs,
														j + jdim * n_ux_dofs) +=
														JxW_face[qp] * beta
																* phi_face[i][qp]
																* normals[qp](
																		idim)
																* phi_face[j][qp]
																* normals[qp](
																		jdim);
											}
										}
									}
								}
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

void Elasticity::solve_system()
{
	LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);

	double tol = 1e-12;
	double max_iter = 2000;
//    system.matrix->print(std::cout);
//    std::cout << "Matrix!\n" << std::endl;
//    system.matrix->print(std::cout);
//    std::cout << "\nRHS!" << std::endl;
//    system.rhs->print(std::cout);

	// set near null space
	{
//        auto petscMatrixPtr = dynamic_cast<libMesh::PetscMatrix<libMesh::Number> *>(system.matrix);
//        auto petscVecPtr = dynamic_cast<libMesh::PetscVector<libMesh::Number> *>(system.request_vector("X"));
//        MatNullSpace matnull;
//        Vec          vec_coords;
//        PetscScalar  *c;

//        VecCreate(MPI_COMM_WORLD,&vec_coords);
//        VecSetBlockSize(vec_coords,3);
//        VecSetSizes(vec_coords,m,PETSC_DECIDE);
//        VecSetUp(vec_coords);
//        VecGetArray(vec_coords,&c);
//        for (i=0; i<m; i++) c[i] = coords[i]; /* Copy since Scalar type might be Complex */
//        VecRestoreArray(vec_coords,&c);
//        MatNullSpaceCreateRigidBody(petscVecPtr->vec(),&matnull);
//        MatSetNearNullSpace(petscMatrixPtr->mat(),matnull);
//        MatNullSpaceDestroy(&matnull);

	}

	std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);
	rval = M_linearSolver->solve(*system.matrix, nullptr,
			system.get_vector("step"), *system.rhs, tol, max_iter);

	std::cout << "num it: " << rval.first << ", res: " << rval.second
			<< std::endl;
//    system.matrix->print(std::cout);
//    system.rhs->print(std::cout);
//    system.get_vector("step").print(std::cout);

}

void Elasticity::project_pressure()
{
	std::cout << "* ELASTICITY: projecting pressure ... " << std::endl;

	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	LinearSystem & system = M_equationSystems.get_system<LinearSystem>(
			M_myName);
	LinearSystem & system_p = M_equationSystems.get_system<LinearSystem>(
			"Pressure_Projection");
	system.update();
	unsigned int ux_var = system.variable_number("displacementx");
	unsigned int uy_var, uz_var, p_var;

	if (dim > 1)
		uy_var = system.variable_number("displacementy");
	if (dim > 2)
		uz_var = system.variable_number("displacementz");
	p_var = system_p.variable_number("pressure");

	const libMesh::DofMap & dof_map = system.get_dof_map();
	libMesh::FEType fe_type = dof_map.variable_type(0);
	libMesh::UniquePtr<libMesh::FEBase> fe(
			libMesh::FEBase::build(dim, fe_type));
	libMesh::QGauss qrule(dim, fe->get_order());
	fe->attach_quadrature_rule(&qrule);

	const libMesh::DofMap & dof_map_p = system_p.get_dof_map();
	libMesh::FEType fe_type_p = dof_map_p.variable_type(0);
	libMesh::UniquePtr<libMesh::FEBase> fe_p(
			libMesh::FEBase::build(dim, fe_type_p));
	libMesh::QGauss qrule_p(dim, libMesh::SECOND);
	fe_p->attach_quadrature_rule(&qrule_p);

	const std::vector<libMesh::Real> & JxW = fe_p->get_JxW();
	const std::vector<std::vector<libMesh::RealGradient> > & dphi =
			fe_p->get_dphi();
	const std::vector<std::vector<libMesh::Real> > & phi = fe_p->get_phi();
	const std::vector<std::vector<libMesh::RealGradient> > & dphi_u =
			fe->get_dphi();
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

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

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
						dUk(idim, jdim) += dphi_u[l][qp](jdim)
								* solution_k[l + idim * n_phi];
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

	rval = M_projectionsLinearSolver->solve(*system_p.matrix, nullptr,
			*system_p.solution, *system_p.rhs, tol, max_iter);
}

void Elasticity::newton(double dt,
		libMesh::NumericVector<libMesh::Number>* activation_ptr)
{
	Timer timer;
	timer.start();
	LinearSystem& system = M_equationSystems.get_system<LinearSystem>(M_myName);
	update_displacements(dt);

	M_currentNewtonIter = 0;
	assemble_residual(dt, activation_ptr);
	auto res_norm = system.rhs->linfty_norm();
	// assume we ask for the same absolute and relative tolerances
	// tol = atol + rtol * res_norm
	double tol = M_newtonData.tol * (1 + res_norm);
	int max_iter = M_newtonData.max_iter;
	int iter = 0;

	std::cout << "* ELASTICITY: Performing Newton solve:  max iterations: "
			<< max_iter << ", target tolerance: " << tol << std::endl;
	std::cout << "\t\t\t  iter: " << iter << ", residual: " << res_norm
			<< std::endl;

	double linear_tol = 1e-12;
	double linear_max_iter = 2000;

	while (res_norm > tol && M_currentNewtonIter < max_iter)
	{
		M_currentNewtonIter++;

		//solve_system();

		std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);
		rval = M_linearSolver->solve(*system.matrix, nullptr,
				system.get_vector("step"), *system.rhs, linear_tol,
				linear_max_iter);
//          std::cout << "RHS!\n" << std::endl;
//          system.rhs->print(std::cout);
//          std::cout << "Step!\n" << std::endl;
// 	     system.get_vector("step").print(std::cout);
//
//         std::cout << "Jac!\n" << std::endl;
// 	    system.matrix->print(std::cout);
//          std::cout << "Solution!\n" << std::endl;
//          system.solution->print(std::cout);

		if (rval.first == 0)
		{
			std::cout
					<< "* ELASTICITY: WARNING: linear solver diverged, try adding -ksp_divtol 1e10 "
					<< std::endl;
			std::cout
					<< "* ELASTICITY: I'm stopping the newton iterations, hopefully you are close enough "
					<< std::endl;
			break;
		}

		(*system.solution) += system.get_vector("step");

		update_displacements(dt);
		assemble_residual(dt, activation_ptr);
		res_norm = system.rhs->linfty_norm();
		std::cout << "\t\t\t  iter: " << M_currentNewtonIter << ", residual: "
				<< res_norm << std::endl;
	}

	std::cout << "* ELASTICITY: Newton solve completed in "
			<< M_currentNewtonIter << " iterations. Final residual: "
			<< res_norm << std::endl;
	timer.stop();
	timer.print(std::cout);
}

void Elasticity::evaluate_nodal_I4f()
{
	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	ParameterSystem& fiber_system =
			M_equationSystems.get_system<ParameterSystem>("fibers");
	ParameterSystem& I4f_system = M_equationSystems.get_system<ParameterSystem>(
			"I4f");
	LinearSystem & system = M_equationSystems.get_system<LinearSystem>(
			M_myName);

	auto& patch_area = I4f_system.get_vector("patch_area");
	auto& F11 = I4f_system.get_vector("F11");
	auto& F12 = I4f_system.get_vector("F12");
	auto& F13 = I4f_system.get_vector("F13");
	auto& F21 = I4f_system.get_vector("F21");
	auto& F22 = I4f_system.get_vector("F22");
	auto& F23 = I4f_system.get_vector("F23");
	auto& F31 = I4f_system.get_vector("F31");
	auto& F32 = I4f_system.get_vector("F32");
	auto& F33 = I4f_system.get_vector("F33");
	auto& nodal_fibersx = I4f_system.get_vector("nodal_fibersx");
	auto& nodal_fibersy = I4f_system.get_vector("nodal_fibersy");
	auto& nodal_fibersz = I4f_system.get_vector("nodal_fibersz");
	patch_area.zero();
	F11.zero();
	F12.zero();
	F13.zero();
	F21.zero();
	F22.zero();
	F23.zero();
	F31.zero();
	F32.zero();
	F33.zero();
	nodal_fibersx.zero();
	nodal_fibersy.zero();
	nodal_fibersz.zero();
	I4f_system.solution->zero();

	const libMesh::DofMap & dof_map_I4f = I4f_system.get_dof_map();
	const libMesh::DofMap & dof_map = system.get_dof_map();
	const libMesh::DofMap & dof_map_fibers = fiber_system.get_dof_map();

	libMesh::FEType fe_type = dof_map_I4f.variable_type(0);
	libMesh::UniquePtr<libMesh::FEBase> fe(
			libMesh::FEBase::build(dim, fe_type));

	// A 5th order Gauss quadrature rule for numerical integration.
	libMesh::QGauss qrule(dim, libMesh::FIRST);

	fe->attach_quadrature_rule(&qrule);
	const std::vector<libMesh::Real> & JxW = fe->get_JxW();
	const std::vector<std::vector<libMesh::Real> > & phi = fe->get_phi();
	const std::vector<std::vector<libMesh::RealGradient> > & dphi =
			fe->get_dphi();

	std::vector<libMesh::dof_id_type> dof_indices;
	std::vector<libMesh::dof_id_type> dof_indices_I4f;
	std::vector<libMesh::dof_id_type> dof_indices_fibers;

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	std::vector<double> disp;
	libMesh::VectorValue<libMesh::Number> f0;
	libMesh::VectorValue<libMesh::Number> f;
	libMesh::VectorValue<libMesh::Number> du1;
	libMesh::VectorValue<libMesh::Number> du2;
	libMesh::VectorValue<libMesh::Number> du3;
	// Grad U
	libMesh::TensorValue<libMesh::Number> dUk;

	double I4f;

	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		const unsigned int elem_id = elem->id();
		dof_map.dof_indices(elem, dof_indices);
		dof_map_I4f.dof_indices(elem, dof_indices_I4f);
		dof_map_fibers.dof_indices(elem, dof_indices_fibers);

		fe->reinit(elem);

		f0(0) = (*fiber_system.solution)(dof_indices_fibers[0]);
		f0(1) = (*fiber_system.solution)(dof_indices_fibers[1]);
		f0(2) = (*fiber_system.solution)(dof_indices_fibers[2]);

		disp.resize(dof_indices.size());
		system.current_local_solution->get(dof_indices, disp);

		const unsigned int n_phi = phi.size();
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
		{
			dUk *= 0.0;
			for (unsigned int l = 0; l < phi.size(); l++)
			{
				for (int idim = 0; idim < dim; idim++)
				{
					for (int jdim = 0; jdim < dim; jdim++)
					{
						dUk(idim, jdim) += dphi[l][qp](jdim)
								* disp[l + idim * n_phi];
					}
				}
			}

		}

		double evol = elem->volume();
		for (auto && index : dof_indices_I4f)
		{
			nodal_fibersx.add(index, evol * f0(0));
			nodal_fibersy.add(index, evol * f0(1));
			nodal_fibersz.add(index, evol * f0(2));
			patch_area.add(index, evol);
			F11.add(index, evol * (1.0 + dUk(0, 0)));
			F12.add(index, evol * (dUk(0, 1)));
			F13.add(index, evol * (dUk(0, 2)));
			F21.add(index, evol * (dUk(1, 0)));
			F22.add(index, evol * (1.0 + dUk(1, 1)));
			F23.add(index, evol * (dUk(1, 2)));
			F31.add(index, evol * (dUk(2, 0)));
			F32.add(index, evol * (dUk(2, 1)));
			F33.add(index, evol * (1.0 + dUk(2, 2)));
		}
	}
	patch_area.close();
	F11.close();
	F12.close();
	F13.close();
	F21.close();
	F22.close();
	F23.close();
	F31.close();
	F32.close();
	F33.close();
	nodal_fibersx.close();
	nodal_fibersy.close();
	nodal_fibersz.close();

	F11 /= patch_area;
	F12 /= patch_area;
	F13 /= patch_area;
	F21 /= patch_area;
	F22 /= patch_area;
	F23 /= patch_area;
	F31 /= patch_area;
	F32 /= patch_area;
	F33 /= patch_area;
	nodal_fibersx /= patch_area;
	nodal_fibersy /= patch_area;
	nodal_fibersz /= patch_area;

	libMesh::TensorValue<libMesh::Number> nodalF;

	auto first_local_index = I4f_system.solution->first_local_index();
	auto last_local_index = I4f_system.solution->last_local_index();

	auto i = first_local_index;
	for (; i < last_local_index; i++)
	{
		nodalF(0, 0) = F11(i);
		nodalF(0, 1) = F12(i);
		nodalF(0, 2) = F13(i);
		nodalF(1, 0) = F21(i);
		nodalF(1, 1) = F22(i);
		nodalF(1, 2) = F23(i);
		nodalF(2, 0) = F31(i);
		nodalF(2, 1) = F32(i);
		nodalF(2, 2) = F33(i);
		f0(0) = nodal_fibersx(i);
		f0(1) = nodal_fibersy(i);
		f0(2) = nodal_fibersz(i);
//        std::cout << std::endl;
//        std::cout << std::endl;
//        nodalF.print(std::cout);
		double normf0 = f0.norm();
		f0 /= normf0;
//        std::cout << std::endl;
//        f0.print(std::cout);
//        std::cout << std::endl;
		f = nodalF * f0;
//        f.print(std::cout);
		I4f = f.contract(f);
//        std::cout << ", I4f = " << I4f << std::endl;

		I4f_system.solution->set(i, I4f);
	}

	I4f_system.solution->close();
	I4f_system.update();

}

void Elasticity::evaluate_L2_J_err()
{
	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	ParameterSystem& I4f_system = M_equationSystems.get_system<ParameterSystem>(
			"I4f");
	LinearSystem & system = M_equationSystems.get_system<LinearSystem>(
			M_myName);

	const libMesh::DofMap & dof_map_I4f = I4f_system.get_dof_map();
	const libMesh::DofMap & dof_map = system.get_dof_map();

	libMesh::FEType fe_type = dof_map_I4f.variable_type(0);
	libMesh::UniquePtr<libMesh::FEBase> fe(
			libMesh::FEBase::build(dim, fe_type));

	// A 5th order Gauss quadrature rule for numerical integration.
	libMesh::QGauss qrule(dim, libMesh::FIRST);

	fe->attach_quadrature_rule(&qrule);
	const std::vector<libMesh::Real> & JxW = fe->get_JxW();
	const std::vector<std::vector<libMesh::Real> > & phi = fe->get_phi();
	const std::vector<std::vector<libMesh::RealGradient> > & dphi =
			fe->get_dphi();

	std::vector<libMesh::dof_id_type> dof_indices;
	std::vector<libMesh::dof_id_type> dof_indices_I4f;

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	std::vector<double> disp;
	libMesh::TensorValue<libMesh::Number> Fk;

	double errJ2 = 0.0;

	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		const unsigned int elem_id = elem->id();
		dof_map.dof_indices(elem, dof_indices);
		dof_map_I4f.dof_indices(elem, dof_indices_I4f);

		fe->reinit(elem);

		disp.resize(dof_indices.size());
		system.current_local_solution->get(dof_indices, disp);

		const unsigned int n_phi = phi.size();
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
		{
			Fk *= 0.0;
			for (unsigned int l = 0; l < phi.size(); l++)
			{
				for (int idim = 0; idim < dim; idim++)
				{
					for (int jdim = 0; jdim < dim; jdim++)
					{
						Fk(idim, jdim) += dphi[l][qp](jdim)
								* disp[l + idim * n_phi];
					}
				}
			}
			Fk(0, 0) += 1.0;
			Fk(1, 1) += 1.0;
			Fk(2, 2) += 1.0;

			double J = Fk.det();
			//std::cout << "J: " << J << std::endl;

			errJ2 += JxW[qp] * (J - 1) * (J - 1);

		}
	}
	mesh.comm().sum(errJ2);

	double errJ = std::sqrt(errJ2);

	std::cout << "\n|| J - 1 || = " << errJ << std::endl;

}

} /* namespace BeatIt */
