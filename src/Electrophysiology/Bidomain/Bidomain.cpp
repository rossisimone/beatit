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

/**
 * \file Bidomain.cpp
 *
 * \class Bidomain
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 11, 2016
 *
 */

// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/type_tensor.h"
#include "PoissonSolver/Poisson.hpp"

// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
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
#include "libmesh/petsc_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

#include "libmesh/exodusII_io.h"
//#include "libmesh/exodusII_io_helper.h"
#include "libmesh/gmv_io.h"

#include "libmesh/perf_log.h"

#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"
#include "libmesh/discontinuity_measure.h"
#include <sys/stat.h>

#include "Electrophysiology/IonicModels/NashPanfilov.hpp"
#include "Electrophysiology/IonicModels/Grandi11.hpp"
#include "Electrophysiology/IonicModels/ORd.hpp"
#include "Electrophysiology/IonicModels/TP06.hpp"
#include "Electrophysiology/Bidomain/Bidomain.hpp"
#include "Util/SpiritFunction.hpp"

#include "libmesh/discontinuity_measure.h"
#include "libmesh/fe_interface.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "Electrophysiology/Pacing/PacingProtocolSpirit.hpp"

namespace BeatIt
{

// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
typedef libMesh::TransientLinearImplicitSystem BidomainSystem;
typedef libMesh::TransientExplicitSystem IonicModelSystem;
typedef libMesh::ExplicitSystem ParameterSystem;

Bidomain::Bidomain(libMesh::EquationSystems& es) :
		M_equationSystems(es), M_bidomainExporter(), M_bidomainExporterNames(), M_ionicModelExporter(), M_ionicModelExporterNames(), M_parametersExporter(), M_parametersExporterNames(), M_outputFolder(), M_datafile(), M_pacing(), M_linearSolver(), M_anisotropy(
				Anisotropy::Orthotropic), M_equationType(
				BidomainEquationType::ParabolicEllipticBidomain), M_artificialDiffusion(
				false), M_penalty(false), M_timeIntegratorType(
				DynamicTimeIntegratorType::Implicit)
{

}

Bidomain::~Bidomain()
{
}

void Bidomain::setup(GetPot& data, std::string section)
{
	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	evaluate_mesh_size();
	std::cout << "* BIDOMAIN: max mesh size:" << M_meshSize << std::endl;

	// Read Input File
	M_datafile = data;

	//Read output folder from datafile
	std::string output_folder = M_datafile(section + "/output_folder",
			"Output");
	M_outputFolder = "./" + output_folder + "/";

	std::cout << "* BIDOMAIN: creating pacing protocol" << std::endl;
	std::string pacing_type = M_datafile(section + "/pacing/type", "S1");
	M_pacing.reset(
			BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_type));
	M_pacing->setup(M_datafile, "bidomain/pacing");

	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	// Starts by creating the equation systems
	// 1) ADR
	std::cout
			<< "* BIDOMAIN: Creating new System for the hyperbolic bidomain equations"
			<< std::endl;
	BidomainSystem& bidomain_system = M_equationSystems.add_system
			< BidomainSystem > ("bidomain");
	// TO DO: Generalize to higher order
	bidomain_system.add_variable("Q", libMesh::FIRST);
	bidomain_system.add_variable("Ve", libMesh::FIRST);

	// Add 3 matrices
	bidomain_system.add_matrix("lumped_mass");
	bidomain_system.add_matrix("high_order_mass");
	// Add lumped mass matrix
	bidomain_system.add_vector("lumped_mass_vector");
	bidomain_system.add_matrix("mass");
	bidomain_system.add_matrix("stiffness");
	bidomain_system.add_vector("ionic_currents");
	bidomain_system.add_vector("old_solution");
	bidomain_system.add_vector("nullspace");
	M_bidomainExporterNames.insert("bidomain");

//    bidomain_system.matrix->init();
	bidomain_system.init();

	// WAVE
	BidomainSystem& wave_system = M_equationSystems.add_system < BidomainSystem
			> ("wave");
	wave_system.add_variable("V", libMesh::FIRST);
	wave_system.add_matrix("Ki");
	wave_system.add_vector("KiV");
	M_bidomainExporterNames.insert("wave");
	wave_system.init();

	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	// 2) ODEs
	std::cout << "* BIDOMAIN: Creating new System for the ionic model "
			<< std::endl;
	IonicModelSystem& ionic_model_system = M_equationSystems.add_system
			< IonicModelSystem > ("ionic_model");
	// Create Ionic Model
	std::string ionic_model = M_datafile(section + "/ionic_model",
			"NashPanfilov");
	M_ionicModelPtr.reset(
			BeatIt::IonicModel::IonicModelFactory::Create(ionic_model));
	M_ionicModelPtr->setup(M_datafile, section);
	int num_vars = M_ionicModelPtr->numVariables();
	// TO DO: Generalize to other conditions
	// We need to exclude the potential
	// therefore we loop up to num_vars-1
	for (int nv = 0; nv < num_vars - 1; ++nv)
	{
		std::string var_name = M_ionicModelPtr->variableName(nv);
		// For the time being we use P1 for the variables
		ionic_model_system.add_variable(&var_name[0], libMesh::FIRST);
	}
	ionic_model_system.init();
	// Add the applied current to this system
	IonicModelSystem& Iion_system = M_equationSystems.add_system
			< IonicModelSystem > ("iion");
	Iion_system.add_variable("iion", libMesh::FIRST);
	Iion_system.add_vector("diion");
	Iion_system.init();
	IonicModelSystem& istim_system = M_equationSystems.add_system
			< IonicModelSystem > ("istim");
	istim_system.add_variable("istim", libMesh::FIRST);
	istim_system.init();
	M_ionicModelExporterNames.insert("ionic_model");
	M_ionicModelExporterNames.insert("iion");
	M_ionicModelExporterNames.insert("istim");

	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	// Distributed Parameters
	std::cout << "* BIDOMAIN: Creating parameters spaces " << std::endl;
	ParameterSystem& activation_times_system = M_equationSystems.add_system
			< ParameterSystem > ("activation_times");
	activation_times_system.add_variable("activation_times", libMesh::FIRST);
	ParameterSystem& fiber_system = M_equationSystems.add_system
			< ParameterSystem > ("fibers");
	fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
	fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
	fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
	ParameterSystem& sheets_system = M_equationSystems.add_system
			< ParameterSystem > ("sheets");
	sheets_system.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
	sheets_system.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
	sheets_system.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
	ParameterSystem& xfiber_system = M_equationSystems.add_system
			< ParameterSystem > ("xfibers");
	xfiber_system.add_variable("xfibersx", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	xfiber_system.add_variable("xfibersy", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	xfiber_system.add_variable("xfibersz", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	ParameterSystem& intra_conductivity_system = M_equationSystems.add_system
			< ParameterSystem > ("intra_conductivity");
	intra_conductivity_system.add_variable("Dffi", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	intra_conductivity_system.add_variable("Dssi", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	intra_conductivity_system.add_variable("Dnni", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	ParameterSystem& extra_conductivity_system = M_equationSystems.add_system
			< ParameterSystem > ("extra_conductivity");
	extra_conductivity_system.add_variable("Dffe", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	extra_conductivity_system.add_variable("Dsse", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	extra_conductivity_system.add_variable("Dnne", libMesh::CONSTANT,
			libMesh::MONOMIAL);
	ParameterSystem& procID_system = M_equationSystems.add_system
			< ParameterSystem > ("ProcID");
	procID_system.add_variable("ProcID", libMesh::CONSTANT, libMesh::MONOMIAL);
	M_parametersExporterNames.insert("activation_times");
	M_parametersExporterNames.insert("fibers");
	M_parametersExporterNames.insert("sheets");
	M_parametersExporterNames.insert("xfibers");
	M_parametersExporterNames.insert("intra_conductivity");
	M_parametersExporterNames.insert("extra_conductivity");
	M_parametersExporterNames.insert("ProcID");
	std::cout << "* BIDOMAIN: Initializing equation systems " << std::endl;
	// Initializing
	activation_times_system.init();
	fiber_system.init();
	sheets_system.init();
	xfiber_system.init();
	intra_conductivity_system.init();
	extra_conductivity_system.init();
	procID_system.init();
//    M_equationSystems.init();
	//M_equationSystems.print_info();

	// Conductivity Tensor in local coordinates
	// Fiber direction
	// Default Value = 1.3342  kOhm^-1 cm^-1
	double Dffi = M_datafile(section + "/Dffi", 1.3342);
	// Sheet direction
	// Default Value = 1.3342  kOhm^-1 cm^-1
	double Dssi = M_datafile(section + "/Dssi", 0.17606);
	// Cross fiber direction
	// Default Value = 1.3342  kOhm^-1 cm^-1
	double Dnni = M_datafile(section + "/Dnni", 0.17606);
	// Default Value = 1.3342  kOhm^-1 cm^-1
	double Dffe = M_datafile(section + "/Dffe", 1.3342);
	// Sheet direction
	// Default Value = 1.3342  kOhm^-1 cm^-1
	double Dsse = M_datafile(section + "/Dsse", 0.17606);
	// Cross fiber direction
	// Default Value = 1.3342  kOhm^-1 cm^-1
	double Dnne = M_datafile(section + "/Dnne", 0.17606);

	// Surface to volume ratio
	double Chi = M_datafile(section + "/Chi", 1400.0);
	// Intracellular and extracellular relaxation times
	double tau_i = M_datafile(section + "/tau_i", 0.0);
	double tau_e = M_datafile(section + "/tau_e", 0.0);
	M_equationSystems.parameters.set<double>("Chi") = Chi;
	M_equationSystems.parameters.set<double>("tau_i") = tau_i;
	M_equationSystems.parameters.set<double>("tau_e") = tau_e;

	std::string anisotropy = M_datafile(section + "/anisotropy", "orhotropic");
	std::map<std::string, Anisotropy> aniso_map;
	aniso_map["orthotropic"] = Anisotropy::Orthotropic;
	aniso_map["isotropic"] = Anisotropy::Isotropic;
	aniso_map["transverse"] = Anisotropy::TransverselyIsotropic;
	std::cout << "* BIDOMAIN: Parameters: " << std::endl;
	std::cout << "              Chi = " << Chi << std::endl;
	std::cout << "              Dffi = " << Dffi << std::endl;
	std::cout << "              Dssi = " << Dssi << std::endl;
	std::cout << "              Dnni = " << Dnni << std::endl;
	std::cout << "              Dffe = " << Dffe << std::endl;
	std::cout << "              Dsse = " << Dsse << std::endl;
	std::cout << "              Dnne = " << Dnne << std::endl;
	std::cout << "              tau_i = " << tau_i << std::endl;
	std::cout << "              tau_e = " << tau_e << std::endl;
	std::cout << "              anisotropy = " << anisotropy << std::endl;

	// ///////////////////////////////////////////////////////////////////////
	// ///////////////////////////////////////////////////////////////////////
	// Setup Exporters
	std::cout << "* BIDOMAIN: Initializing exporters " << std::endl;
	M_bidomainExporter.reset(new Exporter(M_equationSystems.get_mesh()));
	M_ionicModelExporter.reset(new Exporter(M_equationSystems.get_mesh()));
	M_parametersExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
	M_bidomainEXOExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
	M_potentialEXOExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));

	struct stat out_dir;
	if (stat(&M_outputFolder[0], &out_dir) != 0)
	{
		if (bidomain_system.get_mesh().comm().rank() == 0)
		{
			mkdir(M_outputFolder.c_str(), 0777);
		}
	}

	M_equationType = BidomainEquationType::ParabolicParabolicHyperbolic;
	if (tau_e == tau_i)
		M_equationType = BidomainEquationType::ParabolicEllipticHyperbolic;
	if (tau_i == tau_e && 0 == tau_i)
		M_equationType = BidomainEquationType::ParabolicEllipticBidomain;

}

void Bidomain::restart(EXOExporter& importer, int step, bool restart)
{
	if (restart)
	{
		std::cout << "* BIDOMAIN: restarting from step: " << step << std::endl;
		BidomainSystem& bidomain_system  =  M_equationSystems.get_system<BidomainSystem>("bidomain");
		IonicModelSystem& ionic_model_system =  M_equationSystems.get_system<IonicModelSystem>("ionic_model");
		BidomainSystem& wave_system =  M_equationSystems.get_system<BidomainSystem>("wave");
//
	    auto names = importer.get_nodal_var_names();
	    std::cout << "* Bidomain: Importer available nodal fields" << std::endl;
	    for(auto && n : names) std::cout << n << std::endl;
//
		std::cout << "* Bidomain: Importing  V"<< std::endl;
		importer.copy_nodal_solution(wave_system, "V", step);
		std::cout << "* Bidomain: Importing  Q"<< std::endl;
		importer.copy_nodal_solution(bidomain_system, "Q", step);
		std::cout << "* Bidomain: Importing  Ve"<< std::endl;
		importer.copy_nodal_solution(bidomain_system, "Ve", step);
		unsigned int n_vars = ionic_model_system.n_vars();
		std::cout << "* Bidomain: Importing  Ionic Model: num_vars: "<< n_vars <<  std::endl;
		for(unsigned int i = 0; i < n_vars; ++i)
		importer.copy_nodal_solution(ionic_model_system, ionic_model_system.variable_name(i), step);
	}
}

void Bidomain::init(double time)
{
	// Setting initial conditions
	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");
	IonicModelSystem& ionic_model_system = M_equationSystems.get_system
			< IonicModelSystem > ("ionic_model");
	int num_vars = ionic_model_system.n_vars();
	std::vector<double> init_values(num_vars + 1, 0.0);
	M_ionicModelPtr->initialize(init_values);
	// WAVE
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");
	auto first_local_index = wave_system.solution->first_local_index();
	auto last_local_index = wave_system.solution->last_local_index();
	for (auto i = first_local_index; i < last_local_index; ++i)
	{
		wave_system.solution->set(i, init_values[0]);
	}
	first_local_index = ionic_model_system.solution->first_local_index();
	last_local_index = ionic_model_system.solution->last_local_index();
	std::cout << "\n* BIDOMAIN: Setting initial values ... " << std::flush;
	for (auto i = first_local_index; i < last_local_index;)
	{
		for (int nv = 0; nv < num_vars; ++nv)
		{
			ionic_model_system.solution->set(i, init_values[nv + 1]);
			++i;
		}
	}
	std::cout << "done " << std::endl;

	std::string v_ic = M_datafile("bidomain/ic", "");
	if (v_ic != "")
	{
		std::cout << "* BIDOMAIN: Found bidomain initial condition: " << v_ic
				<< std::endl;
		SpiritFunction bidomain_ic;
		bidomain_ic.read(v_ic);
		M_equationSystems.parameters.set < libMesh::Real > ("time") = time;
		bidomain_system.time = time;
		std::cout
				<< "* BIDOMAIN: Projecting initial condition to bidomain system ... "
				<< std::flush;
		wave_system.project_solution(&bidomain_ic);
		std::cout << " done" << std::endl;
	}
	std::cout
			<< "* BIDOMAIN: Copying initial conditions to vectors at n nd at n-1... "
			<< std::flush;

	// Close vectors and update the values in old_local_solution and older_local_solution
	bidomain_system.solution->close();
	bidomain_system.old_local_solution->close();
	bidomain_system.older_local_solution->close();
	ionic_model_system.solution->close();
	ionic_model_system.old_local_solution->close();
	ionic_model_system.older_local_solution->close();
	wave_system.solution->close();
	wave_system.old_local_solution->close();
	wave_system.older_local_solution->close();
	advance();
	std::cout << " done" << std::endl;

	IonicModelSystem& istim_system = M_equationSystems.get_system
			< IonicModelSystem > ("istim");
	std::cout << "* BIDOMAIN: Initializing activation times to -1  ... "
			<< std::flush;
	ParameterSystem& activation_times_system = M_equationSystems.get_system
			< ParameterSystem > ("activation_times");
	first_local_index = activation_times_system.solution->first_local_index();
	last_local_index = activation_times_system.solution->last_local_index();

	for (auto i = first_local_index; i < last_local_index; ++i)
	{
		activation_times_system.solution->set(i, -1.0);
	}
	std::cout << " done" << std::endl;

	std::string fibers_data = M_datafile("bidomain/fibers", "1.0, 0.0, 0.0");
	std::string sheets_data = M_datafile("bidomain/sheets", "0.0, 1.0, 0.0");
	std::string xfibers_data = M_datafile("bidomain/xfibers", "0.0, 0.0, 1.0");

	SpiritFunction fibers_func;
	SpiritFunction sheets_func;
	SpiritFunction xfibers_func;

	fibers_func.read(fibers_data);
	sheets_func.read(sheets_data);
	xfibers_func.read(xfibers_data);

	ParameterSystem& fiber_system = M_equationSystems.get_system
			< ParameterSystem > ("fibers");
	ParameterSystem& sheets_system = M_equationSystems.get_system
			< ParameterSystem > ("sheets");
	ParameterSystem& xfiber_system = M_equationSystems.get_system
			< ParameterSystem > ("xfibers");

	std::cout << "* BIDOMAIN: Creating fibers from function: " << fibers_data
			<< std::flush;
	fiber_system.project_solution(&fibers_func);
	std::cout << " done" << std::endl;
	std::cout << "* BIDOMAIN: Creating sheets from function: " << sheets_data
			<< std::flush;
	sheets_system.project_solution(&sheets_func);
	std::cout << " done" << std::endl;
	std::cout << "* BIDOMAIN: Creating xfibers from function: " << xfibers_data
			<< std::flush;
	xfiber_system.project_solution(&xfibers_func);
	std::cout << " done" << std::endl;

	ParameterSystem& intra_conductivity_system = M_equationSystems.get_system
			< ParameterSystem > ("intra_conductivity");
	std::string Dffi_data = M_datafile("bidomain/Dffi", "2.0");
	std::string Dssi_data = M_datafile("bidomain/Dssi", "0.2");
	std::string Dnni_data = M_datafile("bidomain/Dnni", "0.2");
	SpiritFunction intra_conductivity_func;
	intra_conductivity_func.add_function(Dffi_data);
	intra_conductivity_func.add_function(Dssi_data);
	intra_conductivity_func.add_function(Dnni_data);
	intra_conductivity_system.project_solution(&intra_conductivity_func);
	ParameterSystem& extra_conductivity_system = M_equationSystems.get_system
			< ParameterSystem > ("extra_conductivity");
	std::string Dffe_data = M_datafile("bidomain/Dffe", "1.5");
	std::string Dsse_data = M_datafile("bidomain/Dsse", "1.0");
	std::string Dnne_data = M_datafile("bidomain/Dnne", "1.0");
	SpiritFunction extra_conductivity_func;
	extra_conductivity_func.add_function(Dffe_data);
	extra_conductivity_func.add_function(Dsse_data);
	extra_conductivity_func.add_function(Dnne_data);
	extra_conductivity_system.project_solution(&extra_conductivity_func);

	std::map < std::string, libMesh::SolverType > solver_map;
	solver_map["cg"] = libMesh::CG;
	solver_map["cgs"] = libMesh::CGS;
	solver_map["gmres"] = libMesh::GMRES;
	std::map < std::string, libMesh::PreconditionerType > prec_map;
	prec_map["jacobi"] = libMesh::JACOBI_PRECOND;
	prec_map["sor"] = libMesh::SOR_PRECOND;
	prec_map["ssor"] = libMesh::SSOR_PRECOND;
	prec_map["cholesky"] = libMesh::CHOLESKY_PRECOND;
	prec_map["lu"] = libMesh::LU_PRECOND;
	prec_map["ilu"] = libMesh::ILU_PRECOND;
	prec_map["amg"] = libMesh::AMG_PRECOND;

	ParameterSystem& procID_system = M_equationSystems.get_system
			< ParameterSystem > ("ProcID");
	first_local_index = procID_system.solution->first_local_index();
	last_local_index = procID_system.solution->last_local_index();

	double myrank = static_cast<double>(M_equationSystems.comm().rank());
	for (int el = first_local_index; el < last_local_index; ++el)
	{
		procID_system.solution->set(el, myrank);
	}
	procID_system.solution->close();

}

void Bidomain::update_pacing(double time)
{
	// Add the applied current to this system
	IonicModelSystem& istim_system = M_equationSystems.get_system
			< IonicModelSystem > ("istim");
	auto& istim_ve = istim_system.solution;

	istim_system.time = time;
	M_pacing->update(time);
	istim_system.project_solution(&M_pacing->pacing());
}

void Bidomain::update_activation_time(double time, double threshold)
{
	ParameterSystem& activation_times_system = M_equationSystems.get_system
			< ParameterSystem > ("activation_times");
	// WAVE
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");

	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");

	auto first = wave_system.solution->first_local_index();
	auto last = wave_system.solution->last_local_index();
	double v = 0.0;
	double at = 0.0;
	for (int index = first; index < last; index++)
	{
		v = (*wave_system.solution)(index);
		at = (*activation_times_system.solution)(index);
		if (v > threshold && at < 0.0)
		{
			activation_times_system.solution->set(index, time);
		}
	}
	activation_times_system.solution->close();
}

void Bidomain::generate_fibers(const GetPot& data, const std::string& section)
{
	std::cout << "* BIDOMAIN: Creating fiber fields" << std::endl;

	std::cout << "* BIDOMAIN: Solving Poisson porblem" << std::endl;

	libMesh::Mesh new_mesh(
			dynamic_cast<libMesh::Mesh&>(M_equationSystems.get_mesh()));
	libMesh::EquationSystems es(new_mesh);
	Poisson p(es);
	p.setup(data, section);
	p.assemble_system();
	p.solve_system();
	p.compute_elemental_solution_gradient();
	p.save_exo();

	const auto& sol_ptr = p.get_P0_solution();

	ParameterSystem& fiber_system = M_equationSystems.get_system
			< ParameterSystem > ("fibers");
	ParameterSystem& sheets_system = M_equationSystems.get_system
			< ParameterSystem > ("sheets");
	ParameterSystem& xfiber_system = M_equationSystems.get_system
			< ParameterSystem > ("xfibers");

	*sheets_system.solution = *p.get_gradient();

	auto first_local_index = fiber_system.solution->first_local_index();
	auto last_local_index = fiber_system.solution->last_local_index();
	auto first_local_index_sol = sol_ptr->first_local_index();
	auto last_local_index_sol = sol_ptr->last_local_index();

	double norm = 0.0;
	double fx = 0.0;
	double fy = 0.0;
	double fz = 0.0;
	double cx = data(section + "/centerline_x", 0.0);
	double cy = data(section + "/centerline_y", 0.0);
	double cz = data(section + "/centerline_z", 1.0);
	double cdot = 0.0;
	double sx = 0.0;
	double sy = 0.0;
	double sz = 0.0;
	double xfx = 0.0;
	double xfy = 0.0;
	double xfz = 0.0;
	double epi_angle = data(section + "/epi_angle", -60.0);
	double endo_angle = data(section + "/endo_angle", 60.0);
	double potential = 0.0;
	double W11 = 0.0;
	double W12 = 0.0;
	double W13 = 0.0;
	double W21 = 0.0;
	double W22 = 0.0;
	double W23 = 0.0;
	double W31 = 0.0;
	double W32 = 0.0;
	double W33 = 0.0;
	//
	double R11 = 0.0;
	double R12 = 0.0;
	double R13 = 0.0;
	double R21 = 0.0;
	double R22 = 0.0;
	double R23 = 0.0;
	double R31 = 0.0;
	double R32 = 0.0;
	double R33 = 0.0;
	double sa = 0.0;
	double sa2 = 0.0;
	double teta1 = 0.0;
	double teta2 = 0.0;
	double teta = 0.0;
	double m = 0.0;
	double q = 0.0;
	double f0x = 0.0;
	double f0y = 0.0;
	double f0z = 0.0;

	auto normalize = [](double& x, double& y, double& z,
			double X, double Y, double Z)
	{
		double norm = std::sqrt( x * x + y * y + z * z);
		if(norm >= 1e-12 )
		{
			x /= norm;
			y /= norm;
			z /= norm;
		}
		else
		{
			x = X;
			y = Y;
			z = Z;
		}
	};

	std::cout << "* BIDOMAIN: Computing rule-based fiber fields" << std::endl;

//    sol_ptr->print();

	auto j = first_local_index_sol;

	for (auto i = first_local_index; i < last_local_index;)
	{

//        std::cout << "* MONODOMAIN: getting sol ... " << std::flush;
		potential = (*sol_ptr)(j);
		j++;
//        std::cout << " done" << std::endl;

		sx = (*sheets_system.solution)(i);
		sy = (*sheets_system.solution)(i + 1);
		sz = (*sheets_system.solution)(i + 2);
		normalize(sx, sy, sz, 0.0, 1.0, 0.0);

//        norm = std::sqrt( sx * sx + sy * sy + sz * sz);
//        if(norm >= 1e-12 )
//        {
//            sx /= norm;
//            sy /= norm;
//            sz /= norm;
//        }
//        else
//        {
//            sx = 1.0;
//            sy = 0.0;
//            sz = 0.0;
//        }

		cdot = cx * sx + cy * sy + cz * sz;

		xfx = cx - cdot * sx;
		xfy = cy - cdot * sy;
		xfz = cz - cdot * sz;
		normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

		fx = sy * xfz - sz * xfy;
		fy = sz * xfx - sx * xfz;
		fz = sx * xfy - sy * xfx;
		normalize(fx, fy, fz, 1.0, 0.0, 0.0);

		teta1 = M_PI * epi_angle / 180.0;
		teta2 = M_PI * endo_angle / 180.0;
		m = (teta1 - teta2);
		q = teta2;
		teta = m * potential + q;

		//*************************************************************//
		// The fiber field F is a rotation of the flat fiber field f
		// F = R f
		// where R is the rotation matrix.
		// To compute R we need the sin(teta) and
		// the sin(teta)^2 and the cross-product matrix W (check
		// rodrigues formula on wikipedia :) )
		//*************************************************************//
		sa = std::sin(teta);
		sa2 = 2.0 * std::sin(0.5 * teta) * std::sin(0.5 * teta);

		W11 = 0.0;
		W12 = -sz;
		W13 = sy;
		W21 = sz;
		W22 = 0.0;
		W23 = -sx;
		W31 = -sy;
		W32 = sx;
		W33 = 0.0;
		//
		R11 = 1.0 + sa * W11 + sa2 * (sx * sx - 1.0);
		R12 = 0.0 + sa * W12 + sa2 * (sx * sy);
		R13 = 0.0 + sa * W13 + sa2 * (sx * sz);
		R21 = 0.0 + sa * W21 + sa2 * (sy * sx);
		R22 = 1.0 + sa * W22 + sa2 * (sy * sy - 1.0);
		R23 = 0.0 + sa * W23 + sa2 * (sy * sz);
		R31 = 0.0 + sa * W31 + sa2 * (sz * sx);
		R32 = 0.0 + sa * W32 + sa2 * (sz * sy);
		R33 = 1.0 + sa * W33 + sa2 * (sz * sz - 1.0);

		f0x = R11 * fx + R12 * fy + R13 * fz;
		f0y = R21 * fx + R22 * fy + R23 * fz;
		f0z = R31 * fx + R32 * fy + R33 * fz;
		normalize(f0x, f0y, f0z, 1.0, 0.0, 0.0);

		xfx = f0y * sz - f0z * sy;
		xfy = f0z * sx - f0x * sz;
		xfz = f0x * sy - f0y * sx;
		normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

		sheets_system.solution->set(i, sx);
		sheets_system.solution->set(i + 1, sy);
		sheets_system.solution->set(i + 2, sz);
		xfiber_system.solution->set(i, xfx);
		xfiber_system.solution->set(i + 1, xfy);
		xfiber_system.solution->set(i + 2, xfz);
//        std::cout << "* MONODOMAIN: setting sol ... " << std::flush;

		fiber_system.solution->set(i, f0x);
		fiber_system.solution->set(i + 1, f0y);
		fiber_system.solution->set(i + 2, f0z);
//        std::cout << " done" << std::endl;

		i += 3;
	}

}

double Bidomain::last_activation_time()
{
	ParameterSystem& activation_times_system = M_equationSystems.get_system
			< ParameterSystem > ("activation_times");
	return static_cast<double>(activation_times_system.solution->max());
}

void Bidomain::save(int step)
{
	std::cout << "* BIDOMAIN: GMVIO::Exporting " << step << " in: "
			<< M_outputFolder << " ... " << std::flush;
//    M_equationSystems->
	M_bidomainExporter->write_equation_systems(
			M_outputFolder + "bidomain.gmv." + std::to_string(step),
			M_equationSystems, &M_bidomainExporterNames);
	M_ionicModelExporter->write_equation_systems(
			M_outputFolder + "ionic_model.gmv." + std::to_string(step),
			M_equationSystems, &M_ionicModelExporterNames);
	std::cout << "done " << std::endl;
}
void Bidomain::init_exo_output()
{
	M_bidomainEXOExporter->write_equation_systems(
			M_outputFolder + "bidomain.exo", M_equationSystems);
	M_bidomainEXOExporter->append(true);
	std::vector < std::string > output;
	output.push_back("V");
	output.push_back("Ve");
	output.push_back("Q");
	M_potentialEXOExporter->set_output_variables(output);
	M_potentialEXOExporter->write_equation_systems(
			M_outputFolder + "potential.exo", M_equationSystems);
	M_potentialEXOExporter->append(true);
}

void Bidomain::save_exo(int step, double time)
{
	std::cout << "* BIDOMAIN: EXODUSII::Exporting bidomain.exo at time " << time
			<< " in: " << M_outputFolder << " ... " << std::flush;
	M_bidomainEXOExporter->write_timestep(M_outputFolder + "bidomain.exo",
			M_equationSystems, step, time);
	std::cout << "done " << std::endl;
}

void Bidomain::save_potential(int step, double time)
{
	std::cout << "* BIDOMAIN: EXODUSII::Exporting potential.exo at time "
			<< time << " in: " << M_outputFolder << " ... " << std::flush;

	M_potentialEXOExporter->write_timestep(M_outputFolder + "potential.exo",
			M_equationSystems, step, time);
	std::cout << "done " << std::endl;
}

void Bidomain::amr(libMesh::MeshRefinement& mesh_refinement,
		const std::string& type)
{
	libMesh::ErrorVector error;
	libMesh::ErrorEstimator * p_error_estimator;
	if ("kelly" == type)
		p_error_estimator = new libMesh::KellyErrorEstimator;
	else if ("disc" == type)
		p_error_estimator = new libMesh::DiscontinuityMeasure;
	else
		p_error_estimator = new libMesh::LaplacianErrorEstimator;
	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");
	// WAVE
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");
	p_error_estimator->estimate_error(wave_system, error);
	mesh_refinement.flag_elements_by_error_fraction(error);
	mesh_refinement.refine_and_coarsen_elements();
	M_equationSystems.reinit();
	delete p_error_estimator;
}

void Bidomain::save_parameters()
{
	std::cout << "* BIDOMAIN: EXODUSII::Exporting parameters in: "
			<< M_outputFolder << " ... " << std::flush;
	M_parametersExporter->write_equation_systems(
			M_outputFolder + "parameters.exo", M_equationSystems,
			&M_parametersExporterNames);
	M_parametersExporter->write_element_data(M_equationSystems);
	std::cout << "done " << std::endl;
}

void Bidomain::evaluate_mesh_size()
{
	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	double h = 0;
	double h_max = 0;
	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		h = elem->hmax();
		h_max = std::max(h_max, h);
	}
	mesh.comm().max(h_max);
	M_meshSize = h_max;
}

void Bidomain::assemble_matrices(double dt)
{
	// Implemented only FirstOrderIMEX

	using libMesh::UniquePtr;

	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	const unsigned int max_dim = 3;
	const libMesh::Real Chi = M_equationSystems.parameters.get < libMesh::Real
			> ("Chi");
	const libMesh::Real tau_e = M_equationSystems.parameters.get < libMesh::Real
			> ("tau_e");
	const libMesh::Real tau_i = M_equationSystems.parameters.get < libMesh::Real
			> ("tau_i");
	double Cm = M_ionicModelPtr->membraneCapacitance();

	// Get a reference to the LinearImplicitSystem we are solving
	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");
	IonicModelSystem& ionic_model_system = M_equationSystems.get_system
			< IonicModelSystem > ("ionic_model");
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");

	unsigned int Q_var = bidomain_system.variable_number("Q");
	unsigned int Ve_var = bidomain_system.variable_number("Ve");

	bidomain_system.get_matrix("mass").zero();
	bidomain_system.get_matrix("lumped_mass").zero();
	bidomain_system.get_matrix("high_order_mass").zero();
	bidomain_system.get_matrix("stiffness").zero();
	bidomain_system.get_vector("lumped_mass_vector").zero();

	ParameterSystem& fiber_system = M_equationSystems.get_system
			< ParameterSystem > ("fibers");
	ParameterSystem& sheets_system = M_equationSystems.get_system
			< ParameterSystem > ("sheets");
	ParameterSystem& xfiber_system = M_equationSystems.get_system
			< ParameterSystem > ("xfibers");
	ParameterSystem& intra_conductivity_system = M_equationSystems.get_system
			< ParameterSystem > ("intra_conductivity");
	ParameterSystem& extra_conductivity_system = M_equationSystems.get_system
			< ParameterSystem > ("extra_conductivity");

	// A reference to the  DofMap object for this system.  The  DofMap
	// object handles the index translation from node and element numbers
	// to degree of freedom numbers.  We will talk more about the  DofMap
	// in future examples.
	const libMesh::DofMap & dof_map_bidomain = bidomain_system.get_dof_map();
	const libMesh::DofMap & dof_map_wave = wave_system.get_dof_map();
	const libMesh::DofMap & dof_map_fibers = fiber_system.get_dof_map();

	// Get a constant reference to the Finite Element type
	// for the first (and only) variable in the system.
	libMesh::FEType fe_type_qp1 = dof_map_bidomain.variable_type(Q_var);

	// Build a Finite Element object of the specified type.  Since the
	// FEBase::build() member dynamically creates memory we will
	// store the object as a UniquePtr<FEBase>.  This can be thought
	// of as a pointer that will clean up after itself.  Introduction Example 4
	// describes some advantages of  UniquePtr's in the context of
	// quadrature rules.
	UniquePtr < libMesh::FEBase
			> fe_qp1(libMesh::FEBase::build(dim, fe_type_qp1));

	// A 5th order Gauss quadrature rule for numerical integration.
	libMesh::QGauss qrule(dim, libMesh::THIRD);

	// Tell the finite element object to use our quadrature rule.
	fe_qp1->attach_quadrature_rule(&qrule);

	// Here we define some references to cell-specific data that
	// will be used to assemble the linear system.
	//
	// The element Jacobian * quadrature weight at each integration point.
	const std::vector<libMesh::Real> & JxW_qp1 = fe_qp1->get_JxW();

	// The physical XY locations of the quadrature points on the element.
	// These might be useful for evaluating spatially varying material
	// properties at the quadrature points.
	const std::vector<libMesh::Point> & q_point_qp1 = fe_qp1->get_xyz();

	// The element shape functions evaluated at the quadrature points.
	const std::vector<std::vector<libMesh::Real> > & phi_qp1 =
			fe_qp1->get_phi();

	// The element shape function gradients evaluated at the quadrature
	// points.
	const std::vector<std::vector<libMesh::RealGradient> > & dphi_qp1 =
			fe_qp1->get_dphi();

	const std::vector<std::vector<libMesh::Real> > & dphidx_qp1 =
			fe_qp1->get_dphidx();
	const std::vector<std::vector<libMesh::Real> > & dphidy_qp1 =
			fe_qp1->get_dphidy();
	const std::vector<std::vector<libMesh::Real> > & dphidz_qp1 =
			fe_qp1->get_dphidz();

	// Define data structures to contain the element matrix
	// and right-hand-side vector contribution.  Following
	// basic finite element terminology we will denote these
	// "Ke" and "Fe".  These datatypes are templated on
	//  Number, which allows the same code to work for real
	// or complex numbers.
	libMesh::DenseMatrix<libMesh::Number> Ke;
	libMesh::DenseMatrix<libMesh::Number> Kie;
	libMesh::DenseMatrix<libMesh::Number> Me;
	libMesh::DenseMatrix<libMesh::Number> Mel;
	libMesh::DenseVector<libMesh::Number> Fe;

	// for interior penalty
//     UniquePtr<libMesh::FEBase> fe_elem_face(libMesh::FEBase::build(dim, fe_type_qp1));
//     UniquePtr<libMesh::FEBase> fe_neighbor_face(libMesh::FEBase::build(dim, fe_type_qp1));
	// Tell the finite element object to use our quadrature rule.
	libMesh::QGauss qface(dim - 1, fe_type_qp1.default_quadrature_order());

//     fe_elem_face->attach_quadrature_rule(&qface);
//     fe_neighbor_face->attach_quadrature_rule(&qface);
	// Data for surface integrals on the element boundary
//     const std::vector<std::vector<libMesh::Real> > &  phi_face = fe_elem_face->get_phi();
//     const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_elem_face->get_dphi();
//     const std::vector<libMesh::Real> & JxW_face = fe_elem_face->get_JxW();
//     const std::vector<libMesh::Point> & qface_normals = fe_elem_face->get_normals();
//     const std::vector<libMesh::Point> & qface_points = fe_elem_face->get_xyz();
//     // Data for surface integrals on the neighbor boundary
//     const std::vector<std::vector<libMesh::Real> > &  phi_neighbor_face = fe_neighbor_face->get_phi();
//     const std::vector<std::vector<libMesh::RealGradient> > & dphi_neighbor_face = fe_neighbor_face->get_dphi();
	// Data structures to contain the element and neighbor boundary matrix
	// contribution. This matrices will do the coupling beetwen the dofs of
	// the element and those of his neighbors.
	// Ken: matrix coupling elem and neighbor dofs
//     libMesh::DenseMatrix<libMesh::Number> Kne;
//     libMesh::DenseMatrix<libMesh::Number> Ken;
//     libMesh::DenseMatrix<libMesh::Number> Kee;
//     libMesh::DenseMatrix<libMesh::Number> Knn;

	// This vector will hold the degree of freedom indices for
	// the element.  These define where in the global system
	// the element degrees of freedom get mapped.
	std::vector < libMesh::dof_id_type > dof_indices;
	std::vector < libMesh::dof_id_type > dof_indices_Q;
	std::vector < libMesh::dof_id_type > dof_indices_Ve;
	std::vector < libMesh::dof_id_type > dof_indices_V;

	std::vector < libMesh::dof_id_type > dof_indices_fibers;

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
	libMesh::MeshBase::const_element_iterator el_start =
			mesh.active_local_elements_begin();
	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	// Loop over the elements.  Note that  ++el is preferred to
	// el++ since the latter requires an unnecessary temporary
	// object.
	double f0[3];
	double s0[3];
	double n0[3];
	libMesh::RealGradient DigradV;
	libMesh::RealGradient DigradVe;
	libMesh::RealGradient DiegradVe;
	libMesh::TensorValue<double> D0i;
	libMesh::TensorValue<double> D0e;

	std::cout << "start loop over elements" << std::endl;
	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		const unsigned int elem_id = elem->id();
		dof_map_bidomain.dof_indices(elem, dof_indices);
		dof_map_bidomain.dof_indices(elem, dof_indices_Q, Q_var);
		dof_map_bidomain.dof_indices(elem, dof_indices_Ve, Ve_var);
		dof_map_wave.dof_indices(elem, dof_indices_V, 0);
		dof_map_fibers.dof_indices(elem, dof_indices_fibers);

		// Compute the element-specific data for the current
		// element.  This involves computing the location of the
		// quadrature points (q_point) and the shape functions
		// (phi, dphi) for the current element.
		fe_qp1->reinit(elem);

		// Zero the element matrix and right-hand side before
		// summing them.  We use the resize member here because
		// the number of degrees of freedom might have changed from
		// the last element.  Note that this will be the case if the
		// element type is different (i.e. the last element was a
		// triangle, now we are on a quadrilateral).

		// The  DenseMatrix::resize() and the  DenseVector::resize()
		// members will automatically zero out the matrix  and vector.
		auto n_dofs = dof_indices.size();
		const unsigned int n_Q_dofs = dof_indices_Q.size();
		auto n_dofs_V = dof_indices_V.size();
		Ke.resize(n_dofs, n_dofs);
		Kie.resize(n_dofs_V, n_dofs_V);
		Me.resize(n_dofs, n_dofs);
		Mel.resize(n_dofs, n_dofs);
		Fe.resize(n_dofs);

		//        std::cout << "Fibers" << std::endl;
		// fiber direction
		f0[0] = (*fiber_system.solution)(dof_indices_fibers[0]);
		f0[1] = (*fiber_system.solution)(dof_indices_fibers[1]);
		f0[2] = (*fiber_system.solution)(dof_indices_fibers[2]);
		// sheet direction
		s0[0] = (*sheets_system.solution)(dof_indices_fibers[0]);
		s0[1] = (*sheets_system.solution)(dof_indices_fibers[1]);
		s0[2] = (*sheets_system.solution)(dof_indices_fibers[2]);
		// crossfiber direction
		n0[0] = (*xfiber_system.solution)(dof_indices_fibers[0]);
		n0[1] = (*xfiber_system.solution)(dof_indices_fibers[1]);
		n0[2] = (*xfiber_system.solution)(dof_indices_fibers[2]);
		// Conductivity tensor
		double Dffe = (*extra_conductivity_system.solution)(
				dof_indices_fibers[0]);
		double Dsse = (*extra_conductivity_system.solution)(
				dof_indices_fibers[1]);
		double Dnne = (*extra_conductivity_system.solution)(
				dof_indices_fibers[2]);
		double Dffi = (*intra_conductivity_system.solution)(
				dof_indices_fibers[0]);
		double Dssi = (*intra_conductivity_system.solution)(
				dof_indices_fibers[1]);
		double Dnni = (*intra_conductivity_system.solution)(
				dof_indices_fibers[2]);

		switch (M_anisotropy)
		{
		case Anisotropy::Isotropic:
		{
			for (int idim = 0; idim < max_dim; ++idim)
			{
				D0i(idim, idim) = Dffi;
				D0e(idim, idim) = Dffe;
			}

			break;
		}

		case Anisotropy::TransverselyIsotropic:
		{
			for (int idim = 0; idim < max_dim; ++idim)
			{
				for (int jdim = 0; jdim < max_dim; ++jdim)
				{

					D0i(idim, jdim) = (Dffi - Dssi) * f0[idim] * f0[jdim];
					D0e(idim, jdim) = (Dffe - Dsse) * f0[idim] * f0[jdim];
				}
				D0i(idim, idim) += Dssi;
				D0e(idim, idim) += Dsse;
			}
			break;
		}

		case Anisotropy::Orthotropic:
		default:
		{
			for (int idim = 0; idim < max_dim; ++idim)
			{
				for (int jdim = 0; jdim < max_dim; ++jdim)
				{
					D0i(idim, jdim) = Dffi * f0[idim] * f0[jdim]
							+ Dssi * s0[idim] * s0[jdim]
							+ Dnni * n0[idim] * n0[jdim];

					D0e(idim, jdim) = Dffe * f0[idim] * f0[jdim]
							+ Dsse * s0[idim] * s0[jdim]
							+ Dnne * n0[idim] * n0[jdim];
				}
			}
			break;
		}
		}
		D0i /= Chi;
		D0e /= Chi;

//    	 std::cout << "f  = [" << f0[0] << "," << f0[1]  << ", " << f0[2]  << "]"<< std::endl;
//    	 std::cout << "s = [" << s0[0] << "," << s0[1]  << ", " << s0[2]  << "]"<< std::endl;
//    	 std::cout << "n = [" << n0[0] << "," << n0[1]  << ", " << n0[2]  << "]"<< std::endl;
//	 	 std::cout << "D0 = " << std::endl;
//		 for(int idim = 0; idim < dim; ++idim )
//         {
//        	 std::cout << D0(idim, 0) << "," << D0(idim, 1) << ", " << D0(idim, 2) << std::endl;
//         }

		// Assemble Mass terms
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
		{
			//  Matrix
			for (unsigned int i = 0; i < phi_qp1.size(); i++)
			{
				DigradV = D0i * dphi_qp1[i][qp];
				DigradVe = D0i * dphi_qp1[i][qp];
				DiegradVe = (D0i + D0e) * dphi_qp1[i][qp];

				for (unsigned int j = 0; j < phi_qp1.size(); j++)
				{
					// Q Mass term
					Me(i, j) += JxW_qp1[qp] * (phi_qp1[i][qp] * phi_qp1[j][qp]);
					Mel(i, i) += JxW_qp1[qp]
							* (phi_qp1[i][qp] * phi_qp1[j][qp]);
					Fe(i) += JxW_qp1[qp] * (phi_qp1[i][qp] * phi_qp1[j][qp]);
					// Ve Mass terms
					Me(i + n_Q_dofs, j + n_Q_dofs) += JxW_qp1[qp]
							* (phi_qp1[i][qp] * phi_qp1[j][qp]);
					Mel(i + n_Q_dofs, i + n_Q_dofs) += JxW_qp1[qp]
							* (phi_qp1[i][qp] * phi_qp1[j][qp]);
					Fe(i + n_Q_dofs) += JxW_qp1[qp]
							* (phi_qp1[i][qp] * phi_qp1[j][qp]);

					// Block QQ : ( tau_i + dt ) * Cm * M_L + dt^2 * Ki
					// Using lumped mass matrix
					Ke(i, i) += (tau_i + dt) * Cm * JxW_qp1[qp]
							* (phi_qp1[i][qp] * phi_qp1[j][qp]);
					Ke(i, j) += dt * dt * JxW_qp1[qp] * DigradV
							* dphi_qp1[j][qp];
					// Block QVe : dt * Ki
					// Using lumped mass matrix
					Ke(i, j + n_Q_dofs) += dt * JxW_qp1[qp] * DigradVe
							* dphi_qp1[j][qp];

					// Block VeQ : (tau_i-tau_e) / dt * Cm * M + dt * Ki
					// Using lumped mass matrix
					Ke(i + n_Q_dofs, i) += (tau_i - tau_e) / dt * Cm
							* JxW_qp1[qp] * (phi_qp1[i][qp] * phi_qp1[j][qp]);
					Ke(i + n_Q_dofs, j) += dt * JxW_qp1[qp] * DigradV
							* dphi_qp1[j][qp];

					// Block VeVe : Kie
					// Using lumped mass matrix
					Ke(i + n_Q_dofs, j + n_Q_dofs) += JxW_qp1[qp] * DiegradVe
							* dphi_qp1[j][qp];
					// Ki
					Kie(i, j) += JxW_qp1[qp] * DigradV * dphi_qp1[j][qp];

				}
			}
		}

		if (el == el_start)
		{
//             Ke(n_Q_dofs, n_Q_dofs) +=  1e5;
		}
		bidomain_system.get_matrix("mass").add_matrix(Me, dof_indices);
		bidomain_system.get_matrix("lumped_mass").add_matrix(Mel, dof_indices);
		bidomain_system.get_vector("lumped_mass_vector").add_vector(Fe,
				dof_indices);
		bidomain_system.matrix->add_matrix(Ke, dof_indices);
		wave_system.get_matrix("Ki").add_matrix(Kie, dof_indices_V);
	}
	// closing matrices and vectors
	bidomain_system.get_matrix("mass").close();
	bidomain_system.get_matrix("lumped_mass").close();
	bidomain_system.get_matrix("stiffness").close();
	bidomain_system.get_vector("lumped_mass_vector").close();
	bidomain_system.get_matrix("high_order_mass").close();
	bidomain_system.get_matrix("high_order_mass").add(0.5,
			bidomain_system.get_matrix("mass"));
	bidomain_system.get_matrix("high_order_mass").add(0.5,
			bidomain_system.get_matrix("lumped_mass"));
	bidomain_system.matrix->close();

	{
		typedef libMesh::PetscMatrix<libMesh::Number> Mat;
		typedef libMesh::PetscVector<libMesh::Number> Vec;

		libMesh::DenseVector<libMesh::Number> NSe;
		MatNullSpace   nullspace;

		libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
		const libMesh::MeshBase::const_node_iterator end_node =
				mesh.local_nodes_end();
		for (; node != end_node; ++node)
		{
			const libMesh::Node * nn = *node;
			dof_map_bidomain.dof_indices(nn, dof_indices_Q, Q_var);
			dof_map_bidomain.dof_indices(nn, dof_indices_Ve, Ve_var);
			bidomain_system.get_vector("nullspace").set(dof_indices_Q[0], 0.0);
			bidomain_system.get_vector("nullspace").set(dof_indices_Ve[0], 1.0);
		}
		bidomain_system.get_vector("nullspace").close();
		int size = bidomain_system.get_vector("nullspace").size();
		bidomain_system.get_vector("nullspace") /= std::sqrt(size/2);
		// setting solver type
		std::string solver_type = M_datafile("bidomain/linear_solver/type",
				"gmres");
		std::cout << "* BIDOMAIN: using " << solver_type << std::endl;
		std::string prec_type = M_datafile(
				"bidomain/linear_solver/preconditioner", "amg");
		std::cout << "* BIDOMAIN: using " << prec_type << std::endl;
		//M_linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( M_equationSystems.comm() );
		typedef libMesh::PetscLinearSolver<libMesh::Number> PetscSolver;
		M_linearSolver.reset(new PetscSolver(M_equationSystems.comm()));

		std::cout << "* BIDOMAIN: setting solver  " << prec_type << std::endl;
//        M_linearSolver->set_solver_type(solver_map.find(solver_type)->second);
//        M_linearSolver->set_preconditioner_type(prec_map.find(prec_type)->second);

		std::cout << "* BIDOMAIN: closing matrix  " << prec_type << std::endl;
		Vec& N = (static_cast<Vec&>(bidomain_system.get_vector("nullspace")));
        std::cout << "* BIDOMAIN: nullspace vector  " << prec_type << std::endl;
		auto vec = N.vec();
		std::cout << N.size() <<  ",  " <<
        std::cout << "* BIDOMAIN: nullspace vector 1 " << prec_type << std::endl;
		MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,1,&vec, &nullspace);
        std::cout << "* BIDOMAIN: nullspace vector  2" << prec_type << std::endl;
		Mat * mat = dynamic_cast<Mat *>(bidomain_system.matrix);
        std::cout << "* BIDOMAIN: nullspace created  " << prec_type << std::endl;

		MatSetNullSpace(mat->mat(), nullspace);
		MatNullSpaceDestroy(&nullspace);

		typedef libMesh::PetscMatrix<libMesh::Number> PetscMatrix;
		M_linearSolver->init(
				dynamic_cast<PetscMatrix *>(bidomain_system.matrix));
		std::cout << "* BIDOMAIN: initialized linear solver  " << prec_type
				<< std::endl;
	}
}

void Bidomain::advance()
{
	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");
	IonicModelSystem& ionic_model_system = M_equationSystems.get_system
			< IonicModelSystem > ("ionic_model");

	*bidomain_system.old_local_solution = *bidomain_system.solution;
	*bidomain_system.older_local_solution = *bidomain_system.solution;
	*ionic_model_system.old_local_solution = *ionic_model_system.solution;
	*ionic_model_system.older_local_solution = *ionic_model_system.solution;

	// WAVE
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");
	*wave_system.old_local_solution = *wave_system.solution;
	*wave_system.older_local_solution = *wave_system.solution;

	IonicModelSystem& iion_system = M_equationSystems.get_system
			< IonicModelSystem > ("iion");
	*iion_system.older_local_solution = *iion_system.old_local_solution;
	*iion_system.old_local_solution = *iion_system.solution;

}

void Bidomain::solve_reaction_step(double dt, double time, int step,
		bool useMidpoint, const std::string& mass,
		libMesh::NumericVector<libMesh::Number>* I4f_ptr)
{
	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();

	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain"); //Q
	IonicModelSystem& ionic_model_system = M_equationSystems.get_system
			< IonicModelSystem > ("ionic_model");
	IonicModelSystem& istim_system = M_equationSystems.get_system
			< IonicModelSystem > ("istim");
	// WAVE
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");
	IonicModelSystem& iion_system = M_equationSystems.get_system
			< IonicModelSystem > ("iion");

	bidomain_system.rhs->zero();
	iion_system.solution->zero();
	//iion_system.old_local_solution->zero();
	iion_system.get_vector("diion").zero();

	double Cm = M_ionicModelPtr->membraneCapacitance();
	const libMesh::Real tau_e = M_equationSystems.parameters.get < libMesh::Real
			> ("tau_e");
	const libMesh::Real tau_i = M_equationSystems.parameters.get < libMesh::Real
			> ("tau_i");

	int num_vars = ionic_model_system.n_vars();
	double istim = 0.0;
	double istim_zero = 0.0;
	double dIion = 0.0;
	std::vector<double> values(num_vars + 1, 0.0);
	std::vector<double> old_values(num_vars + 1, 0.0);
	int var_index = 0;

	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	const libMesh::DofMap & dof_map = bidomain_system.get_dof_map();
	const libMesh::DofMap & dof_map_V = wave_system.get_dof_map();
	const libMesh::DofMap & dof_map_gating = ionic_model_system.get_dof_map();

	std::vector < libMesh::dof_id_type > dof_indices;
	std::vector < libMesh::dof_id_type > dof_indices_V;
	std::vector < libMesh::dof_id_type > dof_indices_Ve;
	std::vector < libMesh::dof_id_type > dof_indices_Q;
	std::vector < libMesh::dof_id_type > dof_indices_gating;

	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map.dof_indices(nn, dof_indices_Q, 0);
		dof_map.dof_indices(nn, dof_indices_Ve, 1);
		dof_map_V.dof_indices(nn, dof_indices_V, 0);
		libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));

		double istim = M_pacing->eval(p, time);
		istim_zero = 0.0*istim;

		double Iion_old = 0.0;
		values[0] = (*wave_system.old_local_solution)(dof_indices_V[0]); //V^n
		old_values[0] = (*bidomain_system.old_local_solution)(dof_indices_Q[0]); //Q^n
		Iion_old = (*iion_system.old_local_solution)(dof_indices_V[0]); // gating

		for (int nv = 0; nv < num_vars; ++nv)
		{
			var_index = dof_indices_V[0] * num_vars + nv;
			values[nv + 1] = (*ionic_model_system.old_local_solution)(
					var_index);
			old_values[nv + 1] = values[nv + 1];
		}

		M_ionicModelPtr->updateVariables(values, istim_zero, dt);
		double Iion = M_ionicModelPtr->evaluateIonicCurrent(values, istim_zero,
				dt);
		dIion = M_ionicModelPtr->evaluatedIonicCurrent(values, old_values, dt,
				M_meshSize);
		double Isac = 0.0;
		if (I4f_ptr)
		{
			double I4f = (*I4f_ptr)(dof_indices_V[0]);
			Isac = M_ionicModelPtr->evaluateSAC(values[0], I4f);
		}
		Iion += Isac;

		iion_system.solution->set(dof_indices_V[0], Iion); // contains Istim
		istim_system.solution->set(dof_indices_V[0], istim);
		iion_system.get_vector("diion").set(dof_indices_V[0], dIion);
		for (int nv = 0; nv < num_vars; ++nv)
		{
			var_index = dof_indices_V[0] * num_vars + nv;
			ionic_model_system.solution->set(var_index, values[nv + 1]);
		}


	}

	iion_system.solution->close();
	istim_system.solution->close();
	ionic_model_system.solution->close();
	iion_system.get_vector("diion").close();
	iion_system.update();
	ionic_model_system.update();
	istim_system.update();
}

void Bidomain::form_system_rhs(double dt, bool useMidpoint,
		const std::string& mass)
{
	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");
	// WAVE
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");
	IonicModelSystem& iion_system = M_equationSystems.get_system
			< IonicModelSystem > ("iion");
	IonicModelSystem& istim_system = M_equationSystems.get_system
			< IonicModelSystem > ("istim");
//    istim_system.solution->print();
	bidomain_system.rhs->zero();
	bidomain_system.rhs->close();
	bidomain_system.get_vector("ionic_currents").zero();
	bidomain_system.get_vector("ionic_currents").close();
	bidomain_system.get_vector("old_solution").zero();
	bidomain_system.get_vector("old_solution").close();
	iion_system.get_vector("diion").close();

	double Cm = M_ionicModelPtr->membraneCapacitance();
	const libMesh::Real tau_e = M_equationSystems.parameters.get < libMesh::Real
			> ("tau_e");
	const libMesh::Real tau_i = M_equationSystems.parameters.get < libMesh::Real
			> ("tau_i");

	// Evaluate vector Ki*V^n
	wave_system.get_vector("KiV").zero();
	wave_system.get_vector("KiV").close();
	wave_system.get_matrix("Ki").close();
//	wave_system.get_matrix("Ki").print();
	// NOTE: We use wave_system.solution  because here V has not been updated yet.
	//       So the vector still stores V^n
	wave_system.get_matrix("Ki").vector_mult_add(wave_system.get_vector("KiV"),
			*wave_system.old_local_solution);

	// Evaluate
	// RHS_Q  = M * (tau_i * Cm * Q^n + dt * Iion + dt * tau_i * dIion )
	// RHS_Ve  = M * [ ( tau_i - tau_e ) / dt * Cm * Q^n + ( tau_i - tau_e ) dIion )
	double Qn = 0.0;
	double istim = 0.0;
	double Iion = 0.0;
	double dIion = 0.0;

	double rhsq = 0.0;
	double rhsve = 0.0;
	double rhs_oldq = 0.0;
	double rhs_oldve = 0.0;

	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	const libMesh::DofMap & dof_map = bidomain_system.get_dof_map();
	const libMesh::DofMap & dof_map_V = wave_system.get_dof_map();

	std::vector < libMesh::dof_id_type > dof_indices;
	std::vector < libMesh::dof_id_type > dof_indices_V;
	std::vector < libMesh::dof_id_type > dof_indices_Ve;
	std::vector < libMesh::dof_id_type > dof_indices_Q;

	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map.dof_indices(nn, dof_indices_Q, 0);
		dof_map.dof_indices(nn, dof_indices_Ve, 1);
		dof_map_V.dof_indices(nn, dof_indices_V, 0);

		Qn = (*bidomain_system.old_local_solution)(dof_indices_Q[0]); //Q^n
		istim = (*istim_system.solution)(dof_indices_V[0]); //Istim^n+1
		Iion = (*iion_system.solution)(dof_indices_V[0]); //Iion^* (it contains Istim)
		dIion = iion_system.get_vector("diion")(dof_indices_V[0]); // dIion^*

		rhsq = dt * (Iion + tau_i * dIion + 0.0*istim);
		rhsve = (tau_i - tau_e) * dIion + istim;

		rhs_oldq = tau_i * Cm * Qn;
		rhs_oldve = (tau_i - tau_e) / dt * Cm * Qn;

		bidomain_system.get_vector("old_solution").set(dof_indices_Q[0],
				rhs_oldq);
		bidomain_system.get_vector("old_solution").set(dof_indices_Ve[0],
				rhs_oldve);
		bidomain_system.get_vector("ionic_currents").set(dof_indices_Q[0],
				rhsq);
		bidomain_system.get_vector("ionic_currents").set(dof_indices_Ve[0],
				rhsve);
	}

	bidomain_system.get_vector("old_solution").close();
	bidomain_system.get_vector("ionic_currents").close();
	bidomain_system.rhs->add_vector(
			bidomain_system.get_vector("ionic_currents"),
			bidomain_system.get_matrix("mass"));
	bidomain_system.rhs->add_vector(bidomain_system.get_vector("old_solution"),
			bidomain_system.get_matrix("lumped_mass"));


	// Add KiVn to the rhs
	double KiVn = 0.0;

	node = mesh.local_nodes_begin();
	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map_V.dof_indices(nn, dof_indices_V, 0);
		dof_map.dof_indices(nn, dof_indices_Q, 0);
		dof_map.dof_indices(nn, dof_indices_Ve, 1);
		KiVn = wave_system.get_vector("KiV")(dof_indices_V[0]);
		bidomain_system.rhs->add(dof_indices_Q[0], -dt * KiVn);
		bidomain_system.rhs->add(dof_indices_Ve[0], -KiVn);
	}

	bidomain_system.rhs->close();

}

void Bidomain::solve_diffusion_step(double dt, double time, bool useMidpoint,
		const std::string& mass, bool reassemble)
{
	const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();

	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain"); // Q and Ve
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");// V



//    std::cout << "form rhs" << std::endl;
	form_system_rhs(dt, useMidpoint, mass);
//    std::cout << "done" << std::endl;

	double tol = 1e-12;
	double max_iter = 2000;
//        std::cout << "solving" << std::endl;

	std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);


//	bidomain_system.rhs->print();
//	bidomain_system.matrix->print();
//	double rhs_norm = bidomain_system.rhs->l2_norm();
//	std::cout <<"rhs norm: " << rhs_norm << std::endl;

	rval = M_linearSolver->solve(*bidomain_system.matrix,
			*bidomain_system.solution, *bidomain_system.rhs, tol, max_iter);
//        std::cout << "solving done" << std::endl;
	// WAVE
//	bidomain_system.solution->print();

	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	const libMesh::DofMap & dof_map = bidomain_system.get_dof_map();
	const libMesh::DofMap & dof_map_V = wave_system.get_dof_map();

	std::vector < libMesh::dof_id_type > dof_indices;
	std::vector < libMesh::dof_id_type > dof_indices_V;
	std::vector < libMesh::dof_id_type > dof_indices_Ve;
	std::vector < libMesh::dof_id_type > dof_indices_Q;

	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map.dof_indices(nn, dof_indices_Q, 0);
		dof_map.dof_indices(nn, dof_indices_Ve, 1);
		dof_map_V.dof_indices(nn, dof_indices_V, 0);
		wave_system.solution->set(dof_indices_V[0], (*bidomain_system.solution)(dof_indices_Q[0]) );
	}
	wave_system.solution->close();
	wave_system.solution->scale(dt);
	*wave_system.solution += *wave_system.old_local_solution;
//	double sol_norm = wave_system.solution->l2_norm();
//	std::cout <<"sol norm: " << sol_norm << std::endl;
//    wave_system.solution->print();


}

void Bidomain::reinit_linear_solver()
{
	M_linearSolver->clear();
	M_linearSolver->set_solver_type(libMesh::CG);
	M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
	M_linearSolver->init();
}

std::string Bidomain::get_ionic_model_name() const
{
	return M_ionicModelPtr->ionicModelName();
}

double Bidomain::potential_norm()
{
	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");
	return bidomain_system.solution->l1_norm();
}

const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Bidomain::get_fibers()
{
	ParameterSystem& fiber_system = M_equationSystems.get_system
			< ParameterSystem > ("fibers");
	return fiber_system.solution;
}

const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Bidomain::get_sheets()
{
	ParameterSystem& sheets_system = M_equationSystems.get_system
			< ParameterSystem > ("sheets");
	return sheets_system.solution;
}

const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Bidomain::get_xfibers()
{
	ParameterSystem& xfiber_system = M_equationSystems.get_system
			< ParameterSystem > ("xfibers");
	return xfiber_system.solution;
}

void Bidomain::set_potential_on_boundary(unsigned int boundID, double value)
{
	std::cout
			<< "* BIDOMAIN: WARNING:  set_potential_on_boundary works only for TET4"
			<< std::endl;
	libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	BidomainSystem& bidomain_system = M_equationSystems.get_system
			< BidomainSystem > ("bidomain");
	BidomainSystem& wave_system = M_equationSystems.get_system < BidomainSystem
			> ("wave");

	const libMesh::DofMap & dof_map = wave_system.get_dof_map();
	std::vector < libMesh::dof_id_type > dof_indices;

	libMesh::FEType fe_type = dof_map.variable_type(0);
	// Declare a special finite element object for
	// boundary integration.
	libMesh::UniquePtr < libMesh::FEBase
			> fe_face(libMesh::FEBase::build(dim, fe_type));

	// Boundary integration requires one quadraure rule,
	// with dimensionality one less than the dimensionality
	// of the element.
	libMesh::QGauss qface(dim - 1, libMesh::FIRST);

	// Tell the finite element object to use our
	// quadrature rule.
	fe_face->attach_quadrature_rule(&qface);
	libMesh::DenseVector<libMesh::Number> Fe;

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		dof_map.dof_indices(elem, dof_indices);
		Fe.resize(dof_indices.size());
		for (unsigned int side = 0; side < elem->n_sides(); side++)
		{
			if (elem->neighbor(side) == libmesh_nullptr)
			{
				const unsigned int boundary_id =
						mesh.boundary_info->boundary_id(elem, side);
				if (boundary_id == boundID)
				{
					wave_system.solution->set(dof_indices[0], value);
					wave_system.solution->set(dof_indices[1], value);
					wave_system.solution->set(dof_indices[2], value);
				}
			}

		}
	}
	wave_system.solution->close();

}

} /* namespace BeatIt */
