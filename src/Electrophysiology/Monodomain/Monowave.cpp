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
 * \file Monodomain.cpp
 *
 * \class Monodomain
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
#include "Electrophysiology/Monodomain/Monowave.hpp"
#include "Util/SpiritFunction.hpp"

#include "libmesh/discontinuity_measure.h"
#include "libmesh/fe_interface.h"

#include "Electrophysiology/Pacing/PacingProtocolSpirit.hpp"


namespace BeatIt
{

// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
typedef libMesh::TransientLinearImplicitSystem     MonodomainSystem;
typedef libMesh::TransientExplicitSystem           IonicModelSystem;
typedef libMesh::ExplicitSystem                     ParameterSystem;


//Monowave::Monodomain( libMesh::MeshBase & mesh )
//    : M_equationSystems( libMesh::EquationSystems(mesh) )
//    , M_monodomainExporter()
//    , M_monodomainExporterNames()
//    , M_ionicModelExporter()
//    , M_ionicModelExporterNames()
//    , M_parametersExporter()
//    , M_parametersExporterNames()
//    , M_outputFolder()
//    , M_datafile()
//    , M_pacing()
//    , M_linearSolver()
//{
//
//}

Monowave::Monowave( libMesh::EquationSystems& es )
    : M_equationSystems(es)
    , M_monodomainExporter()
    , M_monodomainExporterNames()
    , M_ionicModelExporter()
    , M_ionicModelExporterNames()
    , M_parametersExporter()
    , M_parametersExporterNames()
    , M_outputFolder()
    , M_datafile()
    , M_pacing()
    , M_linearSolver()
    , M_anisotropy(Anisotropy::Orthotropic)
	, M_equationType(EquationType::ReactionDiffusion)
	, M_artificialDiffusion(false)
	, M_penalty(false)
    , M_timeIntegratorType(DynamicTimeIntegratorType::Implicit)
{

}

Monowave::~Monowave(){}



void Monowave::setup(GetPot& data, std::string section)
{
    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
	evaluate_mesh_size();
    std::cout << "* MONODOMAIN: max mesh size:" << M_meshSize << std::endl;

    // Read Input File
    M_datafile = data;
    // setup time integrator type
    std::string time_integrator = M_datafile(section + "/time_integrator", "implicit");
    if ("implicit" == time_integrator)
        M_timeIntegratorType = DynamicTimeIntegratorType::Implicit;
    else
        M_timeIntegratorType = DynamicTimeIntegratorType::Explicit;

    //Read output folder from datafile
    std::string output_folder = M_datafile(section+"/output_folder",  "Output");
    M_outputFolder = "./" + output_folder + "/";

    std::cout << "* MONODOMAIN: creating pacing protocol" << std::endl;
    std::string pacing_type = M_datafile(section+"/pacing/type", "S1");
    M_pacing.reset( BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_type) );
    M_pacing->setup(M_datafile, "monodomain/pacing");

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // 1) ADR
    std::cout << "* MONODOMAIN: Creating new System for the monodomain diffusion reaction equation" << std::endl;
    MonodomainSystem& monodomain_system  =  M_equationSystems.add_system<MonodomainSystem>("monodomain");
    // TO DO: Generalize to higher order
    monodomain_system.add_variable( "Q", libMesh::FIRST);
    // Add 3 matrices
    monodomain_system.add_matrix("lumped_mass");
    monodomain_system.add_matrix("high_order_mass");
    // Add lumped mass matrix
    monodomain_system.add_vector("lumped_mass_vector");
    monodomain_system.add_matrix("mass");
    monodomain_system.add_matrix("stiffness");
//    monodomain_system.add_vector("ionic_currents");
    M_monodomainExporterNames.insert("monodomain");

//    monodomain_system.matrix->init();
    monodomain_system.init();

    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");
    wave_system.add_variable( "V", libMesh::FIRST);
    M_monodomainExporterNames.insert("wave");
    wave_system.init();

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // 2) ODEs
    std::cout << "* MONODOMAIN: Creating new System for the ionic model " << std::endl;
    IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");
    // Create Ionic Model
    std::string ionic_model = M_datafile(section+"/ionic_model", "NashPanfilov");
    M_ionicModelPtr.reset(  BeatIt::IonicModel::IonicModelFactory::Create(ionic_model) );
    M_ionicModelPtr->setup(M_datafile, section);
   int num_vars = M_ionicModelPtr->numVariables();
    // TO DO: Generalize to other conditions
    // We need to exclude the potential
    // therefore we loop up to num_vars-1
    for (int nv = 0; nv < num_vars-1; ++nv)
    {
        std::string var_name = M_ionicModelPtr->variableName(nv);
        // For the time being we use P1 for the variables
        ionic_model_system.add_variable( &var_name[0], libMesh::FIRST );
    }
    ionic_model_system.init();
    // Add the applied current to this system
    IonicModelSystem& Iion_system = M_equationSystems.add_system<IonicModelSystem>("iion");
    Iion_system.add_variable( "iion", libMesh::FIRST);
    Iion_system.add_vector( "diion");
    Iion_system.init();
    IonicModelSystem& istim_system = M_equationSystems.add_system<IonicModelSystem>("istim");
    istim_system.add_variable( "istim", libMesh::FIRST);
    istim_system.init();
    M_ionicModelExporterNames.insert("ionic_model");
    M_ionicModelExporterNames.insert("iion");
    M_ionicModelExporterNames.insert("istim");
    IonicModelSystem& cut_system = M_equationSystems.add_system<IonicModelSystem>("cut");
    cut_system.add_variable( "cut", libMesh::FIRST);
    cut_system.init();

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Distributed Parameters
    std::cout << "* MONODOMAIN: Creating parameters spaces " << std::endl;
    ParameterSystem& activation_times_system  = M_equationSystems.add_system<ParameterSystem>("activation_times");
    activation_times_system.add_variable( "activation_times", libMesh::FIRST);
    ParameterSystem& fiber_system        = M_equationSystems.add_system<ParameterSystem>("fibers");
    fiber_system.add_variable( "fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable( "fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable( "fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    ParameterSystem& sheets_system       = M_equationSystems.add_system<ParameterSystem>("sheets");
    sheets_system.add_variable( "sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheets_system.add_variable( "sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheets_system.add_variable( "sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    ParameterSystem& xfiber_system       = M_equationSystems.add_system<ParameterSystem>("xfibers");
    xfiber_system.add_variable( "xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.add_variable( "xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.add_variable( "xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    ParameterSystem& conductivity_system = M_equationSystems.add_system<ParameterSystem>("conductivity");
    conductivity_system.add_variable( "Dff", libMesh::CONSTANT, libMesh::MONOMIAL);
    conductivity_system.add_variable( "Dss", libMesh::CONSTANT, libMesh::MONOMIAL);
    conductivity_system.add_variable( "Dnn", libMesh::CONSTANT, libMesh::MONOMIAL);
    ParameterSystem& procID_system = M_equationSystems.add_system<ParameterSystem>("ProcID");
    procID_system.add_variable( "ProcID", libMesh::CONSTANT, libMesh::MONOMIAL);
     M_parametersExporterNames.insert("activation_times");
    M_parametersExporterNames.insert("fibers");
    M_parametersExporterNames.insert("sheets");
    M_parametersExporterNames.insert("xfibers");
    M_parametersExporterNames.insert("conductivity");
    M_parametersExporterNames.insert("ProcID");
    std::cout << "* MONODOMAIN: Initializing equation systems " << std::endl;
    // Initializing
	activation_times_system.init();
	fiber_system.init();
	sheets_system.init();
	xfiber_system.init();
	conductivity_system.init();
	procID_system.init();
//    M_equationSystems.init();
    M_equationSystems.print_info();

    // Conductivity Tensor in local coordinates
    // Fiber direction
    // Default Value = 1.3342  kOhm^-1 cm^-1
    double Dff = M_datafile( section+"/Dff", 1.3342);
    // Sheet direction
    // Default Value = 1.3342  kOhm^-1 cm^-1
    double Dss = M_datafile( section+"/Dss", 0.17606);
    // Cross fiber direction
    // Default Value = 1.3342  kOhm^-1 cm^-1
    double Dnn = M_datafile( section+"/Dnn", 0.17606);
    double Chi = M_datafile( section+"/Chi", 1400.0);
    M_equationSystems.parameters.set<double>("Chi") = Chi;
    double alpha = M_datafile( section+"/alpha", 1.0);
    M_equationSystems.parameters.set<double>("alpha") = alpha;
    std::string anisotropy = M_datafile(section+"/anisotropy", "orhotropic");
    std::map<std::string, Anisotropy> aniso_map;
    aniso_map["orthotropic"] = Anisotropy::Orthotropic;
    aniso_map["isotropic"] = Anisotropy::Isotropic;
    aniso_map["transverse"] = Anisotropy::TransverselyIsotropic;
    aniso_map["user"] = Anisotropy::UserDefined;
    M_anisotropy = aniso_map.find(anisotropy)->second;
    //if( Anisotropy::UserDefined == M_anisotropy)
    {
    	M_conductivity.push_back(Dff);
    	M_conductivity.push_back(Dss);
    	M_conductivity.push_back(Dnn);
        double Dns = M_datafile( section+"/Dns", 0.0);
        double Dnf = M_datafile( section+"/Dnf", 0.0);
        double Dsf = M_datafile( section+"/Dsf", 0.0);
    	M_conductivity.push_back(Dns);
    	M_conductivity.push_back(Dnf);
    	M_conductivity.push_back(Dsf);
	}
    std::cout << "* MONODOMAIN: Parameters: " << std::endl;
    std::cout << "              Chi = " << Chi << std::endl;
    std::cout << "              Dff = " << Dff << std::endl;
    std::cout << "              Dss = " << Dss << std::endl;
    std::cout << "              Dnn = " << Dnn << std::endl;
    std::cout << "              tau = " << alpha << std::endl;
    std::cout << "              anisotropy = " << anisotropy << std::endl;

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Setup Exporters
    std::cout << "* MONODOMAIN: Initializing exporters " << std::endl;
    M_monodomainExporter.reset( new Exporter(M_equationSystems.get_mesh() ) );
    M_ionicModelExporter.reset( new Exporter(M_equationSystems.get_mesh() ) );
    M_parametersExporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );
    M_monodomainEXOExporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );
    M_potentialEXOExporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );


    struct stat out_dir;
    if( stat(&M_outputFolder[0],&out_dir) != 0  )
    {
        if ( monodomain_system.get_mesh().comm().rank() == 0 )
        {
            mkdir ( M_outputFolder.c_str(), 0777 );
        }
    }


    M_equationType = EquationType::Wave;
    std::string eq_type = M_datafile(section+"/type", "wave");
    if(eq_type != "wave")
	{
    	M_equationType = EquationType::ReactionDiffusion;
        std::cout << "* MONODOMAIN type: RD " << std::endl;
	}
    else         std::cout << "* MONODOMAIN type: Wave " << std::endl;

    M_artificialDiffusion = M_datafile(section+"/artificial_diffusion", false);

    M_shellCorrection =  M_datafile(section+"/shell", false);
    std::cout << "* MONODOMAIN: Shell Correction: " << M_shellCorrection << std::endl;

}

void
Monowave::readFibers( EXOExporter& importer, int step)
{
//    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
//    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
//    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");
//
//    auto names = importer.get_elem_var_names();
//    std::cout << "* MONOWAVE: Importer available elem fields" << std::endl;
//    for(auto && n : names) std::cout << n << std::endl;
//    std::cout << "* MONOWAVE: Importing  fibers"<< std::endl;
////    importer.copy_elemental_solution(fiber_system, "fibersx", "fibersx", step);
////    importer.copy_elemental_solution(fiber_system, "fibersy", "fibersy", step);
////    importer.copy_elemental_solution(fiber_system, "fibersz", "fibersz", step);
////    std::cout << "* MONOWAVE: Importing  sheets"<< std::endl;
////    importer.copy_nodal_solution (sheets_system, "sheetsx", "sheetsx", step);
////    importer.copy_nodal_solution (sheets_system, "sheetsy", "sheetsy", step);
////    importer.copy_nodal_solution (sheets_system, "sheetsz", "sheetsz", step);
////    std::cout << "* MONOWAVE: Importing  xfibers"<< std::endl;
////    importer.copy_nodal_solution (xfiber_system, "xfibersx", "xfibersx",step);
////    importer.copy_nodal_solution (xfiber_system, "xfibersy", "xfibersy", step);
////    importer.copy_nodal_solution (xfiber_system, "xfibersz", "xfibersz", step);
//
//    // Map from element ID to elemental variable value.  We need to use
//    // a map here rather than a vector (e.g. elem_var_values) since the
//    // libmesh element numbering can contain "holes".  This is the case
//    // if we are reading elemental var values from an adaptively refined
//    // mesh that has not been sequentially renumbered.
//    std::map<libMesh::dof_id_type, libMesh::Real> elem_var_value_map;
//    importer.exio_helper->read_elemental_var_values("fibersx", step, elem_var_value_map);
//    std::cout << "* MONOWAVE: Importing  fibers 1"<< std::endl;
//
//    const unsigned int var_num = fiber_system.variable_number("fibersx");
//    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
//
//    std::map<libMesh::dof_id_type, libMesh::Real>::iterator
//      it = elem_var_value_map.begin(),
//      end = elem_var_value_map.end();
//    const libMesh::DofMap & dof_map_fibers= fiber_system.get_dof_map();
//    std::vector<libMesh::dof_id_type> dof_indices_fibers;
//
//    std::cout << "* MONOWAVE: Importing  fibers 2"<< std::endl;
//    for (; it!=end; ++it)
//      {
//        const libMesh::Elem * elem = mesh.query_elem_ptr(it->first);
//        std::cout << "el: " << it->first << ", " << elem << std::endl;
//        dof_map_fibers.dof_indices(elem, dof_indices_fibers);
//        std::cout << "* MONOWAVE: getting nc" << fiber_system.number() << ", " << var_num << std::endl;
//        //auto nc = elem->n_comp(fiber_system.number(), var_num) ;
//        std::cout << "* MONOWAVE: Done"<< std::endl;
//        //if (nc> 0)
//          {
//            std::cout << "* MONOWAVE: getting dof index"<< std::endl;
////            libMesh::dof_id_type dof_index = elem->dof_number(fiber_system.number(), var_num, 0);
//            std::cout << "* MONOWAVE: Done"<< std::endl;
//
//            // If the dof_index is local to this processor, set the value
////            if ((dof_index >= fiber_system.solution->first_local_index()) && (dof_index < fiber_system.solution->last_local_index()))
//                std::cout << "dof_index: " << dof_indices_fibers[0] << ", value: " << it->second << std::endl;
//                fiber_system.solution->set (dof_indices_fibers[0], it->second);
//          }
//      }
//    std::cout << "* MONOWAVE: Importing  fibers 3"<< std::endl;
//
//    fiber_system.solution->close();
//    std::cout << "* MONOWAVE: Importing  fibers 4"<< std::endl;
//
//    fiber_system.update();
}

void
Monowave::restart( EXOExporter& importer, int step, bool restart  )
{
	if(restart)
	{
		MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
		IonicModelSystem& ionic_model_system =  M_equationSystems.get_system<IonicModelSystem>("ionic_model");
		IonicModelSystem& wave_system =  M_equationSystems.get_system<IonicModelSystem>("wave");

	    auto names = importer.get_nodal_var_names();
	    std::cout << "* MONOWAVE: Importer available nodal fields" << std::endl;
	    for(auto && n : names) std::cout << n << std::endl;

		std::cout << "* MONOWAVE: Importing  V"<< std::endl;
		importer.copy_nodal_solution(wave_system, "V", step);
		std::cout << "* MONOWAVE: Importing  Q"<< std::endl;
		importer.copy_nodal_solution(monodomain_system, "Q", step);
		unsigned int n_vars = ionic_model_system.n_vars();
		std::cout << "* MONOWAVE: Importing  Ionic Model: num_vars: "<< n_vars <<  std::endl;
		for(unsigned int i = 0; i < n_vars; ++i)
		importer.copy_nodal_solution(ionic_model_system, ionic_model_system.variable_name(i), step);
	}
}

void
Monowave::init(double time)
{
    // Setting initial conditions
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    IonicModelSystem& ionic_model_system =  M_equationSystems.get_system<IonicModelSystem>("ionic_model");
    int num_vars = ionic_model_system.n_vars();
    std::vector<double> init_values(num_vars+1, 0.0);
    M_ionicModelPtr->initialize(init_values);
    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");
    auto first_local_index = wave_system.solution->first_local_index();
    auto last_local_index = wave_system.solution->last_local_index();
    for(auto i = first_local_index; i < last_local_index; ++i)
    {
        wave_system.solution->set(i,init_values[0]);
    }
    first_local_index = ionic_model_system.solution->first_local_index();
    last_local_index = ionic_model_system.solution->last_local_index();
    std::cout << "\n* MONODOMAIN: Setting initial values ... " << std::flush;
    for(auto i = first_local_index; i < last_local_index; )
    {
        for(int nv = 0; nv < num_vars; ++nv)
        {
            ionic_model_system.solution->set(i, init_values[nv+1]);
            ++i;
        }
    }
    std::cout << "done " << std::endl;

    std::string v_ic = M_datafile("monodomain/ic", "");
    if(v_ic != "")
    {
        std::cout << "* MONODOMAIN: Found monodomain initial condition: " << v_ic << std::endl;
        SpiritFunction monodomain_ic;
        monodomain_ic.read(v_ic);
        M_equationSystems.parameters.set<libMesh::Real>("time") =  time;
        monodomain_system.time = time;
        std::cout << "* MONODOMAIN: Projecting initial condition to monodomain system ... " << std::flush;
        wave_system.project_solution(&monodomain_ic);
        std::cout << " done" << std::endl;
    }
    std::cout << "* MONODOMAIN: Copying initial conditions to vectors at n nd at n-1... " << std::flush;

    // Close vectors and update the values in old_local_solution and older_local_solution
    monodomain_system.solution->close();
	monodomain_system.old_local_solution->close();
	monodomain_system.older_local_solution->close();
	ionic_model_system.solution->close();
	ionic_model_system.old_local_solution->close();
	ionic_model_system.older_local_solution->close();
	wave_system.solution->close();
	wave_system.old_local_solution->close();
	wave_system.older_local_solution->close();
    advance();
    std::cout << " done" << std::endl;

    IonicModelSystem& istim_system = M_equationSystems.get_system<IonicModelSystem>("istim");
    std::cout << "* MONODOMAIN: Initializing activation times to -1  ... " << std::flush;
    ParameterSystem& activation_times_system  = M_equationSystems.get_system<ParameterSystem>("activation_times");
    first_local_index = activation_times_system.solution->first_local_index();
    last_local_index = activation_times_system.solution->last_local_index();

    for(auto i = first_local_index; i < last_local_index; ++i)
    {
    	activation_times_system.solution->set(i,-1.0);
    }
    std::cout << " done" << std::endl;

    std::string fibers_data  = M_datafile("monodomain/fibers" , "1.0, 0.0, 0.0");
    std::string sheets_data  = M_datafile("monodomain/sheets" , "0.0, 1.0, 0.0");
    std::string xfibers_data = M_datafile("monodomain/xfibers", "0.0, 0.0, 1.0");

    SpiritFunction fibers_func;
    SpiritFunction sheets_func;
    SpiritFunction xfibers_func;

    fibers_func .read(fibers_data );
    sheets_func .read(sheets_data );
    xfibers_func.read(xfibers_data);

    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");

    std::cout << "* MONODOMAIN: Creating fibers from function: " <<  fibers_data << std::flush;
    fiber_system .project_solution(&fibers_func );
    std::cout << " done" << std::endl;
    std::cout << "* MONODOMAIN: Creating sheets from function: " <<  sheets_data << std::flush;
    sheets_system.project_solution(&sheets_func );
    std::cout << " done" << std::endl;
    std::cout << "* MONODOMAIN: Creating xfibers from function: " <<  xfibers_data << std::flush;
    xfiber_system.project_solution(&xfibers_func);
    std::cout << " done" << std::endl;


    ParameterSystem& conductivity_system = M_equationSystems.get_system<ParameterSystem>("conductivity");
    std::string Dff_data  = M_datafile("monodomain/Dff" , "1.3342");
    std::string Dss_data  = M_datafile("monodomain/Dss" , "0.17606");
    std::string Dnn_data  = M_datafile("monodomain/Dnn" , "0.17606");
    SpiritFunction conductivity_func;
    conductivity_func.add_function(Dff_data);
    conductivity_func.add_function(Dss_data);
    conductivity_func.add_function(Dnn_data);
    conductivity_system .project_solution(&conductivity_func );

    std::map<std::string, libMesh::SolverType> solver_map;
    solver_map["cg"] = libMesh::CG;
    solver_map["cgs"] = libMesh::CGS;
    solver_map["gmres"] = libMesh::GMRES;
    std::map<std::string, libMesh::PreconditionerType> prec_map;
    prec_map["jacobi"] =  libMesh::JACOBI_PRECOND;
    prec_map["sor"] =  libMesh::SOR_PRECOND;
    prec_map["ssor"] =  libMesh::SSOR_PRECOND;
    prec_map["cholesky"] =  libMesh::CHOLESKY_PRECOND;
    prec_map["lu"] =  libMesh::LU_PRECOND;
    prec_map["ilu"] =  libMesh::ILU_PRECOND;
    prec_map["amg"] =  libMesh::AMG_PRECOND;

    std::string solver_type = M_datafile("monodomain/linear_solver/type", "gmres");
    std::cout << "* MONODOMAIN: using " << solver_type << std::endl;
    std::string prec_type = M_datafile("monodomain/linear_solver/preconditioner", "amg");
    std::cout << "* MONODOMAIN: using " << prec_type << std::endl;
     M_linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( M_equationSystems.comm() );
    M_linearSolver->set_solver_type(solver_map.find(solver_type)->second);
    M_linearSolver->set_preconditioner_type(prec_map.find(prec_type)->second);
    M_linearSolver->init();

    ParameterSystem& procID_system = M_equationSystems.get_system<ParameterSystem>("ProcID");
    first_local_index = procID_system.solution->first_local_index();
    last_local_index = procID_system.solution->last_local_index();

    double myrank = static_cast<double>(M_equationSystems.comm().rank());
    for(int el = first_local_index; el < last_local_index; ++el)
    {
    	procID_system.solution->set(el, myrank);
    }
    procID_system.solution->close();

}

void
Monowave::cut(double time, std::string f)
{
    std::cout << "* MONODOMAIN: Cutting potential" << std::endl;

	SpiritFunction func;
	func.read(f);
    // Add the applied current to this system
    IonicModelSystem& cut_system = M_equationSystems.get_system<IonicModelSystem>("cut");
    cut_system.time = time;
    cut_system.project_solution( &func );
    cut_system.solution->close();
//    cut_system.update();
     IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");
    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");
	MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");

  	auto first = cut_system.solution->first_local_index();
	auto last  = cut_system.solution->last_local_index();
	int n_vars = ionic_model_system.n_vars();
	std::vector<double> init_val(n_vars+1, 0.0);
	double cut_value = 0.0;

	M_ionicModelPtr->initialize(init_val);

	ionic_model_system.solution->close();
	wave_system.solution->close();
	monodomain_system.solution->close();
	    cut_system.solution->close();

    for(int index = first; index < last; index++ )
	{
    	cut_value = (*cut_system.solution)(index);
    	if(cut_value < 0.1)
    	{
    		if(M_equationType == EquationType::ReactionDiffusion)
    		{
    			wave_system.solution->set(index, init_val[0]);
    			monodomain_system.solution->set(index, init_val[0]);
    		}
    		else
    		{
    			wave_system.solution->set(index, init_val[0]);
    			monodomain_system.solution->set(index,0.0);
    		}
//       std::cout << "ion value" << std::endl;
			for( int m = 0; m< n_vars; m++ )
			{
				int   var_index =  index  * n_vars + m;
						ionic_model_system.solution->set(var_index, init_val[m+1]);
			}
    	}
	}

    	ionic_model_system.solution->close();
		wave_system.solution->close();
		monodomain_system.solution->close();
	    cut_system.solution->close();

}


void
Monowave::update_pacing(double time)
{
    // Add the applied current to this system
    IonicModelSystem& istim_system = M_equationSystems.get_system<IonicModelSystem>("istim");
    istim_system.time = time;
    M_pacing->update(time);
    istim_system.project_solution( &M_pacing->pacing() );
}

void
Monowave::update_activation_time( double  time, double threshold )
{
	ParameterSystem& activation_times_system  = M_equationSystems.add_system<ParameterSystem>("activation_times");
    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");

	MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");

	auto first = wave_system.solution->first_local_index();
	auto last  = wave_system.solution->last_local_index();
	double v = 0.0;
	double at = 0.0;
	for(int index = first; index < last; index++ )
	{
		v = (*wave_system.solution)(index);
		at = (*activation_times_system.solution)(index);
		if( v > threshold && at < 0.0)
		{
			activation_times_system.solution->set(index, time);
		}
	}
	activation_times_system.solution->close();
}


void
Monowave::generate_fibers(   const GetPot& data,
                               const std::string& section)
{
    std::cout << "* MONODOMAIN: Creating fiber fields" << std::endl;

    std::cout << "* MONODOMAIN: Solving Poisson porblem" << std::endl;

    libMesh::Mesh new_mesh( dynamic_cast< libMesh::Mesh&>(M_equationSystems.get_mesh() ) );
    libMesh::EquationSystems es(new_mesh);
    Poisson p(es);
    p.setup(data, section);
    p.assemble_system();
    p.solve_system();
    p.compute_elemental_solution_gradient();
    p.save_exo();

    const auto& sol_ptr = p.get_P0_solution();


    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");

    *sheets_system.solution = *p.get_gradient();

    auto first_local_index = fiber_system.solution->first_local_index();
    auto last_local_index = fiber_system.solution->last_local_index();
    auto first_local_index_sol = sol_ptr->first_local_index();
    auto last_local_index_sol = sol_ptr->last_local_index();

    double norm = 0.0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    double cx = data (section+"/centerline_x", 0.0);
    double cy = data (section+"/centerline_y", 0.0);
    double cz = data (section+"/centerline_z", 1.0);
    double cdot = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    double xfx = 0.0;
    double xfy = 0.0;
    double xfz = 0.0;
    double epi_angle = data (section+"/epi_angle", -60.0);
    double endo_angle = data (section+"/endo_angle", 60.0);
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
    double sa  = 0.0;
    double sa2 = 0.0;
    double teta1 = 0.0;
    double teta2 = 0.0;
    double teta  = 0.0;
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

    std::cout << "* MONODOMAIN: Computing rule-based fiber fields" << std::endl;

//    sol_ptr->print();

    auto j = first_local_index_sol;

    for(auto i = first_local_index; i < last_local_index; )
    {

//        std::cout << "* MONODOMAIN: getting sol ... " << std::flush;
        potential = (*sol_ptr)(j);
        j++;
//        std::cout << " done" << std::endl;

        sx = (*sheets_system.solution)(i);
        sy = (*sheets_system.solution)(i+1);
        sz = (*sheets_system.solution)(i+2);
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
        m = (teta1 - teta2 );
        q = teta2;
        teta =  m * potential + q;


        //*************************************************************//
        // The fiber field F is a rotation of the flat fiber field f
        // F = R f
        // where R is the rotation matrix.
        // To compute R we need the sin(teta) and
        // the sin(teta)^2 and the cross-product matrix W (check
        // rodrigues formula on wikipedia :) )
        //*************************************************************//
        sa = std::sin (teta);
        sa2 = 2.0 * std::sin (0.5 * teta) * std::sin (0.5 * teta);

        W11 = 0.0; W12 = -sz; W13 =  sy;
        W21 =  sz; W22 = 0.0; W23 = -sx;
        W31 = -sy; W32 =  sx; W33 = 0.0;
        //
        R11 = 1.0 + sa * W11 + sa2 * ( sx * sx - 1.0 );
        R12 = 0.0 + sa * W12 + sa2 * ( sx * sy );
        R13 = 0.0 + sa * W13 + sa2 * ( sx * sz );
        R21 = 0.0 + sa * W21 + sa2 * ( sy * sx );
        R22 = 1.0 + sa * W22 + sa2 * ( sy * sy - 1.0 );
        R23 = 0.0 + sa * W23 + sa2 * ( sy * sz );
        R31 = 0.0 + sa * W31 + sa2 * ( sz * sx );
        R32 = 0.0 + sa * W32 + sa2 * ( sz * sy );
        R33 = 1.0 + sa * W33 + sa2 * ( sz * sz - 1.0 );

        f0x =   R11 * fx + R12 * fy + R13 * fz;
        f0y =   R21 * fx + R22 * fy + R23 * fz;
        f0z =   R31 * fx + R32 * fy + R33 * fz;
        normalize(f0x, f0y, f0z, 1.0, 0.0, 0.0);


        xfx = f0y * sz - f0z * sy;
        xfy = f0z * sx - f0x * sz;
        xfz = f0x * sy - f0y * sx;
        normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

        sheets_system.solution->set(i,   sx);
        sheets_system.solution->set(i+1, sy);
        sheets_system.solution->set(i+2, sz);
        xfiber_system.solution->set(i,   xfx);
        xfiber_system.solution->set(i+1, xfy);
        xfiber_system.solution->set(i+2, xfz);
//        std::cout << "* MONODOMAIN: setting sol ... " << std::flush;

        fiber_system. solution->set(i,   f0x);
        fiber_system. solution->set(i+1, f0y);
        fiber_system. solution->set(i+2, f0z);
//        std::cout << " done" << std::endl;

        i += 3;
    }

}

double
Monowave::last_activation_time()
{
    ParameterSystem& activation_times_system  = M_equationSystems.add_system<ParameterSystem>("activation_times");
    return static_cast<double>(activation_times_system.solution->max());
}


void
Monowave::save(int step)
{
    std::cout << "* MONODOMAIN: GMVIO::Exporting " << step << " in: "  << M_outputFolder << " ... " << std::flush;
//    M_equationSystems->
    M_monodomainExporter->write_equation_systems ( M_outputFolder+"monodomain.gmv."+std::to_string(step),
                                                   M_equationSystems,
                                                   &M_monodomainExporterNames);
    M_ionicModelExporter->write_equation_systems ( M_outputFolder+"ionic_model.gmv."+std::to_string(step),
                                                   M_equationSystems,
                                                   &M_ionicModelExporterNames);
    std::cout << "done " << std::endl;
}
void
Monowave::init_exo_output()
{
    M_monodomainEXOExporter->write_equation_systems (M_outputFolder+"monodomain.exo", M_equationSystems);
    M_monodomainEXOExporter->append(true);
    std::vector<std::string> output;
    output.push_back("V");
    output.push_back("Q");
    M_potentialEXOExporter->set_output_variables(output);
    M_potentialEXOExporter->write_equation_systems (M_outputFolder+"potential.exo", M_equationSystems);
    M_potentialEXOExporter->append(true);
}

void
Monowave::save_exo(int step, double time)
{
    std::cout << "* MONODOMAIN: EXODUSII::Exporting monodomain.exo at time "   << time << " in: "  << M_outputFolder << " ... " << std::flush;
    M_monodomainEXOExporter->write_timestep(  M_outputFolder+"monodomain.exo"
                                         , M_equationSystems
                                         , step, time );
    std::cout << "done " << std::endl;
}

void
Monowave::save_potential(int step, double time)
{
    std::cout << "* MONODOMAIN: EXODUSII::Exporting potential.exo at time "   << time << " in: "  << M_outputFolder << " ... " << std::flush;
    if(M_equationType == EquationType::Wave) std::cout << " wave solver ... " << std::flush;
    else std::cout << " adr solver ... " << std::flush;

    M_potentialEXOExporter->write_timestep(  M_outputFolder+"potential.exo"
                                         , M_equationSystems
                                         , step, time );
//    std::cout << "* MONODOMAIN: GMVIO::Exporting " << step << " in: "  << M_outputFolder << " ... " << std::flush;
//    M_monodomainExporter->write_equation_systems ( M_outputFolder+"monodomain.gmv."+std::to_string(step),
//                                                   M_equationSystems,
//                                                   &M_monodomainExporterNames);
    std::cout << "done " << std::endl;
}

void
Monowave::amr( libMesh:: MeshRefinement& mesh_refinement, const std::string& type)
{
//	std::cout << "* MONODOMAIN: starting AMR ... " << std::flush;
//	BeatIt:Timer timer;
//	std::cout << "Creating Error Vector " << std::endl;
//	timer.start();
	libMesh::ErrorVector error;
//	timer.stop();
//	timer.print(std::cout);
//	std::cout << "Creating Error Estimator " << std::endl;
//	timer.restart();
	libMesh::ErrorEstimator * p_error_estimator;
	if("kelly" == type) p_error_estimator = new  libMesh::KellyErrorEstimator ;
	else if ("disc" == type) p_error_estimator = new libMesh::DiscontinuityMeasure;
	else p_error_estimator = new  libMesh::LaplacianErrorEstimator ;
//	timer.stop();
//	timer.print(std::cout);
//	 libMesh::KellyErrorEstimator error_estimator;
//	 libMesh::LaplacianErrorEstimator error_estimator;
	MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
	    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");

//	std::cout << "Estimating Monodomain Error  " << std::endl;
//	timer.restart();
	p_error_estimator->estimate_error(wave_system, error);
//	timer.stop();
//	timer.print(std::cout);
    // Flag elements to be refined and coarsened
//	std::cout << "Flagging Elements  " << std::endl;
//	timer.restart();
    mesh_refinement.flag_elements_by_error_fraction (error);
//	timer.stop();
//	timer.print(std::cout);
//
//	std::cout << "Refine and Coarsen  " << std::endl;
//	std::cout << " coarsen and refine ...  " << std::flush;
//	timer.restart();
    mesh_refinement.refine_and_coarsen_elements();
//	timer.stop();
//	timer.print(std::cout);
//	std::cout << " reinit ...  " << std::flush;
//	std::cout << "Reinit system  " << std::endl;
//	timer.restart();
    M_equationSystems.reinit();
//	timer.stop();
//	timer.print(std::cout);
//	timer.restart();
//	std::cout << " done  " << std::endl;
    delete p_error_estimator;
}


void
Monowave::save_parameters()
{
    std::cout << "* MONODOMAIN: EXODUSII::Exporting parameters in: "  << M_outputFolder << " ... " << std::flush;
    M_parametersExporter->write_equation_systems (M_outputFolder+"parameters.exo", M_equationSystems, &M_parametersExporterNames);
    M_parametersExporter->write_element_data(M_equationSystems);
    std::cout << "done " << std::endl;
}

void
Monowave::evaluate_mesh_size()
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
void
Monowave::assemble_matrices()
{
//    std::cout << "* MONODOMAIN: Assembling matrices ... " << std::flush;
     using libMesh::UniquePtr;

     const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
     const unsigned int dim = mesh.mesh_dimension();
     const unsigned int max_dim = 3;
     const libMesh::Real Chi = M_equationSystems.parameters.get<libMesh::Real> ("Chi");

     // Get a reference to the LinearImplicitSystem we are solving
     MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
     IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");

     monodomain_system.get_matrix("mass").zero();
     monodomain_system.get_matrix("lumped_mass").zero();
     monodomain_system.get_matrix("high_order_mass").zero();
     monodomain_system.get_matrix("stiffness").zero();
     monodomain_system.get_vector("lumped_mass_vector").zero();



//     MatSetOption( (dynamic_cast<libMesh::PetscMatrix<libMesh::Number> >(monodomain_system.get_matrix("stiffness"))).mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//     MatSetOption( (dynamic_cast<libMesh::PetscMatrix<libMesh::Number> * >(monodomain_system.matrix))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);



     ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
     ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
     ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");
     ParameterSystem& conductivity_system = M_equationSystems.get_system<ParameterSystem>("conductivity");

     // A reference to the  DofMap object for this system.  The  DofMap
     // object handles the index translation from node and element numbers
     // to degree of freedom numbers.  We will talk more about the  DofMap
     // in future examples.
     const libMesh::DofMap & dof_map_monodomain = monodomain_system.get_dof_map();
     const libMesh::DofMap & dof_map_fibers= fiber_system.get_dof_map();

     // Get a constant reference to the Finite Element type
     // for the first (and only) variable in the system.
     libMesh::FEType fe_type_qp1 = dof_map_monodomain.variable_type(0);
     libMesh::FEType fe_type_qp2 = dof_map_monodomain.variable_type(0);

     // Build a Finite Element object of the specified type.  Since the
     // FEBase::build() member dynamically creates memory we will
     // store the object as a UniquePtr<FEBase>.  This can be thought
     // of as a pointer that will clean up after itself.  Introduction Example 4
     // describes some advantages of  UniquePtr's in the context of
     // quadrature rules.
     UniquePtr<libMesh::FEBase> fe_qp1(libMesh::FEBase::build(dim, fe_type_qp1));
     UniquePtr<libMesh::FEBase> fe_qp2(libMesh::FEBase::build(dim, fe_type_qp2));

     // A 5th order Gauss quadrature rule for numerical integration.
     libMesh::QGauss qrule_stiffness(dim, libMesh::SECOND);
     // A 5th order Gauss quadrature rule for numerical integration.
     libMesh::QGauss qrule_mass(dim, libMesh::THIRD);

     // Tell the finite element object to use our quadrature rule.
     fe_qp1->attach_quadrature_rule(&qrule_stiffness);
     fe_qp2->attach_quadrature_rule(&qrule_mass);
     // Here we define some references to cell-specific data that
     // will be used to assemble the linear system.
     //
     // The element Jacobian * quadrature weight at each integration point.
     const std::vector<libMesh::Real> & JxW_qp1 = fe_qp1->get_JxW();
     const std::vector<libMesh::Real> & JxW_qp2 = fe_qp2->get_JxW();

     // The physical XY locations of the quadrature points on the element.
     // These might be useful for evaluating spatially varying material
     // properties at the quadrature points.
     const std::vector<libMesh::Point> & q_point_qp1 = fe_qp1->get_xyz();
     const std::vector<libMesh::Point> & q_point_qp2 = fe_qp2->get_xyz();

     // The element shape functions evaluated at the quadrature points.
     const std::vector<std::vector<libMesh::Real> > & phi_qp1 = fe_qp1->get_phi();
     const std::vector<std::vector<libMesh::Real> > & phi_qp2 = fe_qp2->get_phi();

     // The element shape function gradients evaluated at the quadrature
     // points.
     const std::vector<std::vector<libMesh::RealGradient> > & dphi_qp1 = fe_qp1->get_dphi();
     const std::vector<std::vector<libMesh::RealGradient> > & dphi_qp2 = fe_qp2->get_dphi();

      const std::vector<std::vector<libMesh::Real> > & dphidx_qp1 = fe_qp1->get_dphidx();
      const std::vector<std::vector<libMesh::Real> > & dphidy_qp1 = fe_qp1->get_dphidy();
      const std::vector<std::vector<libMesh::Real> > & dphidz_qp1 = fe_qp1->get_dphidz();




     // Define data structures to contain the element matrix
     // and right-hand-side vector contribution.  Following
     // basic finite element terminology we will denote these
     // "Ke" and "Fe".  These datatypes are templated on
     //  Number, which allows the same code to work for real
     // or complex numbers.
     libMesh::DenseMatrix<libMesh::Number> Ke;
     libMesh::DenseMatrix<libMesh::Number> Me;
     libMesh::DenseMatrix<libMesh::Number> Mel;
     libMesh::DenseVector<libMesh::Number> Fe;



     // for interior penalty
     UniquePtr<libMesh::FEBase> fe_elem_face(libMesh::FEBase::build(dim, fe_type_qp1));
     UniquePtr<libMesh::FEBase> fe_neighbor_face(libMesh::FEBase::build(dim, fe_type_qp1));
     // Tell the finite element object to use our quadrature rule.
     libMesh::QGauss qface(dim-1, fe_type_qp1.default_quadrature_order());

     fe_elem_face->attach_quadrature_rule(&qface);
     fe_neighbor_face->attach_quadrature_rule(&qface);
     // Data for surface integrals on the element boundary
     const std::vector<std::vector<libMesh::Real> > &  phi_face = fe_elem_face->get_phi();
     const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_elem_face->get_dphi();
     const std::vector<libMesh::Real> & JxW_face = fe_elem_face->get_JxW();
     const std::vector<libMesh::Point> & qface_normals = fe_elem_face->get_normals();
     const std::vector<libMesh::Point> & qface_points = fe_elem_face->get_xyz();
     // Data for surface integrals on the neighbor boundary
     const std::vector<std::vector<libMesh::Real> > &  phi_neighbor_face = fe_neighbor_face->get_phi();
     const std::vector<std::vector<libMesh::RealGradient> > & dphi_neighbor_face = fe_neighbor_face->get_dphi();
     // Data structures to contain the element and neighbor boundary matrix
     // contribution. This matrices will do the coupling beetwen the dofs of
     // the element and those of his neighbors.
     // Ken: matrix coupling elem and neighbor dofs
     libMesh::DenseMatrix<libMesh::Number> Kne;
     libMesh::DenseMatrix<libMesh::Number> Ken;
     libMesh::DenseMatrix<libMesh::Number> Kee;
     libMesh::DenseMatrix<libMesh::Number> Knn;

     // This vector will hold the degree of freedom indices for
     // the element.  These define where in the global system
     // the element degrees of freedom get mapped.
     std::vector<libMesh::dof_id_type> dof_indices;
     std::vector<libMesh::dof_id_type> dof_indices_fibers;

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


     // Loop over the elements.  Note that  ++el is preferred to
     // el++ since the latter requires an unnecessary temporary
     // object.
     double f0[3];
     double s0[3];
     double n0[3];
     libMesh::RealGradient DgradV;
     libMesh::TensorValue<double> D0;

     for (; el != end_el; ++el)
     {
         const libMesh::Elem * elem = *el;
         const unsigned int elem_id = elem->id();
         dof_map_monodomain.dof_indices(elem, dof_indices);
         dof_map_fibers.dof_indices(elem, dof_indices_fibers);

         // Compute the element-specific data for the current
         // element.  This involves computing the location of the
         // quadrature points (q_point) and the shape functions
         // (phi, dphi) for the current element.
         fe_qp1->reinit(elem);
         fe_qp2->reinit(elem);


         // Zero the element matrix and right-hand side before
         // summing them.  We use the resize member here because
         // the number of degrees of freedom might have changed from
         // the last element.  Note that this will be the case if the
         // element type is different (i.e. the last element was a
         // triangle, now we are on a quadrilateral).

         // The  DenseMatrix::resize() and the  DenseVector::resize()
         // members will automatically zero out the matrix  and vector.
         auto n_dofs = dof_indices.size();
         Ke.resize(n_dofs, n_dofs);
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
         double Dff = (*conductivity_system.solution)(dof_indices_fibers[0]);
         double Dss = (*conductivity_system.solution)(dof_indices_fibers[1]);
         double Dnn = (*conductivity_system.solution)(dof_indices_fibers[2]);

         if(M_artificialDiffusion)
         {
        	 double delta = 0.0;
             double h = elem->hmax();
             double hf, hs, hn;
             double hfmin, hsmin, hnmin;
             double hfmax, hsmax, hnmax;
             int n_points = elem->n_nodes();
             for(int np = 0; np < n_points; np++)
             {
            	 auto p1 = elem->point(np);
            	 hfmax = std::max(hfmax, p1(0)*f0[0]+p1(1)*f0[1]+p1(2)*f0[2]);
            	 hfmin = std::min(hfmin, p1(0)*f0[0]+p1(1)*f0[1]+p1(2)*f0[2]);
            	 hsmax = std::max(hsmax, p1(0)*s0[0]+p1(1)*s0[1]+p1(2)*s0[2]);
            	 hsmin = std::min(hsmin, p1(0)*s0[0]+p1(1)*s0[1]+p1(2)*s0[2]);
            	 hnmax = std::max(hnmax, p1(0)*n0[0]+p1(1)*n0[1]+p1(2)*n0[2]);
            	 hnmin = std::min(hnmin, p1(0)*n0[0]+p1(1)*n0[1]+p1(2)*n0[2]);
             }
             hf = hfmax - hfmin;
             hs = hsmax - hsmin;
             hn = hnmax - hnmin;
        	 Dff += delta * hf * std::sqrt(Dff/Dff);
        	 Dss += delta * hs * std::sqrt(Dff/Dss);
        	 Dnn += delta * hn * std::sqrt(Dff/Dnn);
         }

		 switch(M_anisotropy)
		 {
			 case Anisotropy::Isotropic:
			 {
				 for(int idim = 0; idim < max_dim; ++idim )
				 {
						  D0(idim,idim) = Dff;
				 }

				 break;
			 }

			 case Anisotropy::TransverselyIsotropic:
			 {
				 for(int idim = 0; idim < max_dim; ++idim )
				 {
					 for(int jdim = 0; jdim < max_dim; ++jdim )
					 {

										  D0(idim,jdim) = ( Dff - Dss ) * f0[idim] * f0[jdim];
					 }
					 D0(idim,idim) += Dss;
				 }
				 break;
			 }
			 case Anisotropy::UserDefined:
			 {
    			 D0(0, 0) = M_conductivity[0]; //Dff
    			 D0(1, 1) = M_conductivity[1]; //Dss
    			 D0(2, 2) = M_conductivity[2]; //Dnn
    			 D0(2, 1) = M_conductivity[3]; //Dns
    			 D0(1, 2) = M_conductivity[3]; //Dns
    			 D0(0, 2) = M_conductivity[4]; //Dnf
    			 D0(2, 0) = M_conductivity[4]; //Dnf
    			 D0(0, 1) = M_conductivity[5]; //Dsf
    			 D0(1, 0) = M_conductivity[5]; //Dsf
				 break;
			 }


			 case Anisotropy::Orthotropic:
			 default:
			 {
				 for(int idim = 0; idim < max_dim; ++idim )
				 {
					 for(int jdim = 0; jdim < max_dim; ++jdim )
					 {
										  D0(idim,jdim) = Dff * f0[idim] * f0[jdim]
									   + Dss * s0[idim] * s0[jdim]
									   + Dnn * n0[idim] * n0[jdim];

					 }
				 }
				 break;
			 }
		 }
		 D0 /= Chi;


//    	 std::cout << "f  = [" << f0[0] << "," << f0[1]  << ", " << f0[2]  << "]"<< std::endl;
//    	 std::cout << "s = [" << s0[0] << "," << s0[1]  << ", " << s0[2]  << "]"<< std::endl;
//    	 std::cout << "n = [" << n0[0] << "," << n0[1]  << ", " << n0[2]  << "]"<< std::endl;
//	 	 std::cout << "D0 = " << std::endl;
//		 for(int idim = 0; idim < dim; ++idim )
//         {
//        	 std::cout << D0(idim, 0) << "," << D0(idim, 1) << ", " << D0(idim, 2) << std::endl;
//         }

         // Assemble Mass terms
         for (unsigned int qp = 0; qp < qrule_mass.n_points(); qp++)
         {
             //  Matrix
             for (unsigned int i = 0; i < phi_qp2.size(); i++)
             {
                 for (unsigned int j = 0; j < phi_qp2.size(); j++)
                 {
                     // Mass term
                     Me(i, j) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                     Mel(i, i) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                     Fe(i) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                 }
             }
         }
         monodomain_system.get_matrix("mass").add_matrix(Me, dof_indices);
         monodomain_system.get_matrix("lumped_mass").add_matrix(Mel, dof_indices);
         monodomain_system.get_vector("lumped_mass_vector").add_vector(Fe, dof_indices);
         for (unsigned int qp = 0; qp < qrule_stiffness.n_points(); qp++)
         {
             for (unsigned int i = 0; i < dphi_qp1.size(); i++)
             {
                 DgradV = D0  * dphi_qp1[i][qp];

//                 std::cout << "DgradV = " << std::endl;
//        		 for(int idim = 0; idim < dim; ++idim )
//                 {
//                	 std::cout << DgradV(idim, 0) << "," << DgradV(idim, 1) << ", " << DgradV(idim, 2) << std::endl;
//                 }

                 for (unsigned int j = 0; j < dphi_qp1.size(); j++)
                 {
                     // stiffness term
                     Ke(i, j) += JxW_qp1[qp] * DgradV * dphi_qp1[j][qp];
                 }
             }
         }
         monodomain_system.get_matrix("stiffness").add_matrix(Ke, dof_indices);


		// If the element has no neighbor on a side then that
		// side MUST live on a boundary of the domain.
		for (unsigned int side=0; side<elem->n_sides(); side++)
		{
			double gamma = 1.0;
			if ( elem->neighbor_ptr(side))
			{
	              // Store a pointer to the neighbor we are currently
	              // working on.
	              const libMesh::Elem * neighbor = elem->neighbor_ptr(side);
	              // Get the global id of the element and the neighbor
	              const unsigned int neighbor_id = neighbor->id();

	              // WARNING!!!! NOTE!!!!
	              // Here I should use some check for amr:
	              // check libmesh test: miscellaneous_ex5
	              if(false)
	              {
	                  // Pointer to the element side
	                  UniquePtr<const libMesh::Elem> elem_side (elem->build_side_ptr(side));

	                  // h dimension to compute the interior penalty penalty parameter
	                  const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
	                  const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face->get_order());
	                  const double side_order = (elem_b_order + neighbor_b_order)/2.;
	                  const double penalty = 0.5*gamma * elem_side->volume() * elem_side->volume();
	                  const double h_elem = (elem->volume()/elem_side->volume())/std::pow(side_order,2.);

	                  // The quadrature point locations on the neighbor side
	                  std::vector<libMesh::Point> qface_neighbor_point;

	                  // The quadrature point locations on the element side
	                  std::vector<libMesh::Point > qface_point;

	                  // Reinitialize shape functions on the element side
	                  fe_elem_face->reinit(elem, side);

	                  // Get the physical locations of the element quadrature points
	                  qface_point = fe_elem_face->get_xyz();

	                  // Find their locations on the neighbor
	                  unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);
	                  libMesh::FEInterface::inverse_map ( elem->dim(),
	                		  	  	  	  	  	  	  	  fe_qp1->get_fe_type(),
														  neighbor,
														  qface_point,
														  qface_neighbor_point);
	                  // Calculate the neighbor element shape functions at those locations
	                  fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);
	                  // Get the degree of freedom indices for the
	                  // neighbor.  These define where in the global
	                  // matrix this neighbor will contribute to.
	                  std::vector<libMesh::dof_id_type> neighbor_dof_indices;
	                  dof_map_monodomain.dof_indices (neighbor, neighbor_dof_indices);
	                  const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
	                  // Zero the element and neighbor side matrix before
	                  // summing them.  We use the resize member here because
	                  // the number of degrees of freedom might have changed from
	                  // the last element or neighbor.
	                  // Note that Kne and Ken are not square matrices if neighbor
	                  // and element have a different p level
	                  Kne.resize (n_neighbor_dofs, n_dofs);
	                  Ken.resize (n_dofs, n_neighbor_dofs);
	                  Kee.resize (n_dofs, n_dofs);
	                  Knn.resize (n_neighbor_dofs, n_neighbor_dofs);


	                  // Now we will build the element and neighbor
	                  // boundary matrices.  This involves
	                  // a double loop to integrate the test funcions
	                  // (i) against the trial functions (j).
	                  for (unsigned int qp=0; qp<qface.n_points(); qp++)
	                    {
	                      // Kee Matrix. Integrate the element test function i
	                      // against the element test function j
	                      for (unsigned int i=0; i<n_dofs; i++)
	                        {
	                          for (unsigned int j=0; j<n_dofs; j++)
	                            {
	                              // consistency
	                              Kee(i,j) += penalty * JxW_face[qp] *
										    ( qface_normals[qp]*dphi_face[i][qp] ) *
											( qface_normals[qp]*dphi_face[j][qp] );
	                              // stability
	                              // Kee(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];
	                            }
	                        }

	                      // Knn Matrix. Integrate the neighbor test function i
	                      // against the neighbor test function j
	                      for (unsigned int i=0; i<n_neighbor_dofs; i++)
	                        {
	                          for (unsigned int j=0; j<n_neighbor_dofs; j++)
	                            {
	                              // consistency
	                              Kee(i,j) += penalty * JxW_face[qp] *
										    ( qface_normals[qp]*dphi_neighbor_face[i][qp] ) *
											( qface_normals[qp]*dphi_neighbor_face[j][qp] );

//	                              Knn(i,j) +=
//	                            		  penalty * JxW_face[qp] *
//	                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) +
//	                                 phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

	                              // stability
	                              // Knn(i,j) +=
	                                //JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp];
	                            }
	                        }

	                      // Kne Matrix. Integrate the neighbor test function i
	                      // against the element test function j
	                      for (unsigned int i=0; i<n_neighbor_dofs; i++)
	                        {
	                          for (unsigned int j=0; j<n_dofs; j++)
	                            {
	                              Kee(i,j) -= penalty * JxW_face[qp] *
										    ( qface_normals[qp]*dphi_neighbor_face[i][qp] ) *
											( qface_normals[qp]*dphi_face[j][qp] );
	                              // consistency
//	                               Kne(i,j) +=
//	                            		   penalty * JxW_face[qp] *
//	                                (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) -
//	                                 phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));

	                              // stability
	                              //Kne(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];
	                            }
	                        }

	                      // Ken Matrix. Integrate the element test function i
	                      // against the neighbor test function j
	                      for (unsigned int i=0; i<n_dofs; i++)
	                        {
	                          for (unsigned int j=0; j<n_neighbor_dofs; j++)
	                            {
	                              // consistency
	                              Kee(i,j) -= penalty * JxW_face[qp] *
										    ( qface_normals[qp]*dphi_face[i][qp] ) *
											( qface_normals[qp]*dphi_neighbor_face[j][qp] );
//	                              Ken(i,j) +=
//	                            		  penalty * JxW_face[qp] *
//	                                (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) -
//	                                 phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

	                              // stability
	                              //Ken(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];
	                            }
	                        }
	                    }

	                  // The element and neighbor boundary matrix are now built
	                  // for this side.  Add them to the global matrix
	                  // The SparseMatrix::add_matrix() members do this for us.
//	                  std::cout << "Assembling jump part: " << std::flush;
//	                  std::cout << "part3,  " << std::flush;
		              monodomain_system.get_matrix("stiffness").add_matrix(Kee, dof_indices);
//	                  std::cout << "part4,  " << std::flush;
		              monodomain_system.get_matrix("stiffness").add_matrix(Knn, neighbor_dof_indices);
//	                  std::cout << "part1,  " << std::flush;
		              monodomain_system.get_matrix("stiffness").add_matrix(Kne, neighbor_dof_indices, dof_indices);
//	                  std::cout << "part2,  " << std::flush;
		              monodomain_system.get_matrix("stiffness").add_matrix(Ken, dof_indices, neighbor_dof_indices);
//	                  std::cout << "Done!\n " << std::flush;

	              }


			}
		}

     }
     // closing matrices and vectors
     monodomain_system.get_matrix("mass").close();
     monodomain_system.get_matrix("lumped_mass").close();
     monodomain_system.get_matrix("stiffness").close();
     monodomain_system.get_vector("lumped_mass_vector").close();
     monodomain_system.get_matrix("high_order_mass").close();
     monodomain_system.get_matrix("high_order_mass").add(0.5, monodomain_system.get_matrix("mass"));
     monodomain_system.get_matrix("high_order_mass").add(0.5, monodomain_system.get_matrix("lumped_mass"));

}

void
Monowave::advance()
{
	MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
	IonicModelSystem& ionic_model_system =  M_equationSystems.get_system<IonicModelSystem>("ionic_model");

    *monodomain_system.old_local_solution    = *monodomain_system.solution;
    *monodomain_system.older_local_solution = *monodomain_system.solution;
    *ionic_model_system.old_local_solution    = *ionic_model_system.solution;
    *ionic_model_system.older_local_solution = *ionic_model_system.solution;

	    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.get_system<IonicModelSystem>("wave");
	*wave_system.old_local_solution    = *wave_system.solution;
    *wave_system.older_local_solution = *wave_system.solution;

    IonicModelSystem& iion_system = M_equationSystems.get_system<IonicModelSystem>("iion");
    *iion_system.older_local_solution = *iion_system.old_local_solution;
    *iion_system.old_local_solution   = *iion_system.solution;


}

void
Monowave::solve_reaction_step( double dt,
                               double time,
                               int    step,
                               bool   useMidpoint,
                               const  std::string& mass,
                               libMesh::NumericVector<libMesh::Number>* I4f_ptr )
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain"); //Q
    IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");
    IonicModelSystem& istim_system = M_equationSystems.get_system<IonicModelSystem>("istim");
	    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");
    IonicModelSystem& iion_system = M_equationSystems.get_system<IonicModelSystem>("iion");
    const libMesh::Real alpha = M_equationSystems.parameters.get<libMesh::Real> ("alpha"); // time constant

    monodomain_system.rhs->zero();
    iion_system.solution->zero();
    //iion_system.old_local_solution->zero();
    iion_system.get_vector("diion").zero();

    // get Dnn from setup method for shell approximation
    double D33 = M_conductivity[2];

    double Cm = M_ionicModelPtr->membraneCapacitance();

    auto first_local_index = wave_system.solution->first_local_index();
    auto last_local_index = wave_system.solution->last_local_index();
    int num_vars = ionic_model_system.n_vars();
    double istim = 0.0;
    double dIion = 0.0;
    std::vector<double> values(num_vars+1, 0.0);
    std::vector<double> old_values(num_vars+1, 0.0);
    auto i = first_local_index;
    int var_index = 0;

    // Solve the first time
    for( ; i < last_local_index; )
    {
        double Iion_old = 0.0;

	   if(M_equationType == EquationType::Wave)
	   {
		   values[0]=(*wave_system.old_local_solution)(i); //V
		   old_values[0]=(*monodomain_system.old_local_solution)(i); //Q
		   Iion_old = (*iion_system.old_local_solution)(i); // gating
	   }
	   else
	   {
		   values[0] = (*monodomain_system.old_local_solution)(i);
	   }
        istim = (*istim_system.solution)(i);

        for(int nv = 0; nv < num_vars; ++nv)
        {
        	var_index =  i * num_vars+ nv;
            values[nv+1] = (*ionic_model_system.old_local_solution)(var_index);
            old_values[nv+1] = values[nv+1];
        }

    	M_ionicModelPtr->updateVariables(values, istim, dt);
        double Iion = M_ionicModelPtr->evaluateIonicCurrent(values, istim, dt);
//        double Iion = M_ionicModelPtr->evaluateIonicCurrentH(values, istim, dt, M_meshSize);

        if(M_equationType == EquationType::Wave)
        {
//            dIion =  (Iion - Iion_old) / dt;
            //std::cout << "Iion: " << Iion << ", Iion_old: " << Iion_old << ", dIion1: " << dIion << std::flush;
//            dIion = M_ionicModelPtr->evaluatedIonicCurrent(values, istim, dt, M_meshSize);
            dIion = M_ionicModelPtr->evaluatedIonicCurrent(values, old_values, dt, M_meshSize);
            //std::cout << ", dIion2: " << dIion << std::endl;
        }
        double Isac = 0.0;
        if(I4f_ptr)
        {
            double I4f = (*I4f_ptr)(i);
            Isac = M_ionicModelPtr->evaluateSAC(values[0], I4f);
        }
        Iion += Isac;
//        monodomain_system.get_vector("ionic_currents").set(i,Iion);
        if(M_shellCorrection)
        {
        	double Lambda = 1.0;
        	Iion += ( old_values[0] > 0 ) ? D33*Lambda*old_values[0] : 0.0;
        }
        iion_system.solution->set(i,Iion);
        //iion_system.old_local_solution->set(i,Iion);
		if(M_equationType == EquationType::Wave)
        {
                	iion_system.get_vector("diion").set(i,dIion);
        }
        for(int nv = 0; nv < num_vars; ++nv)
        {
        	var_index =  i  * num_vars + nv;
            ionic_model_system.solution->set(var_index, values[nv+1]);
        }
        ++i;
    }

//   monodomain_system.get_vector("ionic_currents").close();
	iion_system.solution->close();
    ionic_model_system.solution->close();
    iion_system.get_vector("diion").close();

}

void
Monowave::form_system_matrix(double dt, bool /*useMidpoint */, const std::string& mass)
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
	    // WAVE
    M_systemMass = mass;
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");
    IonicModelSystem& iion_system = M_equationSystems.get_system<IonicModelSystem>("iion");
    double Cm = M_ionicModelPtr->membraneCapacitance();

    const libMesh::Real alpha = M_equationSystems.parameters.get<libMesh::Real> ("alpha");

    monodomain_system.matrix->zero();
    monodomain_system.matrix->close();

    // Matrix part
    if(M_equationType == EquationType::Wave)
    {
    	// S = alpha * Chi * M + dt * Chi * M + dt^2 * K
        // S = alpha * M + dt * M + dt^2/Chi * K
    	//  alpha * Chi * M

        // dt * Chi * M
        //monodomain_system.matrix->add(dt+alpha, monodomain_system.get_matrix(mass) );
        //monodomain_system.matrix->add(dt*dt/Chi, monodomain_system.get_matrix("stiffness" ) );

        monodomain_system.matrix->add(Cm*(dt+alpha), monodomain_system.get_matrix(mass) );
        monodomain_system.matrix->add(dt*dt, monodomain_system.get_matrix("stiffness" ) );
    	//  alpha * M
        // dt^2 * K
    }
    else
    {
    	// S = Chi * M + dt * K
        // S = M + dt K
        monodomain_system.matrix->add(1.0, monodomain_system.get_matrix(mass) );
        monodomain_system.matrix->add(dt, monodomain_system.get_matrix("stiffness" ) );
    }
}





void
Monowave::form_system_rhs(double dt, bool useMidpoint , const std::string& mass)
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
	    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");
    IonicModelSystem& iion_system = M_equationSystems.get_system<IonicModelSystem>("iion");
    double Cm = M_ionicModelPtr->membraneCapacitance();

    const libMesh::Real alpha = M_equationSystems.parameters.get<libMesh::Real> ("alpha"); // time constant

    monodomain_system.rhs->zero();
    monodomain_system.rhs->close();
    // RHS part:
    // RHS WAVE = alpha * M * Q^n - dt / Chi * K * V^n + dt * M * I^n + dt * M * dI^n
    // RHS ADR  =   M * V^n + dt * M * I^n
    // commmon part is: dt * Chi * M * I^n
    //iion_system.solution->scale(dt*Chi);
   monodomain_system.get_matrix(mass).vector_mult_add( *monodomain_system.rhs,
                                                       *iion_system.solution );
   if(M_equationType == EquationType::Wave)
   {
       //iion_system.get_vector("diion").pointwise_mult(*monodomain_system.older_local_solution, iion_system.get_vector("diion"));
       iion_system.get_vector("diion").scale(alpha);
       monodomain_system.get_matrix(mass).vector_mult_add( *monodomain_system.rhs,
                                                           iion_system.get_vector("diion") );
       monodomain_system.rhs->scale(dt);



       *monodomain_system.old_local_solution = *monodomain_system.older_local_solution;
       // + alpha  * M * Q^n
       monodomain_system.older_local_solution->scale(alpha*Cm);
	   monodomain_system.get_matrix(M_systemMass).vector_mult_add( *monodomain_system.rhs,
                                                           *monodomain_system.older_local_solution );
       *monodomain_system.older_local_solution = *monodomain_system.old_local_solution;

	   // - dt  * K * V^n
       *wave_system.older_local_solution = *wave_system.old_local_solution;
	   wave_system.older_local_solution->scale(-dt);
	   monodomain_system.get_matrix("stiffness").vector_mult_add( *monodomain_system.rhs,
                                                                  *wave_system.older_local_solution );
       *wave_system.older_local_solution = *wave_system.old_local_solution;

   }
   else
   {
       monodomain_system.rhs->scale(dt);
	   //    M * V^n
	   monodomain_system.get_matrix(M_systemMass).vector_mult_add( *monodomain_system.rhs,
														   *monodomain_system.older_local_solution );
   }
}





void
Monowave::solve_diffusion_step(double dt, double time,  bool useMidpoint , const std::string& mass, bool reassemble)
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    form_system_rhs(dt, useMidpoint, mass);

    if(M_timeIntegratorType == DynamicTimeIntegratorType::Implicit)
    {
        double tol = 1e-12;
        double max_iter = 2000;

        std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

        rval = M_linearSolver->solve (*monodomain_system.matrix, nullptr,
                                                            *monodomain_system.solution,
                                                            *monodomain_system.rhs, tol, max_iter);
    }
    else
    {
        double Cm = M_ionicModelPtr->membraneCapacitance();
        const libMesh::Real alpha = M_equationSystems.parameters.get<libMesh::Real> ("alpha"); // time constant

        auto first_local_index = monodomain_system.solution->first_local_index();
        auto last_local_index = monodomain_system.solution->last_local_index();
        auto i = first_local_index;
        double rhs = 0.0;
        double ml = 0.0;
        double q_np1 = 0.0;
        for(    i = first_local_index;  i < last_local_index; ++i)
        {
            rhs  =  (*monodomain_system.rhs)(i);
            ml = monodomain_system.get_vector("lumped_mass_vector")(i);
            //v_new = v_old + dt / Cm * Iion / ml;
            q_np1 = rhs / (dt+alpha) / ml / Cm;
            monodomain_system.solution->set(i, q_np1);
        }
        monodomain_system.solution->close();

    }
    // WAVE
	IonicModelSystem& wave_system =  M_equationSystems.add_system<IonicModelSystem>("wave");
    if(M_equationType == EquationType::Wave)
    {
	   *wave_system.solution = *monodomain_system.solution;
	   wave_system.solution->scale(dt);
	   *wave_system.solution += *wave_system.old_local_solution;
    }
    else
    {
 	   *wave_system.solution = *monodomain_system.solution;
    }
}


void
Monowave::reinit_linear_solver()
{
	M_linearSolver->clear();
    M_linearSolver->set_solver_type(libMesh::CG);
    M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
    M_linearSolver->init();
}

std::string
Monowave::get_ionic_model_name() const
{
	return M_ionicModelPtr->ionicModelName();
}

double
Monowave::potential_norm()
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    return monodomain_system.solution->l1_norm();
}


const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Monowave::get_fibers()
{
    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
    return fiber_system.solution;
}

const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Monowave::get_sheets()
{
    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
    return sheets_system.solution;
}

const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Monowave::get_xfibers()
{
    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");
    return xfiber_system.solution;
}


void
Monowave::set_potential_on_boundary(unsigned int boundID, double value )
{
std::cout << "* MONODOMAIN: WARNING:  set_potential_on_boundary works only for TET4" << std::endl;
     libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
	  const unsigned int dim = mesh.mesh_dimension();

    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    const libMesh::DofMap & dof_map = monodomain_system.get_dof_map();
    std::vector<libMesh::dof_id_type> dof_indices;


         libMesh::FEType fe_type = dof_map.variable_type(0);
           // Declare a special finite element object for
  // boundary integration.
  libMesh::UniquePtr<libMesh::FEBase> fe_face (libMesh::FEBase::build(dim, fe_type));

  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  libMesh::QGauss qface(dim-1, libMesh::FIRST);

  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);
 libMesh::DenseVector<libMesh::Number> Fe;


    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();

    for( ; el != end_el; ++el)
    {
    	const libMesh::Elem * elem = *el;
    	dof_map.dof_indices(elem, dof_indices);
    	Fe.resize(dof_indices.size());
        for (unsigned int side=0; side<elem->n_sides(); side++)
        {
			if (elem->neighbor(side) == libmesh_nullptr)
            {
                const unsigned int boundary_id =
                mesh.boundary_info->boundary_id (elem, side);
                if(boundary_id == boundID)
                {
                	monodomain_system.solution->set(dof_indices[0],value);
                	monodomain_system.solution->set(dof_indices[1],value);
                	monodomain_system.solution->set(dof_indices[2],value);
                }
            }

        }
    }
    monodomain_system.solution->close();

}


} /* namespace BeatIt */
