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
#include "libmesh/gmv_io.h"

#include "libmesh/perf_log.h"

#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"

#include <sys/stat.h>

#include "Electrophysiology/IonicModels/NashPanfilov.hpp"
#include "Electrophysiology/IonicModels/ORd.hpp"
#include "Electrophysiology/Monodomain/Monodomain.hpp"
#include "Util/SpiritFunction.hpp"


#include "Electrophysiology/Pacing/PacingProtocolSpirit.hpp"


namespace BeatIt
{

// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
typedef libMesh::TransientLinearImplicitSystem     MonodomainSystem;
typedef libMesh::TransientExplicitSystem           IonicModelSystem;
typedef libMesh::ExplicitSystem                     ParameterSystem;


Monodomain::Monodomain( libMesh::MeshBase & mesh )
    : M_equationSystems(mesh)
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
{

}

Monodomain::~Monodomain(){}



void Monodomain::setup(GetPot& data, std::string section)
{
    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Read Input File
    M_datafile = data;
    //Read output folder from datafile
    std::string output_folder = M_datafile(section+"/output_folder",  "Output");
    M_outputFolder = "./" + output_folder + "/";

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // 1) ADR
    std::cout << "* MONODOMAIN: Creating new System for the monodomain diffusion reaction equation" << std::endl;
    MonodomainSystem& monodomain_system  =  M_equationSystems.add_system<MonodomainSystem>("monodomain");
    // TO DO: Generalize to higher order
    monodomain_system.add_variable( "V", libMesh::FIRST);
    // Add 3 matrices
    monodomain_system.add_matrix("lumped_mass");
    monodomain_system.add_matrix("mass");
    monodomain_system.add_matrix("stiffness");
    monodomain_system.add_vector("ionic_currents");
    M_monodomainExporterNames.insert("monodomain");

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // 2) ODEs
    std::cout << "* MONODOMAIN: Creating new System for the ionic model " << std::endl;
    IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");
    // Create Ionic Model
    std::string ionic_model = M_datafile(section+"/ionic_model", "NashPanfilov");
    M_ionicModelPtr.reset(  BeatIt::IonicModel::IonicModelFactory::Create(ionic_model) );
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
    // Add lumped mass matrix
    ionic_model_system.add_vector("lumped_mass");

    // Add the applied current to this system
    IonicModelSystem& istim_system = M_equationSystems.add_system<IonicModelSystem>("istim");
    istim_system.add_variable( "istim", libMesh::FIRST);

    M_ionicModelExporterNames.insert("ionic_model");
    M_ionicModelExporterNames.insert("istim");

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
    M_parametersExporterNames.insert("activation_times");
    M_parametersExporterNames.insert("fibers");
    M_parametersExporterNames.insert("sheets");
    M_parametersExporterNames.insert("xfibers");
    M_parametersExporterNames.insert("conductivity");

    M_equationSystems.init();
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


    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Setup Exporters
    M_monodomainExporter.reset( new Exporter(M_equationSystems.get_mesh() ) );
    M_ionicModelExporter.reset( new Exporter(M_equationSystems.get_mesh() ) );
    M_parametersExporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );
    M_monodomainEXOExporter.reset( new EXOExporter(M_equationSystems.get_mesh() ) );

    struct stat out_dir;
    if( stat(&M_outputFolder[0],&out_dir) != 0  )
    {
        if ( monodomain_system.get_mesh().comm().rank() == 0 )
        {
            mkdir ( M_outputFolder.c_str(), 0777 );
        }
    }

    std::cout << "* MONODOMAIN: creating pacing protocol" << std::endl;
    M_pacing.reset( new PacingProtocolSpirit() );
    M_pacing->setup(M_datafile, "monodomain/pacing");

}


void
Monodomain::init(double time)
{
    // Setting initial conditions
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    IonicModelSystem& ionic_model_system =  M_equationSystems.get_system<IonicModelSystem>("ionic_model");
    int num_vars = ionic_model_system.n_vars();
    std::vector<double> init_values(num_vars+1, 0.0);
    M_ionicModelPtr->initialize(init_values);
    auto first_local_index = monodomain_system.solution->first_local_index();
    auto last_local_index = monodomain_system.solution->last_local_index();
    for(auto i = first_local_index; i < last_local_index; ++i)
    {
        monodomain_system.solution->set(i,init_values[0]);
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
    if("v_ic" != "")
    {
        std::cout << "* MONODOMAIN: Found monodomain initial condition: " << v_ic << std::endl;
        SpiritFunction monodomain_ic;
        monodomain_ic.read(v_ic);
        M_equationSystems.parameters.set<libMesh::Real>("time") =  time;
        monodomain_system.time = time;
        std::cout << "* MONODOMAIN: Projecting initial condition to monodomain system ... " << std::flush;
        monodomain_system.project_solution(&monodomain_ic);
        std::cout << " done" << std::endl;
    }


    IonicModelSystem& istim_system = M_equationSystems.get_system<IonicModelSystem>("istim");
    ParameterSystem& activation_times_system  = M_equationSystems.get_system<ParameterSystem>("activation_times");
    for(auto i = first_local_index; i < last_local_index; ++i)
    {
    	activation_times_system.solution->set(i,-1.0);
    }

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

    fiber_system .project_solution(&fibers_func );
    sheets_system.project_solution(&sheets_func );
    xfiber_system.project_solution(&xfibers_func);


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

}

void
Monodomain::update_pacing(double time)
{
    // Add the applied current to this system
    IonicModelSystem& istim_system = M_equationSystems.get_system<IonicModelSystem>("istim");
    istim_system.time = time;
    istim_system.project_solution( &M_pacing->pacing() );
}

void
Monodomain::update_activation_time( double  time, double threshold )
{
    ParameterSystem& activation_times_system  = M_equationSystems.add_system<ParameterSystem>("activation_times");
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    auto first = monodomain_system.solution->first_local_index();
    auto last  = monodomain_system.solution->last_local_index();
    double v = 0.0;
    double at = 0.0;
    for(int index = first; index < last; index++ )
    {
    	v = (*monodomain_system.solution)(index);
    	at = (*activation_times_system.solution)(index);
    	if( v > threshold && at < 0.0)
    	{
    		activation_times_system.solution->set(index, time);
    	}
    }
	activation_times_system.solution->close();

}



void
Monodomain::save(int step)
{
    std::cout << "* MONODOMAIN: GMVIO::Exporting " << step << " in: "  << M_outputFolder << " ... " << std::flush;
    M_monodomainExporter->write_equation_systems ( M_outputFolder+"monodomain.gmv."+std::to_string(step),
                                                   M_equationSystems,
                                                   &M_monodomainExporterNames);
    M_ionicModelExporter->write_equation_systems ( M_outputFolder+"ionic_model.gmv."+std::to_string(step),
                                                   M_equationSystems,
                                                   &M_ionicModelExporterNames);
    std::cout << "done " << std::endl;
}
void
Monodomain::init_exo_output()
{
    M_monodomainEXOExporter->write_equation_systems (M_outputFolder+"monodomain.exo", M_equationSystems);
    M_monodomainEXOExporter->append(true);

}

void
Monodomain::save_exo(int step, double time)
{
    std::cout << "* MONODOMAIN: EXODUSII::Exporting time "   << time << " in: "  << M_outputFolder << " ... " << std::flush;
    M_monodomainEXOExporter->write_timestep(  M_outputFolder+"monodomain.exo"
                                         , M_equationSystems
                                         , step, time );
    std::cout << "done " << std::endl;
}

void
Monodomain::save_potential(int step)
{
    std::cout << "* MONODOMAIN: GMVIO::Exporting " << step << " in: "  << M_outputFolder << " ... " << std::flush;
    M_monodomainExporter->write_equation_systems ( M_outputFolder+"monodomain.gmv."+std::to_string(step),
                                                   M_equationSystems,
                                                   &M_monodomainExporterNames);
    std::cout << "done " << std::endl;
}

void
Monodomain::amr( libMesh:: MeshRefinement& mesh_refinement, const std::string& type)
{
//	std::cout << "* MONODOMAIN: starting AMR ... " << std::flush;
	libMesh::ErrorVector error;
	libMesh::ErrorEstimator * p_error_estimator;
	if("kelly" == type) p_error_estimator = new  libMesh::KellyErrorEstimator ;
	else p_error_estimator = new  libMesh::LaplacianErrorEstimator ;
//	 libMesh::KellyErrorEstimator error_estimator;
//	 libMesh::LaplacianErrorEstimator error_estimator;
	MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
	p_error_estimator->estimate_error(monodomain_system, error);
    // Flag elements to be refined and coarsened
    mesh_refinement.flag_elements_by_error_fraction (error);

    // Perform refinement and coarsening
    mesh_refinement.flag_elements_by_error_fraction (error);
//	std::cout << " coarsen and refine ...  " << std::flush;
    mesh_refinement.refine_and_coarsen_elements();
//	std::cout << " reinit ...  " << std::flush;
    M_equationSystems.reinit();
//	std::cout << " done  " << std::endl;

}


void
Monodomain::save_parameters()
{
    std::cout << "* MONODOMAIN: EXODUSII::Exporting parameters in: "  << M_outputFolder << " ... " << std::flush;
    M_parametersExporter->write_equation_systems (M_outputFolder+"parameters.exo", M_equationSystems, &M_parametersExporterNames);
    M_parametersExporter->write_element_data(M_equationSystems);
    std::cout << "done " << std::endl;
}


void
Monodomain::assemble_matrices()
{

//    std::cout << "* MONODOMAIN: Assembling matrices ... " << std::flush;
     using libMesh::UniquePtr;

     const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
     const unsigned int dim = mesh.mesh_dimension();
     // Get a reference to the LinearImplicitSystem we are solving
     MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
     IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");

     monodomain_system.get_matrix("mass").zero();
     monodomain_system.get_matrix("lumped_mass").zero();
     monodomain_system.get_matrix("stiffness").zero();
     ionic_model_system.get_vector("lumped_mass").zero();


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
     libMesh::QGauss qrule_stiffness(dim, libMesh::FIRST);
     // A 5th order Gauss quadrature rule for numerical integration.
     libMesh::QGauss qrule_mass(dim, libMesh::SECOND);

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

     const libMesh::Real Chi = M_equationSystems.parameters.get<libMesh::Real> ("Chi");

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
         Ke.resize(dof_indices.size(), dof_indices.size());
         Me.resize(dof_indices.size(), dof_indices.size());
         Mel.resize(dof_indices.size(), dof_indices.size());
         Fe.resize(dof_indices.size());
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

         for(int idim = 0; idim < dim; ++idim )
         {
             for(int jdim = 0; jdim < dim; ++jdim )
             {
                 D0(idim,jdim) = Dff * f0[idim] * f0[jdim]
                               + Dss * s0[idim] * s0[jdim]
                               + Dnn * n0[idim] * n0[jdim];
             }
         }

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
         for (unsigned int qp = 0; qp < qrule_stiffness.n_points(); qp++)
         {
             for (unsigned int i = 0; i < phi_qp1.size(); i++)
             {
                 DgradV = D0  * dphi_qp1[i][qp];

//                 std::cout << "DgradV = " << std::endl;
//        		 for(int idim = 0; idim < dim; ++idim )
//                 {
//                	 std::cout << DgradV(idim, 0) << "," << DgradV(idim, 1) << ", " << DgradV(idim, 2) << std::endl;
//                 }

                 for (unsigned int j = 0; j < phi_qp1.size(); j++)
                 {
                     // stiffness term
                     Ke(i, j) += JxW_qp1[qp] * DgradV * dphi_qp1[j][qp];
                 }
             }
         }
         monodomain_system.get_matrix("stiffness").add_matrix(Ke, dof_indices);
         ionic_model_system.get_vector("lumped_mass").add_vector(Fe, dof_indices);
     }
//     std::cout << " done" << std::endl;
     monodomain_system.get_matrix("mass").close();
     monodomain_system.get_matrix("lumped_mass").close();
     monodomain_system.get_matrix("stiffness").close();
     ionic_model_system.get_vector("lumped_mass").close();
}

void
Monodomain::advance()
{
	MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
	IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");

    *monodomain_system.old_local_solution    = *monodomain_system.solution;
    *monodomain_system.older_local_solution = *monodomain_system.solution;
    *ionic_model_system.old_local_solution    = *ionic_model_system.solution;
    *ionic_model_system.older_local_solution = *ionic_model_system.solution;
}

void
Monodomain::solve_reaction_step(double dt, double time, int step, bool useMidpoint, const std::string& mass)
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    IonicModelSystem& ionic_model_system =  M_equationSystems.add_system<IonicModelSystem>("ionic_model");
    IonicModelSystem& istim_system = M_equationSystems.get_system<IonicModelSystem>("istim");
    monodomain_system.rhs->zero();

    if(step > 0)
    {
        *monodomain_system.old_local_solution    = *monodomain_system.solution;
        *ionic_model_system.old_local_solution    = *ionic_model_system.solution;
    }
    else
    {
        *monodomain_system.old_local_solution    = *monodomain_system.older_local_solution;
        *ionic_model_system.old_local_solution    = *ionic_model_system.older_local_solution;

    }

    auto first_local_index = monodomain_system.solution->first_local_index();
    auto last_local_index = monodomain_system.solution->last_local_index();
    int num_vars = ionic_model_system.n_vars();
    double istim = 0.0;
    std::vector<double> values(num_vars+1, 0.0);
    auto i = first_local_index;
    int var_index = 0;
    auto& send_list = monodomain_system.get_dof_map().get_send_list();

    // Solve the first time
    for( ; i < last_local_index; )
    {
        values[0] = (*monodomain_system.old_local_solution)(i);
        istim = (*istim_system.solution)(i);

        for(int nv = 0; nv < num_vars; ++nv)
        {
        	var_index =  i * num_vars+ nv;
            values[nv+1] = (*ionic_model_system.old_local_solution)(var_index);
        }

    	M_ionicModelPtr->updateVariables(values, dt);
        double Iion = M_ionicModelPtr->evaluateIonicCurrent(values, istim, dt);
        monodomain_system.get_vector("ionic_currents").set(i,Iion);
//        double new_v = values[0]+dt*Iion;
//        monodomain_system.solution->set(i, new_v);
        for(int nv = 0; nv < num_vars; ++nv)
        {
        	var_index =  i  * num_vars + nv;
            ionic_model_system.solution->set(var_index, values[nv+1]);
        }
        ++i;
    }
//
    monodomain_system.get_vector("ionic_currents").close();
////    double val = 1.0;
	 monodomain_system.get_matrix(mass).vector_mult_add( *monodomain_system.rhs,
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 monodomain_system.get_vector("ionic_currents") );
    i = first_local_index;
    double v_old, Iion, ml, v_new;
//
    for( ; i < last_local_index; ++i)
    {
        v_old = (*monodomain_system.old_local_solution)(i);
        Iion  =  (*monodomain_system.rhs)(i);
//        double Iion2  = monodomain_system. get_vector("ionic_currents")(i);
        ml = ionic_model_system.get_vector("lumped_mass")(i);
        v_new = v_old + dt * Iion / ml;
//        v_new = v_old + dt * Iion2;
        monodomain_system.solution->set(i, v_new);
    }
    monodomain_system.solution->close();
    ionic_model_system.solution->close();
//

    // Solve the second time
    if(useMidpoint)
    {
		std::vector<double> val_n(num_vars+1, 0.0);
		std::vector<double> val_np1(num_vars+1, 0.0);
		i = first_local_index;

		for( ; i < last_local_index; )
		{
			val_np1[0] = (*monodomain_system.solution)(i);
			val_n[0] = (*monodomain_system.old_local_solution)(i);
			istim = (*istim_system.solution)(i);

			for(int nv = 0; nv < num_vars; ++nv)
			{
				var_index =  i * num_vars+ nv;
				val_np1[nv+1] = (*ionic_model_system.solution)(var_index);
				val_n[nv+1] = (*ionic_model_system.old_local_solution)(var_index);
			}

			M_ionicModelPtr->updateVariables(val_n, val_np1, dt);
			double Iion = M_ionicModelPtr->evaluateIonicCurrent(val_n, val_np1, istim, dt);
			monodomain_system.get_vector("ionic_currents").set(i,Iion);
			for(int nv = 0; nv < num_vars; ++nv)
			{
				var_index =  i  * num_vars + nv;
				ionic_model_system.solution->set(var_index, val_np1[nv+1]);
			}
			++i;
		}


		monodomain_system.get_vector("ionic_currents").close();
		monodomain_system.rhs->zero();
		monodomain_system.rhs->close();
		 monodomain_system.get_matrix(mass).vector_mult_add( *monodomain_system.rhs,
				 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 monodomain_system.get_vector("ionic_currents") );
		i = first_local_index;

		for( ; i < last_local_index; ++i)
		{
			v_old = (*monodomain_system.old_local_solution)(i);
			Iion  =  (*monodomain_system.rhs)(i);
			ml = ionic_model_system.get_vector("lumped_mass")(i);
			v_new = v_old + dt * Iion / ml;
			monodomain_system.solution->set(i, v_new);
		}
		monodomain_system.solution->close();
		ionic_model_system.solution->close();
    }

}


void
Monodomain::solve_diffusion_step(double dt, double time,  bool useMidpoint , const std::string& mass)
{
    const libMesh::Real Chi = M_equationSystems.parameters.get<libMesh::Real> ("Chi");
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    monodomain_system.matrix->zero();
    monodomain_system.rhs->zero();
    monodomain_system.matrix->close();
    monodomain_system.matrix->add(1.0, monodomain_system.get_matrix(mass) );

    if(useMidpoint)
    {
    	monodomain_system.matrix->add(dt/(2.0*Chi), monodomain_system.get_matrix("stiffness" ) );
    	monodomain_system.get_matrix("stiffness").vector_mult_add( *monodomain_system.rhs,
      																													   *monodomain_system.solution);
    	monodomain_system.rhs->scale( -dt/(2.0*Chi) );
    }
    else
    {
    	monodomain_system.matrix->add(dt/Chi, monodomain_system.get_matrix("stiffness" ) );
    }

	 monodomain_system.get_matrix(mass).vector_mult_add( *monodomain_system.rhs,
   																				             									   *monodomain_system.solution);

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

    rval = M_linearSolver->solve (*monodomain_system.matrix, nullptr,
														*monodomain_system.solution,
														*monodomain_system.rhs, tol, max_iter);
}


void Monodomain::reinit_linear_solver()
{
	M_linearSolver->clear();
    M_linearSolver->set_solver_type(libMesh::CG);
    M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
    M_linearSolver->init();
}

} /* namespace BeatIt */
