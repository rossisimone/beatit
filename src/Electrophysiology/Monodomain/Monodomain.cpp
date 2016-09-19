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
#include "Electrophysiology/Monodomain/Monodomain.hpp"
#include "Util/SpiritFunction.hpp"


#include "Electrophysiology/Pacing/PacingProtocolSpirit.hpp"

#include "Util/Timer.hpp"

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
    // Add lumped mass matrix
    monodomain_system.add_vector("lumped_mass_vector");
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
    ParameterSystem& procID_system = M_equationSystems.add_system<ParameterSystem>("ProcID");
    procID_system.add_variable( "ProcID", libMesh::CONSTANT, libMesh::MONOMIAL);
     M_parametersExporterNames.insert("activation_times");
    M_parametersExporterNames.insert("fibers");
    M_parametersExporterNames.insert("sheets");
    M_parametersExporterNames.insert("xfibers");
    M_parametersExporterNames.insert("conductivity");
    M_parametersExporterNames.insert("ProcID");

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
    if(v_ic != "")
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
    std::cout << "* MONODOMAIN: Copying initial conditions to vectors at n nd at n-1... " << std::flush;

    // Close vectors and update the values in old_local_solution and older_local_solution
    monodomain_system.solution->close();
	monodomain_system.old_local_solution->close();
	monodomain_system.older_local_solution->close();
	ionic_model_system.solution->close();
	ionic_model_system.old_local_solution->close();
	ionic_model_system.older_local_solution->close();
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
Monodomain::generate_fibers(   const GetPot& data,
                               const std::string& section)
{
    std::cout << "* MONODOMAIN: Creating fiber fields" << std::endl;

    std::cout << "* MONODOMAIN: Solving Poisson porblem" << std::endl;

    libMesh::Mesh new_mesh( dynamic_cast< libMesh::Mesh&>(M_equationSystems.get_mesh() ) );
    Poisson p(new_mesh);
    p.setup(data, section);
    p.assemble_system();
    p.solve_system();
    p.compute_elemental_solution_gradient();
    const auto& sol_ptr = p.get_P0_solution();


    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");

    *sheets_system.solution = *p.get_gradient();

    auto first_local_index = fiber_system.solution->first_local_index();
    auto last_local_index = fiber_system.solution->last_local_index();
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
    auto j = first_local_index;
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
        fiber_system. solution->set(i,   f0x);
        fiber_system. solution->set(i+1, f0y);
        fiber_system. solution->set(i+2, f0z);

        i += 3;
    }

}

double
Monodomain::last_activation_time()
{
    ParameterSystem& activation_times_system  = M_equationSystems.add_system<ParameterSystem>("activation_times");
    return static_cast<double>(activation_times_system.solution->max());
}


void
Monodomain::save(int step)
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
//	std::cout << "Estimating Monodomain Error  " << std::endl;
//	timer.restart();
	p_error_estimator->estimate_error(monodomain_system, error);
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
     monodomain_system.get_vector("lumped_mass_vector").zero();


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
         monodomain_system.get_vector("lumped_mass_vector").add_vector(Fe, dof_indices);
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
     }
     // closing matrices and vectors
     monodomain_system.get_matrix("mass").close();
     monodomain_system.get_matrix("lumped_mass").close();
     monodomain_system.get_matrix("stiffness").close();
     monodomain_system.get_vector("lumped_mass_vector").close();
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
    double Cm = M_ionicModelPtr->membraneCapacitance();

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

    	M_ionicModelPtr->updateVariables(values, istim, dt);
        double Iion = M_ionicModelPtr->evaluateIonicCurrent(values, istim, dt);
        monodomain_system.get_vector("ionic_currents").set(i,Iion);
        for(int nv = 0; nv < num_vars; ++nv)
        {
        	var_index =  i  * num_vars + nv;
            ionic_model_system.solution->set(var_index, values[nv+1]);
        }
        ++i;
    }

    monodomain_system.get_vector("ionic_currents").close();
   monodomain_system.get_matrix(mass).vector_mult_add( *monodomain_system.rhs,
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 monodomain_system.get_vector("ionic_currents") );


    double v_old, Iion, ml, v_new;
    for(    i = first_local_index;  i < last_local_index; ++i)
    {
        v_old = (*monodomain_system.old_local_solution)(i);
        Iion  =  (*monodomain_system.rhs)(i);
        ml = monodomain_system.get_vector("lumped_mass_vector")(i);
        //v_new = v_old + dt / Cm * Iion / ml;
        v_new = v_old + dt * Iion / ml;
        monodomain_system.solution->set(i, v_new);
    }
    monodomain_system.solution->close();
    ionic_model_system.solution->close();

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

			M_ionicModelPtr->updateVariables(val_n, val_np1, istim, dt);
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
			ml = monodomain_system.get_vector("lumped_mass_vector")(i);
			v_new = v_old + dt / Cm * Iion / ml;
			monodomain_system.solution->set(i, v_new);
		}


		monodomain_system.solution->close();
		ionic_model_system.solution->close();
    }

}

void
Monodomain::form_system_matrix(double dt, bool useMidpoint , const std::string& mass)
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    const libMesh::Real Chi = M_equationSystems.parameters.get<libMesh::Real> ("Chi");
    double Cm = M_ionicModelPtr->membraneCapacitance();

    monodomain_system.matrix->zero();
    monodomain_system.matrix->close();
    monodomain_system.matrix->add(1.0, monodomain_system.get_matrix(mass) );
    if(useMidpoint)
    {
        //      monodomain_system.matrix->add(dt/(2.0*Chi*Cm), monodomain_system.get_matrix("stiffness" ) );
        //      monodomain_system.get_matrix("stiffness").vector_mult_add( *monodomain_system.rhs,
        monodomain_system.matrix->add(dt/(2.0*Chi), monodomain_system.get_matrix("stiffness" ) );
    }
    else
    {
//      monodomain_system.matrix->add(dt/Chi/Cm, monodomain_system.get_matrix("stiffness" ) );
        monodomain_system.matrix->add(dt/Chi, monodomain_system.get_matrix("stiffness" ) );
    }

}

void
Monodomain::solve_diffusion_step(double dt, double time,  bool useMidpoint , const std::string& mass, bool reassemble)
{
    const libMesh::Real Chi = M_equationSystems.parameters.get<libMesh::Real> ("Chi");
    double Cm = M_ionicModelPtr->membraneCapacitance();
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    if(reassemble) form_system_matrix(dt, useMidpoint, mass);
    monodomain_system.rhs->zero();

    if(useMidpoint)
    {
        monodomain_system.get_matrix("stiffness").vector_mult_add( *monodomain_system.rhs, *monodomain_system.solution);
        //      monodomain_system.rhs->scale( -dt/(2.0*Chi*Cm) );
    	monodomain_system.rhs->scale( -dt/(2.0*Chi) );
    }


	 monodomain_system.get_matrix(mass).vector_mult_add( *monodomain_system.rhs, *monodomain_system.solution);

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0,0.0);

    rval = M_linearSolver->solve (*monodomain_system.matrix, nullptr,
														*monodomain_system.solution,
														*monodomain_system.rhs, tol, max_iter);

}


void
Monodomain::reinit_linear_solver()
{
	M_linearSolver->clear();
    M_linearSolver->set_solver_type(libMesh::CG);
    M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
    M_linearSolver->init();
}

std::string
Monodomain::get_ionic_model_name() const
{
	return M_ionicModelPtr->ionicModelName();
}

double
Monodomain::potential_norm()
{
    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    return monodomain_system.solution->l1_norm();
}


const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Monodomain::get_fibers()
{
    ParameterSystem& fiber_system        = M_equationSystems.get_system<ParameterSystem>("fibers");
    return fiber_system.solution;
}

const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Monodomain::get_sheets()
{
    ParameterSystem& sheets_system       = M_equationSystems.get_system<ParameterSystem>("sheets");
    return sheets_system.solution;
}

const
libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
Monodomain::get_xfibers()
{
    ParameterSystem& xfiber_system       = M_equationSystems.get_system<ParameterSystem>("xfibers");
    return xfiber_system.solution;
}


void
Monodomain::set_potential_on_boundary(unsigned int boundID, double value )
{
    const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
    libMesh::MeshBase::const_element_iterator el =
            mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el =
            mesh.active_local_elements_end();

    MonodomainSystem& monodomain_system  =  M_equationSystems.get_system<MonodomainSystem>("monodomain");
    std::vector<libMesh::dof_id_type > node_id_list;
    std::vector<libMesh::boundary_id_type > bc_id_list;
    mesh.boundary_info->build_node_list_from_side_list();
    mesh.boundary_info->build_node_list (node_id_list, bc_id_list);
    auto first_local_index = monodomain_system.solution->first_local_index();
    auto last_local_index = monodomain_system.solution->last_local_index();

    for ( int i = 0; i < bc_id_list.size(); ++i )
    {
    	if( bc_id_list[i] == boundID )
    	{
            monodomain_system.solution->set(node_id_list[i], value);
        }
    }

}



} /* namespace BeatIt */
