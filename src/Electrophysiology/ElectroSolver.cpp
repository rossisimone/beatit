/*
 * ElectroSolver.cpp
 *
 *  Created on: Oct 19, 2017
 *      Author: srossi
 */

// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/type_tensor.h"
//#include "PoissonSolver/Poisson.hpp"

// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
//
#include "libmesh/vector_value.h"
#include "libmesh/linear_solver.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/elem.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h"

#include "libmesh/exodusII_io.h"
//#include "libmesh/exodusII_io_helper.h"
#include "libmesh/gmv_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/nemesis_io.h"

#include "libmesh/perf_log.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_preconditioner_type.h"
#include "libmesh/enum_solver_type.h"

#include <sys/stat.h>

#include "Electrophysiology/IonicModels/NashPanfilov.hpp"
#include "Electrophysiology/IonicModels/Grandi11.hpp"
#include "Electrophysiology/IonicModels/ORd.hpp"
#include "Electrophysiology/IonicModels/TP06.hpp"
#include "Electrophysiology/ElectroSolver.hpp"
#include "Util/SpiritFunction.hpp"

#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "Electrophysiology/Pacing/PacingProtocolSpirit.hpp"
#include "Util/IO/io.hpp"

namespace BeatIt
{

// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
    typedef libMesh::TransientLinearImplicitSystem ElectroSystem;
    typedef libMesh::TransientExplicitSystem IonicModelSystem;
    typedef libMesh::ExplicitSystem ParameterSystem;

    ElectroSolver::ElectroSolver(libMesh::EquationSystems& es, std::string model)
            : M_equationSystems(es), M_exporter(), M_exporterNames(), M_ionicModelExporter(), M_ionicModelExporterNames(), M_parametersExporter(), M_parametersExporterNames(), M_outputFolder(), M_datafile(), M_pacing_i(), M_pacing_e(), M_linearSolver(), M_anisotropy(
                    Anisotropy::Orthotropic), M_equationType(EquationType::ParabolicEllipticBidomain), M_timeIntegratorType(DynamicTimeIntegratorType::Implicit), M_useAMR(false), M_assembleMatrix(
                    true), M_systemMass("lumped"), M_intraConductivity(), M_extraConductivity(), M_conductivity(), M_meshSize(1.0), M_model(model), M_ground_ve(Ground::Nullspace), M_timeIntegrator(
                    TimeIntegrator::FirstOrderIMEX), M_timestep_counter(0), M_symmetricOperator(false), M_elapsed_time(), M_num_linear_iters(0), M_order(libMesh::FIRST), M_FEFamily(libMesh::LAGRANGE)
    {
        // TODO Auto-generated constructor stub

    }

    ElectroSolver::~ElectroSolver()
    {
    }

    void ElectroSolver::setup(GetPot& data, std::string section)
    {
        // ///////////////////////////////////////////////////////////////////////
        // ///////////////////////////////////////////////////////////////////////
        // Read Input File
        M_datafile = data;
        M_section = section;

        std::string tissueBlockID = M_datafile(section + "/tissue_blockIDs", "");
        BeatIt::readList(tissueBlockID, M_tissueBlockIDs);
        secret_blockID_key = 0;
        //Read output folder from datafile
        std::string output_folder = M_datafile(M_section + "/output_folder", "Output");
        M_outputFolder = "./" + output_folder + "/";

        std::cout << "* " << M_model << ": creating pacing protocol" << std::endl;

        std::string pacing_type = M_datafile(M_section + "/pacing/type", "NONE");
        std::cout << "* ElectroSolver: pacing protocol type : " << pacing_type << " read from: " << M_section << "/pacing/type" << std::endl;
        if ("NONE" != pacing_type)
        {
            M_pacing.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_type));
            M_pacing->setup(M_datafile, M_section + "/pacing");
            M_pacing->showMe();
        }

        std::string pacing_i_type = M_datafile(M_section + "/pacing_i/type", "NONE");
        if ("NONE" != pacing_i_type)
        {
            M_pacing_i.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_i_type));
            M_pacing_i->setup(M_datafile, M_section + "/pacing_i");
            M_pacing_i->showMe();
        }

        std::string pacing_e_type = M_datafile(M_section + "/pacing_e/type", "NONE");
        if ("NONE" != pacing_e_type)
        {
            M_pacing_e.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_e_type));
            M_pacing_e->setup(M_datafile, M_section + "/pacing_e");
            M_pacing_e->showMe();
        }

        std::string surf_pacing_i_type = M_datafile(M_section + "/surf_pacing_i/type", "NONE");
        if ("NONE" != surf_pacing_i_type)
        {
            M_surf_pacing_i.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(surf_pacing_i_type));
            M_surf_pacing_i->setup(M_datafile, M_section + "/surf_pacing_i");
        }

        std::string surf_pacing_e_type = M_datafile(M_section + "/surf_pacing_e/type", "NONE");
        if ("NONE" != surf_pacing_e_type)
        {
            M_surf_pacing_e.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(surf_pacing_e_type));
            M_surf_pacing_e->setup(M_datafile, M_section + "/surf_pacing_e");
        }

        struct stat out_dir;
        if (stat(&M_outputFolder[0], &out_dir) != 0)
        {
            if (M_equationSystems.get_mesh().comm().rank() == 0)
            {
                //std::system("ls ");
                //std::system("echo $PWD");
                auto m = mkdir(M_outputFolder.c_str(), 0777);
                //std::cout << m << std::endl;
                if (m != 0) perror("mkdir");
            }
        }

        // read order of the basis functions
        std::string order = data(M_section + "/order", "FIRST");
        M_order = libMesh::Utility::string_to_enum < libMesh::Order > (order);
        // read order of finite element family (for DG)
        std::string fefamily = data(M_section + "/fefamily", "LAGRANGE");
        M_FEFamily = libMesh::Utility::string_to_enum < libMesh::FEFamily > (fefamily);
        std::cout << "* ElectroSolver: running with: " << fefamily << " (" << M_FEFamily << "), order: " << order << " (" << M_order << ")" << std::endl;
        // call setup system of the specific class
        setup_systems(M_datafile, M_section);
        // Add ionic current to this system
        IonicModelSystem& Iion_system = M_equationSystems.add_system < IonicModelSystem > ("iion");
        Iion_system.add_variable("iion", M_order, M_FEFamily);
        Iion_system.add_vector("ionic_model_map");
        Iion_system.add_vector("diion");
        Iion_system.add_vector("diion_old");
        Iion_system.add_vector("total_current");
        Iion_system.init();
        M_ionicModelExporterNames.insert("iion");

        // Add the applied current to this system
        IonicModelSystem& istim_system = M_equationSystems.add_system < IonicModelSystem > ("istim");
        istim_system.add_variable("istim", M_order, M_FEFamily);
        istim_system.init();
        istim_system.add_vector("stim_i");
        istim_system.add_vector("stim_e");
        istim_system.add_vector("surf_stim_i");
        istim_system.add_vector("surf_stim_e");

        M_ionicModelExporterNames.insert("istim");
        M_exporterNames.insert("istim");

        std::cout << "* ElectroSolver: Creating parameters spaces " << std::endl;
        ParameterSystem& activation_times_system = M_equationSystems.add_system < ParameterSystem > ("activation_times");
        activation_times_system.add_variable("activation_times", M_order, M_FEFamily);
        activation_times_system.init();
        // CV = Conduction Velocity
        ParameterSystem& CV_system = M_equationSystems.add_system < ParameterSystem > ("CV");
        CV_system.add_variable("cvx", libMesh::CONSTANT, libMesh::MONOMIAL);
        CV_system.add_variable("cvy", libMesh::CONSTANT, libMesh::MONOMIAL);
        CV_system.add_variable("cvz", libMesh::CONSTANT, libMesh::MONOMIAL);
        CV_system.init();

        if (!M_equationSystems.has_system("fibers"))
        {
            ParameterSystem& fiber_system = M_equationSystems.add_system < ParameterSystem > ("fibers");
            fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
            fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
            fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
            fiber_system.init();

        }

        if (!M_equationSystems.has_system("sheets"))
        {
            ParameterSystem& sheets_system = M_equationSystems.add_system < ParameterSystem > ("sheets");
            sheets_system.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
            sheets_system.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
            sheets_system.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
            sheets_system.init();
        }

        if (!M_equationSystems.has_system("xfibers"))
        {
            ParameterSystem& xfiber_system = M_equationSystems.add_system < ParameterSystem > ("xfibers");
            xfiber_system.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
            xfiber_system.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
            xfiber_system.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
            xfiber_system.init();
        }

        ParameterSystem& procID_system = M_equationSystems.add_system < ParameterSystem > ("ProcID");
        procID_system.add_variable("ProcID", libMesh::CONSTANT, libMesh::MONOMIAL);
        procID_system.init();

        M_parametersExporterNames.insert("activation_times");
        M_parametersExporterNames.insert("fibers");
        M_parametersExporterNames.insert("sheets");
        M_parametersExporterNames.insert("xfibers");
        M_parametersExporterNames.insert("ProcID");
        IonicModelSystem& model_map_system = M_equationSystems.add_system < IonicModelSystem > ("model_map");
        model_map_system.add_variable("ionic_model_map", M_order, M_FEFamily);
        model_map_system.init();
        M_parametersExporterNames.insert("model_map");
        ///////////
        // Time integrator
        int time_integrator_order = data(M_section + "/time_integrator_order", 1);
        bool isSecondORderImplemented = true;
        for (auto && map : M_ionicModelPtrMap)
        {
            if (map.second->isSecondOrderImplemented())
            {
                isSecondORderImplemented = false;
                break;
            }
        }
        if (2 == time_integrator_order && isSecondORderImplemented)
        {
            std::cout << "ELECTROSOLVER: using SBDF order 2 " << std::endl;
            M_timeIntegrator = TimeIntegrator::SecondOrderIMEX;
        }
        else
        {
            std::cout << "ELECTROSOLVER: using SBDF order 1 " << std::endl;
            M_timeIntegrator = TimeIntegrator::FirstOrderIMEX;
        }

        // ///////////////////////////////////////////////////////////////////////
        // ///////////////////////////////////////////////////////////////////////
        // Setup Exporters

        std::cout << "* " << M_model << ": Initializing exporters " << std::endl;
        M_exporter.reset(new Exporter(M_equationSystems.get_mesh()));
        //This works only for VTKIO
        M_exporter->set_compression(true);
        M_ionicModelExporter.reset(new Exporter(M_equationSystems.get_mesh()));
        M_ionicModelExporter->set_compression(true);
        M_parametersExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
        M_EXOExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
        M_potentialEXOExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
        M_nemesis_exporter.reset(new NemesisIO(M_equationSystems.get_mesh()));

        M_symmetricOperator = M_datafile(M_section + "/symmetric_operator", false);
        std::cout << "* ElectroSolver: Using Symmetric Operator: " << M_symmetricOperator << std::endl;
    }

    void ElectroSolver::init(double time)
    {
        std::cout << "* ElectroSolver: Init: " << std::endl;

        std::cout << "* ElectroSolver: Setup ionic model map key vector: " << std::endl;
        const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
        libMesh::MeshBase::const_element_iterator el_start = mesh.active_local_elements_begin();
        libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        IonicModelSystem& iion_system = M_equationSystems.get_system < IonicModelSystem > ("iion");
        std::vector < libMesh::dof_id_type > dofs;

        IonicModelSystem& model_map_system = M_equationSystems.get_system < IonicModelSystem > ("model_map");
        M_parametersExporterNames.insert("model_map");

        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;
            const unsigned int bid = elem->subdomain_id();
            if (M_tissueBlockIDs.size() > 0)
            {
                auto it_bid = M_tissueBlockIDs.find(bid);
                if (it_bid != M_tissueBlockIDs.end())
                {
                    iion_system.get_dof_map().dof_indices(elem, dofs);
                    for (auto && d : dofs)
                    {
                        iion_system.get_vector("ionic_model_map").set(d, bid);
                        model_map_system.solution->set(d, bid);
                    }
                }
            }
            else
            {
                for (auto && d : dofs)
                {
                    iion_system.get_vector("ionic_model_map").set(d, secret_blockID_key);
                    model_map_system.solution->set(d, secret_blockID_key);
                }
            }
        }
        if (M_tissueBlockIDs.size() > 0)
        {
            auto first_local_index = iion_system.get_vector("ionic_model_map").first_local_index();
            auto last_local_index = iion_system.get_vector("ionic_model_map").last_local_index();
            for (auto i = first_local_index; i < last_local_index; ++i)
            {
                double value = iion_system.get_vector("ionic_model_map")(i);
                if (value == 0) iion_system.get_vector("ionic_model_map").set(i, -1);
            }

        }

        iion_system.get_vector("ionic_model_map").close();
        std::cout << "* ElectroSolver: Setup ionic model vector: " << std::endl;

        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        // Setting initial conditions
        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);

        libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
        const libMesh::MeshBase::const_node_iterator end_node = mesh.local_nodes_end();
        const libMesh::DofMap & dof_map = system.get_dof_map();
        const libMesh::DofMap & dof_map_V = wave_system.get_dof_map();

        std::vector < libMesh::dof_id_type > dof_indices_V;
        std::vector < libMesh::dof_id_type > dof_indices_Q;
        std::vector < libMesh::dof_id_type > dof_indices_gating;

        int c = 0;
        // loop over nodes
        for (; node != end_node; ++node)
        {
            const libMesh::Node * nn = *node;
            // Are we in the bath?
            auto n_var = nn->n_vars(system.number());
            auto n_dofs = nn->n_dofs(system.number());
            if (n_var == n_dofs)
            {
                dof_map.dof_indices(nn, dof_indices_Q, 0);
                dof_map_V.dof_indices(nn, dof_indices_V, 0);
                int key_iion = iion_system.get_vector("ionic_model_map")(dof_indices_V[0]);
                auto it_ionic_model = M_ionicModelPtrMap.find(key_iion);
                auto it_ionic_model_name = M_ionicModelNameMap.find(key_iion);
                IonicModel * ionicModelPtr = nullptr;
                std::string ionic_model_system_name;
                if (it_ionic_model != M_ionicModelPtrMap.end()) ionicModelPtr = it_ionic_model->second.get();
                if (it_ionic_model_name != M_ionicModelNameMap.end()) ionic_model_system_name = it_ionic_model_name->second;
                if (ionicModelPtr)
                {
                    IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > (ionic_model_system_name);
                    auto& iion_rhs_old = ionic_model_system.get_vector("rhs_old");
                    int num_vars = ionic_model_system.n_vars();
                    std::vector<double> init_values(num_vars + 1, 0.0);
                    const libMesh::DofMap & dof_map_gating = ionic_model_system.get_dof_map();
                    dof_map_gating.dof_indices(nn, dof_indices_gating);
                    ionicModelPtr->initialize(init_values);

                    wave_system.solution->set(dof_indices_V[0], init_values[0]);
                    for (int nv = 0; nv < num_vars; ++nv)
                    {
//                        std::cout << "Setting gatings: " << nv << std::endl;
//                        std::cout << "value: " << std::flush;
//                        std::cout << init_values[nv + 1] << std::flush;
                        ionic_model_system.solution->set(dof_indices_gating[nv], init_values[nv + 1]);
//                        std::cout << " set!" << std::endl;
                    }
                }
                else
                {
                    for(auto && ionicModel : M_ionicModelPtrMap)
                    {
                        std::cout << ionicModel.first << ", " << ionicModel.second->ionicModelName() << std::endl;
                    }
                	std::cout << "NO IONIC MODEL PTR" << std::endl;
                	throw std::runtime_error("No ionic Model Ptr");
                }
            }
        }
        wave_system.solution->close();

        for (auto && name : M_ionic_models_systems_name_vec)
        {
            IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > (name);
            ionic_model_system.solution->close();
        }        //Initialize ionic model
//        for(int c = 0; c < M_ionicModelNameMap.size(); ++c)
//        for (auto && name : M_ionic_models_vec)
//            for(auto && map : M_ionicModelNameMap)
//            {
//            std::vector<std::string>::iterator it = std::find( M_ionic_models_vec.begin(),
//                                                               M_ionic_models_vec.end(),
//                                                               map.second );
//
//            IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > (map.second);
//            int num_vars = ionic_model_system.n_vars();
//            std::vector<double> init_values(num_vars + 1, 0.0);
//            M_ionicModelPtrMap[map.first]->initialize(init_values);
//
//            auto first_local_index = wave_system.solution->first_local_index();
//            auto last_local_index = wave_system.solution->last_local_index();
//            for (auto i = first_local_index; i < last_local_index; ++i)
//            {
//                wave_system.solution->set(i, init_values[0]);
//            }
//            first_local_index = ionic_model_system.solution->first_local_index();
//            last_local_index = ionic_model_system.solution->last_local_index();
//            std::cout << "* " << M_model << ": Setting initial values ... " << std::flush;
//            for (auto i = first_local_index; i < last_local_index;)
//            {
//                for (int nv = 0; nv < num_vars; ++nv)
//                {
//                    ionic_model_system.solution->set(i, init_values[nv + 1]);
//                    ++i;
//                }
//            }
//            ionic_model_system.solution->close();
//            ionic_model_system.old_local_solution->close();
//            ionic_model_system.older_local_solution->close();
//
//        }

//        std::cout << "* Old stuff: " << std::endl;

        system.solution->close();
        system.old_local_solution->close();
        system.older_local_solution->close();

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

        ParameterSystem& procID_system = M_equationSystems.get_system < ParameterSystem > ("ProcID");

        auto first_local_index = procID_system.solution->first_local_index();
        auto last_local_index = procID_system.solution->last_local_index();

        double myrank = static_cast<double>(M_equationSystems.comm().rank());
        for (int el = first_local_index; el < last_local_index; ++el)
        {
            procID_system.solution->set(el, myrank);
        }
        procID_system.solution->close();

        // Call specific system initializations
        init_systems(time);
//    M_linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( M_equationSystems.comm() );
        typedef libMesh::PetscLinearSolver<libMesh::Number> PetscSolver;
        M_linearSolver.reset(new PetscSolver(M_equationSystems.comm()));
        std::string solver_type = M_datafile(M_section + "/linear_solver/type", "gmres");
        std::cout << "* ElectroSolver: using " << solver_type << std::endl;
        std::string prec_type = M_datafile(M_section + "/linear_solver/preconditioner", "amg");
        //M_linearSolver->set_solver_type(solver_map.find(solver_type)->second);
        //M_linearSolver->set_preconditioner_type(prec_map.find(prec_type)->second);
        M_linearSolver->init();

        std::cout << "* ElectroSolver: Init complete " << std::endl;

    }

    void ElectroSolver::init_systems(double time)
    {
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system<ElectroSystem>("wave");

        std::string v_ic = M_datafile(M_section + "/ic", "");
        if (v_ic != "")
        {
            std::cout << "* ElectroSolver: Found bidomain initial condition: " << v_ic << std::endl;
            SpiritFunction bidomain_ic;
            bidomain_ic.read(v_ic);
            setup_ic(bidomain_ic);
        }

        IonicModelSystem& istim_system = M_equationSystems.get_system<IonicModelSystem>("istim");
        std::cout << "* ElectroSolver: Initializing activation times to -1  ... " << std::flush;
        ParameterSystem& activation_times_system = M_equationSystems.get_system<ParameterSystem>("activation_times");
        auto first_local_index = activation_times_system.solution->first_local_index();
        auto last_local_index = activation_times_system.solution->last_local_index();

        for (auto i = first_local_index; i < last_local_index; ++i)
        {
            activation_times_system.solution->set(i, -1.0);
        }
        std::cout << " done" << std::endl;

        bool analytic_fibers = M_datafile(M_section + "/analytic_fibers", true);
        if (analytic_fibers)
        {
            std::string fibers_data = M_datafile(M_section + "/fibers", "1.0, 0.0, 0.0");
            std::string sheets_data = M_datafile(M_section + "/sheets", "0.0, 1.0, 0.0");
            std::string xfibers_data = M_datafile(M_section + "/xfibers", "0.0, 0.0, 1.0");

            SpiritFunction fibers_func;
            SpiritFunction sheets_func;
            SpiritFunction xfibers_func;

            fibers_func.read(fibers_data);
            sheets_func.read(sheets_data);
            xfibers_func.read(xfibers_data);

            ParameterSystem& fiber_system = M_equationSystems.get_system<ParameterSystem>("fibers");
            ParameterSystem& sheets_system = M_equationSystems.get_system<ParameterSystem>("sheets");
            ParameterSystem& xfiber_system = M_equationSystems.get_system<ParameterSystem>("xfibers");

            std::cout << "* ElectroSolver:: Creating fibers from function: " << fibers_data << std::flush;
            fiber_system.project_solution(&fibers_func);
            std::cout << " done" << std::endl;
            std::cout << "* ElectroSolver:: Creating sheets from function: " << sheets_data << std::flush;
            sheets_system.project_solution(&sheets_func);
            std::cout << " done" << std::endl;
            std::cout << "* ElectroSolver:: Creating xfibers from function: " << xfibers_data << std::flush;
            xfiber_system.project_solution(&xfibers_func);
            std::cout << " done" << std::endl;
        }
        std::string conductivity_type = M_datafile(M_section + "/conductivity", "function");

        std::cout << "Conductivity Type: " << conductivity_type << std::endl;

        std::string Dffi_data = M_datafile(M_section + "/Dffi", "2.0");
        std::string Dssi_data = M_datafile(M_section + "/Dssi", "0.2");
        std::string Dnni_data = M_datafile(M_section + "/Dnni", "0.2");
        std::string Dffe_data = M_datafile(M_section + "/Dffe", "1.5");
        std::string Dsse_data = M_datafile(M_section + "/Dsse", "1.0");
        std::string Dnne_data = M_datafile(M_section + "/Dnne", "1.0");
        if(M_model == "monowave")
        {
            ParameterSystem& intra_conductivity_system = M_equationSystems.get_system<ParameterSystem>("conductivity");

            if ("function" == conductivity_type)
            {
                SpiritFunction intra_conductivity_func;
                std::cout << "Dffi_data: " << Dffi_data << std::endl;
                intra_conductivity_func.add_function(Dffi_data);
                intra_conductivity_func.add_function(Dssi_data);
                intra_conductivity_func.add_function(Dnni_data);
                intra_conductivity_system.project_solution(&intra_conductivity_func);
            }
            //list
            else
            {
                std::vector<double> Dffi;
                BeatIt::readList(Dffi_data, Dffi);
                std::vector<double> Dssi;
                BeatIt::readList(Dssi_data, Dssi);
                std::vector<double> Dnni;
                BeatIt::readList(Dnni_data, Dnni);

                std::string i_IDs = M_datafile(M_section + "/i_IDs", "-1");
                std::cout << "i_IDs: " << i_IDs << std::flush;
                std::vector<unsigned int> intracellular_IDs;
                BeatIt::readList(i_IDs, intracellular_IDs);

                const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();

                libMesh::MeshBase::const_element_iterator el_start = mesh.active_local_elements_begin();
                libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
                const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

                const libMesh::DofMap & dof_map = intra_conductivity_system.get_dof_map();
                std::vector<libMesh::dof_id_type> dof_indices;

                for (; el != end_el; ++el)
                {
                    const libMesh::Elem * elem = *el;
                    const unsigned int elem_id = elem->id();
                    int subdomain_ID = elem->subdomain_id();
                    dof_map.dof_indices(elem, dof_indices);

                    double cDffi = 0.;
                    double cDssi = 0.;
                    double cDnni = 0.;
                    for (int k = 0; k < intracellular_IDs.size(); ++k)
                    {
                        if (intracellular_IDs[k] == subdomain_ID)
                        {
                            cDffi = Dffi[k];
                            cDssi = Dssi[k];
                            cDnni = Dnni[k];
                            break;
                        }
                    }
                    intra_conductivity_system.solution->set(dof_indices[0], cDffi);
                    intra_conductivity_system.solution->set(dof_indices[1], cDssi);
                    intra_conductivity_system.solution->set(dof_indices[2], cDnni);
                }

            }
        }
        else
        {
            ParameterSystem& intra_conductivity_system = M_equationSystems.get_system<ParameterSystem>("intra_conductivity");
            ParameterSystem& extra_conductivity_system = M_equationSystems.get_system<ParameterSystem>("extra_conductivity");

            if ("function" == conductivity_type)
            {
                SpiritFunction intra_conductivity_func;
                std::cout << "Dffi_data: " << Dffi_data << std::endl;
                intra_conductivity_func.add_function(Dffi_data);
                intra_conductivity_func.add_function(Dssi_data);
                intra_conductivity_func.add_function(Dnni_data);
                intra_conductivity_system.project_solution(&intra_conductivity_func);
                SpiritFunction extra_conductivity_func;
                extra_conductivity_func.add_function(Dffe_data);
                extra_conductivity_func.add_function(Dsse_data);
                extra_conductivity_func.add_function(Dnne_data);
                extra_conductivity_system.project_solution(&extra_conductivity_func);
            }
            //list
            else
            {
                std::vector<double> Dffi;
                BeatIt::readList(Dffi_data, Dffi);
                std::vector<double> Dssi;
                BeatIt::readList(Dssi_data, Dssi);
                std::vector<double> Dnni;
                BeatIt::readList(Dnni_data, Dnni);

                std::string i_IDs = M_datafile(M_section + "/i_IDs", "-1");
                std::cout << "i_IDs: " << i_IDs << std::flush;
                std::vector<unsigned int> intracellular_IDs;
                BeatIt::readList(i_IDs, intracellular_IDs);
                //for ( auto && v : i_IDs) std::cout << v << std::flush;
                std::vector<double> Dffe;
                BeatIt::readList(Dffe_data, Dffe);
                std::vector<double> Dsse;
                BeatIt::readList(Dsse_data, Dsse);
                std::vector<double> Dnne;
                BeatIt::readList(Dnne_data, Dnne);

                std::string e_IDs = M_datafile(M_section + "/e_IDs", "-1");
                std::cout << "\ne_IDs: " << e_IDs << std::endl;

                std::vector<unsigned int> extracellular_IDs;
                BeatIt::readList(e_IDs, extracellular_IDs);
                //for ( auto && v : e_IDs) std::cout << v << std::endl;

                const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();

                libMesh::MeshBase::const_element_iterator el_start = mesh.active_local_elements_begin();
                libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
                const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

                const libMesh::DofMap & dof_map = extra_conductivity_system.get_dof_map();
                std::vector<libMesh::dof_id_type> dof_indices;

                for (; el != end_el; ++el)
                {
                    const libMesh::Elem * elem = *el;
                    const unsigned int elem_id = elem->id();
                    int subdomain_ID = elem->subdomain_id();
                    dof_map.dof_indices(elem, dof_indices);

                    double cDffi = 0.;
                    double cDssi = 0.;
                    double cDnni = 0.;
                    for (int k = 0; k < intracellular_IDs.size(); ++k)
                    {
                        if (intracellular_IDs[k] == subdomain_ID)
                        {
                            cDffi = Dffi[k];
                            cDssi = Dssi[k];
                            cDnni = Dnni[k];
                            break;
                        }
                    }
                    double cDffe = 0.;
                    double cDsse = 0.;
                    double cDnne = 0.;
                    for (int k = 0; k < extracellular_IDs.size(); ++k)
                    {
                        if (extracellular_IDs[k] == subdomain_ID)
                        {
                            cDffe = Dffe[k];
                            cDsse = Dsse[k];
                            cDnne = Dnne[k];
                            break;
                        }
                    }
                    extra_conductivity_system.solution->set(dof_indices[0], cDffe);
                    extra_conductivity_system.solution->set(dof_indices[1], cDsse);
                    extra_conductivity_system.solution->set(dof_indices[2], cDnne);
                    intra_conductivity_system.solution->set(dof_indices[0], cDffi);
                    intra_conductivity_system.solution->set(dof_indices[1], cDssi);
                    intra_conductivity_system.solution->set(dof_indices[2], cDnni);
                }

            }
        }

    }

    void ElectroSolver::init_endocardial_ve(std::set<libMesh::boundary_id_type>& IDs, std::set<unsigned short>& subdomainIDs)
    {
        std::cout << "* ElectroSolver: Initializing endocardial mesh " << std::endl;

        M_boundary_ve.M_endocardium.reset(new libMesh::BoundaryMesh(M_equationSystems.comm()));
        M_equationSystems.get_mesh().boundary_info->sync(IDs, *(M_boundary_ve.M_endocardium), subdomainIDs);
        M_boundary_ve.M_endocardium->allow_renumbering(false);
        M_boundary_ve.M_endocardium->prepare_for_use();
        //std::map< libMesh::dof_id_type, libMesh::dof_id_type > node_id_map;
        std::map<libMesh::dof_id_type, unsigned char> side_id_map;
        M_equationSystems.get_mesh().boundary_info->get_side_and_node_maps(*M_boundary_ve.M_endocardium, M_boundary_ve.M_node_id_map, side_id_map);
        for (auto it = M_boundary_ve.M_node_id_map.begin(); it != M_boundary_ve.M_node_id_map.end(); ++it)
        {
            M_boundary_ve.M_reverse_node_id_map[it->second] = it->first;
        }

        std::cout << "* ElectroSolver: creating endocardial equation systems " << std::endl;
        M_boundary_ve.M_boundary_es.reset(new libMesh::EquationSystems(*M_boundary_ve.M_endocardium));
        std::cout << "* ElectroSolver: creating new endocardial system " << std::endl;
        auto& ve_sys = M_boundary_ve.M_boundary_es->add_system < libMesh::ExplicitSystem > ("ve");
        ve_sys.add_variable("phie", libMesh::FIRST);
        std::cout << "* ElectroSolver: initialize endocardial equation systems " << std::endl;
        M_boundary_ve.M_boundary_es->init();
        std::cout << "* ElectroSolver: create new exporter only for the endocardium ... " << std::flush;
        M_boundary_ve.M_EXOExporter.reset(new libMesh::ExodusII_IO(*M_boundary_ve.M_endocardium));
        std::cout << " done." << std::endl;

    }

    void ElectroSolver::setup_ic(libMesh::FunctionBase<libMesh::Number>& ic, double time) // setup initial conditions
    {
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        M_equationSystems.parameters.set < libMesh::Real > ("time") = time;
        wave_system.time = time;
        std::cout << "* BIDOMAIN+BATH: Projecting initial condition to bidomain system ... " << std::flush;
        wave_system.project_solution(&ic);
        std::cout << "* BIDOMAIN+BATH: Copying initial conditions to vectors at n nd at n-1... " << std::flush;
        // Close vectors and update the values in old_local_solution and older_local_solution
        advance();
        std::cout << " done" << std::endl;
    }

    void ElectroSolver::restart(EXOExporter& importer, int step, bool restart)
    {
        if (restart)
        {
            auto& nodal_var_names = importer.get_nodal_var_names();
            auto nodal_first = nodal_var_names.begin();
            auto nodal_end = nodal_var_names.end();
            auto& elem_var_names = importer.get_elem_var_names();
            auto elem_first = elem_var_names.begin();
            auto elem_end = elem_var_names.end();

            int num_systems = M_equationSystems.n_systems();
            for (int k = 0; k < num_systems; ++k)
            {
                libMesh::System& system = M_equationSystems.get_system(k);
                std::string name = system.name();
                std::cout << "Importing System: " << name << std::endl;
                int n_vars = system.n_vars();
                for (int l = 0; l < n_vars; ++l)
                {
                    std::string var_name = system.variable_name(l);

                    auto elem_it = std::find(elem_first, elem_end, var_name);
                    if (elem_it != elem_end)
                    {
                        std::cout << "\t elemental variable: " << *elem_it << std::endl;
                        importer.copy_elemental_solution(system, *elem_it, *elem_it, step);
                    }
                    else
                    {
                        auto nodal_it = std::find(nodal_first, nodal_end, var_name);
                        if (nodal_it != nodal_end)
                        {
                            std::cout << "\t nodal variable: " << *nodal_it << std::endl;
                            importer.copy_nodal_solution(system, *nodal_it, *nodal_it, step);
                        }
                    }
                }
            }
        }
    }

    void ElectroSolver::read_fibers(EXOExporter& importer, int step)
    {
        const int num_fiber_systems = 3;
        std::vector < std::string > fibers(num_fiber_systems);
        fibers[0] = "fibers";
        fibers[1] = "sheets";
        fibers[2] = "xfibers";

        auto& elem_var_names = importer.get_elem_var_names();
        auto elem_first = elem_var_names.begin();
        auto elem_end = elem_var_names.end();

        for (int k = 0; k < num_fiber_systems; ++k)
        {
            libMesh::System& system = M_equationSystems.get_system(fibers[k]);
            std::string name = system.name();
            std::cout << "Importing System: " << name << std::endl;
            int n_vars = system.n_vars();
            for (int l = 0; l < n_vars; ++l)
            {
                std::string var_name = system.variable_name(l);
                auto elem_it = std::find(elem_first, elem_end, var_name);
                if (elem_it != elem_end)
                {
                    std::cout << "\t elemental variable: " << *elem_it << std::endl;
                    importer.copy_elemental_solution(system, *elem_it, *elem_it, step);
                }
            }
        }
    }

    void ElectroSolver::init_exo_output()
    {
        std::cout << "* " << M_model << ": EXODUSII::Exporting " << M_model << ".exo at time 0.0" << " in: " << M_outputFolder << " ... " << std::flush;

        M_EXOExporter->write_equation_systems(M_outputFolder + M_model + ".exo", M_equationSystems);
        M_EXOExporter->append(true);
        M_EXOExporter->write_element_data(M_equationSystems);
        std::vector < std::string > output;
        output.push_back("V");
        output.push_back("Q");
        output.push_back("Ve");
        M_potentialEXOExporter->set_output_variables(output);
        M_potentialEXOExporter->write_equation_systems(M_outputFolder + "potential.exo", M_equationSystems);
        M_potentialEXOExporter->append(true);
        std::cout << " done!" << std::endl;
    }

    void ElectroSolver::save_ve_timestep(int step, double time)
    {
        std::cout << "* ElectroSolver: copying data to endocardial surface " << std::endl;

        std::vector < libMesh::dof_id_type > dof_indices;
        std::vector < libMesh::dof_id_type > dof_indices_endo;

        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);
        auto& boundary_sys = M_boundary_ve.M_boundary_es->get_system < libMesh::ExplicitSystem > ("ve");
        const libMesh::DofMap & dof_map = system.get_dof_map();
        const libMesh::DofMap & dof_map_endo = boundary_sys.get_dof_map();

//    libMesh::MeshBase::const_node_iterator node = M_boundary_ve.M_endocardium->local_nodes_begin();
//    const libMesh::MeshBase::const_node_iterator end_node =
//            M_boundary_ve.M_endocardium->local_nodes_end();
//    for( ; node  != end_node; ++node)
//    {
//        const libMesh::Node * endo_nn = *node;
//        dof_map_endo.dof_indices(endo_nn, dof_indices_endo, 0);
//
//        int endo_node_id = endo_nn->id();
//        auto it = M_boundary_ve.M_reverse_node_id_map.find(endo_node_id);
//        if(it != M_boundary_ve.M_reverse_node_id_map.end( ))
//        {
//            int id = it->second;
//            const libMesh::Node * nn = M_equationSystems.get_mesh().node_ptr(  id );
//            dof_map.dof_indices(nn, dof_indices, 1);
//
//            //std::cout << "get_ve: " << dof_indices.size() << std::endl;
//            double ve = (*system.solution)(dof_indices[0]);
//            //std::cout << "set_ve" << std::endl;
//            boundary_sys.solution->set(dof_indices_endo[0], ve);
//        }
//
//
//    }

        //std::cout << "* ElectroSolver: Loop over the map " << std::endl;
        // Loop over node of the
        for (auto it = M_boundary_ve.M_node_id_map.begin(); it != M_boundary_ve.M_node_id_map.end(); ++it)
        {
            //std::cout << "nn" << std::endl;
            const libMesh::Node * nn = M_equationSystems.get_mesh().node_ptr(it->first);
            dof_map.dof_indices(nn, dof_indices, 1);
            if (M_equationSystems.comm().rank() == nn->processor_id())
            {

                //std::cout << "endo_nn" << std::endl;
                const libMesh::Node * endo_nn = M_boundary_ve.M_endocardium->node_ptr(it->second);
                dof_map_endo.dof_indices(endo_nn, dof_indices_endo, 0);

                //std::cout << "get_ve: " << dof_indices.size() << std::endl;
                double ve = (*system.solution)(dof_indices[0]);
                //std::cout << "set_ve" << std::endl;
                boundary_sys.solution->set(dof_indices_endo[0], ve);
            }
        }
        std::cout << "* ElectroSolver: Closing solution " << std::endl;

        boundary_sys.solution->close();

        std::cout << "* " << M_model << ": EXODUSII::Exporting " << M_model << "_endo.exo at time " << time << " in: " << M_outputFolder << " ... " << std::flush;
        if (1 == step)
        {
            M_boundary_ve.M_EXOExporter->write_equation_systems(M_outputFolder + M_model + "_endo.exo", *M_boundary_ve.M_boundary_es);
            M_boundary_ve.M_EXOExporter->append(true);
        }
        M_boundary_ve.M_EXOExporter->write_timestep(M_outputFolder + M_model + "_endo.exo", *M_boundary_ve.M_boundary_es, step, time);

        std::cout << "done " << std::endl;

    }

    void ElectroSolver::save_exo_timestep(int step, double time)
    {
        std::cout << "* " << M_model << ": EXODUSII::Exporting " << M_model << ".exo at time " << time << " in: " << M_outputFolder << " ... " << std::flush;
        M_EXOExporter->write_timestep(M_outputFolder + M_model + ".exo", M_equationSystems, step, time);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save_parameters()
    {
//        std::cout << "* " << M_model << ": VTKIO::Exporting parameters in: " << M_outputFolder << " ... " << std::flush;
//        Exporter vtk(M_equationSystems.get_mesh());
//        vtk.write_equation_systems(M_outputFolder + "parameters.pvtu", M_equationSystems, &M_parametersExporterNames);
//        std::cout << "done " << std::endl;
        std::cout << "* " << M_model << ": EXODUSII::Exporting parameters in: " << M_outputFolder << " ... " << std::flush;
        EXOExporter exo(M_equationSystems.get_mesh());
        exo.write_equation_systems(M_outputFolder+"parameters.exo", M_equationSystems, &M_parametersExporterNames);
//        exo.write(M_outputFolder + "parameters.exo");
        exo.write_element_data(M_equationSystems);
        std::cout << "done " << std::endl;

    }

    void ElectroSolver::save_potential(int step, double time)
    {
        std::cout << "* " << M_model << ": VTKIO::Exporting potential*.pvtu at step " << step << " for time: " << time << " in: " << M_outputFolder << " ... " << std::flush;

        //M_potentialEXOExporter->write_timestep(M_outputFolder + "potential.exo", M_equationSystems, step, time);
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << step;
        std::string step_str = ss.str();
        M_exporter->write_equation_systems(M_outputFolder + "potential_" + step_str + ".pvtu", M_equationSystems, &M_exporterNames);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save_potential_nemesis(int step, double time)
    {

        //M_potentialEXOExporter->write_timestep(M_outputFolder + "potential.exo", M_equationSystems, step, time);
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << step;
        std::string step_str = ss.str();
        std::string filename = M_outputFolder + "potential_" + step_str + ".e";
        std::cout << "* " << M_model << ": NEMESIS_IO::Exporting potential*.e at step " << step << " for time: " << time << " in: " << M_outputFolder << " ... " << std::flush;
        std::cout << "\n" << filename << std::endl;
        M_nemesis_exporter.reset(new NemesisIO(M_equationSystems.get_mesh()));
        M_nemesis_exporter->write_equation_systems(filename, M_equationSystems, &M_exporterNames);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save_activation_times(int step)
    {
        std::cout << "* " << M_model << ": VTKIO::Exporting activaton times: " << M_outputFolder << " ... " << std::flush;
        Exporter vtk(M_equationSystems.get_mesh());
        std::set < std::string > output;
        output.insert("activation_times");

        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << step;
        std::string step_str = ss.str();
        vtk.write_equation_systems(M_outputFolder + "activation_times_" + step_str + ".pvtu", M_equationSystems, &output);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save_conduction_velocity(int step)
    {
        std::cout << "* " << M_model << ": VTKIO::Exporting Conduction Velocity: " << M_outputFolder << " ... " << std::flush;
        Exporter vtk(M_equationSystems.get_mesh());
        std::set < std::string > output;
        output.insert("CV");
        vtk.write_equation_systems(M_outputFolder + "CV" + std::to_string(step) + ".pvtu", M_equationSystems, &output);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save(int step)
    {
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << step;
        std::string step_str = ss.str();

        //save in subfolder

        std::cout << "* " << M_model << ": VTKIO::Exporting " << step << " in: " << M_outputFolder << " ... " << std::flush;
        M_exporter->write_equation_systems(M_outputFolder + M_model + "_" + step_str + ".pvtu", M_equationSystems, &M_exporterNames);
        M_ionicModelExporter->write_equation_systems(M_outputFolder + "ionic_model_" + step_str + ".pvtu", M_equationSystems, &M_ionicModelExporterNames);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::reinit_linear_solver()
    {
        M_linearSolver->clear();
        M_linearSolver->set_solver_type(libMesh::CG);
        M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
        M_linearSolver->init();
    }

    void ElectroSolver::evaluate_conduction_velocity()
    {
        ParameterSystem& activation_times_system = M_equationSystems.get_system < ParameterSystem > ("activation_times");
        activation_times_system.update();
        // WAVE
        ParameterSystem& CV_system = M_equationSystems.get_system < ParameterSystem > ("CV");

        const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

        const libMesh::DofMap & dof_map_CV = CV_system.get_dof_map();
        const libMesh::DofMap & dof_map_at = activation_times_system.get_dof_map();

        std::vector < libMesh::dof_id_type > dof_indices_CV;
        std::vector < libMesh::dof_id_type > dof_indices_at;

        std::vector<double> at;

        libMesh::FEType fe_type = dof_map_at.variable_type(0);
        std::unique_ptr < libMesh::FEBase > fe_qp1(libMesh::FEBase::build(dim, fe_type));
        // A 5th order Gauss quadrature rule for numerical integration.
        libMesh::QGauss qrule(dim, libMesh::FIRST);
        // Tell the finite element object to use our quadrature rule.
        fe_qp1->attach_quadrature_rule(&qrule);
        // The element shape functions evaluated at the quadrature points.
        const std::vector<std::vector<libMesh::Real> > & phi = fe_qp1->get_phi();

        // The element shape function gradients evaluated at the quadrature
        // points.
        const std::vector<std::vector<libMesh::RealGradient> > & dphi = fe_qp1->get_dphi();

        libMesh::VectorValue < libMesh::Number > CV;
        libMesh::VectorValue < libMesh::Number > CV_direction;
        libMesh::VectorValue < libMesh::Number > grad_at;

        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;
            dof_map_at.dof_indices(elem, dof_indices_at);
            dof_map_CV.dof_indices(elem, dof_indices_CV);

            const unsigned int n_dofs = dof_indices_at.size();
            at.resize(n_dofs);
            activation_times_system.current_local_solution->get(dof_indices_at, at);
            fe_qp1->reinit(elem);
            CV *= 0.0;
            CV_direction *= 0.0;
            grad_at *= 0.0;
            const unsigned int n_phi = phi.size();
            // use single quadrature point;
            int qp = 0;
            for (unsigned int l = 0; l < n_phi; l++)
            {
                grad_at += dphi[l][qp] * at[l];
            }
            double grad_at_mag2 = grad_at.norm_sq();
            double grad_at_mag = std::sqrt(grad_at_mag2);

            double cv0 = 1000 * grad_at(0) / grad_at_mag2;
            double cv1 = 1000 * grad_at(1) / grad_at_mag2;
            double cv2 = 1000 * grad_at(2) / grad_at_mag2;
            // std::cout << "grad_at_mag2: " << grad_at_mag2 << ", cv = " << cv0 << ", " << cv1 << ", " << cv2 << std::endl;

            CV_system.solution->set(dof_indices_CV[0], cv0);
            CV_system.solution->set(dof_indices_CV[1], cv1);
            CV_system.solution->set(dof_indices_CV[2], cv2);
        }
        // WAVE

    }

    void ElectroSolver::update_activation_time(double time, double threshold)
    {
        ParameterSystem& activation_times_system = M_equationSystems.get_system < ParameterSystem > ("activation_times");
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        // Setting initial conditions
        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);

        const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();

        libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
        const libMesh::MeshBase::const_node_iterator end_node = mesh.local_nodes_end();

        const libMesh::DofMap & dof_map = wave_system.get_dof_map();
        const libMesh::DofMap & dof_map_q = system.get_dof_map();
        const libMesh::DofMap & dof_map_at = activation_times_system.get_dof_map();

        std::vector < libMesh::dof_id_type > dof_indices_Q;
        std::vector < libMesh::dof_id_type > dof_indices_V;
        std::vector < libMesh::dof_id_type > dof_indices_at;

        double v = 0.0;
        double q = 0.0;
        double at = 0.0;

        for (; node != end_node; ++node)
        {
            const libMesh::Node * nn = *node;
            // Are we in the bath?

            dof_map.dof_indices(nn, dof_indices_V, 0);
            dof_map_at.dof_indices(nn, dof_indices_at, 0);
            dof_map_q.dof_indices(nn, dof_indices_Q, 0);

            if (dof_indices_V.size() > 0)
            {
                v = (*wave_system.solution)(dof_indices_V[0]);
                at = (*activation_times_system.solution)(dof_indices_at[0]);
                q = (*system.solution)(dof_indices_Q[0]);
                if (v > threshold && at < 0.0)
                {
                    activation_times_system.solution->set(dof_indices_at[0], time);
                }
                else if( v <= -40.0 && at > 0.0)
                {
                    activation_times_system.solution->set(dof_indices_at[0], -1.0);
                }
            }
        }

        activation_times_system.solution->close();
        activation_times_system.update();
    }

    void ElectroSolver::advance()
    {
        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);

        system.solution->close();
        system.old_local_solution->close();
        system.older_local_solution->close();
        *system.older_local_solution = *system.old_local_solution;
        *system.old_local_solution = *system.solution;
        system.update();

        for (auto && name : M_ionic_models_systems_name_vec)
        {
            IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > (name);
            ionic_model_system.solution->close();
            ionic_model_system.old_local_solution->close();
            ionic_model_system.older_local_solution->close();
            *ionic_model_system.older_local_solution = *ionic_model_system.old_local_solution;
            *ionic_model_system.old_local_solution = *ionic_model_system.solution;
            ionic_model_system.get_vector("rhs_old") = *ionic_model_system.rhs;
            ionic_model_system.update();
        }
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        wave_system.solution->close();
        wave_system.old_local_solution->close();
        wave_system.older_local_solution->close();
        *wave_system.older_local_solution = *wave_system.old_local_solution;
        *wave_system.old_local_solution = *wave_system.solution;
        wave_system.update();

        IonicModelSystem& iion_system = M_equationSystems.get_system < IonicModelSystem > ("iion");
        iion_system.solution->close();
        iion_system.old_local_solution->close();
        iion_system.older_local_solution->close();
        *iion_system.older_local_solution = *iion_system.old_local_solution;
        *iion_system.old_local_solution = *iion_system.solution;
        iion_system.get_vector("diion_old") = iion_system.get_vector("diion");
        iion_system.update();
    }

    std::string ElectroSolver::get_ionic_model_name(unsigned int key) const
    {
        auto it = M_ionicModelPtrMap.find(key);
        if (it != M_ionicModelPtrMap.end()) return it->second->ionicModelName();
        else throw std::runtime_error("ElectroSolver::get_ionic_model_name");
    }

    double ElectroSolver::last_activation_time()
    {
        ParameterSystem& activation_times_system = M_equationSystems.get_system < ParameterSystem > ("activation_times");
        return static_cast<double>(activation_times_system.solution->max());
    }

    double ElectroSolver::potential_norm()
    {
        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);
        return system.solution->l1_norm();
    }

    void ElectroSolver::set_potential_on_boundary(unsigned int boundID, double value, int subdomain)
    {
        std::cout << "* " << M_model << ": WARNING:  set_potential_on_boundary works only for TET4" << std::flush;
        libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");

        const libMesh::DofMap & dof_map = wave_system.get_dof_map();
        std::vector < libMesh::dof_id_type > dof_indices;

        libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

        if (subdomain >= 0)
        {
            el = mesh.active_local_subdomain_elements_begin(subdomain);
            end_el = mesh.active_local_subdomain_elements_end(subdomain);
        }

        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;

            {
                dof_map.dof_indices(elem, dof_indices);
                for (unsigned int side = 0; side < elem->n_sides(); side++)
                {
                    if (elem->neighbor_ptr(side) == libmesh_nullptr)
                    {
                        unsigned int n_boundary_ids=mesh.boundary_info->n_boundary_ids(elem,side);
                        std::vector<short int> boundary_ids_vec(n_boundary_ids);
                        mesh.boundary_info->boundary_ids(elem,side, boundary_ids_vec);
                        const unsigned int boundary_id = boundary_ids_vec[0]; 
                        if (boundary_id == boundID)
                        {
                            wave_system.solution->set(dof_indices[0], value);
                            wave_system.solution->set(dof_indices[1], value);
                            wave_system.solution->set(dof_indices[2], value);
                        }
                    }
                }
            }
        }
        std::cout << "done" << std::endl;
        wave_system.solution->close();
        std::cout << "done" << std::endl;
    }

    const std::unique_ptr<libMesh::NumericVector<libMesh::Number> >&
    ElectroSolver::get_fibers()
    {
        ParameterSystem& fiber_system = M_equationSystems.get_system < ParameterSystem > ("fibers");
        return fiber_system.solution;
    }

    const std::unique_ptr<libMesh::NumericVector<libMesh::Number> >&
    ElectroSolver::get_sheets()
    {
        ParameterSystem& sheets_system = M_equationSystems.get_system < ParameterSystem > ("sheets");
        return sheets_system.solution;
    }

    const std::unique_ptr<libMesh::NumericVector<libMesh::Number> >&
    ElectroSolver::get_xfibers()
    {
        ParameterSystem& xfiber_system = M_equationSystems.get_system < ParameterSystem > ("xfibers");
        return xfiber_system.solution;
    }

    void ElectroSolver::solve_reaction_step(double dt, double time, int step, bool useMidpoint, const std::string& mass, libMesh::NumericVector<libMesh::Number>* I4f_ptr)
    {
        if (M_FEFamily == libMesh::MONOMIAL || M_FEFamily == libMesh::L2_LAGRANGE)
        {
            solve_reaction_step_dg(dt, time, step, useMidpoint, mass, I4f_ptr);
        }
        else solve_reaction_step_cg(dt, time, step, useMidpoint, mass, I4f_ptr);

    }

    void ElectroSolver::solve_reaction_step_cg(double dt, double time, int step, bool useMidpoint, const std::string& mass, libMesh::NumericVector<libMesh::Number>* I4f_ptr)
    {
        const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
//    BidomainSystem& bidomain_system = M_equationSystems.get_system
//            < BidomainSystem > ("bidomain"); //Q
        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);
        IonicModelSystem& istim_system = M_equationSystems.get_system < IonicModelSystem > ("istim");
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        IonicModelSystem& iion_system = M_equationSystems.get_system < IonicModelSystem > ("iion");

        system.rhs->zero();
        istim_system.solution->zero();
        istim_system.get_vector("stim_i").zero();
        istim_system.get_vector("stim_e").zero();
        istim_system.get_vector("surf_stim_i").zero();
        istim_system.get_vector("surf_stim_e").zero();
        iion_system.solution->zero();
        iion_system.get_vector("diion").zero();

        //iion_system.old_local_solution->zero();

        double Cm = 1.0;        //M_ionicModelPtr->membraneCapacitance();

        int var_index = 0;

        libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
        const libMesh::MeshBase::const_node_iterator end_node = mesh.local_nodes_end();

        const libMesh::DofMap & dof_map = system.get_dof_map();
        const libMesh::DofMap & dof_map_V = wave_system.get_dof_map();
        const libMesh::DofMap & dof_map_istim = istim_system.get_dof_map();

        std::vector < libMesh::dof_id_type > dof_indices;
        std::vector < libMesh::dof_id_type > dof_indices_V;
        std::vector < libMesh::dof_id_type > dof_indices_Q;
        std::vector < libMesh::dof_id_type > dof_indices_istim;
        std::vector < libMesh::dof_id_type > dof_indices_gating;

        if (M_pacing) M_pacing->update(time);
        if (M_pacing_i) M_pacing_i->update(time);
        if (M_pacing_e) M_pacing_e->update(time);
        if (M_surf_pacing_i) M_surf_pacing_i->update(time);
        if (M_surf_pacing_e) M_surf_pacing_e->update(time);

        int c = 0;
        for (; node != end_node; ++node)
        {
            double istim = 0.0;
            double istim_zero = 0.0;
            double stim_i = 0.0;
            double surf_stim_i = 0.0;
            double stim_e = 0.0;
            double surf_stim_e = 0.0;

            double dIion = 0.0;

            const libMesh::Node * nn = *node;
            // Are we in the bath?
            auto n_var = nn->n_vars(system.number());
            auto n_dofs = nn->n_dofs(system.number());
            // std::cout << "n_dofs per node: " << n_dofs << ", n_var: " << n_var << std::endl;
            libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));
            if (n_var == n_dofs)
            {
                dof_map.dof_indices(nn, dof_indices_Q, 0);
                dof_map_V.dof_indices(nn, dof_indices_V, 0);
                dof_map_istim.dof_indices(nn, dof_indices_istim, 0);

                if (M_pacing_i) stim_i = M_pacing_i->eval(p, time);
                if (M_pacing_e) stim_e = M_pacing_e->eval(p, time);
                if (M_surf_pacing_i) surf_stim_i = M_surf_pacing_i->eval(p, time);
                if (M_surf_pacing_e) surf_stim_e = M_surf_pacing_e->eval(p, time);
                if (M_pacing) istim = M_pacing->eval(p, time);

                double Iion_old = 0.0;
                double v = (*wave_system.old_local_solution)(dof_indices_V[0]); //V^n
                double gating_rhs_ = (*system.old_local_solution)(dof_indices_Q[0]); //Q^n
                Iion_old = (*iion_system.old_local_solution)(dof_indices_istim[0]); // gating
                int key = iion_system.get_vector("ionic_model_map")(dof_indices_istim[0]);
                auto it_ionic_model = M_ionicModelPtrMap.find(key);
                auto it_ionic_model_name = M_ionicModelNameMap.find(key);
//                std::cout << "key: " << key << std::endl;

                IonicModel * ionicModelPtr = nullptr;
                if (it_ionic_model != M_ionicModelPtrMap.end()) ionicModelPtr = it_ionic_model->second.get();
                else
                {
                   throw std::runtime_error("node without ionicModelPtr!!!");
                }
//                if (it_ionic_model_name != M_ionicModelNameMap.end()) ionic_model_system_name = it_ionic_model_name->second;
                std::string ionic_model_system_name = M_ionicModelNameMap[it_ionic_model->first];
//                std::cout << "node ID: " << nn->id() << std::endl;
//                std::cout  << "ionicModelPtr: " << ionicModelPtr << std::end;
//                std::cout  << "ionic_model_system_name: " << ionic_model_system_name << std::end;
               // std::cout << "ptr: " << ionicModelPtr << std::endl;
                if (ionicModelPtr)
                {
                    IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > (ionic_model_system_name);
                    auto& iion_rhs_old = ionic_model_system.get_vector("rhs_old");
                    int num_vars = ionic_model_system.n_vars();
                   // std::cout << "ionic_model_system_name: " << ionicModelPtr->ionicModelName() << ", num_vars: " << num_vars << std::endl;
                    std::vector<double> values(num_vars + 1, 0.0);
                    values[0] = v;
                    std::vector<double> old_values(num_vars + 1, 0.0);
                    std::vector<double> gating_rhs(num_vars + 1, 0.0); // First entry is reserved to Q^n
                    gating_rhs[0] = gating_rhs_;
                    std::vector<double> gating_rhs_old(num_vars + 1, 0.0); // First entry is reserved to Q^n-1
                    const libMesh::DofMap & dof_map_gating = ionic_model_system.get_dof_map();
                    dof_map_gating.dof_indices(nn, dof_indices_gating);

                    for (int nv = 0; nv < num_vars; ++nv)
                    {
                        var_index = dof_indices_gating[nv];
                        values[nv + 1] = (*ionic_model_system.old_local_solution)(var_index);
                        old_values[nv + 1] = values[nv + 1];

                    }

                    if (TimeIntegrator::FirstOrderIMEX == M_timeIntegrator)
                    {
                        ionicModelPtr->updateVariables(values, istim, dt);

                    }
                    else // using SBDF2
                    {
                        if (M_timestep_counter >= 0)
                        {
                            bool overwrite = true;
                            // Recall: gating_rhs[0] = Q^n
                            ionicModelPtr->updateVariables(values, gating_rhs, istim, dt, overwrite);
                        }
                        else
                        {
                            bool overwrite = false;
                            // Recall: gating_rhs[0] = Q^n
                            ionicModelPtr->updateVariables(values, gating_rhs, istim, dt, overwrite);
                            for (int nv = 0; nv < values.size(); ++nv)
                            {
                                var_index = dof_indices_gating[nv];
                                ionic_model_system.rhs->set(var_index, gating_rhs[nv + 1]);

                                old_values[nv + 1] = (*ionic_model_system.older_local_solution)(var_index);
                                double f_nm1 = iion_rhs_old(var_index);
                                double f_n = gating_rhs[nv + 1];
                                // Update using SBDF2
                                // w^n+1 = 4/3 * w^n - 1/3 * w^n-1 + 2/3*dt * (2*f^n - f^n-1)
                                // w^n+1 = ( 4 * w^n - w^n-1 + 2 * dt * (2*f^n - f^n-1) ) / 3
                                values[nv + 1] = (4.0 * values[nv + 1] - old_values[nv + 1] + 2.0 * dt * (2 * f_n - f_nm1)) / 3;
                            }
                        }

                    }
                    double Iion = ionicModelPtr->current_scaling() * ionicModelPtr->evaluateIonicCurrent(values, istim, dt);
                    // Recall: gating_rhs[0] = Q^n
                    // HACK: For now, as I've implemented the second order scheme only for a dew ionic models
                    //       I keep everything as it was before I started the implementation of SBDF2
                    if (ionicModelPtr->isSecondOrderImplemented()) dIion = ionicModelPtr->evaluateIonicCurrentTimeDerivative(values, gating_rhs, dt, M_meshSize);
                    else dIion = ionicModelPtr->evaluateIonicCurrentTimeDerivative(values, old_values, dt, M_meshSize);

                    double Isac = 0.0;
                    if (I4f_ptr)
                    {
                        double I4f = (*I4f_ptr)(dof_indices_V[0]);
                        Isac = ionicModelPtr->evaluateSAC(values[0], I4f);
                    }
                    Iion += Isac;

                    iion_system.solution->set(dof_indices_istim[0], Iion); // contains Istim
                    istim_system.solution->set(dof_indices_istim[0], istim);
                    istim_system.get_vector("stim_i").set(dof_indices_istim[0], stim_i); //Istim^n+1
                    istim_system.get_vector("surf_stim_i").set(dof_indices_istim[0], surf_stim_i); //Istim^n+1
                    istim_system.get_vector("stim_e").set(dof_indices_istim[0], stim_e); //Istim^n+1
                    istim_system.get_vector("surf_stim_e").set(dof_indices_istim[0], surf_stim_e); //Istim^n+1

                    iion_system.get_vector("diion").set(dof_indices_istim[0], dIion);
                    for (int nv = 0; nv < num_vars; ++nv)
                    {
                        var_index = dof_indices_gating[nv];

                        ionic_model_system.solution->set(var_index, values[nv + 1]);
                    }
                }
            }
            c++;

        }
        iion_system.solution->close();

        istim_system.solution->close();
//    std::cout << "Iion sol: " << std::endl;
//    iion_system.solution->print();
//    istim_system.solution->print();
//    std::cout << "Iion sol: " << std::endl;

        istim_system.get_vector("stim_i").close();
        istim_system.get_vector("stim_e").close();
        istim_system.get_vector("surf_stim_i").close();
        istim_system.get_vector("surf_stim_e").close();
        iion_system.get_vector("diion").close();
        iion_system.get_vector("diion_old").close();

        for (auto && name : M_ionic_models_systems_name_vec)
        {
            IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > (name);
            ionic_model_system.solution->close();
            ionic_model_system.get_vector("rhs_old").close();
            ionic_model_system.rhs->close();
            ionic_model_system.update();
        }

        iion_system.update();
        istim_system.update();
    }

    void ElectroSolver::solve_reaction_step_dg(double dt, double time, int step, bool useMidpoint, const std::string& mass, libMesh::NumericVector<libMesh::Number>* I4f_ptr)
    {
        throw std::runtime_error("DG NOT CODED!");
    }

    void ElectroSolver::setup_ODE_systems(GetPot& data, std::string section)
    {
        // ///////////////////////////////////////////////////////////////////////
        // ///////////////////////////////////////////////////////////////////////
        // 2) ODEs
        std::cout << "* ElectroSolver: Creating new System for the ionic model. Data section: " << section << std::endl;
        // Create Ionic Model
        std::string ionic_models_list = data(section + "/ionic_models_list", "NONE");
        BeatIt::readList(ionic_models_list, M_ionic_models_vec);
        int num_ionic_models = M_ionic_models_vec.size();

        std::string ionic_models_IDs_list = data(section + "/ionic_models_IDs_list", "NONE");
        BeatIt::readList(ionic_models_IDs_list, M_ionic_models_IDs_vec);
        if (M_ionic_models_IDs_vec.size() < 1)
        {
            for (auto && m : M_ionic_models_vec)
                M_ionic_models_IDs_vec.push_back(secret_blockID_key);
        }

        for (int k = 0; k < num_ionic_models; ++k)
        {
        	std::string ionic_model_system_name = M_ionic_models_vec[k] + "_" + std::to_string(k);
        	M_ionic_models_systems_name_vec.push_back(ionic_model_system_name);
            //        std::string ionic_model = M_datafile(section + "/ionic_model", "NashPanfilov");
            IonicModelSystem& ionic_model_system = M_equationSystems.add_system < IonicModelSystem > (M_ionic_models_systems_name_vec[k]);
            M_ionicModelExporterNames.insert(M_ionic_models_systems_name_vec[k]);

            std::shared_ptr<IonicModel> ionicModelPtr(BeatIt::IonicModel::IonicModelFactory::Create(M_ionic_models_vec[k]));
            M_ionicModelPtrMap[M_ionic_models_IDs_vec[k]] = ionicModelPtr;
            M_ionicModelNameMap[M_ionic_models_IDs_vec[k]] = M_ionic_models_systems_name_vec[k];
            ionicModelPtr->setup(M_datafile, section);
            int num_vars = ionicModelPtr->numVariables();
            // TO DO: Generalize to other conditions
            // We need to exclude the potential
            // therefore we loop up to num_vars-1
            std::set<short unsigned int> gating_block_id;
            gating_block_id.insert(M_ionic_models_IDs_vec[k]);
            for (int nv = 0; nv < num_vars - 1; ++nv)
            {
                std::string var_name = ionicModelPtr->variableName(nv);
                // For the time being we use P1 for the variables
                //if (M_tissueBlockID >= 0)
                ionic_model_system.add_variable(&var_name[0], M_order, M_FEFamily, &gating_block_id);
                //        else
                //            ionic_model_system.add_variable(&var_name[0], M_order);
            }
            ionic_model_system.add_vector("rhs_old");
            //ionic_model_system.add_vector("rhs_older");
            ionic_model_system.init();
        }
        std::cout << "Ionic Model Name Map: " << std::endl;
        for (auto && m : M_ionicModelNameMap)
            std::cout << m.first << ", " << m.second << std::endl;
        std::cout << "Ionic Model Ptr Map: " << std::endl;
        for (auto && m : M_ionicModelPtrMap)
            std::cout << m.first << ", " << m.second << std::endl;

    }
} /* namespace BeatIt */
