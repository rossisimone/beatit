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

#include "libmesh/exodusII_io.h"
//#include "libmesh/exodusII_io_helper.h"
#include "libmesh/gmv_io.h"

#include "libmesh/perf_log.h"

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

namespace BeatIt
{

// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
    typedef libMesh::TransientLinearImplicitSystem ElectroSystem;
    typedef libMesh::TransientExplicitSystem IonicModelSystem;
    typedef libMesh::ExplicitSystem ParameterSystem;

    ElectroSolver::ElectroSolver(libMesh::EquationSystems& es, std::string model) :
            M_equationSystems(es), M_exporter(), M_exporterNames(), M_ionicModelExporter(), M_ionicModelExporterNames(), M_parametersExporter(), M_parametersExporterNames(), M_outputFolder(), M_datafile(), M_pacing_i(), M_pacing_e(), M_linearSolver(), M_anisotropy(
                    Anisotropy::Orthotropic), M_equationType(EquationType::ParabolicEllipticBidomain), M_timeIntegratorType(
                    DynamicTimeIntegratorType::Implicit), M_useAMR(false), M_assembleMatrix(true), M_systemMass("lumped"), M_intraConductivity(), M_extraConductivity(), M_conductivity(), M_meshSize(
                    1.0), M_model(model), M_ground_ve(false), M_timeIntegrator(TimeIntegrator::FirstOrderIMEX), M_timestep_counter(0)
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

        //Read output folder from datafile
        std::string output_folder = M_datafile(M_section + "/output_folder", "Output");
        M_outputFolder = "./" + output_folder + "/";

        std::cout << "* " << M_model << ": creating pacing protocol" << std::endl;

        std::string pacing_type = M_datafile(M_section + "/pacing/type", "NONE");
        if ("NONE" != pacing_type)
        {
            M_pacing.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_type));
            M_pacing->setup(M_datafile, M_section + "/pacing");
        }

        std::string pacing_i_type = M_datafile(M_section + "/pacing_i/type", "NONE");
        if ("NONE" != pacing_i_type)
        {
            M_pacing_i.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_i_type));
            M_pacing_i->setup(M_datafile, M_section + "/pacing_i");
        }

        std::string pacing_e_type = M_datafile(M_section + "/pacing_e/type", "NONE");
        if ("NONE" != pacing_e_type)
        {
            M_pacing_e.reset(BeatIt::PacingProtocol::PacingProtocolFactory::Create(pacing_e_type));
            M_pacing_e->setup(M_datafile, M_section + "/pacing_e");
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
                mkdir(M_outputFolder.c_str(), 0777);
            }
        }

        // call setup system of the specific class
        setupSystems(M_datafile, M_section);
        // Add ionic current to this system
        IonicModelSystem& Iion_system = M_equationSystems.add_system < IonicModelSystem > ("iion");
        Iion_system.add_variable("iion", libMesh::FIRST);
        Iion_system.add_vector("diion");
        Iion_system.add_vector("diion_old");
        Iion_system.add_vector("total_current");
        Iion_system.init();
        M_ionicModelExporterNames.insert("iion");

        // Add the applied current to this system
        IonicModelSystem& istim_system = M_equationSystems.add_system < IonicModelSystem > ("istim");
        istim_system.add_variable("istim", libMesh::FIRST);
        istim_system.init();
        istim_system.add_vector("stim_i");
        istim_system.add_vector("stim_e");
        istim_system.add_vector("surf_stim_i");
        istim_system.add_vector("surf_stim_e");

        M_ionicModelExporterNames.insert("istim");

        std::cout << "* ElectroSolver: Creating parameters spaces " << std::endl;
        ParameterSystem& activation_times_system = M_equationSystems.add_system < ParameterSystem > ("activation_times");
        activation_times_system.add_variable("activation_times", libMesh::FIRST);
        activation_times_system.init();

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

        ///////////
        // Time integrator
        int time_integrator_order = data(M_section + "/time_integrator_order", 1);
        if (2 == time_integrator_order && M_ionicModelPtr->isSecondOrderImplemented()) M_timeIntegrator = TimeIntegrator::SecondOrderIMEX;

        // ///////////////////////////////////////////////////////////////////////
        // ///////////////////////////////////////////////////////////////////////
        // Setup Exporters

        std::cout << "* " << M_model << ": Initializing exporters " << std::endl;
        M_exporter.reset(new Exporter(M_equationSystems.get_mesh()));
        M_ionicModelExporter.reset(new Exporter(M_equationSystems.get_mesh()));
        M_parametersExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
        M_EXOExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));
        M_potentialEXOExporter.reset(new EXOExporter(M_equationSystems.get_mesh()));

    }

    void ElectroSolver::init(double time)
    {
        std::cout << "* ElectroSolver: Init: " << std::endl;

        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        // Setting initial conditions
        ElectroSystem& electro_system = M_equationSystems.get_system < ElectroSystem > (M_model);
        IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > ("ionic_model");
        int num_vars = ionic_model_system.n_vars();
        std::vector<double> init_values(num_vars + 1, 0.0);
        M_ionicModelPtr->initialize(init_values);

        auto first_local_index = wave_system.solution->first_local_index();
        auto last_local_index = wave_system.solution->last_local_index();
        for (auto i = first_local_index; i < last_local_index; ++i)
        {
            wave_system.solution->set(i, init_values[0]);
        }
        first_local_index = ionic_model_system.solution->first_local_index();
        last_local_index = ionic_model_system.solution->last_local_index();
        std::cout << "* " << M_model << ": Setting initial values ... " << std::flush;
        for (auto i = first_local_index; i < last_local_index;)
        {
            for (int nv = 0; nv < num_vars; ++nv)
            {
                ionic_model_system.solution->set(i, init_values[nv + 1]);
                ++i;
            }
        }
        electro_system.solution->close();
        electro_system.old_local_solution->close();
        electro_system.older_local_solution->close();
        ionic_model_system.solution->close();
        ionic_model_system.old_local_solution->close();
        ionic_model_system.older_local_solution->close();

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

        first_local_index = procID_system.solution->first_local_index();
        last_local_index = procID_system.solution->last_local_index();

        double myrank = static_cast<double>(M_equationSystems.comm().rank());
        for (int el = first_local_index; el < last_local_index; ++el)
        {
            procID_system.solution->set(el, myrank);
        }
        procID_system.solution->close();

        // Call specific system initializations
        initSystems(time);
//    M_linearSolver =  libMesh::LinearSolver<libMesh::Number>::build( M_equationSystems.comm() );
        typedef libMesh::PetscLinearSolver<libMesh::Number> PetscSolver;
        M_linearSolver.reset(new PetscSolver(M_equationSystems.comm()));
        std::string solver_type = M_datafile(M_section + "/linear_solver/type", "gmres");
        std::cout << "* ElectroSolver: using " << solver_type << std::endl;
        std::string prec_type = M_datafile(M_section + "/linear_solver/preconditioner", "amg");
        M_linearSolver->set_solver_type(solver_map.find(solver_type)->second);
        M_linearSolver->set_preconditioner_type(prec_map.find(prec_type)->second);
        M_linearSolver->init();

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

    void ElectroSolver::readFibers(EXOExporter& importer, int step)
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
        std::cout << "* " << M_model << ": EXODUSII::Exporting " << M_model << ".exo at time 0.0" << " in: " << M_outputFolder << " ... "
                << std::flush;

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

        std::cout << "* " << M_model << ": EXODUSII::Exporting " << M_model << "_endo.exo at time " << time << " in: " << M_outputFolder << " ... "
                << std::flush;
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
        std::cout << "* " << M_model << ": EXODUSII::Exporting " << M_model << ".exo at time " << time << " in: " << M_outputFolder << " ... "
                << std::flush;
        M_EXOExporter->write_timestep(M_outputFolder + M_model + ".exo", M_equationSystems, step, time);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save_parameters()
    {
        std::cout << "* " << M_model << ": EXODUSII::Exporting parameters in: " << M_outputFolder << " ... " << std::flush;
        EXOExporter exo(M_equationSystems.get_mesh());
        exo.write("parameters.exo");
        exo.write_element_data(M_equationSystems);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save_potential(int step, double time)
    {
        std::cout << "* " << M_model << ": EXODUSII::Exporting potential.exo at time " << time << " in: " << M_outputFolder << " ... " << std::flush;

        M_potentialEXOExporter->write_timestep(M_outputFolder + "potential.exo", M_equationSystems, step, time);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save_activation_times(int step)
    {
        std::cout << "* " << M_model << ": GMVIO::Exporting activaton times: " << M_outputFolder << " ... " << std::flush;
        Exporter gmv(M_equationSystems.get_mesh());
        std::set < std::string > output;
        output.insert("activation_times");
        gmv.write_equation_systems(M_outputFolder + "activation_times.gmv." + std::to_string(step), M_equationSystems, &output);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::save(int step)
    {
        std::cout << "* " << M_model << ": GMVIO::Exporting " << step << " in: " << M_outputFolder << " ... " << std::flush;
        M_exporter->write_equation_systems(M_outputFolder + M_model + ".gmv." + std::to_string(step), M_equationSystems, &M_exporterNames);
        M_ionicModelExporter->write_equation_systems(M_outputFolder + "ionic_model.gmv." + std::to_string(step), M_equationSystems,
                &M_ionicModelExporterNames);
        std::cout << "done " << std::endl;
    }

    void ElectroSolver::reinit_linear_solver()
    {
        M_linearSolver->clear();
        M_linearSolver->set_solver_type(libMesh::CG);
        M_linearSolver->set_preconditioner_type(libMesh::AMG_PRECOND);
        M_linearSolver->init();
    }

    void ElectroSolver::update_activation_time(double time, double threshold)
    {
        ParameterSystem& activation_times_system = M_equationSystems.get_system < ParameterSystem > ("activation_times");
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");

        const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();

        libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
        const libMesh::MeshBase::const_node_iterator end_node = mesh.local_nodes_end();

        const libMesh::DofMap & dof_map = wave_system.get_dof_map();
        const libMesh::DofMap & dof_map_at = activation_times_system.get_dof_map();

        std::vector < libMesh::dof_id_type > dof_indices_V;
        std::vector < libMesh::dof_id_type > dof_indices_at;

        double v = 0.0;
        double at = 0.0;

        for (; node != end_node; ++node)
        {
            const libMesh::Node * nn = *node;
            // Are we in the bath?

            dof_map.dof_indices(nn, dof_indices_V, 0);
            dof_map_at.dof_indices(nn, dof_indices_at, 0);

            if (dof_indices_V.size() > 0)
            {
                v = (*wave_system.solution)(dof_indices_V[0]);
                at = (*activation_times_system.solution)(dof_indices_at[0]);
                if (v > threshold && at < 0.0)
                {
                    activation_times_system.solution->set(dof_indices_at[0], time);
                }

            }
        }

        activation_times_system.solution->close();
    }

    void ElectroSolver::advance()
    {
        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);
        IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > ("ionic_model");

        *system.older_local_solution = *system.old_local_solution;
        *system.old_local_solution = *system.solution;
        *ionic_model_system.older_local_solution = *ionic_model_system.old_local_solution;
        *ionic_model_system.old_local_solution = *ionic_model_system.solution;
        ionic_model_system.get_vector("rhs_old") = *ionic_model_system.rhs;
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        *wave_system.older_local_solution = *wave_system.old_local_solution;
        *wave_system.old_local_solution = *wave_system.solution;

        IonicModelSystem& iion_system = M_equationSystems.get_system < IonicModelSystem > ("iion");
        *iion_system.older_local_solution = *iion_system.old_local_solution;
        *iion_system.old_local_solution = *iion_system.solution;
        iion_system.get_vector("diion_old") = iion_system.get_vector("diion");

    }

    std::string ElectroSolver::get_ionic_model_name() const
    {
        return M_ionicModelPtr->ionicModelName();
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

    void ElectroSolver::set_potential_on_boundary(unsigned int boundID, double value)
    {
        std::cout << "* " << M_model << ": WARNING:  set_potential_on_boundary works only for TET4" << std::endl;
        libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");

        const libMesh::DofMap & dof_map = wave_system.get_dof_map();
        std::vector < libMesh::dof_id_type > dof_indices;

        libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
        const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

        for (; el != end_el; ++el)
        {
            const libMesh::Elem * elem = *el;
            dof_map.dof_indices(elem, dof_indices);
            for (unsigned int side = 0; side < elem->n_sides(); side++)
            {
                if (elem->neighbor(side) == libmesh_nullptr)
                {
                    const unsigned int boundary_id = mesh.boundary_info->boundary_id(elem, side);
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

    const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
    ElectroSolver::get_fibers()
    {
        ParameterSystem& fiber_system = M_equationSystems.get_system < ParameterSystem > ("fibers");
        return fiber_system.solution;
    }

    const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
    ElectroSolver::get_sheets()
    {
        ParameterSystem& sheets_system = M_equationSystems.get_system < ParameterSystem > ("sheets");
        return sheets_system.solution;
    }

    const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
    ElectroSolver::get_xfibers()
    {
        ParameterSystem& xfiber_system = M_equationSystems.get_system < ParameterSystem > ("xfibers");
        return xfiber_system.solution;
    }

    void ElectroSolver::solve_reaction_step(double dt, double time, int step, bool useMidpoint, const std::string& mass,
            libMesh::NumericVector<libMesh::Number>* I4f_ptr)
    {
        const libMesh::MeshBase & mesh = M_equationSystems.get_mesh();
//    BidomainSystem& bidomain_system = M_equationSystems.get_system
//            < BidomainSystem > ("bidomain"); //Q
        ElectroSystem& system = M_equationSystems.get_system < ElectroSystem > (M_model);
        IonicModelSystem& ionic_model_system = M_equationSystems.get_system < IonicModelSystem > ("ionic_model");
        IonicModelSystem& istim_system = M_equationSystems.get_system < IonicModelSystem > ("istim");
        // WAVE
        ElectroSystem& wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
        IonicModelSystem& iion_system = M_equationSystems.get_system < IonicModelSystem > ("iion");

        system.rhs->zero();
        iion_system.solution->zero();
        istim_system.solution->zero();
        auto& iion_rhs_old = ionic_model_system.get_vector("rhs_old");
        istim_system.get_vector("stim_i").zero();
        istim_system.get_vector("stim_e").zero();
        istim_system.get_vector("surf_stim_i").zero();
        istim_system.get_vector("surf_stim_e").zero();
        //iion_system.old_local_solution->zero();
        iion_system.get_vector("diion").zero();

        double Cm = M_ionicModelPtr->membraneCapacitance();

        int num_vars = ionic_model_system.n_vars();
        std::vector<double> values(num_vars + 1, 0.0);
        std::vector<double> old_values(num_vars + 1, 0.0);

        std::vector<double> gating_rhs(num_vars + 1, 0.0); // First entry is reserved to Q^n
        std::vector<double> gating_rhs_old(num_vars + 1, 0.0); // First entry is reserved to Q^n-1

        int var_index = 0;

        libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
        const libMesh::MeshBase::const_node_iterator end_node = mesh.local_nodes_end();

        const libMesh::DofMap & dof_map = system.get_dof_map();
        const libMesh::DofMap & dof_map_V = wave_system.get_dof_map();
        const libMesh::DofMap & dof_map_gating = ionic_model_system.get_dof_map();
        const libMesh::DofMap & dof_map_istim = istim_system.get_dof_map();

        std::vector < libMesh::dof_id_type > dof_indices;
        std::vector < libMesh::dof_id_type > dof_indices_V;
        std::vector < libMesh::dof_id_type > dof_indices_Q;
        std::vector < libMesh::dof_id_type > dof_indices_gating;
        std::vector < libMesh::dof_id_type > dof_indices_istim;

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
            libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));
            if (n_var == n_dofs)
            {
                dof_map.dof_indices(nn, dof_indices_Q, 0);
                dof_map_V.dof_indices(nn, dof_indices_V, 0);
                dof_map_istim.dof_indices(nn, dof_indices_istim, 0);
                dof_map_gating.dof_indices(nn, dof_indices_gating);

                if (M_pacing_i) stim_i = M_pacing_i->eval(p, time);
                if (M_pacing_e) stim_e = M_pacing_e->eval(p, time);
                if (M_surf_pacing_i) surf_stim_i = M_surf_pacing_i->eval(p, time);
                if (M_surf_pacing_e) surf_stim_e = M_surf_pacing_e->eval(p, time);
                if (M_pacing) istim = M_pacing->eval(p, time);

                double Iion_old = 0.0;
                values[0] = (*wave_system.old_local_solution)(dof_indices_V[0]); //V^n
                gating_rhs[0] = (*system.old_local_solution)(dof_indices_Q[0]); //Q^n
                Iion_old = (*iion_system.old_local_solution)(dof_indices_istim[0]); // gating

                for (int nv = 0; nv < num_vars; ++nv)
                {
                    var_index = dof_indices_gating[nv];
                    values[nv + 1] = (*ionic_model_system.old_local_solution)(var_index);
                    old_values[nv + 1] = values[nv + 1];

                }

                if (TimeIntegrator::FirstOrderIMEX == M_timeIntegrator)
                {
                    M_ionicModelPtr->updateVariables(values, istim_zero, dt);

                }
                else // using SBDF2
                {
                    if (M_timestep_counter == 0)
                    {
                        bool overwrite = true;
                        // Recall: gating_rhs[0] = Q^n
                        M_ionicModelPtr->updateVariables(values, gating_rhs, istim_zero, dt, overwrite);

                    }
                    else
                    {
                        bool overwrite = false;
                        // Recall: gating_rhs[0] = Q^n
                        M_ionicModelPtr->updateVariables(values, gating_rhs, istim_zero, dt, overwrite);
                        for (int nv = 0; nv < num_vars; ++nv)
                        {
                            var_index = dof_indices_gating[nv];
                            ionic_model_system.rhs->set(var_index, gating_rhs[nv + 1]);

                            old_values[nv + 1] = (*ionic_model_system.older_local_solution)(var_index);
                            double f_nm1 = iion_rhs_old(var_index);
                            double f_n = gating_rhs[nv + 1];
                            // Update using SBDF2
                            // w^n+1 = 4/3 * w^n - 1/3 * w^n-1 + 2/3*dt * (2*f^n - f^n-1)
                            // w^n+1 = ( 4 * w^n - w^n-1 + 2 * dt * (2*f^n - f^n-1) ) / 3
                            values[nv + 1] += (4.0 * values[nv + 1] - old_values[nv + 1] + 2.0 * dt * (2 * f_n - f_nm1)) / 3;
                        }
                    }

                }
                double Iion = M_ionicModelPtr->evaluateIonicCurrent(values, istim_zero, dt);
                // Recall: gating_rhs[0] = Q^n
                // HACK: For now, as I've implemented the second order scheme only for a dew ionic models
                //       I keep everything as it was before I started the implementation of SBDF2
                if (M_ionicModelPtr->isSecondOrderImplemented()) dIion = M_ionicModelPtr->evaluatedIonicCurrent(values, gating_rhs, dt, M_meshSize);
                else dIion = M_ionicModelPtr->evaluatedIonicCurrent(values, old_values, dt, M_meshSize);

                double Isac = 0.0;
                if (I4f_ptr)
                {
                    double I4f = (*I4f_ptr)(dof_indices_V[0]);
                    Isac = M_ionicModelPtr->evaluateSAC(values[0], I4f);
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
        ionic_model_system.solution->close();
        iion_system.get_vector("diion").close();
        iion_system.update();
        ionic_model_system.update();
        istim_system.update();
    }

} /* namespace BeatIt */
