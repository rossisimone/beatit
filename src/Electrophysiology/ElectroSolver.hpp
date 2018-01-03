/*
 * ElectroSolver.hpp
 *
 *  Created on: Oct 19, 2017
 *      Author: srossi
 */

#ifndef SRC_ELECTROPHYSIOLOGY_ELECTROSOLVER_HPP_
#define SRC_ELECTROPHYSIOLOGY_ELECTROSOLVER_HPP_


#include "Util/TimeData.hpp"
// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_linear_solver.h"

#include <memory>
#include "Util/SpiritFunction.hpp"
#include "Util/Enums.hpp"
#include "Util/Factory.hpp"

#include "libmesh/id_types.h"

// Forward Definition
namespace libMesh
{
class Mesh;
class ExodusII_IO;
class VTKIO;
class GMVIO;
class MeshRefinement;
class TimeData;
class ExplicitSystem;
class LinearImplicitSystem;
class PacingProtocol;
class ErrorVector;
class  MeshRefinement;
class BoundaryMesh;
}

namespace BeatIt
{

/*!
 *  forward declarations
 */
class IonicModel;
class PacingProtocol;

enum class Anisotropy { Isotropic,
                        TransverselyIsotropic,
                        Orthotropic,
                        UserDefined };

enum class EquationType { ReactionDiffusion,
                          Wave,
                          ParabolicEllipticBidomain,
                          ParabolicEllipticHyperbolic,
                          ParabolicParabolicHyperbolic   };
enum class ModelType { Monodomain, Bidomain, BidomainWithBath };
enum class TimeIntegrator { FirstOrderIMEX,     // FORWARD-BACKWARD EULER
                            SecondOrderIMEX  }; // SBDF2


class ElectroSolver
{
public:
    typedef FactoryArg<ElectroSolver, std::string, libMesh::EquationSystems >     ElectroFactory;
    typedef libMesh::GMVIO Exporter;
    // Another alternative when not using AMR
    typedef libMesh::ExodusII_IO EXOExporter;

    ElectroSolver( libMesh::EquationSystems& es, std::string model = "bidomain" );
    virtual ~ElectroSolver();


    void setup(GetPot& data, std::string section);
    virtual void setupSystems(GetPot& data, std::string section) = 0;
    void restart( EXOExporter& importer, int step = 0, bool restart = true );
    void readFibers( EXOExporter& importer, int step = 0);

    void init(double time);
    virtual void initSystems(double time) = 0;
    void save(int step);
    void save_exo_timestep(int step, double time);
    void save_ve_timestep(int step, double time);
    void init_exo_output();
    void save_potential(int step, double time = 0.0);
    void save_parameters();
    void save_activation_times(int step = 1);

    virtual void amr( libMesh:: MeshRefinement& mesh_refinement, const std::string& type = "kelly" ) {}
    void reinit_linear_solver();
    //void update_pacing(double time);
    void update_activation_time(double time, double threshold = 0.8);


    //void cut(double time, std::string f);
    virtual void assemble_matrices(double dt = 1.0)  = 0;
    virtual void form_system_matrix(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass") = 0;
    virtual void form_system_rhs(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass") = 0;
    void advance();
    virtual void solve_reaction_step( double dt,
                              double time,
                              int step = 0,
                              bool useMidpoint = true,
                              const std::string& mass = "mass",
                              libMesh::NumericVector<libMesh::Number>* I4f_ptr = nullptr);

    virtual void solve_diffusion_step(double dt, double time,  bool useMidpoint = true, const std::string& mass = "lumped_mass", bool reassemble = true) = 0;
    virtual void generate_fibers(   const GetPot& data,
                            const std::string& section = "rule_based_fibers" ) {}
   double last_activation_time();
   double potential_norm();

   void set_potential_on_boundary(unsigned int boundID, double value = 1.0);

   std::string  get_ionic_model_name() const;


   const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
   get_fibers();
   const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
   get_sheets();
   const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
   get_xfibers();
    //protected:


    /// input file
    GetPot                     M_datafile;
    /// Store pointer to the ionic model
    std::unique_ptr<IonicModel> M_ionicModelPtr;
    /// Equation Systems: One for the potential and one for the other variables
    /*!
     *  Use separate systems to avoid saving in all the variables
     */
    libMesh::EquationSystems&  M_equationSystems;

    std::unique_ptr<Exporter> M_exporter;
    std::set<std::string> M_exporterNames;
    std::unique_ptr<Exporter> M_ionicModelExporter;
    std::set<std::string> M_ionicModelExporterNames;
    std::unique_ptr<EXOExporter> M_parametersExporter;
    std::set<std::string> M_parametersExporterNames;
    std::unique_ptr<EXOExporter> M_EXOExporter;
    std::unique_ptr<EXOExporter> M_potentialEXOExporter;

    std::string  M_outputFolder;
    bool M_assembleMatrix;
    bool M_useAMR;
    std::string  M_systemMass;

    std::unique_ptr<PacingProtocol> M_pacing; // intracellular
    std::unique_ptr<PacingProtocol> M_pacing_i; // intracellular
    std::unique_ptr<PacingProtocol> M_pacing_e; // extracellular
    std::unique_ptr<PacingProtocol> M_surf_pacing_i; // intracellular
    std::unique_ptr<PacingProtocol> M_surf_pacing_e; // extracellular


    std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> > M_linearSolver;
    std::vector<double>  M_intraConductivity;
    std::vector<double>  M_extraConductivity;

    Anisotropy M_anisotropy;
    std::vector<double>  M_conductivity;
    EquationType M_equationType;
    DynamicTimeIntegratorType M_timeIntegratorType;
    TimeIntegrator M_timeIntegrator;
protected:
    long int M_timestep_counter;
public:


    std::string M_model;

    double M_meshSize;
    ModelType M_modelType;
    std::string M_section;

    std::set<unsigned int> M_transmembranePotentialActiveSubdomains;
    int M_constraint_dof_id;
    bool M_ground_ve;


    struct EndocardialVe
    {
        std::unique_ptr<libMesh::BoundaryMesh> M_endocardium;
        std::unique_ptr<libMesh::EquationSystems>  M_boundary_es;
        std::map<libMesh::dof_id_type, libMesh::dof_id_type> M_reverse_node_id_map;
        std::map<libMesh::dof_id_type, libMesh::dof_id_type> M_node_id_map;

        std::unique_ptr<EXOExporter> M_EXOExporter;
    };
    EndocardialVe M_boundary_ve;
    void init_endocardial_ve(std::set<libMesh::boundary_id_type>& IDs, std::set<unsigned short>& subdomainIDs);
};

} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_ELECTROSOLVER_HPP_ */
