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

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/elem.h"

#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/quadrature_gauss.h"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"

#include <sys/stat.h>
#include "libmesh/getpot.h"

#include "libmesh/vtk_io.h"
#include <iomanip>      // std::setw

// Domain numbering
enum class Subdomain : unsigned int
{
    TISSUE = 1, BATH = 2
};
// Time integrators
enum class TimeIntegrator : unsigned int
{
    SBDF1 = 1, SBDF2 = 2, SBDF3 = 3, EXPLICIT_EXTRACELLULAR = 4, EXPLICIT_INTRACELLULAR = 5, SEMI_IMPLICIT = 6, SEMI_IMPLICIT_HALF_STEP = 7
};
enum class EquationType : unsigned int
{
    PARABOLIC = 1, ELLIPTIC = 2
};

// Store here details about time stepping
struct TimeData
{
    // constructor
    TimeData(GetPot &data);
    // show parameters
    void print();
    // parameters
    double time;      // Time
    double end_time;      // Final Time
    double dt;            // Timestep
    int timestep;         // timestep coutner
    int export_timesteps; // interval at which export the solution
};


struct IonicModel
{
    void setup(GetPot& data)
    {
        v0 = data("v0", 0.0);
        v1 = data("v1", 0.1);
        v2 = data("v2", 1.0);
        k = data("k", 8.0);
    }

    double iion(double vn)
    {
        return k * (vn-v0) * (vn-v1) * (vn-v2);
    }

    void print()
    {
        std::cout << "**************************************** " << std::endl;
        std::cout << "Ionic model parameters:" << std::endl;
        std::cout << "k: " << k << std::endl;
        std::cout << "v0: " << v0 << std::endl;
        std::cout << "v1: " << v1 << std::endl;
        std::cout << "v2: " << v2 << std::endl;
        std::cout << "**************************************** " << std::endl;
    }

    double k;
    double v0, v1, v2;
};

struct Pacing
{
    Pacing(GetPot& data)
    {
        duration = data("stimulus_duration", 2.0);
        amplitude = data("stimulus_amplitude", 50.0);
        start_time = data("stimulus_start_time", 0.0);
        double maxx = data("stimulus_maxx", 1.0);
        double maxy = data("stimulus_maxy", 1.0);
        double maxz = data("stimulus_maxz", 1.0);
        double minx = data("stimulus_minx", 0.0);
        double miny = data("stimulus_miny", 0.0);
        double minz = data("stimulus_minz", 0.0);
        box = libMesh::BoundingBox(libMesh::Point(minx, miny, minz), libMesh::Point(maxx, maxy, maxz));
    }

    double istim(const libMesh::Point& p, double time)
    {
        double I = 0.0;
        if(time > start_time && time <= start_time+ duration)
        {
            if(box.contains_point(p))
            {
                I = amplitude;
            }
        }
        return I;
    }

    void print()
    {
        std::cout << "**************************************** " << std::endl;
        std::cout << "Pacing parameters:" << std::endl;
        std::cout << "start_time: " << start_time << std::endl;
        std::cout << "duration: " << duration << std::endl;
        std::cout << "amplitude: " << amplitude << std::endl;
        std::cout << "min: " << std::endl;
        box.min().print();
        std::cout << "\nmax: " << std::endl;
        box.max().print();
        std::cout << "\n**************************************** " << std::endl;
        std::cout << std::endl;
    }

    double duration;
    double amplitude;
    double start_time;
    libMesh::BoundingBox box;
};


// Read mesh as specified in the input file
void read_mesh(const GetPot &data, libMesh::ParallelMesh &mesh);
// read BC sidesets from string: e.g. bc = "1 2 3 5", or bc = "1, 55, 2, 33"
void read_bc_list(std::string &bc, std::set<int> &bc_sidesets);
// Assemble matrices
void assemble_matrices(libMesh::EquationSystems &es, const TimeData &timedata, TimeIntegrator time_integrator, libMesh::Order p_order, const GetPot &data);
// Solve ionic model / evaluate ionic currents
void solve_ionic_model_evaluate_ionic_currents(libMesh::EquationSystems &es,
        IonicModel& ionic_model,
        Pacing& pacing,
        TimeData& datatime,
        TimeIntegrator& time_integrator);
// void assemble RHS
void assemble_rhs(  libMesh::EquationSystems &es,
                    const TimeData &timedata,
                    TimeIntegrator time_integrator,
                    libMesh::Order p_order,
                    const GetPot &data,
                    EquationType type = EquationType::PARABOLIC );


int main(int argc, char **argv)
{
    // Read input file
    GetPot cl(argc, argv);
    std::string datafile_name = cl.follow("data.beat", 2, "-i", "--input");
    GetPot data(datafile_name);
    // Use libMesh stuff without writing libMesh:: everytime
    using namespace libMesh;
    // Initialize libMesh
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // Create folde in which we save the output of the solution
    std::string output_folder = data("output_folder", "Output");
    struct stat out_dir;
    if (stat(&output_folder[0], &out_dir) != 0)
    {
        if (init.comm().rank() == 0)
        {
            mkdir(output_folder.c_str(), 0777);
        }
    }

    // Create empty mesh
    ParallelMesh mesh(init.comm());
    read_mesh(data, mesh);
    // output the details about the mesh
    mesh.print_info();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Time stuff
    TimeData datatime(data);
    datatime.print();

    // Pick time integrator
    // This will define if we use 1 or 2 systems
    std::map<std::string, TimeIntegrator> time_integrator_map =
    {
    { "SBDF1", TimeIntegrator::SBDF1 },
    { "SBDF2", TimeIntegrator::SBDF2 },
    { "SBDF3", TimeIntegrator::SBDF3 },
    { "EXPLICIT_EXTRACELLULAR", TimeIntegrator::EXPLICIT_EXTRACELLULAR },
    { "EXPLICIT_INTRACELLULAR", TimeIntegrator::EXPLICIT_INTRACELLULAR },
    { "SEMI_IMPLICIT", TimeIntegrator::SEMI_IMPLICIT },
    { "SEMI_IMPLICIT_HALF_STEP", TimeIntegrator::SEMI_IMPLICIT_HALF_STEP } };

    std::string integrator = data("integrator", "SBDF1");
    auto it = time_integrator_map.find(integrator);
    TimeIntegrator time_integrator = it->second;
    bool using_implicit_time_integrator = false;
    if (time_integrator == TimeIntegrator::SBDF1 || time_integrator == TimeIntegrator::SBDF2 || time_integrator == TimeIntegrator::SBDF3)
    {
        using_implicit_time_integrator = true;
    }

    // Create libMesh Equations systems
    // This will hold the mesh and create the corresponding
    // finite element spaces
    libMesh::EquationSystems es(mesh);

    // Define tissue active subdomain for subdomain restricted variables
    std::set < libMesh::subdomain_id_type > tissue_subdomains;
    tissue_subdomains.insert(static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE));
    // Define finite element ORDER
    int p_order = data("p_order", 1);
    std::cout << "P order: " << p_order << std::endl;
    Order order = FIRST;
    if (p_order == 2)
        order = SECOND;


    // Bidomain System
    switch (time_integrator)
    {
        case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
        case TimeIntegrator::EXPLICIT_INTRACELLULAR:
        case TimeIntegrator::SEMI_IMPLICIT:
        case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
        {
            libMesh::LinearImplicitSystem &elliptic_system = es.add_system < libMesh::LinearImplicitSystem > ("elliptic");
            libMesh::TransientLinearImplicitSystem & parabolic_system = es.add_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            std::cout << "Using element of order p = " << order << std::endl;
            parabolic_system.add_variable("V", order, LAGRANGE, &tissue_subdomains);
            elliptic_system.add_variable("Ve", order, LAGRANGE);
            parabolic_system.add_matrix("Ki");
            parabolic_system.add_matrix("Ke");
            parabolic_system.add_matrix("Kg");
            parabolic_system.add_matrix("M");
            parabolic_system.add_vector("ML");
            parabolic_system.add_vector("In");
            parabolic_system.add_vector("aux1"); // auxilliary vector for assembling the RHS
            break;
        }
        case TimeIntegrator::SBDF1:
        case TimeIntegrator::SBDF2:
        case TimeIntegrator::SBDF3:
        default:
        {
            libMesh::TransientLinearImplicitSystem &system = es.add_system < libMesh::TransientLinearImplicitSystem > ("bidomain");
            std::cout << "Using element of order p = " << order << std::endl;
            system.add_variable("V", order, LAGRANGE, &tissue_subdomains);
            system.add_variable("Ve", order, LAGRANGE);
            system.add_matrix("K");
            system.add_matrix("M");
            system.add_vector("ML");
            system.add_vector("In");
            system.add_vector("Inm1");
            system.add_vector("Inm2");
            system.add_vector("Vnm2");
            system.add_vector("aux1"); // auxilliary vector for assembling the RHS
            break;
        }
    }

    es.init();
    es.print_info();

    // setup pacing protocol
    Pacing pacing(data);
    pacing.print();
    IonicModel ionic_model;
    ionic_model.setup(data);
    ionic_model.print();


    int save_iter = 0;
    {
        std::cout << "Time: " << datatime.time << ", " << std::flush;
        std::cout << "Exporting  ... " << std::flush;
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << save_iter;
        std::string step_str = ss.str();
        libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/bath_" + step_str + ".pvtu", es);
        std::cout << "done " << std::endl;
    }


    // Start loop in time
    for (; datatime.time < datatime.end_time; )
    {
        // advance time
        datatime.time += datatime.dt;
        datatime.timestep++;
        //std::cout << "Time: " << datatime.time << std::endl;
        // advance vectors and solve
        // Bidomain System
        switch (time_integrator)
        {
            case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
            case TimeIntegrator::EXPLICIT_INTRACELLULAR:
            case TimeIntegrator::SEMI_IMPLICIT:
            case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
            {
                if(datatime.timestep == 1)
                {
                    assemble_matrices(es, datatime, time_integrator, order, data);
                }
                libMesh::LinearImplicitSystem &elliptic_system = es.get_system < libMesh::LinearImplicitSystem > ("elliptic");
                libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");

                *parabolic_system.old_local_solution = *parabolic_system.solution;
                solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                // Solve Parabolic
                assemble_rhs(es, datatime, time_integrator, order, data, EquationType::PARABOLIC);
                if(time_integrator == TimeIntegrator::SEMI_IMPLICIT ||
                   time_integrator == TimeIntegrator::SEMI_IMPLICIT_HALF_STEP)
                {
                    parabolic_system.solve();
                }
                else
                {
                    *parabolic_system.solution += *parabolic_system.rhs;
                }
                parabolic_system.update();
                // Solve Elliptic
                assemble_rhs(es, datatime, time_integrator, order, data, EquationType::ELLIPTIC);
                elliptic_system.solve();
                elliptic_system.update();

                if(time_integrator == TimeIntegrator::SEMI_IMPLICIT_HALF_STEP)
                {
                     // Solve Parabolic
                     assemble_rhs(es, datatime, time_integrator, order, data, EquationType::PARABOLIC);
                     parabolic_system.solve();
                     parabolic_system.update();
                }
                break;
            }
            case TimeIntegrator::SBDF3:
            {
                if(datatime.timestep == 1)
                {
                    assemble_matrices(es, datatime, TimeIntegrator::SBDF1, order, data);
                }
                if(datatime.timestep == 2)
                {
                    assemble_matrices(es, datatime, TimeIntegrator::SBDF2, order, data);
                }
                if(datatime.timestep == 3)
                {
                    assemble_matrices(es, datatime, time_integrator, order, data);
                }

                libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("bidomain");

                system.get_vector("Inm2") = system.get_vector("Inm1");
                system.get_vector("Inm1") = system.get_vector("In");
                system.get_vector("Vnm2") = *system.older_local_solution;
                *system.older_local_solution = *system.old_local_solution;
                *system.old_local_solution = *system.solution;
                 solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                if(datatime.timestep == 1) assemble_rhs(es, datatime, TimeIntegrator::SBDF1, order, data);
                else if(datatime.timestep == 2) assemble_rhs(es, datatime, TimeIntegrator::SBDF2, order, data);
                else assemble_rhs(es, datatime, time_integrator, order, data);
                system.solve();
                system.update();
                break;
            }
            case TimeIntegrator::SBDF2:
            {
                if(datatime.timestep == 1)
                {
                    assemble_matrices(es, datatime, TimeIntegrator::SBDF1, order, data);
                }
                if(datatime.timestep == 2)
                {
                    assemble_matrices(es, datatime, time_integrator, order, data);
                }

                libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("bidomain");

                system.get_vector("Inm2") = system.get_vector("Inm1");
                system.get_vector("Inm1") = system.get_vector("In");
                system.get_vector("Vnm2") = *system.older_local_solution;
                *system.older_local_solution = *system.old_local_solution;
                *system.old_local_solution = *system.solution;
                 solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                if(datatime.timestep == 1) assemble_rhs(es, datatime, TimeIntegrator::SBDF1, order, data);
                else assemble_rhs(es, datatime, time_integrator, order, data);
                //system.rhs->print();
                //system.get_vector("ML").print();
                //system.get_matrix("M").print();
                //system.matrix->print();
                system.solve();
                //system.solution->print();
                system.update();
                break;
            }
            case TimeIntegrator::SBDF1:
            default:
            {
                if(datatime.timestep == 1) assemble_matrices(es, datatime, time_integrator, order, data);

                libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("bidomain");
                system.get_vector("Inm2") = system.get_vector("Inm1");
                system.get_vector("Inm1") = system.get_vector("In");
                system.get_vector("Vnm2") = *system.older_local_solution;
                *system.older_local_solution = *system.old_local_solution;
                *system.old_local_solution = *system.solution;
                solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                assemble_rhs(es, datatime, time_integrator, order, data);
                system.solve();
                system.update();
                break;
            }
        }

        //Export the solution if at the right timestep
        if (0 == datatime.timestep % datatime.export_timesteps)
        {
            std::cout << "Time: " << datatime.time << ", " << std::flush;
            std::cout << "Exporting  ... " << std::flush;
            // update output file counter
            save_iter++;//
            std::ostringstream ss;
            ss << std::setw(4) << std::setfill('0') << save_iter;
            std::string step_str = ss.str();
            libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/bath_" + step_str + ".pvtu", es);
            std::cout << "done " << std::endl;
        }
    }

    return 0;
}

void read_mesh(const GetPot &data, libMesh::ParallelMesh &mesh)
{
    using namespace libMesh;
    // read mesh from file?
    std::string mesh_filename = data("mesh_filename", "NONE");
    if (mesh_filename.compare("NONE") != 0)
    {
        // READ MESH
        std::cout << "I will read the mesh here: " << mesh_filename << std::endl;
    }
    else
    {
        // Create Mesh
        // number of elements in the x,y and z direction
        int nelx = data("nelx", 10);
        int nely = data("nely", 10);
        // if nelz = 0 we run in 2D
        int nelz = data("nelz", 0);

        // the cube dimensions are defined as [minx, maxx] x [miny, maxy] x [minz, maxz]
        double maxx = data("maxx", 1.0);
        double maxy = data("maxy", 1.0);
        double maxz = data("maxz", 1.0);

        double minx = data("minx", 0.0);
        double miny = data("miny", 0.0);
        double minz = data("minz", 0.0);

        // Get mesh parameters
        std::string eltype = data("eltype", "simplex");
        std::map < std::string, ElemType > elem_type_map =
        {
        { "TRI3", TRI3 },
        { "TRI6", TRI6 },
        { "QUAD4", QUAD4 },
        { "QUAD9", QUAD9 },
        { "TET4", TET4 },
        { "HEX8", HEX8 } };
        auto elType = TRI3;
        auto elem_type_it = elem_type_map.find(eltype);
        if (elem_type_it != elem_type_map.end())
            elType = elem_type_it->second;

        std::cout << "Creating the cube [" << minx << ", " << maxx << "] x [" << miny << ", " << maxy << "] x [" << minx << ", " << maxx << "] " << std::endl;
        std::cout << "Using " << nelx << " x " << nely << " x " << nelz << " elements." << std::endl;
        std::cout << "Element type " << elem_type_it->first << std::endl;

        // Create mesh
        MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, minx, maxx, miny, maxy, minz, maxz, elType);

        // setup subdomains
        // tissue domain box
        double tissue_maxx = data("tissue_maxx", 1.0);
        double tissue_maxy = data("tissue_maxy", 1.0);
        double tissue_maxz = data("tissue_maxz", 1.0);
        double tissue_minx = data("tissue_minx", 0.0);
        double tissue_miny = data("tissue_miny", 0.0);
        double tissue_minz = data("tissue_minz", 0.0);
        libMesh::BoundingBox box(libMesh::Point(tissue_minx, tissue_miny, tissue_minz), libMesh::Point(tissue_maxx, tissue_maxy, tissue_maxz));

        // right_interface
        for (auto el : mesh.element_ptr_range())
        {
            auto c = el->centroid();
            // Are we in the tissue region?
            if (box.contains_point(c))
            {
                el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
            }
            // we are not
            else
            {
                el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::BATH);
            }
        }
    }
}

// Initialize TimeData from input
TimeData::TimeData(GetPot &data) :
        time(0.0),
        end_time(data("end_time", 1.0)),
        dt(data("dt", 0.125)),
        timestep(0),
        export_timesteps(data("export_timesteps", 1))
{
}
// Initialize TimeData from input
void TimeData::print()
{
    std::cout << "**************************************** " << std::endl;
    std::cout << "TimeData Parameters: " << std::endl;
    std::cout << "End time: " << end_time << std::endl;
    std::cout << "Dt: " << dt << std::endl;
    std::cout << "Current Timestep: " << timestep << std::endl;
    std::cout << "Export Timesteps: " << export_timesteps << std::endl;
    std::cout << "**************************************** " << std::endl;
}

void assemble_matrices(libMesh::EquationSystems &es, const TimeData &timedata, TimeIntegrator time_integrator, libMesh::Order p_order, const GetPot &data)
{
    using namespace libMesh;
    // Create vector of BC sidesets ids
    std::set<int> bc_sidesets;
    std::string bcs = data("bcs", "");
    read_bc_list(bcs, bc_sidesets);


    const MeshBase &mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // volume element
    FEType fe_type(p_order);
    std::unique_ptr < FEBase > fe(FEBase::build(dim, fe_type));
    Order qp_order = THIRD;
    if (p_order > 1)
        qp_order = FIFTH;
    libMesh::QGauss qrule(dim, qp_order);
    fe->attach_quadrature_rule(&qrule);
    // surface element
    std::unique_ptr < FEBase > fe_face(FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, qp_order);
    fe_face->attach_quadrature_rule(&qface);

    // quantities for volume integration
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<Point> &q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // quantities for surface integration
    const std::vector<Real> &JxW_face = fe_face->get_JxW();
    const std::vector<Point> &q_point_face = fe_face->get_xyz();
    const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi_face = fe_face->get_dphi();
    const std::vector<Point> &normal = fe_face->get_normals();

    // define fiber field
    double fx = data("fx", 1.0), fy = data("fy", 0.0), fz = data("fz", 0.0);
    double sx = data("sx", 0.0), sy = data("sy", 1.0), sz = data("sz", 0.0);
    double nx = data("nx", 0.0), ny = data("ny", 0.0), nz = data("nz", 1.0);
    VectorValue<Number>  f0(fx, fy, fz);
    VectorValue<Number>  s0(sx, sy, sz);
    VectorValue<Number>  n0(nx, ny, nz);
    // setup conductivities:
    // Default parameters from
    // Cardiac excitation mechanisms, wavefront dynamics and strength�interval
    // curves predicted by 3D orthotropic bidomain simulations

    // read parameters
    double sigma_f_i = data("sigma_f_i", 2.3172);
    double sigma_s_i = data("sigma_s_i", 0.2435); // sigma_t in the paper
    double sigma_n_i = data("sigma_n_i", 0.0569);
    double sigma_f_e = data("sigma_f_e", 1.5448);
    double sigma_s_e = data("sigma_s_e", 1.0438);  // sigma_t in the paper
    double sigma_n_e = data("sigma_n_e", 0.37222);
    double sigma_b_ie = data("sigma_b", 6.0);
    double chi = data("chi", 1e3);
    double Cm = data("Cm", 1.5);
    double penalty = data("penalty", 1e8);
    double interface_penalty = data("interface_penalty", 1e4);

    // setup tensors parameters
    // f0 \otimes f0
    TensorValue<Number> fof(fx * fx, fx * fy, fx * fz, fy * fx, fy * fy, fy * fz, fz * fx, fz * fy, fz * fz);
    // s0 \otimes s0
    TensorValue<Number> sos(sx * sx, sx * sy, sx * sz, sy * sx, sy * sy, sy * sz, sz * sx, sz * sy, sz * sz);
    // n0 \otimes n0
    TensorValue<Number> non(nx * nx, nx * ny, nx * nz, ny * nx, ny * ny, ny * nz, nz * nx, nz * ny, nz * nz);

    TensorValue<Number> sigma_i = ( sigma_f_i * fof + sigma_s_i * sos + sigma_n_i * non ) /chi;
    TensorValue<Number> sigma_e = ( sigma_f_e * fof + sigma_s_e * sos + sigma_n_e * non) / chi;
    TensorValue<Number> sigma_b(sigma_b_ie/chi, 0.0, 0.0, 0.0, sigma_b_ie/chi, 0.0, 0.0, 0.0, sigma_b_ie/chi);

    DenseMatrix < Number > K;
    DenseMatrix < Number > M;
    DenseVector < Number > M_lumped;

    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > parabolic_dof_indices;
    std::vector < dof_id_type > elliptic_dof_indices;

    // Assemble matrices differently based on time integrator
    switch (time_integrator)
    {
        case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
        case TimeIntegrator::EXPLICIT_INTRACELLULAR:
        case TimeIntegrator::SEMI_IMPLICIT:
        case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
        {
            std::cout << "Assembling matrices for EXPLICIT EXTRACELLULAR: " << static_cast<int>(time_integrator) << std::endl;
            libMesh::LinearImplicitSystem &elliptic_system = es.get_system < libMesh::LinearImplicitSystem > ("elliptic");
            libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            elliptic_system.zero_out_matrix_and_rhs = false;
            elliptic_system.assemble_before_solve = false;
            parabolic_system.zero_out_matrix_and_rhs = false;
            parabolic_system.assemble_before_solve = false;
            elliptic_system.matrix->zero();
            parabolic_system.matrix->zero();
            parabolic_system.get_matrix("Ki").zero();
            parabolic_system.get_matrix("Ke").zero();
            parabolic_system.get_matrix("Kg").zero();
            parabolic_system.get_matrix("M").zero();
            parabolic_system.get_vector("ML").zero();
            DenseMatrix < Number > Ki;
            DenseMatrix < Number > Ke;
            DenseMatrix < Number > Kg;
            DenseMatrix < Number > Ke_interface;
            DenseMatrix < Number > S;

            const DofMap & elliptic_dof_map = elliptic_system.get_dof_map();
            const DofMap & parabolic_dof_map = parabolic_system.get_dof_map();


            for (auto elem : mesh.element_ptr_range())
            {
                parabolic_dof_map.dof_indices(elem, parabolic_dof_indices);
                elliptic_dof_map.dof_indices(elem, elliptic_dof_indices);
                fe->reinit(elem);
                int elliptic_ndofs = elliptic_dof_indices.size();
                int parabolic_ndofs = parabolic_dof_indices.size();

                // resize local elemental matrices
                K.resize(elliptic_ndofs, elliptic_ndofs);
                Ki.resize(parabolic_ndofs, parabolic_ndofs);
                Ke.resize(parabolic_ndofs, parabolic_ndofs);
                Kg.resize(parabolic_ndofs, parabolic_ndofs);
                M.resize(parabolic_ndofs, parabolic_ndofs);
                S.resize(parabolic_ndofs, parabolic_ndofs);
                M_lumped.resize(parabolic_ndofs);

                auto subdomain_id = elem->subdomain_id();

                //std::cout << "Loop over volume: " << std::flush;

                // if we are in the bath
                if (subdomain_id == static_cast<short unsigned int>(Subdomain::BATH))
                {
                    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i != elliptic_ndofs; i++)
                        {
                            for (unsigned int j = 0; j != elliptic_ndofs; j++)
                            {
                                K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_b * dphi[j][qp]);
                            }
                        }
                    }
                }
                else
                {
                    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i != elliptic_ndofs; i++)
                        {
                            for (unsigned int j = 0; j != elliptic_ndofs; j++)
                            {
                                // Elliptic equation matrix
                                K(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i + sigma_e) * dphi[j][qp]);

                                // Parabolic matrices
                                // add mass matrix
                                M_lumped(i) += JxW[qp] * Cm / timedata.dt * phi[i][qp] * phi[j][qp];
                                // stiffness matrix
                                Ki(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i) * dphi[j][qp]);
                                Ke(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_e) * dphi[j][qp]);
                                M(i,  j) += JxW[qp] * phi[i][qp] * phi[j][qp];
                                S(i, j)  += JxW[qp] * dphi[i][qp] * ((sigma_i) * dphi[j][qp]);
                                S(i, i)  += JxW[qp] * Cm / timedata.dt * phi[i][qp] * phi[j][qp];
                            }
                        }
                    }
                }
                //std::cout << ", Loop over interface: " << std::flush;

                // interface
                if( time_integrator == TimeIntegrator::EXPLICIT_EXTRACELLULAR )
                {
                    if (subdomain_id == static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE))
                    {
                        for (auto side : elem->side_index_range())
                        {
                            if (elem->neighbor_ptr(side) != nullptr)
                            {
                                auto neighbor_subdomain_id = elem->neighbor_ptr(side)->subdomain_id();
                                if(subdomain_id != neighbor_subdomain_id)
                                {
                                    //std::cout << "Assembling interface" << std::endl;
                                    fe_face->reinit(elem, side);
                                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                                    {
                                        for (unsigned int i = 0; i != parabolic_ndofs; i++)
                                        {
                                            for (unsigned int j = 0; j != parabolic_ndofs; j++)
                                            {
                                                Ke(i, j) -= JxW_face[qp] * ( ( sigma_e * dphi_face[j][qp] ) * normal[qp] ) * phi_face[i][qp];
                                                Kg(i, j) -= interface_penalty * JxW_face[qp] * ( ( sigma_i * dphi_face[j][qp] ) * normal[qp] ) * phi_face[i][qp];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // boundaries
                //std::cout << ", Loop over boundaries: " << std::flush;
                if (subdomain_id == static_cast<libMesh::subdomain_id_type>(Subdomain::BATH))
                {
                    for (auto side : elem->side_index_range())
                    {
                        if (elem->neighbor_ptr(side) == nullptr)
                        {
                            auto boundary_id = mesh.get_boundary_info().boundary_id(elem, side);
                            auto it = bc_sidesets.find(static_cast<int>(boundary_id));
                            if (it != bc_sidesets.end())
                            {
                                //std::cout << "found_boundary " << std::endl;
                                fe_face->reinit(elem, side);
                                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                                {
                                    for (unsigned int i = 0; i != elliptic_ndofs; i++)
                                    {
                                        for (unsigned int j = 0; j != elliptic_ndofs; j++)
                                        {
                                            K(i, j) += JxW_face[qp] * penalty * phi_face[j][qp] * phi_face[i][qp];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                //std::cout << "Add matrix " << std::endl;

                elliptic_system.matrix->add_matrix (K, elliptic_dof_indices, elliptic_dof_indices);
                parabolic_system.matrix->add_matrix(S, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("M").add_matrix(M, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_vector("ML").add_vector(M_lumped, parabolic_dof_indices);
                parabolic_system.get_matrix("Ki").add_matrix(Ki, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Ke").add_matrix(Ke, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Kg").add_matrix(Kg, parabolic_dof_indices, parabolic_dof_indices);
            }
            std::cout << "closing" << std::endl;
            elliptic_system.matrix->close();
            parabolic_system.matrix->close();
            parabolic_system.get_matrix("M").close();
            parabolic_system.get_vector("ML").close();
            parabolic_system.get_matrix("Ki").close();
            parabolic_system.get_matrix("Ke").close();
            parabolic_system.get_matrix("Kg").close();
            std::cout << "done" << std::endl;
            break;
        }
        case TimeIntegrator::SBDF1:
            // Cm * M * ( V^n+1 - Vn ) / dt + Ki * V^n+1 + Ki * Ve^n+1 = -I^n
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //         -         -
            //        | Cm /dt * M + Ki,   Ki    |
            //   K  = | Ki,                Ki+Ke |;
            //         -         -
        case TimeIntegrator::SBDF2:
            // Cm * M * ( 3 * V^n+1 - 4 * Vn + V^n-1 ) / (2*dt) + Ki * V^n+1 + Ki * Ve^n+1 = -2*I^n + I^n-1
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //         -         -
            //        | 3/2 * Cm /dt * M + Ki,   Ki    |
            //   K  = | Ki,                      Ki+Ke |;
            //         -         -
        case TimeIntegrator::SBDF3:
            // Cm * M * ( 11/6 * V^n+1 - 3 * Vn + 3/2 V^n-1 -1/3 V^n-2) / (dt) + Ki * V^n+1 + Ki * Ve^n+1 = -3*I^n + 3*I^n-1 - I^n-2
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //         -         -
            //        | 11/6 * Cm /dt * M + Ki,   Ki    |
            //   K  = | Ki,                      Ki+Ke |;
            //         -         -
        default:
        {
            std::cout << "Assembling matrices for SBDF" << static_cast<int>(time_integrator) << std::endl;

            TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("bidomain");
            system.zero_out_matrix_and_rhs = false;
            system.assemble_before_solve = false;
            system.matrix->zero();
            system.get_matrix("M").zero();
            system.get_vector("ML").zero();
            const DofMap &dof_map = system.get_dof_map();
            int V_var_number = system.variable_number("V");
            int Ve_var_number = system.variable_number("Ve");

            DenseSubMatrix<Number> KVV(K), KVVe(K), KVeV(K), KVeVe(K);
            DenseSubMatrix < Number > MVV(M);
            DenseSubVector < Number > ML(M_lumped);

            double coefficient = Cm / timedata.dt;
            if (time_integrator == TimeIntegrator::SBDF2)
                coefficient *= 1.5;
            else if (time_integrator == TimeIntegrator::SBDF3)
                coefficient *= 11 / 6.0;
            std::cout << "SBDF coefficient: " << coefficient << std::endl;

            for (auto elem : mesh.element_ptr_range())
            {
                dof_map.dof_indices(elem, dof_indices);
                dof_map.dof_indices(elem, parabolic_dof_indices, V_var_number);
                dof_map.dof_indices(elem, elliptic_dof_indices, Ve_var_number);

                fe->reinit(elem);
                int ndofs = dof_indices.size();
                int elliptic_ndofs = elliptic_dof_indices.size();
                int parabolic_ndofs = parabolic_dof_indices.size();

                // resize local elemental matrices
                K.resize(ndofs, ndofs);
                M.resize(ndofs, ndofs);
                M_lumped.resize(ndofs);

                // reposition submatrices
                // Reposition the submatrices...  The idea is this:
                //
                //         -         -
                //        | KVV  KVVe  |
                //   K  = | KVeV KVeVe |;
                //         -         -

                // if we are in the tissue we are using the submatrices
                KVV.reposition(V_var_number * parabolic_ndofs, V_var_number * parabolic_ndofs, parabolic_ndofs, parabolic_ndofs);
                KVVe.reposition(V_var_number * parabolic_ndofs, Ve_var_number * parabolic_ndofs, parabolic_ndofs, elliptic_ndofs);
                KVeV.reposition(Ve_var_number * elliptic_ndofs, V_var_number * elliptic_ndofs, elliptic_ndofs, parabolic_ndofs);
                KVeVe.reposition(Ve_var_number * elliptic_ndofs, Ve_var_number * elliptic_ndofs, elliptic_ndofs, elliptic_ndofs);
                MVV.reposition(V_var_number * parabolic_ndofs, V_var_number * parabolic_ndofs, parabolic_ndofs, parabolic_ndofs);
                ML.reposition(V_var_number * parabolic_ndofs, parabolic_ndofs);

                auto subdomain_id = elem->subdomain_id();
                // if we are in the bath
                if (subdomain_id == static_cast<short unsigned int>(Subdomain::BATH))
                {
                    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i != elliptic_ndofs; i++)
                        {
                            for (unsigned int j = 0; j != elliptic_ndofs; j++)
                            {
                                K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_b * dphi[j][qp]);
                            }
                        }
                    }
                }
                else
                {
                    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i != elliptic_ndofs; i++)
                        {
                            for (unsigned int j = 0; j != elliptic_ndofs; j++)
                            {
                                // add mass matrix
                                if (p_order == SECOND) // we cannot do the lumping
                                    KVV(i, j) += JxW[qp] * coefficient * phi[i][qp] * phi[j][qp];
                                else
                                    // we can do the lumping
                                    KVV(i, i) += JxW[qp] * coefficient * phi[i][qp] * phi[j][qp];

                                // stiffness matrix
                                KVV(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i) * dphi[j][qp]);
                                KVVe(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i) * dphi[j][qp]);
                                KVeV(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i) * dphi[j][qp]);
                                KVeVe(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i + sigma_e) * dphi[j][qp]);

                                // These are useful for the RHS
                                MVV(i, j) += JxW[qp] * phi[i][qp] * phi[j][qp];
                                ML(i) += JxW[qp] * phi[i][qp] * phi[j][qp];
                            }
                        }
                    }
                }
                // boundaries
                if (subdomain_id == static_cast<libMesh::subdomain_id_type>(Subdomain::BATH))
                {
                    for (auto side : elem->side_index_range())
                    {
                        if (elem->neighbor_ptr(side) == nullptr)
                        {
                            auto boundary_id = mesh.get_boundary_info().boundary_id(elem, side);
                            auto it = bc_sidesets.find(static_cast<int>(boundary_id));
                            if (it != bc_sidesets.end())
                            {
                                fe_face->reinit(elem, side);
                                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                                {
                                    for (unsigned int i = 0; i != elliptic_ndofs; i++)
                                    {
                                        for (unsigned int j = 0; j != elliptic_ndofs; j++)
                                        {
                                            K(i, j) += JxW_face[qp] * penalty * phi_face[j][qp] * phi_face[i][qp];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                system.matrix->add_matrix (K, dof_indices, dof_indices);
                system.get_matrix("M").add_matrix(M, dof_indices, dof_indices);
                system.get_vector("ML").add_vector(M_lumped, dof_indices);
            }
            system.matrix->close();
            system.get_matrix("M").close();
            system.get_vector("ML").close();
            break;
        }
    }
}

void assemble_rhs(  libMesh::EquationSystems &es,
                    const TimeData &timedata,
                    TimeIntegrator time_integrator,
                    libMesh::Order p_order,
                    const GetPot &data,
                    EquationType type)
{
    using namespace libMesh;
    double Cm = data("Cm", 1.0);
    //std::cout << " assemble_rhs " << std::endl;
    switch (time_integrator)
    {
        case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
        case TimeIntegrator::EXPLICIT_INTRACELLULAR:
        case TimeIntegrator::SEMI_IMPLICIT:
        case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
         {
            const MeshBase& mesh = es.get_mesh();
            libMesh::LinearImplicitSystem &elliptic_system = es.get_system < libMesh::LinearImplicitSystem > ("elliptic");
            libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            const DofMap & elliptic_dof_map = elliptic_system.get_dof_map();
            const DofMap & parabolic_dof_map = parabolic_system.get_dof_map();
            std::vector < dof_id_type > parabolic_dof_indices;
            std::vector < dof_id_type > elliptic_dof_indices;

            if(type == EquationType::PARABOLIC)
            {
                parabolic_system.get_vector("aux1").zero();
                parabolic_system.rhs->zero();

                // Transfer Ve from elliptic to parabolic
                for (auto node : mesh.node_ptr_range())
                {
                    parabolic_dof_map.dof_indices(node, parabolic_dof_indices);
                    elliptic_dof_map.dof_indices(node, elliptic_dof_indices);
                    int elliptic_ndofs = elliptic_dof_indices.size();
                    int parabolic_ndofs = parabolic_dof_indices.size();

                    if (elliptic_ndofs == parabolic_ndofs)
                    {
                        double ven = (*elliptic_system.solution)(elliptic_dof_indices[0]);
                        parabolic_system.get_vector("aux1").set(parabolic_dof_indices[0], ven);
                    }
                }
                parabolic_system.get_vector("aux1").close();

                if(time_integrator == TimeIntegrator::EXPLICIT_EXTRACELLULAR)
                {
                    // Assemble RHS:
                    parabolic_system.rhs->add_vector(parabolic_system.get_vector("aux1"), parabolic_system.get_matrix("Ke"));
                    parabolic_system.rhs->add_vector(parabolic_system.get_vector("aux1"), parabolic_system.get_matrix("Kg"));
                    parabolic_system.rhs->add_vector(*parabolic_system.solution, parabolic_system.get_matrix("Kg"));
                }
                else if(time_integrator == TimeIntegrator::EXPLICIT_INTRACELLULAR)
                {
                    // Assemble RHS:
                    parabolic_system.rhs->add_vector(parabolic_system.get_vector("aux1"), parabolic_system.get_matrix("Ki"));
                    parabolic_system.rhs->add_vector(*parabolic_system.solution, parabolic_system.get_matrix("Ki"));
                    parabolic_system.rhs->scale(-1.0);
                }
                else // SEMI IMPLICIT
                {
                    // Assemble RHS:
                    parabolic_system.rhs->add_vector(parabolic_system.get_vector("aux1"), parabolic_system.get_matrix("Ki"));
                    parabolic_system.rhs->scale(-1.0);
                }
                parabolic_system.get_vector("aux1").zero();
                parabolic_system.get_vector("aux1").add(-1.0, parabolic_system.get_vector("In"));
                // add  M * In
                parabolic_system.rhs->add_vector(parabolic_system.get_vector("aux1"), parabolic_system.get_matrix("M"));

                if(time_integrator == TimeIntegrator::SEMI_IMPLICIT ||
                   time_integrator == TimeIntegrator::SEMI_IMPLICIT_HALF_STEP)
                {
                    parabolic_system.get_vector("aux1").zero();
                    parabolic_system.get_vector("aux1").pointwise_mult(*parabolic_system.solution, parabolic_system.get_vector("ML"));
                    *parabolic_system.rhs += parabolic_system.get_vector("aux1");
                }
                else
                {
                    (*parabolic_system.rhs) /= parabolic_system.get_vector("ML");
                }

            }
            else
            {
                parabolic_system.get_vector("aux1").zero();
                elliptic_system.rhs->zero();

                parabolic_system.get_vector("aux1").add_vector(*parabolic_system.solution, parabolic_system.get_matrix("Ki"));
                // Transfer KiV from parabolic to elliptic
                for (auto node : mesh.node_ptr_range())
                {
                    parabolic_dof_map.dof_indices(node, parabolic_dof_indices);
                    elliptic_dof_map.dof_indices(node, elliptic_dof_indices);
                    int elliptic_ndofs = elliptic_dof_indices.size();
                    int parabolic_ndofs = parabolic_dof_indices.size();

                    if (elliptic_ndofs == parabolic_ndofs)
                    {
                        double KiV = parabolic_system.get_vector("aux1")(parabolic_dof_indices[0]);
                        elliptic_system.rhs->set(elliptic_dof_indices[0], -KiV);
                    }
                }
            }
            break;
        }
        case TimeIntegrator::SBDF3:
            // Cm * M * ( 11/6 * V^n+1 - 3 * Vn + 3/2 V^n-1 -1/3 V^n-2) / (dt) + Ki * V^n+1 + Ki * Ve^n+1 = -3*I^n + 3*I^n-1 - I^n-2
            // Ki * V^n+1 + Kie * Ve^n+1 = 0            //
            //   RHS = 3 Cm /dt * M * Vn - 3/2 Cm /dt * M * Vnm1 +Cm/dt/3 M Vnm2 - 3 * M In + 3 Inm1 - Inm2

        {
            TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("bidomain");
            //system.get_vector("In").print();
            system.rhs->zero();
            system.get_vector("aux1").zero();
            // eval: 3 Cm /dt * M * Vn
            system.get_vector("aux1").add(3.0*Cm/timedata.dt, *system.old_local_solution);
            // add: -1.5 Cm /dt * M * Vnm1
            system.get_vector("aux1").add(-1.5*Cm/timedata.dt, *system.older_local_solution);
            // add: +1/3 Cm /dt * M * Vnm2
            system.get_vector("aux1").add(Cm/(3*timedata.dt), system.get_vector("Vnm2"));
            //system.get_vector("aux1").print();

            if(p_order == SECOND)
            {
                system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            }
            else
            {
                system.rhs->pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            }
            // add  M * (2*In-Inm1)
            system.get_vector("aux1").zero();
            system.get_vector("aux1").add(-3.0, system.get_vector("In"));
            system.get_vector("aux1").add(3.0, system.get_vector("Inm1"));
            system.get_vector("aux1").add(-1.0, system.get_vector("Inm2"));
            //system.get_vector("aux1").print();

            // add  M * In
            system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            break;
        }

        case TimeIntegrator::SBDF2:
            // Cm * M * ( 3 * V^n+1 - 4 * Vn + V^n-1 ) / (2*dt) + Ki * V^n+1 + Ki * Ve^n+1 = -2*I^n + I^n-1
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //
            //   RHS = 2 Cm /dt * M * Vn - 1/2 Cm /dt * M * Vnm1 + 2 * M In - Inm1
            // Cm /dt M ( 2 * Vn - 1/2  Vnm1)
        {
            TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("bidomain");
            system.rhs->zero();
            system.get_vector("aux1").zero();
            //system.get_vector("In").print();
            // eval: 2 Cm /dt * M * Vn
            system.get_vector("aux1").add(2*Cm/timedata.dt, *system.old_local_solution);
            // add: -0.5 Cm /dt * M * Vnm1
            system.get_vector("aux1").add(-.5*Cm/timedata.dt, *system.older_local_solution);
            //system.get_vector("aux1").print();
            if(p_order == SECOND)
            {
                system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            }
            else
            {
                system.rhs->pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            }
            // add  M * (2*In-Inm1)
            system.get_vector("aux1").zero();
            system.get_vector("aux1").add(-2.0, system.get_vector("In"));
            system.get_vector("aux1").add(1.0, system.get_vector("Inm1"));
            //system.get_vector("aux1").print();
            // add  M * In
            system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            //system.get_vector("aux1").pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            //*system.rhs += system.get_vector("aux1");
            break;
        }
        case TimeIntegrator::SBDF1:
            // Cm * M * ( V^n+1 - Vn ) / dt + Ki * V^n+1 + Ki * Ve^n+1 = -I^n
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //
            //   RHS = Cm /dt * M * Vn + M In
        default:
        {

            TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("bidomain");
            system.rhs->zero();
            system.get_vector("aux1").zero();
            // eval: Cm /dt * M * Vn
            system.get_vector("aux1").add(Cm/timedata.dt, *system.old_local_solution);
            if(p_order == SECOND)
            {
                system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            }
            else
            {
                system.rhs->pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            }
            // add  M * In
            system.get_vector("aux1").zero();
            system.get_vector("aux1").add(-1.0, system.get_vector("In"));
            // add  M * In
            system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            break;
        }
    }
}




// TODO: update for SBDF2 and SBDF3
void solve_ionic_model_evaluate_ionic_currents(libMesh::EquationSystems &es,
                                               IonicModel& ionic_model,
                                               Pacing& pacing,
                                               TimeData& datatime,
                                               TimeIntegrator &time_integrator)
{
    using namespace libMesh;
    const MeshBase& mesh = es.get_mesh();
    if( time_integrator == TimeIntegrator::EXPLICIT_EXTRACELLULAR ||
        time_integrator == TimeIntegrator::EXPLICIT_INTRACELLULAR ||
        time_integrator == TimeIntegrator::SEMI_IMPLICIT ||
        time_integrator == TimeIntegrator::SEMI_IMPLICIT_HALF_STEP )
    {
        libMesh::LinearImplicitSystem &elliptic_system = es.get_system < libMesh::LinearImplicitSystem > ("elliptic");
        libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
        const DofMap & elliptic_dof_map = elliptic_system.get_dof_map();
        const DofMap & parabolic_dof_map = parabolic_system.get_dof_map();
        std::vector < dof_id_type > parabolic_dof_indices;
        std::vector < dof_id_type > elliptic_dof_indices;

        for (auto node : mesh.node_ptr_range())
        {
            parabolic_dof_map.dof_indices(node, parabolic_dof_indices);
            elliptic_dof_map.dof_indices(node, elliptic_dof_indices);
            int elliptic_ndofs = elliptic_dof_indices.size();
            int parabolic_ndofs = parabolic_dof_indices.size();

            if (elliptic_ndofs == parabolic_ndofs)
            {
                double stimulus = pacing.istim(*node, datatime.time);
                double vn = (*parabolic_system.solution)(parabolic_dof_indices[0]);
                double i_ion = ionic_model.iion(vn);
                double i_tot = i_ion + stimulus;
                parabolic_system.get_vector("In").set(parabolic_dof_indices[0], i_tot);
            }
        }
        parabolic_system.get_vector("In").close();
    }
    else
    {
        TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("bidomain");
        const DofMap &dof_map = system.get_dof_map();
        int V_var_number = system.variable_number("V");
        int Ve_var_number = system.variable_number("Ve");
        std::vector < dof_id_type > dof_indices;
        std::vector < dof_id_type > parabolic_dof_indices;
        std::vector < dof_id_type > elliptic_dof_indices;

        for (auto node : mesh.node_ptr_range())
        {
            dof_map.dof_indices(node, parabolic_dof_indices, V_var_number);
            dof_map.dof_indices(node, elliptic_dof_indices, Ve_var_number);
            int elliptic_ndofs = elliptic_dof_indices.size();
            int parabolic_ndofs = parabolic_dof_indices.size();

            if (elliptic_ndofs == parabolic_ndofs)
            {
                double stimulus = pacing.istim(*node, datatime.time);
                double vn = (*system.solution)(parabolic_dof_indices[0]);
                double i_ion = ionic_model.iion(vn);
                double i_tot = i_ion + stimulus;
                system.get_vector("In").set(parabolic_dof_indices[0], i_tot);
            }
        }
        system.get_vector("In").close();
    }
}


// read BC sidesets from string: e.g. bc = "1 2 3 5", or bc = "1, 55, 2, 33"
void read_bc_list(std::string &bc, std::set<int> &bc_sidesets)
{
    std::string number;
    for (auto it = bc.begin(); it != bc.end(); it++)
    {
        auto character = *it;
        if (std::isdigit(character))
            number += character;
        else
        {
            if (number.size() > 0)
            {
                bc_sidesets.insert(std::stoi(number));
                number.clear();
            }
        }
        if (it + 1 == bc.end())
        {
            bc_sidesets.insert(std::stoi(number));
        }
    }
    std::cout << "BC sideset list: " << std::flush;
    for (auto &&sideset : bc_sidesets)
        std::cout << sideset << ", " << std::flush;
    std::cout << std::endl;
}
