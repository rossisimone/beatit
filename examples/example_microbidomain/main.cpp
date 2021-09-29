// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/gmv_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "petsctime.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petsc/private/kspimpl.h"
#include "libmesh/getpot.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
// This example will solve a linear transient system,
// so we need to include the TransientLinearImplicitSystem definition.
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"
#include "libmesh/fe_interface.h"
// The definition of a geometric element
#include "libmesh/elem.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/mesh_modification.h"
#include "Electrophysiology/IonicModels/IonicModel.hpp"
#include "Electrophysiology/IonicModels/Fabbri17.hpp"
#include "Util/IO/io.hpp"
// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side at each time step.  Note that
// since the system is linear we technically do not need to
// assemble the matrix at each time step, but we will anyway.
// In subsequent examples we will employ adaptive mesh refinement,
// and with a changing mesh it will be necessary to rebuild the
// system matrix.
void assemble_matrix(EquationSystems & es);
void assemble_matrix_separate_qrs(EquationSystems & es);
void assemble_rhs(EquationSystems & es);
void advance_ionic_model(EquationSystems & es);
void assemble_matrix_lagrange_multiplier(EquationSystems & es);
//void assemble_rhs_lagrange_multiplier(EquationSystems & es);

// Function prototype.  This function will initialize the system.
// Initialization functions are optional for systems.  They allow
// you to specify the initial values of the solution.  If an
// initialization function is not provided then the default (0)
// solution is provided.
void init_phi(EquationSystems & es);


enum class Formulation { PENALTY,
                         LAGRANGE_MULTIPLIER };
// Exact solution function prototype.  This gives the exact
// solution as a function of space and time.  In this case the
// initial condition will be taken as the exact solution at time 0,
// as will the Dirichlet boundary conditions at time t.
namespace Methods
{
    static Formulation formulation;
}

namespace MeshUtility
{
    // For each element on the interface we store
    // id of the interface element,
    // pointer to the element on the extracellular space,
    // pointer to the element on the intracellular space,
    static std::unordered_map< unsigned int, std::pair<Elem *, Elem *> > interface_elements_connectivity;

    static std::unordered_map< unsigned int, Elem * > extracellular_elements_connectivity;

    static std::unordered_map< unsigned int, Elem * > intracellular_elements_connectivity;

    // cells are rectangles
    struct Cell
    {
        Cell( double a, double b, double c, double d) : xmin(a), ymin(b), xmax(c), ymax(d) {}
        double xmin;
        double ymin;
        double xmax;
        double ymax;
    };

    struct Cells
    {
        void setup(GetPot& data)
        {
            std::string xmins = data("xmins", "0.0");
            std::string ymins = data("ymins", "0.0");
            std::string xmaxs = data("xmaxs", "1.0");
            std::string ymaxs = data("ymaxs", "1.0");

            std::vector<double> xmin, ymin, xmax, ymax;
            BeatIt::readList(xmins, xmin);
            BeatIt::readList(ymins, ymin);
            BeatIt::readList(xmaxs, xmax);
            BeatIt::readList(ymaxs, ymax);
            // settin up number of cells
            N = xmin.size();

            for ( int k = 0; k < N; ++k)
            {
                cells.emplace_back(xmin[k], ymin[k], xmax[k], ymax[k]);
            }

        }
        std::vector<Cell> cells;
        int N;
    }
}



namespace IonicModel
{
    static std::unique_ptr<BeatIt::IonicModel> ionic_model;
    static std::vector<double> variables;
    static int num_vars;
}
// We can now begin the main program.  Note that this
// example will fail if you are using complex numbers
// since it was designed to be run only with real numbers.
int main(int argc, char ** argv) {
    // Initialize libMesh.
    LibMeshInit init(argc, argv);
    
    // Read the mesh from file.  This is the coarse mesh that will be used
    // in example 10 to demonstrate adaptive mesh refinement.  Here we will
    // simply read it in and uniformly refine it 5 times before we compute
    // with it.
    //
    // Create a mesh object, with dimension to be overridden later,
    // distributed across the default MPI communicator.
    ReplicatedMesh mesh(init.comm());
    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data.beat", 2, "-i", "--input");
    GetPot data(datafile_name);
    
    
    std::string meshname = data("input_mesh", "NONE");
    
    if(meshname == "NONE")
    {
        std::cout << "Creating mesh: " << std::endl;
                  double xmin = data("xmin", 0.0);
        double ymin = data("ymin", 0.0);
        double xmax = data("xmax", 1.0);
        double ymax = data("ymax", 1.0);
        int elx = data("elx", 10);
        int ely = data("ely", 10);

        MeshUtility::Cells cells;
        cells.setup(data);

        libMesh::MeshTools::build_square(mesh, elx, ely, xmin, xmax, ymin, ymax, QUAD4);
        for(auto & elem : mesh.local_active_element_ptr_range)
        {
            libMesh::Point X = elem->centroid();
            int extracellullar_blockID = 1;
            int intracellullar_blockID = 2;

            elem->subdomain_id() = extracellullar_blockID;
            for( int k = 0; k < cells.N ; ++k)
            {
                if( X(0) > cells.cells[k].xmin &&
                    X(0) < cells.cells[k].xmax &&
                    X(1) > cells.cells[k].ymin &&
                    X(1) < cells.cells[k].xmax )
                {
                    elem->subdomain_id() = intracellullar_blockID;
                }
            }
        }
    }
    else
    {
        std::cout << "Reading mesh: " << std::endl;
        mesh.read(meshname);
    }
    double scale = data("scale", 1.0);
    MeshTools::Modification::scale(mesh, scale, scale, scale);
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    // Preprocess mesh
    std::cout << "Preprocessing interface: " << std::endl;
    // change sideset ID
    int counter = 0;
    int nn_index = mesh.max_elem_id() + 1;
    const libMesh::MeshBase::const_element_iterator end_el = mesh.elements_end();
    for (libMesh::MeshBase::const_element_iterator el = mesh.elements_begin(); el != end_el; ++el)
    {
        Elem *const elem = *el;
        unsigned int blockID = elem->subdomain_id();

        ElemType elType = elem->type();

        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            if (elem->neighbor_ptr(side) != libmesh_nullptr)
            {
                Elem *neighbor = elem->neighbor_ptr(side);
                unsigned int neighborID = neighbor->subdomain_id();
                // If we are at the interface let's add a sideset
                if (1 == blockID && 2 <= neighborID)
                {
                    counter ++ ;
                    libMesh::BoundaryInfo *boundary_info = mesh.boundary_info.get();
                    boundary_info->add_side(elem, side, 3);
                    Elem * side_el;
                    if( elType == TRI3 || elType == QUAD4 ) side_el = new Edge2;
                    else if( elType == TRI6) side_el = new Edge3;
                    side_el->set_id(nn_index);
                    nn_index++;
                    side_el = mesh.add_elem(side_el);
                    side_el->subdomain_id() = 3;
                    auto side_elem = elem->build_side_ptr(side);
                    side_el->set_node(0) = side_elem->node_ptr(0);
                    side_el->set_node(1) = side_elem->node_ptr(1);
                    if( elType == TRI6) side_el->set_node(2) = side_elem->node_ptr(2);
                    MeshUtility::interface_elements_connectivity[ side_el->id() ] = std::pair<Elem *, Elem *> ( elem, elem->neighbor_ptr(side) );
                    MeshUtility::extracellular_elements_connectivity[ elem->id() ] = side_el;
                    MeshUtility::intracellular_elements_connectivity[ elem->neighbor_ptr(side)->id() ] = side_el;
                }
                if (1 == neighborID && 2 <= blockID)
                {
                    libMesh::BoundaryInfo *boundary_info = mesh.boundary_info.get();
                    boundary_info->add_side(elem, side, 3);
                }
            }
        }
    }
//    mesh.all_second_order();
    mesh.prepare_for_use(true);
    mesh.print_info();
    std::set< subdomain_id_type > ids;
    mesh.subdomain_ids(ids);
    for(auto && id : ids) std::cout << id << std::endl;

    std::cout << "Done: " << counter << std::endl;


    double dt = data("dt",0.0125);
    double Cm = data("Cm",0.01); // pF / um^2
    double sigma_i = data("sigma_i",0.3); // uS / um
    double sigma_e = data("sigma_e",2.0); // uS / um
    double penalty = data("penalty",1e8); // pF / um^2
    std::string order = data("order", "SECOND");

    Methods::formulation = Formulation::PENALTY;
    if(penalty <= 0.0) Methods::formulation = Formulation::LAGRANGE_MULTIPLIER;
    Order phi_order = FIRST;
    if(order == "SECOND") phi_order = SECOND;

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);
    


    // Add a transient system to the EquationSystems
    // object named "Convection-Diffusion".
    TransientLinearImplicitSystem & system = equation_systems.add_system
    < TransientLinearImplicitSystem > ("MicroBidomain");
    

    std::set<unsigned short> extracellular_space;
    std::set<unsigned short> intracellular_space;
    std::set<unsigned short> cell_membrane;
    unsigned int extracellular_space_blockID = 1;
    unsigned int intracellular_space_blockID = 2;
    unsigned int cell_membrane_blockID = 3;
    extracellular_space.insert(extracellular_space_blockID);
    intracellular_space.insert(intracellular_space_blockID);
    cell_membrane.insert(cell_membrane_blockID);
    int n_subdomains = mesh.n_subdomains();
    if( 3 < n_subdomains )
    {
        for(int blockID = 4; blockID <= n_subdomains; blockID++)  intracellular_space.insert(blockID);
    }


    system.add_variable("phi_e", phi_order, LAGRANGE, &extracellular_space);
    system.add_variable("phi_i", phi_order, LAGRANGE, &intracellular_space);
    system.add_variable("phi_m", phi_order, LAGRANGE, &cell_membrane);
    // if using lagrange multiplier add it to the system variables
    if(Methods::formulation == Formulation::LAGRANGE_MULTIPLIER)
    system.add_variable("lambda", FIRST, LAGRANGE, &cell_membrane);
    
    TransientExplicitSystem & system_iion = equation_systems.add_system
    < TransientExplicitSystem > ("Iion");
    system_iion.add_variable("iion", phi_order, LAGRANGE, &cell_membrane);

    IonicModel::ionic_model.reset( BeatIt::IonicModel::IonicModelFactory::Create("Fabbri17") );
    IonicModel::num_vars = IonicModel::ionic_model->numVariables();
    IonicModel::variables.resize(IonicModel::num_vars);
    auto gating_variables_names = IonicModel::ionic_model->variablesNames();
    for (auto && g_var : gating_variables_names)
        system_iion.add_vector(g_var);

    IonicModel:: ionic_model->setup(data,"ionic_model");

//    // Give the system a pointer to the matrix assembly
//    // and initialization functions.
//    //system.attach_assemble_function(assemble_cd);
//    system.attach_init_function(init_cd);
    
    // Initialize the data structures for the equation system.
    std::cout << "Initialize Equation systems" << std::endl;
    equation_systems.init();
    
    // Prints information about the system to the screen.
    std::cout << "Print info about the equation systems" << std::endl;
    equation_systems.print_info();

    // set parameters
    std::cout << "Set Parameters" << std::endl;

    equation_systems.parameters.set < Real > ("dt") = dt;
    equation_systems.parameters.set < Real > ("sigma_i") = sigma_i;
    equation_systems.parameters.set < Real > ("sigma_e") = sigma_e;
    equation_systems.parameters.set < Real > ("Cm") = Cm;
    equation_systems.parameters.set < Real > ("penalty") = penalty;

    std::cout << "Set initial conditions" << std::endl;
    init_phi(equation_systems);


    std::cout << "Create Output" << std::endl;
    ExodusII_IO exo(mesh);
    std::string exodus_filename = data("output", "output.e");
    exo.write_equation_systems(exodus_filename, equation_systems);
    exo.append(true);
    int timestep = 1;
    double time = 0.0;
    exo.write_timestep(exodus_filename, equation_systems, timestep, time);

    std::cout << "Assemble Matrix" << std::endl;
    if(Methods::formulation == Formulation::LAGRANGE_MULTIPLIER)
    {
        assemble_matrix_lagrange_multiplier(equation_systems);
        std::cout << "Assemble Matrix Done" << std::endl;
    }
    else assemble_matrix_separate_qrs(equation_systems);

    system.assemble_before_solve = false;

    std::cout << "Start Timeloop" << std::endl;
    int n_steps = data("n_steps",10);
    int save_steps = data("save_steps",10);
    int save_step = timestep;
    for (unsigned int t_step = 0; t_step < n_steps; t_step++)
    {
        timestep++;
        time += dt;
        std::cout << "Time: " << time << std::flush;
        system.update();
        *system.old_local_solution = *system.current_local_solution;
        std::cout << " Advancing Ionic Model ... " << std::flush;
        advance_ionic_model(equation_systems);
        std::cout << " Assembling RHS ... " << std::flush;
        assemble_rhs(equation_systems);
        std::cout << " Solving ... " << std::flush;
        system.solve();
        if(t_step % save_steps == 0)
        {
            save_step ++;
            std::cout << " Exporting ... " << std::flush;
             exo.write_timestep(exodus_filename, equation_systems, save_step, time);
        }
         std::cout << std::endl;
    }
    // All done.
    return 0;
}




// We now define the function which provides the
// initialization routines for the "Convection-Diffusion"
// system.  This handles things like setting initial
// conditions and boundary conditions.
void init_phi(EquationSystems & es)
{
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("MicroBidomain");
    TransientExplicitSystem & system_iion = es.get_system
    < TransientExplicitSystem > ("Iion");

    const DofMap & dof_map = system.get_dof_map();
    const DofMap & dof_map_iion = system_iion.get_dof_map();

    
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_phi_i;
    std::vector < dof_id_type > dof_indices_phi_e;
    std::vector < dof_id_type > dof_indices_phi_m;
    std::vector < dof_id_type > dof_indices_iion;

    IonicModel::ionic_model->initialize(IonicModel::variables);
    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  Since the mesh
    // will be refined we want to only consider the ACTIVE elements,
    // hence we use a variant of the active_elem_iterator.
    libMesh::BoundaryInfo *boundary_info = mesh.boundary_info.get();

    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices(elem, dof_indices);

        unsigned int blockID = elem->subdomain_id();
        unsigned int elID = elem->id();
        double init_value = 0.0;
        switch(blockID)
        {
            default:
            case 1:
            {
                init_value = 0.0;
                for (auto && index : dof_indices) system.solution->set(index, init_value);
                break;
            }
            case 2:
            case 4:
            {
                init_value = 0.0;
                for (auto && index : dof_indices) system.solution->set(index, init_value);
                break;
            }
            case 3:
            {
                double offset = 0.0;
                auto it =MeshUtility::interface_elements_connectivity.find(elID);
                if(it != MeshUtility::interface_elements_connectivity.end())
                {
                    if(4 == it->second.second->subdomain_id() ) offset = 60.0;
                }
                int num_nodes = elem->n_nodes();
                for(unsigned int n = 0; n < num_nodes; ++n)
                {
                    Node * node = elem->node_ptr(n);
                    dof_map.dof_indices(node, dof_indices_phi_m, 2);
                    dof_map_iion.dof_indices(node, dof_indices_iion);
                    for (auto && index : dof_indices_phi_m) system.solution->set(index, IonicModel::variables[0]+offset);
                    for (int nv = 0; nv < IonicModel::num_vars-1; nv++)
                    {
                        for (auto && index : dof_indices_iion) system_iion.get_vector(nv).set(index, IonicModel::variables[nv+1]);
                    }
                }
                break;
            }
        }

    }
    system.solution->close();
    for (int nv = 0; nv < IonicModel::num_vars-1; nv++)
    {
        system_iion.get_vector(nv).close();
    }

}

// Now we define the assemble function which will be used
// by the EquationSystems object at each timestep to assemble
// the linear system for solution.
void assemble_matrix(EquationSystems & es)
{
    std::cout << "Assmble Matrix with penalty method" << std::endl;
   // Ignore unused parameter warnings when !LIBMESH_ENABLE_AMR.
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("MicroBidomain");

    {
        system.matrix->zero();
        typedef libMesh::PetscMatrix<libMesh::Number> Mat;
        Mat * mat = dynamic_cast<Mat *>(system.matrix);
        MatSetOption(mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }

    
    unsigned int phi_e_var = 0;
    unsigned int phi_i_var = 1;
    unsigned int phi_m_var = 2;
    
    Real dt = es.parameters.get<Real>("dt");
    Real Cm = es.parameters.set < Real > ("Cm"); // uF mm^-2
    Real sigma_i = es.parameters.set < Real > ("sigma_i"); // mS mm^-1
    Real sigma_e = es.parameters.set < Real > ("sigma_e"); // mS mm^-1
    Real penalty = es.parameters.set < Real > ("penalty");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_phi_e_type = system.variable_type(phi_e_var);
    FEType fe_phi_i_type = system.variable_type(phi_i_var);
    FEType fe_phi_m_type = system.variable_type(phi_m_var);
    
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr < FEBase > fe_phi_e(FEBase::build(dim, fe_phi_e_type));
    std::unique_ptr < FEBase > fe_phi_i(FEBase::build(dim, fe_phi_i_type));
    std::unique_ptr < FEBase > fe_phi_m(FEBase::build(dim-1, fe_phi_m_type));
    std::unique_ptr < FEBase > fe_phi_m_SRI(FEBase::build(dim-1, fe_phi_m_type));

    
    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    QGauss qrule_phi_e(dim,   libMesh::SECOND);
    QGauss qrule_phi_i(dim,   libMesh::SECOND);
    QGauss qrule_phi_m(dim-1, libMesh::SECOND);
    QGauss qrule_phi_m_SRI(dim-1, libMesh::FIRST);
    
    // Tell the finite element object to use our quadrature rule.
    fe_phi_e->attach_quadrature_rule(&qrule_phi_e);
    fe_phi_i->attach_quadrature_rule(&qrule_phi_i);
    fe_phi_m->attach_quadrature_rule(&qrule_phi_m);
    fe_phi_m_SRI->attach_quadrature_rule(&qrule_phi_m_SRI);
    
    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.  We will start
    // with the element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW_e = fe_phi_e->get_JxW();
    const std::vector<Real> & JxW_i = fe_phi_i->get_JxW();
    const std::vector<Real> & JxW_m = fe_phi_m->get_JxW();
    const std::vector<Real> & JxW_m_SRI = fe_phi_m_SRI->get_JxW();
    
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> & phi_i = fe_phi_i->get_phi();
    const std::vector<std::vector<Real>> & phi_e = fe_phi_e->get_phi();
    const std::vector<std::vector<Real>> & phi_m = fe_phi_m->get_phi();
    const std::vector<std::vector<Real>> & phi_m_SRI = fe_phi_m_SRI->get_phi();
    
    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>> & dphi_i = fe_phi_i->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_e = fe_phi_e->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_m = fe_phi_m->get_dphi();
    
    const std::vector<Point >& qps_phi_m =fe_phi_m->get_xyz();
    const std::vector<Point >& qps_phi_m_SRI =fe_phi_m_SRI->get_xyz();


    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix < Number > Kee, Kei, Kem,
                           Kie, Kii, Kim,
                           Kme, Kmi, Kmm;
    
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_phi_i;
    std::vector < dof_id_type > dof_indices_phi_e;
    std::vector < dof_id_type > dof_indices_phi_m;

    // Boundary elements
    // Boundary integration requires one quadraure rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    std::unique_ptr<libMesh::FEBase> fe_phi_e_face (libMesh::FEBase::build(dim, fe_phi_e_type));
    libMesh::QGauss qface_phi_e(dim-1, libMesh::SECOND);

    fe_phi_e_face->attach_quadrature_rule (&qface_phi_e);
    const std::vector< libMesh::Real> & JxW_e_face = fe_phi_e_face->get_JxW();
    const std::vector<std::vector< libMesh::Real> > & phi_e_face = fe_phi_e_face->get_phi();
    const std::vector<Point> & qface_e_normals = fe_phi_e_face->get_normals();
    const std::vector<Point >& qps_phi_e_face =fe_phi_e_face->get_xyz();

    std::unique_ptr<libMesh::FEBase> fe_phi_i_face (libMesh::FEBase::build(dim, fe_phi_i_type));
    libMesh::QGauss qface_phi_i(dim-1, libMesh::SECOND);
    fe_phi_i_face->attach_quadrature_rule (&qface_phi_i);
    const std::vector<Point >& qps_phi_i_face =fe_phi_i_face->get_xyz();
    const std::vector<libMesh::Real> &JxW_i_face = fe_phi_i_face->get_JxW();
    const std::vector<std::vector<libMesh::Real> > &phi_i_face = fe_phi_i_face->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi_i_face = fe_phi_i_face->get_dphi();
    const std::vector<Point> & qface_i_normals = fe_phi_i_face->get_normals();



    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  Since the mesh
    // will be refined we want to only consider the ACTIVE elements,
    // hence we use a variant of the active_elem_iterator.
    for (const auto & elem : mesh.active_local_element_ptr_range()) {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        unsigned int blockID = elem->subdomain_id();
        unsigned int elID = elem->id();
        
        
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_phi_e, phi_e_var);
        dof_map.dof_indices(elem, dof_indices_phi_i, phi_i_var);
        dof_map.dof_indices(elem, dof_indices_phi_m, phi_m_var);
        
        const unsigned int n_dofs = dof_indices.size();
        unsigned int n_phi_i_dofs = dof_indices_phi_i.size();
        unsigned int n_phi_e_dofs = dof_indices_phi_e.size();
        unsigned int n_phi_m_dofs = dof_indices_phi_m.size();
        

        Kee.resize(n_phi_e_dofs, n_phi_e_dofs);
        Kii.resize(n_phi_i_dofs, n_phi_i_dofs);
        Kmm.resize(n_phi_m_dofs, n_phi_m_dofs);

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        switch(blockID)
        {
            // Extracellular space
            case 1:
            {
                fe_phi_e->reinit(elem);
                for (unsigned int qp = 0; qp < qrule_phi_e.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_e.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_e.size(); j++)
                        {
                            Kee(i,j) += JxW_e[qp] * sigma_e * dphi_e[i][qp] * dphi_e[j][qp];
                        }
                    }
                }

                for (unsigned int side=0; side<elem->n_sides(); side++)
                {
                    // Homogeneuous Dirichlet BC
                    if (elem->neighbor_ptr(side) == libmesh_nullptr)
                    {
                        // The value of the shape functions at the quadrature
                        // points.
                        const std::vector<std::vector< libMesh::Real> > & phi_face = fe_phi_e_face->get_phi();

                        // The Jacobian * Quadrature Weight at the quadrature
                        // points on the face.
                        const std::vector< libMesh::Real> & JxW_face = fe_phi_e_face->get_JxW();

                        // Compute the shape function values on the element
                        // face.
                        fe_phi_e_face->reinit(elem, side);
                        // Loop over the face quadrature points for integration.
                        for (unsigned int qp=0; qp<qface_phi_e.n_points(); qp++)
                        {
                            for (unsigned int i=0; i<phi_face.size(); i++)
                                for (unsigned int j=0; j<phi_face.size(); j++)
                                    Kee(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
                        }

                    }
                    else
                    {
                        // Flux conditions at the cell membrane interface
                        if( mesh.boundary_info->has_boundary_id(elem, side, 3) )
                        {
                            const Elem * neighbor = elem->neighbor_ptr(side);

                            // The quadrature point locations on the neighbor side
                            std::vector<Point> qface_neighbor_point;

                            // Compute the shape function values on the element
                            // face.
                            fe_phi_e_face->reinit(elem, side);

                            FEInterface::inverse_map (elem->dim(),
                                                      fe_phi_e->get_fe_type(),
                                                      neighbor,
                                                      qps_phi_e_face,
                                                      qface_neighbor_point);
                            // Calculate the neighbor element shape functions at those locations
                            fe_phi_i_face->reinit(neighbor, &qface_neighbor_point);

                            dof_map.dof_indices(neighbor, dof_indices_phi_i, phi_i_var);
                            n_phi_i_dofs = dof_indices_phi_i.size();

                            Kei.resize(n_phi_e_dofs, n_phi_i_dofs);

                            for (unsigned int qp=0; qp<qface_phi_e.n_points(); qp++)
                            {
                                for (unsigned int i=0; i<phi_e_face.size(); i++)
                                    for (unsigned int j=0; j<dphi_i_face.size(); j++)
                                        Kei(i,j) -= JxW_e_face[qp] * sigma_i * ( dphi_i_face[j][qp]*qface_e_normals[qp] ) * phi_e_face[i][qp];

                            }

                            system.matrix->add_matrix(Kei, dof_indices_phi_e, dof_indices_phi_i);
                        }
                    } // if on boundary
                } // loop on sides
                system.matrix->add_matrix(Kee, dof_indices_phi_e, dof_indices_phi_e);

                break;
            }
            // Intracellular space
            default:
            case 2:
            {
                fe_phi_i->reinit(elem);
                // Laplacian
                for (unsigned int qp = 0; qp < qrule_phi_i.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_i.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_i.size(); j++)
                        {
                            Kii(i, j) += JxW_i[qp] * sigma_i * dphi_i[i][qp] * dphi_i[j][qp];
                        }
                    }
                }
                // Dirichlet at the interface
                for (unsigned int side = 0; side < elem->n_sides(); side++)
                {
                    // Penalty Dirichlet BC
                    if (elem->neighbor_ptr(side) != libmesh_nullptr)
                    {
                        if (mesh.boundary_info->has_boundary_id(elem, side, 3))
                        {
                            // The value of the shape functions at the quadrature
                            // points.

                            // The Jacobian * Quadrature Weight at the quadrature
                            // points on the face.

                            // Compute the shape function values on the element
                            // face.
                            fe_phi_i_face->reinit(elem, side);
                            // Loop over the face quadrature points for integration.
                            for (unsigned int qp = 0; qp < qface_phi_i.n_points(); qp++)
                            {
                                for (unsigned int i = 0; i < phi_i_face.size(); i++)
                                    for (unsigned int j = 0; j < phi_i_face.size(); j++)
                                        Kii(i, j) += JxW_i_face[qp] * penalty * phi_i_face[i][qp] * phi_i_face[j][qp];
                            }

                            const Elem * neighbor = elem->neighbor_ptr(side);

                            // The quadrature point locations on the neighbor side
                            std::vector<Point> qface_neighbor_point;
                            // The quadrature point locations on the element side

                            FEInterface::inverse_map (elem->dim(),
                                                      fe_phi_i->get_fe_type(),
                                                      neighbor,
                                                      qps_phi_i_face,
                                                      qface_neighbor_point);
                            // Calculate the neighbor element shape functions at those locations
                            fe_phi_e_face->reinit(neighbor, &qface_neighbor_point);

                            dof_map.dof_indices(neighbor, dof_indices_phi_e, phi_e_var);
                            n_phi_e_dofs = dof_indices_phi_e.size();

                            Kie.resize(n_phi_i_dofs, n_phi_e_dofs);

                            for (unsigned int qp=0; qp<qface_phi_i.n_points(); qp++)
                            {
                                for (unsigned int i=0; i<phi_i_face.size(); i++)
                                    for (unsigned int j=0; j<phi_e_face.size(); j++)
                                        Kie(i,j) -= JxW_i_face[qp] * penalty * phi_e_face[j][qp] * phi_i_face[i][qp];
                            }

                            system.matrix->add_matrix(Kie, dof_indices_phi_i, dof_indices_phi_e);
                        }
                    }
                }
                system.matrix->add_matrix(Kii, dof_indices_phi_i, dof_indices_phi_i);

                break;
            }
            // cell membrane
            case 3:
            {
                std::cout << "Assembling block 3: " << std::flush;
                fe_phi_m->reinit(elem);
                fe_phi_m_SRI->reinit(elem);
                // Mass Matrix
                std::cout << "Kmm: " << std::flush;

                for (unsigned int qp = 0; qp < qrule_phi_m.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_m.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_m.size(); j++)
                        {
                            Kmm(i, j) += JxW_m[qp] * phi_m[i][qp] * phi_m[j][qp];
                        }
                    }
                }
                system.matrix->add_matrix(Kmm, dof_indices_phi_m, dof_indices_phi_m);

                auto it = MeshUtility::interface_elements_connectivity.find(elID);
                if(it != MeshUtility::interface_elements_connectivity.end() )
                {
                    // Calculate the neighbor element shape functions at those locations
                    Elem * intracellular_elem = it->second.second;
                    Elem * extracellular_elem = it->second.first;

                    dof_map.dof_indices(intracellular_elem, dof_indices_phi_i, phi_i_var);
                    n_phi_i_dofs = dof_indices_phi_i.size();


                    std::vector<Point> qface_intracellular_point;

                    FEInterface::inverse_map (intracellular_elem->dim(),
                                              fe_phi_i->get_fe_type(),
                                              intracellular_elem,
                                              qps_phi_m,
                                              qface_intracellular_point);

                    fe_phi_i_face->reinit(intracellular_elem, &qface_intracellular_point);

                    // Assemble Kmi

                    // Note:
                    // I'm using the extracellular element to evaluate the normal
                    // This is because reinitializing the intracellular element on the qface_intracellular_point
                    // does not update the normal
                    // I will assume that the normal is constant over the element side
                    // and I will use the only the normal at the first quadrature point
                    int side = extracellular_elem->which_neighbor_am_i(intracellular_elem);
                    fe_phi_e_face->reinit(extracellular_elem, side);

                    Kmi.resize(n_phi_m_dofs, n_phi_i_dofs);
//                    std::cout << "Kmi: " << std::flush;
//                    std::cout << n_phi_m_dofs << ", " << n_phi_i_dofs << ", " <<  phi_m.size() << ", " << phi_i_face.size() << std::endl;
//                    std::cout << qrule_phi_m.n_points() << ", " << JxW_m.size() << ", " <<  phi_m[0].size() << ", " <<  dphi_i_face[0].size() << ", " << qface_e_normals.size() << std::endl;
                    for (unsigned int qp = 0; qp < qrule_phi_m.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i < phi_m.size(); i++)
                        {
                            for (unsigned int j = 0; j < phi_i_face.size(); j++)
                            {
                                Kmi(i, j) -= JxW_m[qp] * dt / Cm * phi_m[i][qp] * ( sigma_i * dphi_i_face[j][qp] * qface_e_normals[qp]);
                            }
                        }
                    }

                    system.matrix->add_matrix(Kmi, dof_indices_phi_m, dof_indices_phi_i);


                    // Assemble Kim penalty


                    // To test selective reduced integration
                    // this integral will be evaluated using the SRI shape functions
                    // I will have to reinitialize the intracellular shape functions
                    qface_intracellular_point.clear();
                    FEInterface::inverse_map (intracellular_elem->dim(),
                                              fe_phi_i->get_fe_type(),
                                              intracellular_elem,
                                              qps_phi_m_SRI,
                                              qface_intracellular_point);
                    fe_phi_i_face->reinit(intracellular_elem, &qface_intracellular_point);

                    Kim.resize(n_phi_i_dofs, n_phi_m_dofs);

                    std::cout << "Kim: " << std::flush;
                   for (unsigned int qp = 0; qp < qrule_phi_m_SRI.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i < phi_i_face.size(); i++)
                        {
                            for (unsigned int j = 0; j < phi_m_SRI.size(); j++)
                            {
                                Kim(i, j) -= JxW_m_SRI[qp] * penalty * phi_m_SRI[j][qp] * phi_i_face[i][qp];
                            }
                        }
                    }
                    system.matrix->add_matrix(Kim, dof_indices_phi_i, dof_indices_phi_m);
                    std::cout << "done! " << std::endl;




                }
                break;
            }
//            default:
//            {
//                std::cout << "We should not be here!" << std::endl;
//                throw std::runtime_error("");
//                break;
//            }

        }

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The SparseMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
    }
    // That concludes the system matrix assembly routine.
    system.matrix->close();
    std::cout <<"Matrix assembled"  << std::endl;
   // system.matrix->print();
}


// Now we define the assemble function which will be used
// by the EquationSystems object at each timestep to assemble
// the linear system for solution.
void assemble_matrix_separate_qrs(EquationSystems & es)
{
    std::cout << "Assmble Matrix with penalty method" << std::endl;
   // Ignore unused parameter warnings when !LIBMESH_ENABLE_AMR.
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("MicroBidomain");

    {
        system.matrix->zero();
        typedef libMesh::PetscMatrix<libMesh::Number> Mat;
        Mat * mat = dynamic_cast<Mat *>(system.matrix);
        MatSetOption(mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }


    unsigned int phi_e_var = 0;
    unsigned int phi_i_var = 1;
    unsigned int phi_m_var = 2;

    Real dt = es.parameters.get<Real>("dt");
    Real Cm = es.parameters.set < Real > ("Cm"); // uF mm^-2
    Real sigma_i = es.parameters.set < Real > ("sigma_i"); // mS mm^-1
    Real sigma_e = es.parameters.set < Real > ("sigma_e"); // mS mm^-1
    Real penalty = es.parameters.set < Real > ("penalty");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_phi_e_type = system.variable_type(phi_e_var);
    FEType fe_phi_i_type = system.variable_type(phi_i_var);
    FEType fe_phi_m_type = system.variable_type(phi_m_var);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr < FEBase > fe_phi_e(FEBase::build(dim, fe_phi_e_type));
    std::unique_ptr < FEBase > fe_phi_i(FEBase::build(dim, fe_phi_i_type));
    std::unique_ptr < FEBase > fe_phi_m(FEBase::build(dim-1, fe_phi_m_type));

    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    QGauss qrule_first(dim,   libMesh::FIRST);
    QGauss qrule_second(dim,  libMesh::SECOND);
    QGauss qface_first(dim-1,   libMesh::FIRST);
    QGauss qface_second(dim-1,  libMesh::SECOND);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.  We will start
    // with the element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW_e = fe_phi_e->get_JxW();
    const std::vector<Real> & JxW_i = fe_phi_i->get_JxW();
    const std::vector<Real> & JxW_m = fe_phi_m->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> & phi_i = fe_phi_i->get_phi();
    const std::vector<std::vector<Real>> & phi_e = fe_phi_e->get_phi();
    const std::vector<std::vector<Real>> & phi_m = fe_phi_m->get_phi();
    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>> & dphi_i = fe_phi_i->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_e = fe_phi_e->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_m = fe_phi_m->get_dphi();
    const std::vector<Point >& qps_phi_m =fe_phi_m->get_xyz();


    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix < Number > Kee, Kei, Kem,
                           Kie, Kii, Kim,
                           Kme, Kmi, Kmm;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_phi_i;
    std::vector < dof_id_type > dof_indices_phi_e;
    std::vector < dof_id_type > dof_indices_phi_m;

    // Boundary elements
    // Boundary integration requires one quadraure rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    std::unique_ptr<libMesh::FEBase> fe_phi_e_face (libMesh::FEBase::build(dim, fe_phi_e_type));
    libMesh::QGauss qface1(dim-1, libMesh::SECOND);
    const std::vector< libMesh::Real> & JxW_e_face = fe_phi_e_face->get_JxW();
    const std::vector<std::vector< libMesh::Real> > & phi_e_face = fe_phi_e_face->get_phi();
    const std::vector<Point> & qface_e_normals = fe_phi_e_face->get_normals();
    const std::vector<Point >& qps_phi_e_face =fe_phi_e_face->get_xyz();
    std::unique_ptr<libMesh::FEBase> fe_phi_i_face (libMesh::FEBase::build(dim, fe_phi_i_type));
    libMesh::QGauss qface_phi_i(dim-1, libMesh::SECOND);
    const std::vector<Point >& qps_phi_i_face =fe_phi_i_face->get_xyz();
    const std::vector<libMesh::Real> &JxW_i_face = fe_phi_i_face->get_JxW();
    const std::vector<std::vector<libMesh::Real> > &phi_i_face = fe_phi_i_face->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi_i_face = fe_phi_i_face->get_dphi();
    const std::vector<Point> & qface_i_normals = fe_phi_i_face->get_normals();



    bool assemble_outer_dirichlet = true;
    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  Since the mesh
    // will be refined we want to only consider the ACTIVE elements,
    // hence we use a variant of the active_elem_iterator.
    for (const auto & elem : mesh.active_local_element_ptr_range()) {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        unsigned int blockID = elem->subdomain_id();
        unsigned int elID = elem->id();


        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_phi_e, phi_e_var);
        dof_map.dof_indices(elem, dof_indices_phi_i, phi_i_var);
        dof_map.dof_indices(elem, dof_indices_phi_m, phi_m_var);

        const unsigned int n_dofs = dof_indices.size();
        unsigned int n_phi_i_dofs = dof_indices_phi_i.size();
        unsigned int n_phi_e_dofs = dof_indices_phi_e.size();
        unsigned int n_phi_m_dofs = dof_indices_phi_m.size();


        Kee.resize(n_phi_e_dofs, n_phi_e_dofs);
        Kii.resize(n_phi_i_dofs, n_phi_i_dofs);
        Kmm.resize(n_phi_m_dofs, n_phi_m_dofs);

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        switch(blockID)
        {
            // Extracellular space
            case 1:
            {

                fe_phi_e->attach_quadrature_rule(&qrule_first);
                fe_phi_e->reinit(elem);
                for (unsigned int qp = 0; qp < qrule_first.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_e.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_e.size(); j++)
                        {
                            Kee(i,j) += JxW_e[qp] * sigma_e * dphi_e[i][qp] * dphi_e[j][qp];
                        }
                    }
                }

                for (unsigned int side=0; side<elem->n_sides(); side++)
                {
                    // Homogeneuous Dirichlet BC
                    if (elem->neighbor_ptr(side) == libmesh_nullptr && true == assemble_outer_dirichlet)
                    {
                        // The value of the shape functions at the quadrature
                        // points.
                        const std::vector<std::vector< libMesh::Real> > & phi_face = fe_phi_e_face->get_phi();

                        // The Jacobian * Quadrature Weight at the quadrature
                        // points on the face.
                        const std::vector< libMesh::Real> & JxW_face = fe_phi_e_face->get_JxW();

                        // Compute the shape function values on the element
                        // face.
                        fe_phi_e_face->attach_quadrature_rule (&qface_second);

                        fe_phi_e_face->reinit(elem, side);
                        // Loop over the face quadrature points for integration.
                        for (unsigned int qp=0; qp<qface_second.n_points(); qp++)
                        {
                            for (unsigned int i=0; i<phi_face.size(); i++)
                                for (unsigned int j=0; j<phi_face.size(); j++)
                                    Kee(i,j) += 0*JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
                        }
                        assemble_outer_dirichlet = false;
                    }
                    else
                    {
                        // Flux conditions at the cell membrane interface
                        if( mesh.boundary_info->has_boundary_id(elem, side, 3) )
                        {
                            const Elem * neighbor = elem->neighbor_ptr(side);

                            // The quadrature point locations on the neighbor side
                            std::vector<Point> qface_neighbor_point;

                            // Compute the shape function values on the element
                            // face.
                            fe_phi_e_face->attach_quadrature_rule (&qface_first);
                            fe_phi_e_face->reinit(elem, side);

                            FEInterface::inverse_map (elem->dim(),
                                                      fe_phi_e->get_fe_type(),
                                                      neighbor,
                                                      qps_phi_e_face,
                                                      qface_neighbor_point);
                            // Calculate the neighbor element shape functions at those locations
                            fe_phi_i_face->reinit(neighbor, &qface_neighbor_point);

                            dof_map.dof_indices(neighbor, dof_indices_phi_i, phi_i_var);
                            n_phi_i_dofs = dof_indices_phi_i.size();

                            Kei.resize(n_phi_e_dofs, n_phi_i_dofs);

                            for (unsigned int qp=0; qp<qface_first.n_points(); qp++)
                            {
                                for (unsigned int i=0; i<phi_e_face.size(); i++)
                                    for (unsigned int j=0; j<dphi_i_face.size(); j++)
                                        Kei(i,j) -= JxW_e_face[qp] * sigma_i * ( dphi_i_face[j][qp]*qface_e_normals[qp] ) * phi_e_face[i][qp];

                            }

                            system.matrix->add_matrix(Kei, dof_indices_phi_e, dof_indices_phi_i);
                        }
                    } // if on boundary
                } // loop on sides
                system.matrix->add_matrix(Kee, dof_indices_phi_e, dof_indices_phi_e);

                break;
            }
            // Intracellular space
            case 2:
            default:
            {
                fe_phi_i->attach_quadrature_rule (&qrule_first);

                fe_phi_i->reinit(elem);
                // Laplacian
                for (unsigned int qp = 0; qp < qrule_first.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_i.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_i.size(); j++)
                        {
                            Kii(i, j) += JxW_i[qp] * sigma_i * dphi_i[i][qp] * dphi_i[j][qp];
                        }
                    }
                }
                // Dirichlet at the interface
                for (unsigned int side = 0; side < elem->n_sides(); side++)
                {
                    // Penalty Dirichlet BC
                    if (elem->neighbor_ptr(side) != libmesh_nullptr)
                    {
                        if (mesh.boundary_info->has_boundary_id(elem, side, 3))
                        {
                            // The value of the shape functions at the quadrature
                            // points.

                            // The Jacobian * Quadrature Weight at the quadrature
                            // points on the face.

                            // Compute the shape function values on the element
                            // face.
                            fe_phi_i_face->attach_quadrature_rule (&qface_second);

                            fe_phi_i_face->reinit(elem, side);
                            // Loop over the face quadrature points for integration.
                            for (unsigned int qp = 0; qp < qface_second.n_points(); qp++)
                            {
                                for (unsigned int i = 0; i < phi_i_face.size(); i++)
                                    for (unsigned int j = 0; j < phi_i_face.size(); j++)
                                        Kii(i, j) += JxW_i_face[qp] * penalty * phi_i_face[i][qp] * phi_i_face[j][qp];
                            }

                            const Elem * neighbor = elem->neighbor_ptr(side);

                            // The quadrature point locations on the neighbor side
                            std::vector<Point> qface_neighbor_point;
                            // The quadrature point locations on the element side

                            FEInterface::inverse_map (elem->dim(),
                                                      fe_phi_i->get_fe_type(),
                                                      neighbor,
                                                      qps_phi_i_face,
                                                      qface_neighbor_point);
                            // Calculate the neighbor element shape functions at those locations
                            fe_phi_e_face->reinit(neighbor, &qface_neighbor_point);

                            dof_map.dof_indices(neighbor, dof_indices_phi_e, phi_e_var);
                            n_phi_e_dofs = dof_indices_phi_e.size();

                            Kie.resize(n_phi_i_dofs, n_phi_e_dofs);

                            for (unsigned int qp=0; qp<qface_second.n_points(); qp++)
                            {
                                for (unsigned int i=0; i<phi_i_face.size(); i++)
                                    for (unsigned int j=0; j<phi_e_face.size(); j++)
                                        Kie(i,j) -= JxW_i_face[qp] * penalty * phi_e_face[j][qp] * phi_i_face[i][qp];
                            }

                            system.matrix->add_matrix(Kie, dof_indices_phi_i, dof_indices_phi_e);

                            auto it = MeshUtility::intracellular_elements_connectivity.find(elID);
                            if(it != MeshUtility::intracellular_elements_connectivity.end() )
                            {
                                Elem * interface_elem = it->second;
                                dof_map.dof_indices(interface_elem, dof_indices_phi_m, phi_m_var);
                                n_phi_m_dofs = dof_indices_phi_m.size();

                                // The quadrature point locations on the neighbor side
                                std::vector<Point> qface_point;
                                FEInterface::inverse_map (elem->dim()-1,
                                                          fe_phi_m->get_fe_type(),
                                                          interface_elem,
                                                          qps_phi_i_face,
                                                          qface_point);
                                // Calculate the neighbor element shape functions at those locations
                                fe_phi_m->reinit(interface_elem, &qface_point);

                                Kim.resize(n_phi_i_dofs, n_phi_m_dofs);
                                for (unsigned int qp=0; qp<qface_second.n_points(); qp++)
                                {
                                    for (unsigned int i=0; i<phi_i_face.size(); i++)
                                        for (unsigned int j=0; j<phi_m.size(); j++)
                                            Kim(i,j) -= JxW_i_face[qp] * penalty * phi_m[j][qp] * phi_i_face[i][qp];
                                }
                                system.matrix->add_matrix(Kim, dof_indices_phi_i, dof_indices_phi_m);


                            }
                        }
                    }
                }
                system.matrix->add_matrix(Kii, dof_indices_phi_i, dof_indices_phi_i);

                break;
            }
            // cell membrane
            case 3:
            {
                fe_phi_m->attach_quadrature_rule (&qface_second);
                fe_phi_m->reinit(elem);
                // Mass Matrix
                for (unsigned int qp = 0; qp < qface_second.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_m.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_m.size(); j++)
                        {
                            Kmm(i, j) += JxW_m[qp] * phi_m[i][qp] * phi_m[j][qp];
                        }
                    }
                }
                system.matrix->add_matrix(Kmm, dof_indices_phi_m, dof_indices_phi_m);

                auto it = MeshUtility::interface_elements_connectivity.find(elID);
                if(it != MeshUtility::interface_elements_connectivity.end() )
                {
                    // Calculate the neighbor element shape functions at those locations
                    std::cout << "Prepare the elements: " << std::flush;
                    Elem * intracellular_elem = it->second.second;
                    Elem * extracellular_elem = it->second.first;

                    dof_map.dof_indices(intracellular_elem, dof_indices_phi_i, phi_i_var);
                    n_phi_i_dofs = dof_indices_phi_i.size();

                    fe_phi_i_face->attach_quadrature_rule (&qface_first);
                    fe_phi_e_face->attach_quadrature_rule (&qface_first);
                    fe_phi_m->attach_quadrature_rule (&qface_first);
                    fe_phi_m->reinit(elem);

                    std::vector<Point> qface_intracellular_point;

                    FEInterface::inverse_map (intracellular_elem->dim(),
                                              fe_phi_i->get_fe_type(),
                                              intracellular_elem,
                                              qps_phi_m,
                                              qface_intracellular_point);

                    fe_phi_i_face->reinit(intracellular_elem, &qface_intracellular_point);

                    // Assemble Kmi

                    // Note:
                    // I'm using the extracellular element to evaluate the normal
                    // This is because reinitializing the intracellular element on the qface_intracellular_point
                    // does not update the normal
                    // I will assume that the normal is constant over the element side
                    // and I will use the only the normal at the first quadrature point
                    int side = extracellular_elem->which_neighbor_am_i(intracellular_elem);
                    fe_phi_e_face->reinit(extracellular_elem, side);

                    Kmi.resize(n_phi_m_dofs, n_phi_i_dofs);
                    for (unsigned int qp = 0; qp < qface_first.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i < phi_m.size(); i++)
                        {
                            for (unsigned int j = 0; j < phi_i_face.size(); j++)
                            {
                                Kmi(i, j) -= JxW_m[qp] * dt / Cm * phi_m[i][qp] * ( sigma_i * dphi_i_face[j][qp] * qface_e_normals[qp]);
                            }
                        }
                    }

                    system.matrix->add_matrix(Kmi, dof_indices_phi_m, dof_indices_phi_i);


                    // Assemble Kim penalty

                    std::cout << "Prepare the elements 2: " << std::flush;

//                    // To test selective reduced integration
//                    // this integral will be evaluated using the SRI shape functions
//                    // I will have to reinitialize the intracellular shape functions
//                    fe_phi_i_face->attach_quadrature_rule (&qface_second);
//                    fe_phi_m->attach_quadrature_rule (&qface_second);
//                    fe_phi_m->reinit(elem);
//
//                    qface_intracellular_point.clear();
//                    FEInterface::inverse_map (intracellular_elem->dim(),
//                                              fe_phi_i->get_fe_type(),
//                                              intracellular_elem,
//                                              qps_phi_m,
//                                              qface_intracellular_point);
//                    fe_phi_i_face->reinit(intracellular_elem, &qface_intracellular_point);
//
//                    Kim.resize(n_phi_i_dofs, n_phi_m_dofs);
//
//                    std::cout << "Kim: " << std::flush;
//                   for (unsigned int qp = 0; qp < qface_second.n_points(); qp++)
//                    {
//                        for (unsigned int i = 0; i < phi_i_face.size(); i++)
//                        {
//                            for (unsigned int j = 0; j < phi_m.size(); j++)
//                            {
//                                Kim(i, j) -= JxW_m[qp] * penalty * phi_m[j][qp] * phi_i_face[i][qp];
//                            }
//                        }
//                    }
//                    system.matrix->add_matrix(Kim, dof_indices_phi_i, dof_indices_phi_m);
                    std::cout << "done! " << std::endl;




                }
                break;
            }
//            default:
//            {
//                std::cout << "We should not be here!" << std::endl;
//                throw std::runtime_error("");
//                break;
//            }

        }

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The SparseMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
    }
    // That concludes the system matrix assembly routine.
    system.matrix->close();
    std::cout <<"Matrix assembled"  << std::endl;
   // system.matrix->print();
}



// Now we define the assemble function which will be used
// by the EquationSystems object at each timestep to assemble
// the linear system for solution.
void assemble_matrix_lagrange_multiplier(EquationSystems & es)
{
    std::cout << "Assmble Matrix with Lagrange Multiplier" << std::endl;

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("MicroBidomain");

    // Increase matrix sparsity
    {
        system.matrix->zero();
        typedef libMesh::PetscMatrix<libMesh::Number> Mat;
        Mat * mat = dynamic_cast<Mat *>(system.matrix);
        MatSetOption(mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }


    // extracellular potential
    unsigned int phi_e_var = 0;
    // intracellular potential
    unsigned int phi_i_var = 1;
    // transmembrane potential
    unsigned int phi_m_var = 2;
    // Lagrange Multiplier
    unsigned int lambda_var = 3;

    // Get parameters
    Real dt = es.parameters.get<Real>("dt");
    Real Cm = es.parameters.set < Real > ("Cm"); // uF mm^-2
    Real sigma_i = es.parameters.set < Real > ("sigma_i"); // mS mm^-1
    Real sigma_e = es.parameters.set < Real > ("sigma_e"); // mS mm^-1
    Real penalty = es.parameters.set < Real > ("penalty");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_phi_e_type = system.variable_type(phi_e_var);
    FEType fe_phi_i_type = system.variable_type(phi_i_var);
    FEType fe_phi_m_type = system.variable_type(phi_m_var);
    FEType fe_lambda_type = system.variable_type(lambda_var);

    // Build a Finite Element object of the specified type.
    // extracellular FE
    std::unique_ptr < FEBase > fe_phi_e(FEBase::build(dim, fe_phi_e_type));
    std::unique_ptr<libMesh::FEBase> fe_phi_e_face (libMesh::FEBase::build(dim, fe_phi_e_type));

    // extracellular FE = extracellular FE
    std::unique_ptr < FEBase > fe_phi_i(FEBase::build(dim, fe_phi_i_type));
    std::unique_ptr<libMesh::FEBase> fe_phi_i_face (libMesh::FEBase::build(dim, fe_phi_e_type));

    // Transmembrane FE on the surface
    std::unique_ptr < FEBase > fe_phi_m(FEBase::build(dim-1, fe_phi_m_type));
    // Lagrange Multiplier FE on the surface
    std::unique_ptr < FEBase > fe_lambda(FEBase::build(dim-1, fe_lambda_type));


    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    QGauss qrule_first(dim,   libMesh::FIRST);
    QGauss qrule_second(dim,   libMesh::SECOND);
    QGauss qface_first(dim-1,   libMesh::FIRST);
    QGauss qface_second(dim-1,   libMesh::SECOND);
    QGauss qface_fifth(dim-1,   libMesh::FIFTH);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.  We will start
    // with the element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW_e = fe_phi_e->get_JxW();
    const std::vector<Real> & JxW_e_face = fe_phi_e_face->get_JxW();
    const std::vector<Real> & JxW_i = fe_phi_i->get_JxW();
    const std::vector<Real> & JxW_i_face = fe_phi_i_face->get_JxW();
    const std::vector<Real> & JxW_m = fe_phi_m->get_JxW();
    const std::vector<Real> & JxW_lambda = fe_lambda->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> & phi_i = fe_phi_i->get_phi();
    const std::vector<std::vector<Real>> & phi_i_face = fe_phi_i_face->get_phi();
    const std::vector<std::vector<Real>> & phi_e = fe_phi_e->get_phi();
    const std::vector<std::vector<Real>> & phi_e_face = fe_phi_e_face->get_phi();
    const std::vector<std::vector<Real>> & phi_m = fe_phi_m->get_phi();
    const std::vector<std::vector<Real>> & phi_lambda = fe_lambda->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>> & dphi_i = fe_phi_i->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_e = fe_phi_e->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_i_face = fe_phi_i_face->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_e_face = fe_phi_e_face->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_lambda = fe_lambda->get_dphi();

    // Quadrature Points on surface
    const std::vector<Point >& qps_phi_m      = fe_phi_m->get_xyz();
    const std::vector<Point >& qps_phi_lambda = fe_lambda->get_xyz();
    const std::vector<Point >& qps_phi_i_face = fe_phi_i_face->get_xyz();
    const std::vector<Point >& qps_phi_e_face = fe_phi_e_face->get_xyz();

    // normals
    const std::vector<Point> & qface_phi_i_normals = fe_phi_i_face->get_normals();
    const std::vector<Point> & qface_phi_e_normals = fe_phi_e_face->get_normals();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    //DenseMatrix < Number > Kee, Kii, Kmm, Kei, Kie, Kmi, Kim, Kme, Kem;
    DenseMatrix < Number > Kee, Kei, Kem, Kel,
                           Kie, Kii, Kim, Kil,
                           Kme, Kmi, Kmm, Kml,
                           Kle, Kli, Klm, Kll;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_phi_i;
    std::vector < dof_id_type > dof_indices_phi_e;
    std::vector < dof_id_type > dof_indices_phi_m;
    std::vector < dof_id_type > dof_indices_lambda;


    int counter = 0;

    std::cout << "Loop over elements" << std::endl;
    for (const auto & elem : mesh.active_local_element_ptr_range()) {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        unsigned int blockID = elem->subdomain_id();
        unsigned int elID = elem->id();


        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_phi_e, phi_e_var);
        dof_map.dof_indices(elem, dof_indices_phi_i, phi_i_var);
        dof_map.dof_indices(elem, dof_indices_phi_m, phi_m_var);
        dof_map.dof_indices(elem, dof_indices_lambda, lambda_var);

        const unsigned int n_dofs = dof_indices.size();
        unsigned int n_phi_i_dofs = dof_indices_phi_i.size();
        unsigned int n_phi_e_dofs = dof_indices_phi_e.size();
        unsigned int n_phi_m_dofs = dof_indices_phi_m.size();
        unsigned int n_lambda_dofs = dof_indices_lambda.size();



        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        switch(blockID)
        {
            // Extracellular space
            case 1:
            {
                // Here we need only a single point quadrature rule
                Kee.resize(n_phi_e_dofs, n_phi_e_dofs);
                fe_phi_e->attach_quadrature_rule(&qrule_second);
                fe_phi_e->reinit(elem);
                for (unsigned int qp = 0; qp < qrule_second.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_e.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_e.size(); j++)
                        {
                            Kee(i,j) += JxW_e[qp] * sigma_e * dphi_e[i][qp] * dphi_e[j][qp];
                        }
                    }
                }

                for (unsigned int side=0; side<elem->n_sides(); side++)
                {
                    // Homogeneuous Dirichlet BC on the physical boundary
                    if (elem->neighbor_ptr(side) == libmesh_nullptr)
                    {
                        // Compute the shape function values on the element
                        // face.
                        fe_phi_e_face->attach_quadrature_rule(&qface_second);
                        fe_phi_e_face->reinit(elem, side);
                        // Loop over the face quadrature points for integration.
                        for (unsigned int qp=0; qp<qface_second.n_points(); qp++)
                        {
                            for (unsigned int i=0; i<phi_e_face.size(); i++)
                                for (unsigned int j=0; j<phi_e_face.size(); j++)
                                    Kee(i,j) += JxW_e_face[qp]*1e6*phi_e_face[i][qp]*phi_e_face[j][qp];
                        }
                    }
                    else
                    {
                        // Flux conditions at the cell membrane interface
                        if( mesh.boundary_info->has_boundary_id(elem, side, 3) )
                        {
                            auto it = MeshUtility::extracellular_elements_connectivity.find(elID);
                            if(it != MeshUtility::extracellular_elements_connectivity.end() )
                            {

                                Elem * interface_elem = it->second;
                                dof_map.dof_indices(interface_elem, dof_indices_lambda, lambda_var);
                                n_lambda_dofs = dof_indices_lambda.size();

                                // The quadrature point locations on the neighbor side
                                std::vector<Point> qface_point;

                                // Compute the shape function values on the element
                                // face.
                                fe_phi_e_face->attach_quadrature_rule(&qface_second);
                                fe_phi_e_face->reinit(elem, side);

                                FEInterface::inverse_map (elem->dim()-1,
                                                          fe_lambda->get_fe_type(),
                                                          interface_elem,
                                                          qps_phi_e_face,
                                                          qface_point);
                                // Calculate the neighbor element shape functions at those locations
                                fe_lambda->reinit(interface_elem, &qface_point);


                                Kel.resize(n_phi_e_dofs, n_lambda_dofs);
                                Kle.resize(n_lambda_dofs, n_phi_e_dofs);

                                for (unsigned int qp=0; qp<qface_second.n_points(); qp++)
                                {
                                    for (unsigned int i=0; i<phi_e_face.size(); i++)
                                    {
                                        for (unsigned int j=0; j< phi_lambda.size(); j++)
                                        {
                                            Kel(i,j) -= JxW_e_face[qp] * phi_e_face[i][qp] * phi_lambda[j][qp];
                                        }
                                    }
                                    for (unsigned int i=0; i<phi_lambda.size(); i++)
                                    {
                                        for (unsigned int j=0; j< phi_e_face.size(); j++)
                                        {
                                            Kle(i,j) -= JxW_e_face[qp] * phi_e_face[j][qp] * phi_lambda[i][qp];
                                        }
                                    }
                                }
                                system.matrix->add_matrix(Kel, dof_indices_phi_e, dof_indices_lambda);
                                system.matrix->add_matrix(Kle, dof_indices_lambda, dof_indices_phi_e);
                            }
                        }

                    } // if on boundary
                } // loop on sides
                system.matrix->add_matrix(Kee, dof_indices_phi_e, dof_indices_phi_e);

                break;
            }
            // Intracellular space
            case 2:
            {
                Kii.resize(n_phi_i_dofs, n_phi_i_dofs);
                fe_phi_i->attach_quadrature_rule(&qrule_second);
                fe_phi_i->reinit(elem);
                // Laplacian
                for (unsigned int qp = 0; qp < qrule_second.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_i.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_i.size(); j++)
                        {
                            Kii(i, j) += JxW_i[qp] * sigma_i * dphi_i[i][qp] * dphi_i[j][qp];
                        }
                    }
                }
                system.matrix->add_matrix(Kii, dof_indices_phi_i, dof_indices_phi_i);

                for (unsigned int side=0; side<elem->n_sides(); side++)
                {
                    // Flux conditions at the cell membrane interface
                    if( mesh.boundary_info->has_boundary_id(elem, side, 3) )
                    {
                        auto it = MeshUtility::intracellular_elements_connectivity.find(elID);
                        if(it != MeshUtility::intracellular_elements_connectivity.end() )
                        {

                            Elem * interface_elem = it->second;
                            dof_map.dof_indices(interface_elem, dof_indices_lambda, lambda_var);
                            n_lambda_dofs = dof_indices_lambda.size();

                            // The quadrature point locations on the neighbor side
                            std::vector<Point> qface_point;

                            // Compute the shape function values on the element
                            // face.
                            fe_phi_i_face->attach_quadrature_rule(&qface_second);
                            fe_phi_i_face->reinit(elem, side);

                            FEInterface::inverse_map (elem->dim()-1,
                                                      fe_lambda->get_fe_type(),
                                                      interface_elem,
                                                      qps_phi_i_face,
                                                      qface_point);
                            // Calculate the neighbor element shape functions at those locations
                            fe_lambda->reinit(interface_elem, &qface_point);


                            Kil.resize(n_phi_i_dofs, n_lambda_dofs);
                            Kli.resize(n_lambda_dofs, n_phi_i_dofs);

                            for (unsigned int qp=0; qp<qface_second.n_points(); qp++)
                            {
                                for (unsigned int i=0; i<phi_i_face.size(); i++)
                                {
                                    for (unsigned int j=0; j< phi_lambda.size(); j++)
                                    {
                                        Kil(i,j) += JxW_i_face[qp] * phi_i_face[i][qp] * phi_lambda[j][qp];
                                    }
                                }
                                for (unsigned int i=0; i<phi_lambda.size(); i++)
                                {
                                    for (unsigned int j=0; j< phi_i_face.size(); j++)
                                    {
                                        Kli(i,j) += JxW_i_face[qp] * phi_i_face[j][qp] * phi_lambda[i][qp];
                                    }
                                }

                            }

                            system.matrix->add_matrix(Kil, dof_indices_phi_i, dof_indices_lambda);
                            system.matrix->add_matrix(Kli, dof_indices_lambda, dof_indices_phi_i);
                        }
                    }

                }

                break;
            }
            // cell membrane
            case 3:
            {
                std::cout << "Kmm: " << std::endl;
                Kmm.resize(n_phi_m_dofs, n_phi_m_dofs);
                Kll.resize(n_lambda_dofs, n_lambda_dofs);

                fe_phi_m->attach_quadrature_rule(&qface_second);
                fe_phi_m->reinit(elem);
                // Mass Matrix
                for (unsigned int qp = 0; qp < qface_second.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_m.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_m.size(); j++)
                        {
                            Kmm(i, j) += JxW_m[qp] * phi_m[i][qp] * phi_m[j][qp];
                        }
                    }
                }
                system.matrix->add_matrix(Kmm, dof_indices_phi_m, dof_indices_phi_m);
                // Stabilization Mass Matrix
                double tau = elem->hmax() * elem->hmax() / std::sqrt( sigma_i * sigma_e/ (sigma_i + sigma_e ) );
                // double tau =0* elem->hmax() / std::sqrt( sigma_i * sigma_e/ (sigma_i + sigma_e ) );
                for (unsigned int qp = 0; qp < qface_second.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_lambda.size(); i++)
                    {
                        for (unsigned int j = 0; j < phi_lambda.size(); j++)
                        {
                            Kll(i, j) += tau * JxW_lambda[qp] * dphi_lambda[i][qp] * dphi_lambda[j][qp];
                            //Kll(i, j) += tau * JxW_lambda[qp] * phi_lambda[i][qp] * phi_lambda[j][qp];
                        }
                    }
                }
                system.matrix->add_matrix(Kll, dof_indices_lambda, dof_indices_lambda);

                auto it = MeshUtility::interface_elements_connectivity.find(elID);
                if(it != MeshUtility::interface_elements_connectivity.end() )
                {
                    counter ++;
                    // Assemble Kim
                    // dt / Cm  * n_i cdot \sigma_i \nabla \phi_i
                    //  = - dt / Cm  * n_e cdot \sigma_i \nabla \phi_i
                    Elem * intracellular_elem = it->second.second;
                    Elem * extracellular_elem = it->second.first;
                    int side = extracellular_elem->which_neighbor_am_i(intracellular_elem);

                    dof_map.dof_indices(extracellular_elem, dof_indices_phi_e, phi_e_var);
                    dof_map.dof_indices(intracellular_elem, dof_indices_phi_i, phi_i_var);
                    dof_map.dof_indices(elem, dof_indices_phi_m, phi_m_var);
                    dof_map.dof_indices(elem, dof_indices_lambda, lambda_var);


                    n_phi_i_dofs = dof_indices_phi_i.size();
                    n_phi_e_dofs = dof_indices_phi_e.size();
                    Kmi.resize(n_phi_m_dofs, n_phi_i_dofs);

                    std::vector<Point> qface_point;

                    // We need only 1 quadrature point:
                    fe_phi_m->attach_quadrature_rule(&qface_second);
                    fe_phi_m->reinit(elem);

                    FEInterface::inverse_map (intracellular_elem->dim(),
                                              fe_phi_i->get_fe_type(),
                                              intracellular_elem,
                                              qps_phi_m,
                                              qface_point);

                    fe_phi_i_face->reinit(intracellular_elem, &qface_point);
                    // use this to get the normal working OK
                    fe_phi_e_face->reinit(extracellular_elem, side);

                    for (unsigned int qp = 0; qp < qface_second.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i < phi_m.size(); i++)
                        {
                            for (unsigned int j = 0; j < phi_i_face.size(); j++)
                            {
                                Kmi(i, j) -= JxW_m[qp] * dt / Cm * phi_m[i][qp] * ( sigma_i * dphi_i_face[j][qp] * qface_phi_e_normals[0]);
                            }
                        }
                    }

                    system.matrix->add_matrix(Kmi, dof_indices_phi_m, dof_indices_phi_i);
                    qface_point.clear();

                    // Assemble Lagrange Multiplier terms

                    // Klm = - mu * phi_m
                    fe_lambda->attach_quadrature_rule(&qface_second);
                    fe_lambda->reinit(elem);
                    fe_phi_m->attach_quadrature_rule(&qface_second);
                    fe_phi_m->reinit(elem);
                    Klm.resize(n_lambda_dofs, n_phi_m_dofs);
                    for (unsigned int qp = 0; qp < qface_second.n_points(); qp++)
                    {
                        // test
                        for (unsigned int i = 0; i < phi_lambda.size(); i++)
                        {
                            // trial
                            for (unsigned int j = 0; j < phi_m.size(); j++)
                            {
                                Klm(i, j) -= 0*JxW_lambda[qp] * phi_lambda[i][qp] * phi_m[j][qp];
                            }
                        }
                    }
                    std::cout << "add matrix: " << std::endl;
                    system.matrix->add_matrix(Klm, dof_indices_lambda, dof_indices_phi_m);

                }
                break;
            }
            default:
            {
                std::cout << "We should not be here!" << std::endl;
                throw std::runtime_error("");
                break;
            }

        }

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The SparseMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
    }
    std::cout <<" counter matrix: " << counter  << std::endl;

    // That concludes the system matrix assembly routine.
    system.matrix->close();
    std::cout <<"Matrix assembled"  << std::endl;
    // system.matrix->print();
    std::cout <<"Matrix assembled Are we sure?"  << std::endl;
}



// Now we define the assemble function which will be used
// by the EquationSystems object at each timestep to assemble
// the linear system for solution.
void assemble_rhs(EquationSystems & es)
{

    // Ignore unused parameter warnings when !LIBMESH_ENABLE_AMR.
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("MicroBidomain");
    TransientExplicitSystem & system_iion = es.get_system
    < TransientExplicitSystem > ("Iion");

    system.rhs->zero();


    unsigned int phi_e_var = 0;
    unsigned int phi_i_var = 1;
    unsigned int phi_m_var = 2;
    unsigned int lambda_var = 3;
    
    Real dt = es.parameters.get<Real>("dt");


    Real Cm = es.parameters.set < Real > ("Cm"); // uF mm^-2
    Real Chi = es.parameters.set < Real > ("Chi"); //  mm^-1
    Real sigma_i = es.parameters.set < Real > ("sigma_i"); // mS mm^-1
    Real sigma_e = es.parameters.set < Real > ("sigma_e"); // mS mm^-1

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_phi_e_type = system.variable_type(phi_e_var);
    FEType fe_phi_i_type = system.variable_type(phi_i_var);
    FEType fe_phi_m_type = system.variable_type(phi_m_var);
    
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr < FEBase > fe_phi_e(FEBase::build(dim, fe_phi_e_type));
    std::unique_ptr < FEBase > fe_phi_i(FEBase::build(dim, fe_phi_i_type));
    std::unique_ptr < FEBase > fe_phi_m(FEBase::build(dim-1, fe_phi_m_type));

    
    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    QGauss qrule_phi_e(dim,   libMesh::SECOND);
    QGauss qrule_phi_i(dim,   libMesh::SECOND);
    QGauss qrule_phi_m(dim-1, libMesh::SECOND);
    
    // Tell the finite element object to use our quadrature rule.
    fe_phi_e->attach_quadrature_rule(&qrule_phi_e);
    fe_phi_i->attach_quadrature_rule(&qrule_phi_i);
    fe_phi_m->attach_quadrature_rule(&qrule_phi_m);
    
    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.  We will start
    // with the element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW_e = fe_phi_e->get_JxW();
    const std::vector<Real> & JxW_i = fe_phi_i->get_JxW();
    const std::vector<Real> & JxW_m = fe_phi_m->get_JxW();
    
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> & phi_i = fe_phi_i->get_phi();
    const std::vector<std::vector<Real>> & phi_e = fe_phi_e->get_phi();
    const std::vector<std::vector<Real>> & phi_m = fe_phi_m->get_phi();
    
    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>> & dphi_i = fe_phi_i->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_e = fe_phi_e->get_dphi();
    const std::vector<std::vector<RealGradient>> & dphi_m = fe_phi_m->get_dphi();
    

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    const libMesh::DofMap & dof_map_iion = system_iion.get_dof_map();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseVector < Number > Fm;
    DenseVector < Number > Fl;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_phi_i;
    std::vector < dof_id_type > dof_indices_phi_e;
    std::vector < dof_id_type > dof_indices_phi_m;
    std::vector < dof_id_type > dof_indices_lambda;
    std::vector < libMesh::dof_id_type > dof_indices_iion;

    // Boundary elements
    // Boundary integration requires one quadraure rule,
    // with dimensionality one less than the dimensionality
    // of the element.

    // Now we will loop over all the elements in the mesh that
    // live on the local processor. We will compute the element
    // matrix and right-hand-side contribution.  Since the mesh
    // will be refined we want to only consider the ACTIVE elements,
    // hence we use a variant of the active_elem_iterator.
    int counter = 0;
    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        unsigned int blockID = elem->subdomain_id();
        unsigned int elID = elem->id();
        switch(blockID)
        {
            case 3:
            {
                counter ++;
                dof_map.dof_indices(elem, dof_indices_phi_m, phi_m_var);
                //dof_map.dof_indices(elem, dof_indices_lambda, lambda_var);
                dof_map_iion.dof_indices(elem, dof_indices_iion);
                unsigned int n_phi_m_dofs = dof_indices_phi_m.size();
                //unsigned int n_lambda_dofs = dof_indices_lambda.size();
                Fm.resize(n_phi_m_dofs);
                //Fl.resize(n_lambda_dofs);
                 fe_phi_m->reinit(elem);
                // Mass Matrix
                for (unsigned int qp = 0; qp < qrule_phi_m.n_points(); qp++)
                {
                    double phi_m_old = 0.;
                    double I_ion = 0;
                    for (int l = 0; l < phi_m.size(); l++)
                    {
                        phi_m_old += phi_m[l][qp] * (*system.old_local_solution)(dof_indices_phi_m[l]);
                        I_ion += phi_m[l][qp] * (*system_iion.current_local_solution)(dof_indices_iion[l]);
                    }
                    for (int i = 0; i < phi_m.size(); i++)
                    {
                        Fm(i) += JxW_m[qp] * phi_m[i][qp] * ( phi_m_old - dt / Cm * I_ion);
                  //      Fl(i) += JxW_m[qp] * phi_m[i][qp] * phi_m_old;
                    }
                }
                system.rhs->add_vector(Fm, dof_indices_phi_m);
                //system.rhs->add_vector(Fl, dof_indices_lambda);
                break;
            }
            default:
            {
                break;
            }
        }

    }
    system.rhs->close();
    std::cout <<" counter rhs: " << counter  << std::endl;
    //system.rhs->print();
    
    
    // That concludes the system matrix assembly routine.
}

void advance_ionic_model(EquationSystems & es)
{
    std::cout << "Advance Ionic Model: " << std::flush;
    Real dt = es.parameters.get<Real>("dt");

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("MicroBidomain");
    TransientExplicitSystem & system_iion = es.get_system
    < TransientExplicitSystem > ("Iion");

    system_iion.solution->zero();
    libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
    const libMesh::MeshBase::const_node_iterator end_node = mesh.local_nodes_end();

    const libMesh::DofMap & dof_map_phi_m = system.get_dof_map();
    const libMesh::DofMap & dof_map_iion = system_iion.get_dof_map();

    std::vector < libMesh::dof_id_type > dof_indices_phi_m;
    std::vector < libMesh::dof_id_type > dof_indices_iion;

    int counter= 0;
    for (; node != end_node; ++node)
    {
        const libMesh::Node * nn = *node;
        auto n_var = nn->n_vars(system_iion.number());
        auto n_dofs = nn->n_dofs(system_iion.number());
//        std::cout << "Doing stuff " << n_var << " " << n_dofs << std::endl;
//        nn->print_info();
        // If we are on the cell membrane n_var > 0
        if(n_var ==  n_dofs)
        {
            counter++;
//            std::cout << "Doing stuff " << std::endl;
//            nn->print_dof_info();
            dof_map_phi_m.dof_indices(nn, dof_indices_phi_m, 2);
            dof_map_iion.dof_indices(nn, dof_indices_iion);

//            std::cout << "Accessing phi_m" << std::endl;

            double phi_m_n = (*system.old_local_solution)(dof_indices_phi_m[0]); //V^n
//            std::cout << "Done" << std::endl;
            std::vector<double> values(IonicModel::num_vars, 0.0);
            values[0] = phi_m_n;
//            std::cout << "Accessing gatings: " << IonicModel::num_vars << ", " << values.size() << std::endl;
            for (int nv = 0; nv < IonicModel::num_vars-1; ++nv)
            {
//                std::cout << nv << " " << std::flush;
                values[nv+1] = system_iion.get_vector(nv)(dof_indices_iion[0]);
            }
//            std::cout << "\nUpdating Ionic Model" << std::endl;
            IonicModel::ionic_model->updateVariables(values, 0.0, dt);
//            std::cout << "Iion" << std::endl;
            double Iion = IonicModel::ionic_model->current_scaling() * IonicModel::ionic_model->evaluateIonicCurrent(values, 0.0, dt);
//            std::cout << "set Iion" << std::endl;
            system_iion.solution->set(dof_indices_iion[0], Iion);
//            std::cout << "set gatings" << std::endl;
            for (int nv = 0; nv < IonicModel::num_vars-1; ++nv)
            {
//                std::cout << nv << " " << std::flush;
                system_iion.get_vector(nv).set(dof_indices_iion[0], values[nv+1]);
            }
//            std::cout <<  std::endl;
        }
    }


    system_iion.solution->close();
    system_iion.update();
    std::cout << "Done! " << std::endl;

}

