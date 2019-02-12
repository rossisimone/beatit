// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
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
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "petsctime.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petsc/private/kspimpl.h"
#include "libmesh/getpot.h"


// This example will solve a linear transient system,
// so we need to include the TransientLinearImplicitSystem definition.
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side at each time step.  Note that
// since the system is linear we technically do not need to
// assemble the matrix at each time step, but we will anyway.
// In subsequent examples we will employ adaptive mesh refinement,
// and with a changing mesh it will be necessary to rebuild the
// system matrix.
void assemble_matrix(EquationSystems & es, const std::string & system_name);
void assemble_rhs(EquationSystems & es, const std::string & system_name);

// Function prototype.  This function will initialize the system.
// Initialization functions are optional for systems.  They allow
// you to specify the initial values of the solution.  If an
// initialization function is not provided then the default (0)
// solution is provided.
void init_cd(EquationSystems & es, const std::string & system_name);

// Exact solution function prototype.  This gives the exact
// solution as a function of space and time.  In this case the
// initial condition will be taken as the exact solution at time 0,
// as will the Dirichlet boundary conditions at time t.




Number exact_value(const Point & p, const Parameters & parameters,
                   const std::string & /*system_name*/, const std::string & variable_name )
{
    //double time = parameters.get < Real > ("time");
    double x = p(0);
    double y = p(1);
    double z = p(2);
    double x1=pow(x,2);
    double y1=pow(y,2);
    
    
    if( "u" == variable_name)
    {
        
        
        //      if(y < -0.5)return 1.0;
        //       else return 0.0;
        if(x1+y1 < 1) return 1.0;
        else return 0.0;
    }
    else if( "v" == variable_name)
    {
        return 0.0;
        //return 1.0;
    }
    return 0.0;
    
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
    Mesh mesh(init.comm());
    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data.beat", 2, "-i", "--input");
    GetPot data(datafile_name);
    
    
    std::string meshname = data("input_mesh_name", "NONE");
    
    if(meshname == "NONE") throw std::runtime_error("Specify Mesh Name in input file!");
    mesh.read(meshname);

    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    // Create an equation systems object.
    EquationSystems equation_systems(mesh);
    
    // Add a transient system to the EquationSystems
    // object named "Convection-Diffusion".
    TransientLinearImplicitSystem & system = equation_systems.add_system
    < TransientLinearImplicitSystem > ("Convection-Diffusion");
    
    // Adds the variable "u" to "Convection-Diffusion".  "u"
    // will be approximated using first-order approximation.
    int tissueid = data("tissueid",1);
    int bathid = data("bathid",2);
    equation_systems.parameters.set <unsigned int > ("tissueid") = tissueid;
    equation_systems.parameters.set <unsigned int > ("bathid") = bathid;
    std::set<unsigned short> active_subdomains;
    active_subdomains.insert(tissueid);
    system.add_variable("u", FIRST, LAGRANGE, &active_subdomains);
    system.add_variable("v", FIRST);
    
    system.add_vector("nullspace");
    
    // Give the system a pointer to the matrix assembly
    // and initialization functions.
    //system.attach_assemble_function(assemble_cd);
    system.attach_init_function(init_cd);
    
    // Initialize the data structures for the equation system.
    equation_systems.init();
    
    // Prints information about the system to the screen.
    equation_systems.print_info();
    
    double Cm = data("Cm",0.01);
    double Chi = data("Chi",140);
    double sigma_i_l = data("sigma_i_l",0.17);
    double sigma_i_t = data("sigma_i_t",0.019);
    double sigma_e_l = data("sigma_e_l",0.62);
    double sigma_e_t = data("sigma_e_t",0.24);
    double sigma_b_l = data("sigma_b_l",2.0);
    equation_systems.parameters.set < Real > ("Cm") = Cm; // uF mm^-2
    equation_systems.parameters.set < Real > ("Chi") = Chi; //  mm^-1
    equation_systems.parameters.set < Real > ("sigma_i_l") = sigma_i_l; // mS mm^-1
    equation_systems.parameters.set < Real > ("sigma_i_t") = sigma_i_t; // mS mm^-1
    equation_systems.parameters.set < Real > ("sigma_e_l") = sigma_e_l; // mS mm^-1
    equation_systems.parameters.set < Real > ("sigma_e_t") = sigma_e_t; // mS mm^-1
    equation_systems.parameters.set < Real > ("sigma_b_l") = sigma_b_l;
    
    
    
    // The Convection-Diffusion system requires that we specify
    // the flow velocity.  We will specify it as a RealVectorValue
    // data type and then use the Parameters object to pass it to
    // the assemble function.
    equation_systems.parameters.set < RealVectorValue > ("velocity") =
    RealVectorValue(0.8, 0.8);
    
    // Solve the system "Convection-Diffusion".  This will be done by
    // looping over the specified time interval and calling the
    // solve() member at each time step.  This will assemble the
    // system and call the linear solver.
    double dt = data("dt",0.0125);
    equation_systems.parameters.set < Real > ("dt") = dt;
    
    system.time = 0.;
    
    // Write out the initial conditions.
    // If Exodus is available, we'll write all timesteps to the same file
    // rather than one file per timestep.
    std::string exodus_filename = "bath.e";
    ExodusII_IO exo(mesh);
    exo.write_equation_systems(exodus_filename, equation_systems);
    exo.append(true);
    //exo.write_timestep(exodus_filename, equation_systems, 0, system.time);
    int save_it = 0;
    
    
    
    assemble_matrix(equation_systems, "Convection-Diffusion");
    
    system.get_linear_solver()->init();
    
    libMesh::LinearSolver<Number>* my_solver = system.get_linear_solver();
    
    libMesh::PetscLinearSolver<Number> petsc_solver(init.comm());
    
    { // set BLOCKS
        
        IS is_u_local;
        IS is_u_global;
        IS is_v_local;
        const DofMap & dof_map = system.get_dof_map();
        
        std::vector<libMesh::dof_id_type> u_indices;
        dof_map.local_variable_indices(u_indices, mesh, 0);
        
        ISCreateGeneral(PETSC_COMM_SELF, u_indices.size(), reinterpret_cast<int*>(&u_indices[0]),PETSC_COPY_VALUES,&is_u_local);
        
        int nmin =  dof_map.first_dof();
        int nmax = dof_map.end_dof();
        
        ISComplement(is_u_local, nmin, nmax, &is_v_local);
        
        typedef libMesh::PetscMatrix<libMesh::Number> PetscMatrix;
        petsc_solver.init(dynamic_cast<PetscMatrix *>(system.matrix));
        
        PCFieldSplitSetIS(petsc_solver.ksp()->pc,"u",is_u_local);
        PCFieldSplitSetIS(petsc_solver.ksp()->pc,"v",is_v_local);
    }
    
    int n_steps = data("n_steps",100);
    for (unsigned int t_step = 0; t_step < n_steps; t_step++) {
        // Increment the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        system.time += dt;
        
        equation_systems.parameters.set < Real > ("time") = system.time;
        
        // A pretty update message
        libMesh::out << " Solving time step ";
        
        // Do fancy zero-padded formatting of the current time.
        {
            std::ostringstream out;
            
            out << std::setw(2) << std::right << t_step << ", time="
            << std::fixed << std::setw(6) << std::setprecision(3)
            << std::setfill('0') << std::left << system.time << "...";
            
            libMesh::out << out.str() << std::endl;
        }
        
        // At this point we need to update the old
        // solution vector.  The old solution vector
        // will be the current solution vector from the
        // previous time step.  We will do this by extracting the
        // system from the EquationSystems object and using
        // vector assignment.  Since only TransientSystems
        // (and systems derived from them) contain old solutions
        // we need to specify the system type when we ask for it.
        system.update();
        *system.old_local_solution = *system.current_local_solution;
        
        // Assemble & solve the linear system
        //equation_systems.get_system("Convection-Diffusion").solve();
        
        // Assemble the RHS
        assemble_rhs(equation_systems, "Convection-Diffusion");
        
        // TODO:measure time it takes to solve the linear system
        // Solve Linear System
        int maxits = 2000000;
        double tol=1e-12;
        
        std::pair<unsigned int, double> rval = std::make_pair(0,0.0);
        //system.rhs->print(std::cout);
        // system.solution->zero();
        //TODO ADD TIMESTART
        double v1,v2,solvetime;
        //PetscErrorCode ierr;
        //ierr =
        PetscTime(&v1);
        //CHKERRQ(ierr);
        
        
        //        rval = my_solver->solve (*system.matrix, system.request_matrix("Preconditioner"), *system.solution, *system.rhs, tol, maxits);
        rval = petsc_solver.solve (*system.matrix, *system.solution, *system.rhs, tol, maxits);
        //TODO ADD TIME END
        //ierr =
        PetscTime(&v2);
        //CHKERRQ(ierr);
        
        //TODO ELAPSED TIME = TIME END -TIME Start
        solvetime = v2 - v1;
        
        std::cout << "time for solving: " << solvetime << " s. Iterations: " << rval.first << " Residual: " << rval.second <<std::endl;
        
        int n_linear_iter   = rval.first;
        double final_linear_residual = rval.second;
        // Output every 10 timesteps to file.
        if ((t_step + 1) % 1 == 0) {
            save_it++;
            std::cout << "Exporting the solution at step " << t_step
            << std::endl;
            exo.write_timestep(exodus_filename, equation_systems, save_it,
                               system.time);
            
        }
    }
    
    
    // All done.
    return 0;
}




// We now define the function which provides the
// initialization routines for the "Convection-Diffusion"
// system.  This handles things like setting initial
// conditions and boundary conditions.
void init_cd(EquationSystems & es,
             const std::string & libmesh_dbg_var(system_name)) {
    // It is a good idea to make sure we are initializing
    // the proper system.
    libmesh_assert_equal_to(system_name, "Convection-Diffusion");
    
    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("Convection-Diffusion");
    
    // Project initial conditions at time 0
    es.parameters.set < Real > ("time") = system.time = 0;
    
    system.project_solution(exact_value, libmesh_nullptr, es.parameters);
}

// Now we define the assemble function which will be used
// by the EquationSystems object at each timestep to assemble
// the linear system for solution.
void assemble_matrix(EquationSystems & es, const std::string & system_name) {
    // Ignore unused parameter warnings when !LIBMESH_ENABLE_AMR.
    libmesh_ignore(es);
    libmesh_ignore(system_name);
    
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Convection-Diffusion");
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("Convection-Diffusion");
    
    
    unsigned int vn_u = 0;
    unsigned int vn_v = 1;
    
    
    //Real time = es.parameters.get<Real>("time");
    Real dt = es.parameters.get<Real>("dt");
    
    
    Real Cm = es.parameters.set < Real > ("Cm"); // uF mm^-2
    Real Chi = es.parameters.set < Real > ("Chi"); //  mm^-1
    Real sigma_i_l = es.parameters.set < Real > ("sigma_i_l"); // mS mm^-1
    Real sigma_i_t = es.parameters.set < Real > ("sigma_i_t"); // mS mm^-1
    Real sigma_e_l = es.parameters.set < Real > ("sigma_e_l"); // mS mm^-1
    Real sigma_e_t = es.parameters.set < Real > ("sigma_e_t"); // mS mm^-1
    unsigned int tissueid=es.parameters.set <unsigned int> ("tissueid");
    unsigned int bathid=es.parameters.set <unsigned int> ("bathid");
    Real sigma_b_l = es.parameters.set < Real > ("sigma_b_l");
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = system.variable_type(0);
    
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr < FEBase > fe(FEBase::build(dim, fe_type));
    
    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    QGauss qrule(dim, fe_type.default_quadrature_order());
    
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);
    
    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.  We will start
    // with the element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();
    
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    
    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
    
    // The XY locations of the quadrature points used for face integration
    
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseMatrix < Number > Ke;
    DenseSubMatrix<Number> ma11(Ke),ma12(Ke),ma21(Ke),ma22(Ke);
    
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_u;
    std::vector < dof_id_type > dof_indices_v;
    
    // Here we extract the velocity & parameters that we put in the
    // EquationSystems object.
    const RealVectorValue velocity = es.parameters.get < RealVectorValue
    > ("velocity");
    
    TensorValue<Number> sigma_i;
    sigma_i(0,0) = sigma_i_l;
    sigma_i(1,1) = sigma_i_t;//0.1
    
    TensorValue<Number> sigma_e;
    sigma_e(0,0) = sigma_e_l;
    sigma_e(1,1) = sigma_e_t;//0.1
    
    TensorValue<Number> sigma_b;
    sigma_b(0,0) = sigma_b_l;
    sigma_b(1,1) = sigma_b_l;//0.1
    
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
        unsigned int elID = elem->subdomain_id();
        
        
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_u, vn_u);
        dof_map.dof_indices(elem, dof_indices_v, vn_v);
        
        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        
        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);
        
        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize(n_dofs, n_dofs);
        ma11.reposition(vn_u*n_u_dofs, vn_u*n_u_dofs, n_u_dofs, n_u_dofs);
        
        ma12.reposition(vn_u*n_u_dofs, vn_v*n_u_dofs, n_u_dofs, n_v_dofs);
        
        ma21.reposition(vn_v*n_u_dofs, vn_u*n_u_dofs, n_v_dofs, n_u_dofs);
        
        ma22.reposition(vn_v*n_u_dofs, vn_v*n_u_dofs, n_v_dofs, n_v_dofs);
        
        
        // Now we will build the element matrix and right-hand-side.
        // Constructing the RHS requires the solution and its
        // gradient from the previous timestep.  This myst be
        // calculated at each quadrature point by summing the
        // solution degree-of-freedom values by the appropriate
        // weight functions.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
            // Values to hold the old solution & its gradient.
            
            // if(time<1)
            // {
            //   if(x1+y1<1) I_stim= 10.0;
            //  else    I_stim=0.0;
            // }
            
            if(elID==tissueid){
                
                for (std::size_t i = 0; i < phi.size(); i++) {
                    for (std::size_t j = 0; j < phi.size(); j++) {
                        // The matrix contribution
                        ma11(i, j) +=
                        JxW[qp]
                        * ( (1 / dt) * phi[j][qp] * phi[i][qp] +
                           (1.0 / (Chi*Cm)) * ( dphi[i][qp] * ( sigma_i * dphi[j][qp] ) )
                           );
                        
                        ma12(i,j)+=JxW[qp]*((1.0 / (Chi*Cm))
                                            * (dphi[i][qp] * ( sigma_i * dphi[j][qp] ) ));
                        ma21(i,j)+=JxW[qp]*((1.0 / (Chi*Cm))
                                            * (dphi[i][qp] * ( sigma_i * dphi[j][qp] ) ));
                        ma22(i,j)+=JxW[qp] * (
                                              (1.0 / (Chi*Cm) )
                                              * (dphi[j][qp] * ( ( sigma_i  + sigma_e ) * dphi[i][qp] ) ) );
                    }
                }
            }
            else
            {
                for (std::size_t i = 0; i < phi.size(); i++) {
                    for (std::size_t j = 0; j < phi.size(); j++) {
                        ma22(i,j)+=JxW[qp]*( (1.0 / (Chi*Cm)) * (dphi[i][qp] * ( sigma_b * dphi[j][qp] ) ));
                    }
                }
            }
            
        }
        
        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The SparseMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
    }
    // That concludes the system matrix assembly routine.
    system.matrix->close();
    std::cout <<"Matrix assembled"  << std::endl;
    
    {
        // CREATE NULLSPACE
        typedef libMesh::PetscMatrix<libMesh::Number> Mat;
        typedef libMesh::PetscVector<libMesh::Number> Vec;
        libMesh::DenseVector<libMesh::Number> NSe;
        
        {
            
            libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
            const libMesh::MeshBase::const_node_iterator end_node = mesh.local_nodes_end();
            
            for (; node != end_node; ++node)
            {
                
                const libMesh::Node * nn = *node;
                
                dof_map.dof_indices(nn, dof_indices_u, vn_u);
                dof_map.dof_indices(nn, dof_indices_v, vn_v);
                
                if(dof_indices_u.size() > 0) system.get_vector("nullspace").set(dof_indices_u[0], 0.0);
                system.get_vector("nullspace").set(dof_indices_v[0], 1.0);
                
            }
            
            system.get_vector("nullspace").close();
            // Normalization
            system.get_vector("nullspace") /= system.get_vector("nullspace").l2_norm();
            
        }
        
        Vec& N = (static_cast<Vec&>(system.get_vector("nullspace")));
        auto vec = N.vec();
        // Create NullSPace
        MatNullSpace nullspace;
        MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, &vec, &nullspace);
        // Set the NullSpace
        Mat * mat = dynamic_cast<Mat *>(system.matrix);
        MatSetNullSpace(mat->mat(), nullspace);
        // Destroy NullSpace
        MatNullSpaceDestroy(&nullspace);
        
    }
    
    
    
    
}


// Now we define the assemble function which will be used
// by the EquationSystems object at each timestep to assemble
// the linear system for solution.
void assemble_rhs(EquationSystems & es, const std::string & system_name) {
    // Ignore unused parameter warnings when !LIBMESH_ENABLE_AMR.
    libmesh_ignore(es);
    libmesh_ignore(system_name);
    
    std::cout <<" assembling RHS"  << std::endl;
    
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Convection-Diffusion");
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = es.get_mesh();
    
    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();
    
    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & system = es.get_system
    < TransientLinearImplicitSystem > ("Convection-Diffusion");
    system.rhs->zero();
    unsigned int vn_u = 0;
    unsigned int vn_v = 1;
    
    //Real time = es.parameters.get<Real>("time");
    Real dt = es.parameters.get<Real>("dt");
    Real Cm = es.parameters.set < Real > ("Cm"); // uF mm^-2
    Real Chi = es.parameters.set < Real > ("Chi"); //  mm^-1
    Real sigma_i_l = es.parameters.set < Real > ("sigma_i_l"); // mS mm^-1
    Real sigma_i_t = es.parameters.set < Real > ("sigma_i_t"); // mS mm^-1
    Real sigma_e_l = es.parameters.set < Real > ("sigma_e_l"); // mS mm^-1
    Real sigma_e_t = es.parameters.set < Real > ("sigma_e_t"); // mS mm^-1
    
    unsigned int tissueid=es.parameters.set <unsigned int> ("tissueid");
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = system.variable_type(0);
    
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr < FEBase > fe(FEBase::build(dim, fe_type));
    
    // A Gauss quadrature rule for numerical integration.
    // Let the FEType object decide what order rule is appropriate.
    QGauss qrule(dim, fe_type.default_quadrature_order());
    
    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);
    
    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.  We will start
    // with the element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();
    
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<Point> & qp_xyz = fe->get_xyz();
    
    // The element shape function gradients evaluated at the quadrature
    // points.
    //const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
    
    
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap & dof_map = system.get_dof_map();
    
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".
    DenseVector < Number > Fe;
    DenseSubVector < Number > vec1(Fe),vec2(Fe);
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > dof_indices_u;
    std::vector < dof_id_type > dof_indices_v;
    
    // Here we extract the velocity & parameters that we put in the
    // EquationSystems object.
    const RealVectorValue velocity = es.parameters.get < RealVectorValue
    > ("velocity");
    
    //  TensorValue<Number> sigma;
    // sigma(0,0) = 1.0;
    // sigma(1,1) = 0.1;
    
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
        unsigned int elID = elem->subdomain_id();
        if(elID==tissueid)
        {
            
            
            dof_map.dof_indices(elem, dof_indices);
            dof_map.dof_indices(elem, dof_indices_u, vn_u);
            dof_map.dof_indices(elem, dof_indices_v, vn_v);
            const unsigned int n_dofs = dof_indices.size();
            const unsigned int n_u_dofs = dof_indices_u.size();
            const unsigned int n_v_dofs = dof_indices_v.size();
            // Compute the element-specific data for the current
            // element.  This involves computing the location of the
            // quadrature points (q_point) and the shape functions
            // (phi, dphi) for the current element.
            fe->reinit(elem);
            
            // Zero the element matrix and right-hand side before
            // summing them.  We use the resize member here because
            // the number of degrees of freedom might have changed from
            // the last element.  Note that this will be the case if the
            // element type is different (i.e. the last element was a
            // triangle, now we are on a quadrilateral).
            Fe.resize(dof_indices.size());
            vec1.reposition(vn_u * n_u_dofs, n_u_dofs);
            vec2.reposition(vn_v * n_v_dofs, n_v_dofs);
            
            
            // Now we will build the element matrix and right-hand-side.
            // Constructing the RHS requires the solution and its
            // gradient from the previous timestep.  This myst be
            // calculated at each quadrature point by summing the
            // solution degree-of-freedom values by the appropriate
            // weight functions.
            for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
                // Values to hold the old solution & its gradient.
                
                Point p = qp_xyz[qp];
                double x = p(0);
                double y = p(1);
                double z = p(2);
                
                Number u_old = 0.;
                Gradient grad_u_old;
                Number I_ion = 0.;
                
                // double x1=pow(x,2);
                // double y1=pow(y,2);
                // double I_stim;
                // Number Q_t = 0.;
                // if(time<1)
                // {
                //   if(x1+y1<1) I_stim= 10.0;
                //  else    I_stim=0.0;
                // }
                // Compute the old solution & its gradient.

                for (std::size_t l = 0; l < phi.size(); l++) {
                    u_old += phi[l][qp] * (*system.old_local_solution)(dof_indices_u[l]);
                    double v_old = (*system.old_local_solution)(dof_indices_u[l]);
                    double I_ion_l = v_old * (v_old - 1) * (v_old - 0.1);
                    I_ion += phi[l][qp] * I_ion_l;
                    // This will work,
                    // grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
                    // but we can do it without creating a temporary like this:
                    //grad_u_old.add_scaled (dphi[l][qp], system.old_solution (dof_indices[l]));
                }
                
                //I_ion = u_old * (u_old - 1) * (u_old - 0.1);
                 double rhs = 1 / dt * u_old - I_ion/Cm;
                
                // Now compute the element matrix and RHS contributions.
                //   double rhs = 1 / dt * u_old - I_ion  + I_stim;
                //vec1.print(std::cout);
                for (std::size_t i = 0; i < phi.size(); i++) {
                    // The RHS contribution
                    vec1(i) += JxW[qp] * rhs * phi[i][qp];
                    vec2(i) = 0.0;
                }
                
            }
            
            // The element matrix and right-hand-side are now built
            // for this element.  Add them to the global matrix and
            // right-hand-side vector.  The SparseMatrix::add_matrix()
            // and NumericVector::add_vector() members do this for us.
            system.rhs->add_vector(Fe, dof_indices);
        }
    }
    system.rhs->close();
    std::cout <<" done"  << std::endl;
    
    
    // That concludes the system matrix assembly routine.
}
