// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"

#include "libmesh/quadrature_gauss.h"
#include "libmesh/quadrature_trap.h"
#include "libmesh/quadrature_simpson.h"
#include "libmesh/quadrature_grid.h"

#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/patch.h"

#include "libmesh/point_locator_tree.h"
#include "libmesh/tree.h"
#include "libmesh/tree_base.h"
#include "libmesh/tree_node.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/fe_interface.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/exact_solution.h"
#include "fstream"
#include "libmesh/distributed_vector.h"
#include <algorithm>
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/eigen_sparse_matrix.h"
#include <cstring>
#include "libmesh/patch.h"
#include "libmesh/timestamp.h"
#include <time.h>
#include <sys/time.h>


#define BOUNDARY_ID_MIN_Y 0
#define BOUNDARY_ID_MAX_X 1
#define BOUNDARY_ID_MAX_Y 2
#define BOUNDARY_ID_MIN_X 3

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Matrix and right-hand side assemble
void preassemble_elasticity(EquationSystems & es,
                    const std::string & libmesh_dbg_var(system_name));
void compute_mass(EquationSystems & es,
                    const std::string & libmesh_dbg_var(system_name));
void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name);

RealTensor kernel2( Point xq , Point xq_nl);
// Define the elasticity tensor, which is a fourth-order tensor
// i.e. it has four indices i, j, k, l
Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l);

void compute_nonlocal_info(EquationSystems & es,
                       const std::string & libmesh_dbg_var(system_name));

Real PD_stretch( Point xi , Point eta);

Point PD_BB_force( Point x_deformed, Real stretch, const std::vector<Real> & PD_params);


//add to separate function file?
Real norm_difference( UniquePtr<NumericVector<Real>> &v1, UniquePtr<NumericVector<Real>> &v2, int pnorm)
{
  UniquePtr<NumericVector<Real>> v3 = v1->clone();
  v3->add( -1.0, *v2);
  if(pnorm==1)
    return v3->l1_norm();
  else if(pnorm==2)
    return v3->l2_norm();
  else if(pnorm>2)
    return v3->linfty_norm();
}

Real ratio_norm_difference( UniquePtr<NumericVector<Real>> &v1, UniquePtr<NumericVector<Real>> &v2, UniquePtr<NumericVector<Real>> &v3, int pnorm)
{
  return norm_difference( v1, v2, pnorm) / norm_difference( v2, v3, pnorm);
}

void get_simname(std::string & simname, int argc, char** argv)
{
  char str1[2] = "_";
  simname = std::strcat(argv[1],str1);
  for(int kk=2;kk< argc; kk++)
  {
    simname += argv[kk];
    simname += "_";
  }
    simname+=".e";
}

double get_wall_time(){
struct timeval time;
if (gettimeofday(&time,NULL)){
//  Handle error
return 0;
}
return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
return (double)clock() / CLOCKS_PER_SEC;
}

//GLOBAL
Real poisson_ratio = 0.3333;
Real bulk_modulus = 1.0;
Real young_modulus = 1.0;
int PD_model_type = 0;
Real applied_force_x = 1e2;
Real applied_force_y = 0;
int nelemx = 32;
Real dt = 1.0;
Real PD_horizon_factor = 2.0;
Real PD_horizon = 1.0;
Real PD_horizon2 = 0.;
Real PD_horizon3 = 0.;
Real PD_horizon4 = 0.;
Real PD_horizon_sq = 1.0;
Real force_increment = 1.1;
const Real convergence_tol = 1e-6;
Real applied_damping = 1e2;
unsigned int target_patch_size = 26;
const Real max_extension = 1.0e-4;
Real lambda_1=0., lambda_2=0., poisson_ratio_factor=0.;
Real bc_width_right = 0.05;
Real bc_width_left = 0.05;
Real bc_width_top = 0.05;
Real bc_width_bottom = 0.05;
const Real set_penalty = 1;
const bool flag_implicit_damping = 1;

/*
0: PD
1: CG
*/
const int use_method = 1;

Real force_attenuated(Real f, Real xcoord)
{
  Real xmin = 0.5 * bc_width_right;
  Real xmax = 0.5;
  return f * (xcoord-xmin) / (xmax-xmin);
  // return f;

}


// Begin the main program.
int main (int argc, char ** argv)
{

  // Initialize libMesh and any dependent libraries
  LibMeshInit init (argc, argv);

  unsigned int iarg = 1;
  for(int ii=0;ii<argc;ii++)
    std::cout << argv[ii] << std::endl;

  Real bulk_modulus =  static_cast<Real>(std::atof(argv[iarg]));iarg+=1;
  Real poisson_ratio =  static_cast<Real>(std::atof(argv[iarg]));iarg+=1;
  Real applied_damping = static_cast<int>(std::atoi(argv[iarg]));iarg+=1;
  Real applied_force_x = static_cast<Real>(std::atof(argv[iarg]));iarg+=1;
  Real applied_force_y = static_cast<Real>(std::atof(argv[iarg]));iarg+=1;
  unsigned int nelem = static_cast<int>(std::atoi(argv[iarg]));iarg+=1;
  Real dt = static_cast<Real>(std::atof(argv[iarg]));iarg+=1;
  Real PD_horizon_factor = static_cast<Real>(std::atof(argv[iarg]));iarg+=1;
  Real force_increment = static_cast<Real>(std::atof(argv[iarg]));iarg+=1;
  const FEFamily fe_family = static_cast<FEFamily>(std::atoi(argv[iarg]));iarg+=1;
  const Order fe_order = static_cast<Order>(std::atoi(argv[iarg]));iarg+=1;
  const ElemType elem_type = static_cast<ElemType>(std::atoi(argv[iarg]));iarg+=1;
  Real applied_force_ramp_up_duration = static_cast<Real>(std::atof(argv[iarg]));iarg+=1;

  //reset params
  poisson_ratio_factor = (3 * poisson_ratio - 1) / (1 - poisson_ratio);
  lambda_1 = 3 * bulk_modulus * poisson_ratio / ((1. + poisson_ratio) ); //lambda
  lambda_2 = 3 * bulk_modulus * (1. - 2. * poisson_ratio) / ((1. + poisson_ratio) ); //G
  young_modulus = 3 * bulk_modulus * ( 1 - 2 * poisson_ratio );
  std::cout << "poisson_ratio_factor " << poisson_ratio_factor << "Lame 1" << lambda_1 << "Lame 2" << lambda_2 << "young_modulus" << young_modulus << std::endl;
  std::string simname(80,' ');
  get_simname(simname, argc, argv);
  std::cout << simname << std::endl;


  std::cout << simname << std::endl;

  const unsigned int nt = static_cast<unsigned int>(1e6);
  const unsigned int inc_write = 100;

  const unsigned int dim = 2;


  ReplicatedMesh mesh(init.comm(), dim);

  // Real Fx_applied_max = applied_force_x;
  // Real Fy_applied_max = applied_force_y; 
  Point applied_force = Point(applied_force_x, applied_force_y, 0.0);
  applied_force *= nelem;
  applied_force *= dt / applied_force_ramp_up_duration;

  Real elem_length_avg = 1.0 / ((double) nelem ) * 1;
  if(PD_horizon_factor > 0){
    if(elem_type==QUAD4 || elem_type==QUAD8 || elem_type==QUAD9)
    {
      PD_horizon = 1.0 / ((double) nelem ) * PD_horizon_factor;
      target_patch_size = (2.0 * ceil(PD_horizon_factor) + 1.0); 
      target_patch_size = ceil( target_patch_size * target_patch_size);      
    }
    if(elem_type==TRI3 || elem_type==TRI6)
    {
      PD_horizon = 1.0 / ((double) nelem ) * PD_horizon_factor;
      target_patch_size = (2.0 * ceil(PD_horizon_factor) + 1.0); 
      target_patch_size = ceil( target_patch_size * target_patch_size);  
    }
    }
  else
  {
    PD_horizon = fabs(PD_horizon_factor);
    PD_horizon_sq = PD_horizon * PD_horizon;
    target_patch_size = (2.0 * ceil(PD_horizon * ((double) nelem) ) + 1.0); 
    target_patch_size = ceil( target_patch_size * target_patch_size);    
  }
  PD_horizon2 = PD_horizon * PD_horizon;
  PD_horizon3 = PD_horizon2 * PD_horizon;
  PD_horizon4 = PD_horizon2 * PD_horizon2;

  std::cout << "PD_horizon " << PD_horizon << std::endl;
  std::cout << "target_patch_size " << target_patch_size << std::endl;
  std::cout << "  "  << std::endl;


  //convergence

  const Real tol_ratio_1norm_increment_displacement = convergence_tol, tol_ratio_2norm_increment_displacement = convergence_tol, tol_ratio_infnorm_increment_displacement = convergence_tol;
  const Real tol_ratio_1norm_increment_velocity = convergence_tol, tol_ratio_2norm_increment_velocity = convergence_tol, tol_ratio_infnorm_increment_velocity = convergence_tol;


  //domain length
  const Real L = 0.5;
  //bounding boxes
  const Point &ptmin = Point(-0.5,-0.5,0.0);
  const Point &ptmax = Point(0.5,0.5,0.0);
  const Point &ptmin_physical = Point(-0.5 + PD_horizon, -0.5 + PD_horizon, 0.0);
  const Point &ptmax_physical = Point(0.5 - PD_horizon , 0.5 - PD_horizon , 0.0);
  BoundingBox  bbox(ptmin,ptmax);
  BoundingBox  bbox_physical(ptmin_physical, ptmax_physical);
  BoundingBox  bbox0( Point( -L, -L, 0.) , Point( L, -L+PD_horizon, 0.) ); //bottom
  BoundingBox  bbox1( Point( L-PD_horizon, -L, 0.) , Point( L, L, 0.) ); //right
  BoundingBox  bbox2( Point( -L, L-PD_horizon, 0.) , Point( L, L, 0.) ); //top
  BoundingBox  bbox3( Point( -L, -L, 0.) , Point( -L+PD_horizon, L, 0.) ); //left

  MeshTools::Generation::build_square (mesh, nelem, nelem, ptmin(0), ptmax(0), ptmin(1), ptmax(1), elem_type);

  // Print information about the mesh to the screen.
  mesh.print_info();

  bc_width_right = elem_length_avg;
  bc_width_top = elem_length_avg;
  bc_width_bottom = elem_length_avg;
  bc_width_left = elem_length_avg;
  unsigned int nelem_bc_right=0, nnode_bc_right=0;
  MeshBase::element_iterator       elem_it  = mesh.elements_begin();
  const MeshBase::element_iterator elem_end = mesh.elements_end();
  for (; elem_it != elem_end; ++elem_it)
    {
      Elem * elem = *elem_it;
      if ( bbox.contains_point( elem->centroid() ) )
          elem->subdomain_id() = 5; //the physical domain
      if ( elem->centroid()(1) < -(L-bc_width_bottom) )
          elem->subdomain_id() = 2; //bottom volume boundary
      if ( elem->centroid()(1) > (L-bc_width_top) )
          elem->subdomain_id() = 4; //top volume boundary
      if ( elem->centroid()(0) < -(L-bc_width_left) )
          elem->subdomain_id() = 3; //left volume boundary
      if ( elem->centroid()(0) > (L - bc_width_right) )
      {
          elem->subdomain_id() = 1; //right volume boundary
          nelem_bc_right +=1;
      }
      // if ( elem->centroid()(1) < -(L-elem_length_avg))
    }
  



  // Create an equation systems object.
  EquationSystems equation_systems (mesh);



  // Declare the system and its variables.
  // Create a system named "Elasticity" governing displacements and "dElasticity" governing velocities
  TransientLinearImplicitSystem & system = equation_systems.add_system<TransientLinearImplicitSystem> ("Elasticity");
  unsigned int u_var = system.add_variable("displacementsx", fe_order, fe_family);
  unsigned int v_var = system.add_variable("displacementsy", fe_order, fe_family);
  TransientLinearImplicitSystem & dsystem = equation_systems.add_system<TransientLinearImplicitSystem> ("dElasticity");
  unsigned int du_var = dsystem.add_variable("velocitiesx", fe_order, fe_family);
  unsigned int dv_var = dsystem.add_variable("velocitiesy", fe_order, fe_family);
  ExplicitSystem & fsystem = equation_systems.add_system<ExplicitSystem> ("fElasticity");
  unsigned int fu_var = fsystem.add_variable("forcesx", fe_order, fe_family);
  unsigned int fv_var = fsystem.add_variable("forcesy", fe_order, fe_family);



  /*
  PD parameters
  */
  unsigned int nelem2 = mesh.n_elem();
  // PD_horizon = 1.0 / ( sqrt((double) nelem2) ) * 2.0;
  std::cout<< "PD_horizon " << PD_horizon << std::endl;
  std::vector<Real> PD_params(2,0.);
  PD_params[0] = 18.0 * bulk_modulus / (5. * PD_horizon*PD_horizon);





  // Assembly and Initialization for displacement, velocity systems
  dsystem.attach_assemble_function (assemble_elasticity);
  // dsystem.attach_init_function (init_elasticity);

  Real density = 1.;
  Real shear_wave_speed = sqrt( PD_horizon /  (lambda_1 + 2* lambda_2 / density) );
  Real elem_diameter = 1.0/((double)nelem);
  Real dt_max = elem_diameter / shear_wave_speed;
  // dt = max(dt, 0.9*dt_max); 
  Real alpha_CFL = 0.1;
  // applied_force *= shear_wave_speed;

  equation_systems.parameters.set<Real> ("time") = 0;
  equation_systems.parameters.set<Real> ("dt")   = dt;
  equation_systems.parameters.set<VectorValue<Real>> ("applied_force")   = applied_force;
  equation_systems.parameters.set<Real> ("applied_damping")   = applied_damping;
  equation_systems.parameters.set<Real> ("PD_horizon") = PD_horizon;
  equation_systems.parameters.set<std::vector<Real>> ("PD_params") = PD_params;
  equation_systems.parameters.set<bool> ("es_compute_matrix") = 1;

  dsystem.add_matrix("mass");
  dsystem.add_matrix("deformation_gradient");
  dsystem.add_matrix("elem_adjacency");
  dsystem.add_matrix("quad_adjacency");
  // Initialize the data structures for the equation system.
  equation_systems.init();

  //preassble
  std::cout << " preassemble" << std::endl;
  preassemble_elasticity(equation_systems, "Elasticity");


{
  compute_mass(equation_systems, "Elasticity");
}


  // cache_assemble(equation_systems, "Elasticity");
  std::cout << "Matrix assembled" << std::endl;
  // Print information about the system to the screen.
  equation_systems.print_info();


  //timesteps
  std::string exodus_filename = simname;
  std::string exodus_filename_eq =  "EQ_"+simname;;
  
  ExodusII_IO exo(mesh);
  exo.write_equation_systems (exodus_filename, equation_systems);
  exo.append(true);

  ExodusII_IO exoeq(mesh);
  exoeq.write_equation_systems (exodus_filename_eq, equation_systems);
  exoeq.append(true);


  system.time = 0.;
  unsigned int force_stage=0;
  double wc0 = get_wall_time();
  Real max_speed = 1.0;
  std::cout << shear_wave_speed << std::endl;
  int dt_stable_inc =0;
  bool first_t_step=1;

  int save_counter = 0;
  for (unsigned int t_step = 0; t_step < nt; t_step++)
    {
      dsystem.solve();


      equation_systems.parameters.set<bool> ("es_compute_matrix") = 1;
      *dsystem.solution = *dsystem.current_local_solution;
      system.current_local_solution->add(  dt ,*dsystem.current_local_solution );
      *system.solution = *system.current_local_solution;

      Real ratio_increment_displacement = ratio_norm_difference(system.current_local_solution , system.old_local_solution , system.older_local_solution,1);
      Real ratio_increment_velocity = ratio_norm_difference(dsystem.current_local_solution , dsystem.old_local_solution , dsystem.older_local_solution,1);
      if (t_step==0 || (t_step+1)%inc_write == 0)
      {
        libMesh::out << " Time= " << system.time << " dt " << dt;
        libMesh::out << " Max Displ " << " " << system.solution->max() << " Vel = " << dsystem.solution->max() << " CFL " << max_speed << " " << shear_wave_speed;
        libMesh::out << " Convg u:"  <<  ratio_increment_displacement << " du: " << ratio_increment_velocity << " wct " << get_wall_time()-wc0 << std::endl;
        wc0 = get_wall_time();
      }

      system.time += dt;
      dsystem.time += dt;
      *system.old_local_solution = *system.current_local_solution;
      *dsystem.old_local_solution = *dsystem.current_local_solution;
      system.solution->close();
      dsystem.solution->close();
      equation_systems.update();

      //ramp up forcing
      applied_force *= std::min(1.0 , dsystem.time / applied_force_ramp_up_duration );

      if ((t_step+1)%inc_write == 0)
      {
          save_counter++;
          exo.write_timestep (exodus_filename, equation_systems, save_counter, system.time);
      }

      //update force after EQ reached
      if ( ratio_increment_displacement < tol_ratio_1norm_increment_displacement )
      {
        std::cout << "EQ reached, Forcing Stage: " << force_stage << std::endl;
        exoeq.write_timestep (exodus_filename_eq, equation_systems, t_step+1, system.time);


        if( system.solution->max() > max_extension || dt < TOLERANCE)
          break;

        //update force
        applied_force(0) *= force_increment;
        applied_force(1) *= force_increment;
        // applied_damping *= 10 * force_increment;
        compute_mass(equation_systems, "Elasticity");

        force_stage++;
        dt_stable_inc=0;
      }
      dt_stable_inc++;
      
  // std::cout << system.old_local_solution.size() << std::endl;
    }
ExodusII_IO(mesh).write_equation_systems ("final.e", equation_systems);

double wctotal = get_wall_time();
std::cout << "Total Simulation Time: " << wctotal << std::endl;
return 0;
}























void preassemble_elasticity(EquationSystems & es,
                    const std::string & libmesh_dbg_var(system_name))
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & system = es.get_system<TransientLinearImplicitSystem>("Elasticity");
  TransientLinearImplicitSystem & dsystem = es.get_system<TransientLinearImplicitSystem>("dElasticity");
  const unsigned int u_var = system.variable_number ("displacementsx");
  const unsigned int v_var = system.variable_number ("displacementsy");
  const unsigned int du_var = dsystem.variable_number ("velocitiesx");
  const unsigned int dv_var = dsystem.variable_number ("velocitiesy");
  const DofMap & dof_map = system.get_dof_map();
  Real dt = es.parameters.get<Real>   ("dt");
  const Real applied_damping = es.parameters.get<Real> ("applied_damping");
  VectorValue<Real> applied_force = es.parameters.get<VectorValue<Real>> ("applied_force");
  const Real PD_horizon = es.parameters.get<Real> ("PD_horizon");
  const std::vector<Real> PD_params = es.parameters.get<std::vector<Real>> ("PD_params");
  bool es_compute_matrix = es.parameters.get<bool>   ("es_compute_matrix");

  //  std::vector<Real> nonlocal_quadrature_weights = es.parameters.get<std::vector<Real>> ("nonlocal_quadrature_weights");
  //  std::vector<Point> nonlocal_quadrature_points =  es.parameters.get<std::vector<Point>> ("nonlocal_quadrature_points");
  // const unsigned int num_nonlocal_quadrature = es.parameters.get<unsigned int> ("num_nonlocal_quadrature");


  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order()); //QUAD_ORDER
  fe->attach_quadrature_rule (&qrule);
  DenseMatrix<Number> Me;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Muu(Me), Muv(Me),
    Mvu(Me), Mvv(Me);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe);

  const unsigned int nelem = mesh.n_elem();
  const unsigned int nquad = qrule.n_points();
  const unsigned int nreserve = nelem * nquad * 1;
  std::map<dof_id_type,unsigned int> local_elem_id;

  const std::vector<Real> & JxW = fe->get_JxW();
  std::vector<std::vector<Real>> JxW_elem;
  const std::vector<Point> & q_point = fe->get_xyz();
  std::vector<std::vector<Point>> q_point_elem;
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  std::vector<std::vector<std::vector<Real>>> phi_elem;
  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_elem;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<std::vector<dof_id_type>> dof_indices_elem_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<std::vector<dof_id_type>> dof_indices_elem_v;
  std::vector<Point> q_point_bc_right;

  JxW_elem.reserve(nreserve);
  q_point_elem.reserve(nreserve);
  phi_elem.reserve(nreserve);
  dof_indices_elem.reserve(nreserve);
  dof_indices_elem_u.reserve(nreserve);
  dof_indices_elem_v.reserve(nreserve);
  q_point_bc_right.reserve(nreserve);

  // std::vector<

  // const std::vector<Point> & q_point = fe->get_xyz();
  const Patch patch(mesh.processor_id());
  // std::vector< const Patch > patch_elem;
  std::vector<std::vector< const Elem *>> patch_elem;
  patch_elem.reserve(nelem*100);


  // Pre-compute.
  unsigned int local_elem_index = 0;
  for (const auto & elem : mesh.active_element_ptr_range())
    {
      // local_elem_id[elem->id()] = local_elem_index;
      local_elem_index = int(elem->id());

      //DOF INDICES
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);

      fe->reinit (elem);

      local_elem_id[elem->id()] = local_elem_index;
      JxW_elem.push_back(JxW);
      q_point_elem.push_back(q_point);
      if(elem->subdomain_id()==1)
        for(int qp=0; qp<q_point.size();qp++)
        {
          q_point_bc_right.push_back(q_point[qp]);        
        }


      phi_elem.push_back(phi);
      dof_indices_elem.push_back(dof_indices);
      dof_indices_elem_u.push_back(dof_indices_u);
      dof_indices_elem_v.push_back(dof_indices_v);


      //patch
      libMesh::Patch patch(mesh.processor_id());
      libMesh::Patch::PMF patchtype = &Patch::add_face_neighbors;
      patch.build_around_element(elem, target_patch_size, patchtype);
      std::vector<const Elem *> patchvec;
        patchvec.reserve(patch.size());
      Patch::const_iterator        patch_it  = patch.begin();
      const Patch::const_iterator  patch_end = patch.end();
      for (; patch_it != patch_end; ++patch_it)
      {
        const Elem * elem2 = *patch_it;
        patchvec.push_back(elem2);
      }
      patch_elem.push_back(patchvec);

      ++local_elem_index;
    }

  es.parameters.set<std::vector<std::vector<Real>> > ("JxW_elem") = JxW_elem;
  es.parameters.set<std::vector<std::vector<Point>>> ("q_point_elem") = q_point_elem;;
  es.parameters.set<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem") = dof_indices_elem;;
  es.parameters.set<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem_u") = dof_indices_elem_u;;
  es.parameters.set<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem_v") = dof_indices_elem_v;;
  es.parameters.set<std::vector<std::vector< const Elem *>> > ("patch_elem") = patch_elem;;
  es.parameters.set<std::vector<Point>> ("q_point_bc_right") = q_point_bc_right;

}






void compute_mass(EquationSystems & es,
                    const std::string & libmesh_dbg_var(system_name))
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & system = es.get_system<TransientLinearImplicitSystem>("Elasticity");
  TransientLinearImplicitSystem & dsystem = es.get_system<TransientLinearImplicitSystem>("dElasticity");
  const unsigned int u_var = system.variable_number ("displacementsx");
  const unsigned int v_var = system.variable_number ("displacementsy");
  const unsigned int du_var = dsystem.variable_number ("velocitiesx");
  const unsigned int dv_var = dsystem.variable_number ("velocitiesy");
  const DofMap & dof_map = system.get_dof_map();
  Real dt = es.parameters.get<Real>   ("dt");
  const Real applied_damping = es.parameters.get<Real> ("applied_damping");
  VectorValue<Real> applied_force = es.parameters.get<VectorValue<Real>> ("applied_force");
  const Real PD_horizon = es.parameters.get<Real> ("PD_horizon");
  const std::vector<Real> PD_params = es.parameters.get<std::vector<Real>> ("PD_params");
  bool es_compute_matrix = es.parameters.get<bool>   ("es_compute_matrix");
  SparseMatrix<Number> & mass = dsystem.get_matrix("mass");
  mass.zero();

  //preassembled values
  auto & JxW_elem = es.parameters.get<std::vector<std::vector<Real>> > ("JxW_elem") ;
  auto & q_point_elem = es.parameters.get<std::vector<std::vector<Point>>> ("q_point_elem") ;
  auto & dof_indices_elem = es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem") ;
  auto & dof_indices_elem_u = es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem_u") ;
  auto & dof_indices_elem_v = es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem_v") ;
  auto & patch_elem = es.parameters.get<std::vector<std::vector< const Elem *>> > ("patch_elem");

  // es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem") = dof_indices_elem;;
  //  std::vector<Real> nonlocal_quadrature_weights = es.parameters.get<std::vector<Real>> ("nonlocal_quadrature_weights");
  //  std::vector<Point> nonlocal_quadrature_points =  es.parameters.get<std::vector<Point>> ("nonlocal_quadrature_points");
  // const unsigned int num_nonlocal_quadrature = es.parameters.get<unsigned int> ("num_nonlocal_quadrature");


  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order()); //QUAD_ORDER
  fe->attach_quadrature_rule (&qrule);
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Real> & JxW = fe->get_JxW();
  // const std::vector<Point> & q_point = fe->get_xyz();
  // const std::vector<std::vector<Real>> & phi = fe->get_phi();
  DenseMatrix<Number> Me;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Muu(Me), Muv(Me),
    Mvu(Me), Mvv(Me);

    std::vector<dof_id_type> dof_indices_v,dof_indices_u,dof_indices;
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {

      dof_map.dof_indices(elem, dof_indices);
      dof_map.dof_indices(elem, dof_indices_v,u_var);
      dof_map.dof_indices(elem, dof_indices_u,v_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();

      fe->reinit(elem);

      Me.resize (n_dofs, n_dofs);
      Muu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Mvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Real penalty = 1.0;


      if(elem->subdomain_id()==1)
        penalty +=set_penalty;
      for (unsigned int qp=0; qp<qrule.n_points(); qp++){
        for (unsigned int i=0; i<n_u_dofs; i++){
          for (unsigned int j=0; j<n_u_dofs; j++){
              Muu(i,j) += JxW[qp] * phi[i][qp] * phi[j][qp] * (1.* penalty +  dt * applied_damping * flag_implicit_damping) ;
              Mvv(i,j) += JxW[qp] * phi[i][qp] * phi[j][qp] * (1.* penalty +  dt * applied_damping * flag_implicit_damping);
            }
        }
      }

      
      mass.add_matrix(Me,dof_indices);
    }//elem
    mass.close();

}













void assemble_elasticity(EquationSystems & es,
                    const std::string & libmesh_dbg_var(system_name))
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & system = es.get_system<TransientLinearImplicitSystem>("Elasticity");
  const unsigned int u_var = system.variable_number ("displacementsx");
  const unsigned int v_var = system.variable_number ("displacementsy");
  TransientLinearImplicitSystem & dsystem = es.get_system<TransientLinearImplicitSystem>("dElasticity");
  const unsigned int du_var = dsystem.variable_number ("velocitiesx");
  const unsigned int dv_var = dsystem.variable_number ("velocitiesy");
  ExplicitSystem & fsystem = es.get_system<ExplicitSystem>("fElasticity");
  const unsigned int fu_var = fsystem.variable_number ("forcesx");
  const unsigned int fv_var = fsystem.variable_number ("forcesy");
  const DofMap & dof_map = system.get_dof_map();
  Real dt = es.parameters.get<Real>   ("dt");
  const Real applied_damping = es.parameters.get<Real> ("applied_damping");
  VectorValue<Real> applied_force = es.parameters.get<VectorValue<Real>> ("applied_force");
  const Real PD_horizon = es.parameters.get<Real> ("PD_horizon");
  const std::vector<Real> PD_params = es.parameters.get<std::vector<Real>> ("PD_params");
  bool es_compute_matrix = es.parameters.get<bool>   ("es_compute_matrix");
  SparseMatrix<Number> & mass = dsystem.get_matrix("mass");


  //preassembled values
  auto & JxW_elem = es.parameters.get<std::vector<std::vector<Real>> > ("JxW_elem") ;
  auto & q_point_elem = es.parameters.get<std::vector<std::vector<Point>>> ("q_point_elem") ;
  auto & dof_indices_elem = es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem") ;
  auto & dof_indices_elem_u = es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem_u") ;
  auto & dof_indices_elem_v = es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem_v") ;
  auto & patch_elem = es.parameters.get<std::vector<std::vector< const Elem *>> > ("patch_elem");
  auto & q_point_bc_right = es.parameters.get<std::vector<Point>> ("q_point_bc_right");
  unsigned int num_quad_points_forced = q_point_bc_right.size();
  // es.parameters.get<std::vector<std::vector<dof_id_type>>> ("dof_indices_elem") = dof_indices_elem;;
  //  std::vector<Real> nonlocal_quadrature_weights = es.parameters.get<std::vector<Real>> ("nonlocal_quadrature_weights");
  //  std::vector<Point> nonlocal_quadrature_points =  es.parameters.get<std::vector<Point>> ("nonlocal_quadrature_points");
  // const unsigned int num_nonlocal_quadrature = es.parameters.get<unsigned int> ("num_nonlocal_quadrature");


  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order()); //QUAD_ORDER
  fe->attach_quadrature_rule (&qrule);
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
  const std::vector<std::vector<RealTensor> >& d2phi = fe->get_d2phi();
  DenseMatrix<Number> Me;
  DenseVector<Number> Fe;

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe);

  const unsigned int nelem = mesh.n_elem();
  const unsigned int nquad = qrule.n_points();
  const unsigned int nreserve = nelem * nquad * 10000;
  std::map<dof_id_type,unsigned int> local_elem_id;

  std::vector<double>  u_dof,v_dof, du_dof, dv_dof;
  system.old_local_solution->localize(u_dof);//, &dof_indices_u);
  dsystem.old_local_solution->localize(du_dof);//, dof_indices_v);

  // Compute.
  RealVectorValue u_old(2), u_old2(2),xi(2),xihat(2),eta(2),etahat(2),xieta(2),xietahat(2),pd_force(2);
  Real xinorm, etanorm, xietanorm, xietahat_dot_xihat, dilatation;

  Real PD_mass = pi * PD_horizon4;
  Real PD_C = 0.75 * young_modulus * 2 * 2 / PD_mass; 
  Real kern;
  Real Lambda;
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // unsigned int local_elem_index = local_elem_id[elem->id()];
      unsigned int local_elem_index = int(elem->id());
      auto & JxW = JxW_elem[local_elem_index];
      const auto & q_point = q_point_elem[local_elem_index];
      // const auto & phi = phi_elem[local_elem_index];
      fe->reinit(elem);
      // fe0->reinit(elem);

      const auto & dof_indices = dof_indices_elem[local_elem_index]; //non-const?
      const auto & dof_indices_u = dof_indices_elem_u[local_elem_index];
      const auto & dof_indices_v = dof_indices_elem_v[local_elem_index];
      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();

      Me.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);
      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_v_dofs, n_v_dofs);



      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {

        //Either PD or CG
        Number du_old = 0., dv_old = 0.;
        for (std::size_t l=0; l<phi.size(); l++){
            du_old += phi[l][qp] * du_dof[dof_indices_u[l]];
            dv_old += phi[l][qp] * du_dof[dof_indices_v[l]];
        }

        for (unsigned int i=0; i < n_u_dofs; i++){
          Fu(i) += JxW[qp] * du_old * phi[i][qp] * (1 - dt * applied_damping * (1-flag_implicit_damping) );
          Fv(i) += JxW[qp] * dv_old * phi[i][qp] * (1 - dt * applied_damping * (1-flag_implicit_damping) );
        }


        //If using CG-FE
        if(use_method==0)
        {
          DenseMatrix<Number> gradu_old(2, 2);
          gradu_old.zero();

          for (unsigned int var_j=0; var_j<2; var_j++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              gradu_old(0,var_j) += dphi[j][qp](var_j) * u_dof[dof_indices_u[j]];

          for (unsigned int var_j=0; var_j<2; var_j++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              gradu_old(1,var_j) += dphi[j][qp](var_j) * u_dof[dof_indices_v[j]];

          for (unsigned int i=0; i < n_u_dofs; i++){
            for (unsigned int var_j=0; var_j<2; var_j++){
              Fu(i) += -JxW[qp] * gradu_old(0,var_j) * dphi[i][qp](var_j) * dt;
              Fv(i) += -JxW[qp] * gradu_old(1,var_j) * dphi[i][qp](var_j) * dt;            
            }
          }          

          // for (unsigned int i=0; i < n_u_dofs; i++){
          //   for (unsigned int var_j=0; var_j<2; var_j++){
          //     Fu(i) += JxW[qp] * gradu_old(0,var_j) * phi[i][qp];
          //     Fv(i) += JxW[qp] * gradu_old(1,var_j) * phi[i][qp];            
          //   }
          // }

          // DenseVector<Number> divgradu_old(2);
          // // for (unsigned int j=0; j<n_u_dofs; j++)
          // // {
          // //     divgradu_old(0) += d2phi[j][qp]( 0,0) * u_dof[dof_indices_u[j]] + d2phi[j][qp]( 0,1) * u_dof[dof_indices_v[j]];
          // //     divgradu_old(1) += d2phi[j][qp]( 1,0) * u_dof[dof_indices_u[j]] + d2phi[j][qp]( 1,1) * u_dof[dof_indices_v[j]];
          // // }
          // // for (unsigned int j=0; j<n_u_dofs; j++)
          // // {
          // //     divgradu_old(0) += d2phi[j][qp]( 0,0) * u_dof[dof_indices_u[j]] + d2phi[j][qp]( 1,0) * u_dof[dof_indices_v[j]];
          // //     divgradu_old(1) += d2phi[j][qp]( 0,1) * u_dof[dof_indices_u[j]] + d2phi[j][qp]( 1,1) * u_dof[dof_indices_v[j]];
          // // }
          //   for (unsigned int j=0; j<n_u_dofs; j++){
          //     for (unsigned int var_j=0; var_j<2; var_j++){
          //       divgradu_old(0) += d2phi[j][qp]( var_j,0) * u_dof[dof_indices_u[j]];
          //       divgradu_old(1) += d2phi[j][qp]( var_j,1) * u_dof[dof_indices_v[j]];
          //     }
          //   }

          // for (unsigned int i=0; i < n_u_dofs; i++){
          //   Fu(i) += 2 * 0.75 * young_modulus * JxW[qp]  * phi[i][qp] * divgradu_old(0) * dt;
          //   Fv(i) += 2 * 0.75 * young_modulus * JxW[qp]  * phi[i][qp] * divgradu_old(1) * dt;
          // }

        }

      }


      //If using PD-FE
      if(use_method==1)
      {
        const auto & patch = patch_elem[local_elem_index];
        for (auto elem2 : patch)
        {
          unsigned int local_elem_index2 = int(elem2->id());
          const std::vector<Real> & JxW2 = JxW_elem[local_elem_index2];
          const std::vector<Point> & q_point2 = q_point_elem[local_elem_index2];
          // const std::vector<std::vector<Real> > & phi2 = phi_elem[local_elem_index2];
          const std::vector<dof_id_type> & dof_indices2 = dof_indices_elem[local_elem_index2];
          const std::vector<dof_id_type> & dof_indices_u2 = dof_indices_elem_u[local_elem_index2];
          const std::vector<dof_id_type> & dof_indices_v2 = dof_indices_elem_v[local_elem_index2];

          const unsigned int n_dofs2   = dof_indices2.size();
          const unsigned int n_u_dofs2 = dof_indices_u2.size();
          const unsigned int n_v_dofs2 = dof_indices_v2.size();

          for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
            xi.assign(RealVectorValue( q_point[qp](0), q_point[qp](1) ));
            u_old.zero();
            for (std::size_t l=0; l<phi.size(); l++){
              u_old += RealVectorValue(phi[l][qp] * u_dof[dof_indices_u[l]] , phi[l][qp] * u_dof[dof_indices_v[l]]);
            }

            for (unsigned int qp2=0; qp2<qrule.n_points(); qp2++){
              xi.assign(RealVectorValue( q_point[qp](0) - q_point2[qp2](0)  ,  q_point[qp](1) - q_point2[qp2](1) ));
              xinorm = xi.norm();
              if(xinorm > TOLERANCE)
              {
                u_old2.zero();
                for (std::size_t l=0; l<phi.size(); l++){
                    u_old2.slice(0) += phi[l][qp2] * u_dof[dof_indices_u2[l]]; //dof_indices_u2
                    u_old2.slice(1) += phi[l][qp2] * u_dof[dof_indices_v2[l]];
                }
                eta.assign( u_old-u_old2);
                xieta = xi + eta;
                Lambda = (PD_horizon2 / xi.norm_sq())  * ( 1. - (xieta * xi) / xieta.norm_sq() );
                pd_force = 2 * PD_C * Lambda * xieta;

                // pd_force.print();
                Real penalty = 1.0;
                if(elem->subdomain_id()==1)
                  penalty+=set_penalty;
              
                for (unsigned int i=0; i<n_u_dofs; i++)
                {
                  Real tmpF = dt * (-1.0) * JxW[qp]* JxW2[qp2]*phi[i][qp];
                  Fu(i) += tmpF * pd_force(0)  * 1 ;
                  Fv(i) += tmpF * pd_force(1)  * 1 ;
                }             
              }//if(xi.norm() > TOLERANCE)
            }//qp2
          }//qp
        }//elem2        
      }


      
      Real penalty = 1.0;
      if(elem->subdomain_id()==1){
        penalty+=set_penalty;
        for (unsigned int qp=0; qp<qrule.n_points(); qp++){
          for (unsigned int i=0; i<n_u_dofs; i++){
            Fu(i) +=  penalty * JxW[qp]  * phi[i][qp] * force_attenuated( applied_force(0) , q_point[qp](0) );
            Fv(i) +=  penalty * JxW[qp]  * phi[i][qp] * force_attenuated( applied_force(1) , q_point[qp](0) );
          }
        }
      }
      if(elem->subdomain_id()==3){
        penalty+=set_penalty;
        for (unsigned int qp=0; qp<qrule.n_points(); qp++){
          for (unsigned int i=0; i<n_u_dofs; i++){
            Fu(i) =  0.0;
            Fv(i) =  0.0;
          }
        }
      }

      dsystem.rhs->add_vector    (Fe, dof_indices);
      dsystem.matrix->add_matrix    (Me, dof_indices);
      fsystem.solution->add_vector(Fe, dof_indices);
    }//elem
    dsystem.matrix->close();
      dsystem.matrix->add    (1.,mass);
      dsystem.matrix->close();
      fsystem.solution->close();
      fsystem.update();
}







Real PD_stretch( Point xi , Point eta)
{
  Real norm_eta = (xi + eta).norm();
  Real norm_xi = xi.norm();
  return (norm_eta / norm_xi - 1.0);
}

Point PD_BB_force( Point x_deformed, Real stretch, const std::vector<Real> & PD_params )
{
  Real norm = x_deformed.norm();
  if(norm > 0)
    return (PD_params[0] * stretch) * x_deformed / norm;
  else
    return Point(0.,0.,0.);
}


Real eval_elasticity_tensor(unsigned int i,
                            unsigned int j,
                            unsigned int k,
                            unsigned int l)
{
  // Define the Poisson ratio and Young's Modulus
  const Real nu = 0.3;
  const Real youngs_mod = 1.0;

  // Define the Lame constants (lambda_1 and lambda_2) based on Poisson ratio
  // Lame 1
  const Real lambda_1 = youngs_mod * nu / ((1. + nu) * (1. - 2.*nu));
  // Lame 2 = Shear modulus
  const Real lambda_2 = youngs_mod * 0.5 / (1 + nu);

  // Define the Kronecker delta functions that we need here
  Real delta_ij = (i == j) ? 1. : 0.;
  Real delta_il = (i == l) ? 1. : 0.;
  Real delta_ik = (i == k) ? 1. : 0.;
  Real delta_jl = (j == l) ? 1. : 0.;
  Real delta_jk = (j == k) ? 1. : 0.;
  Real delta_kl = (k == l) ? 1. : 0.;

  return lambda_1 * delta_ij * delta_kl + lambda_2 * (delta_ik * delta_jl + delta_il * delta_jk);
}





RealTensor kernel2( Point xq , Point xq_nl)
{
  Point diff = xq_nl - xq;

  if(diff.norm_sq() < TOLERANCE)
    return TensorValue<Real>(0.,0.,0.,0.,0.,0.,0.,0.,0.);

  Real knorm = 1.0 / diff.norm();
  knorm *= knorm;
  knorm *= knorm;
  // std::cout << knorm << std::endl;
  return TensorValue<Real>( knorm*diff(0)*diff(0), knorm*diff(0)*diff(1), 0., knorm*diff(0)*diff(1), knorm*diff(1)*diff(1), 0., 0.);
}
