//
//  dynamicEq.C
//  EquilibriumSolver
//
//  Created by Ben Vadala-Roth on 11/16/15.
//  Copyright © 2015 ben. All rights reserved.
//
#include <algorithm>
#include <cmath>
#include <iostream>
#include "libmesh/analytic_function.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/fe.h"
#include "libmesh/fe_map.h"
#include "libmesh/getpot.h"
#include "libmesh/gmv_io.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
//#include "dynamicEq.h"
using namespace libMesh;
using namespace std;
static int dim;
static double rho;
static double eta;
static double mu;
static double kappa;
static double nu;
static double T1x, T1y, T2x, T2y, T3x, T3y, T4x, T4y;
static double fx, fy;
static double gx, gy;
static Real integral_J;
static Real max_int_J;
RealVectorValue evaluate_body_force(const Real x, const Real y, const Real z, const Real t);
RealTensor evaluate_first_PK_stress(RealTensor FF);
void assemble_equilibrium(EquationSystems& es, const string& system_name);
void assemble_position(EquationSystems& es, const string& system_name);
void apply_initial_position(EquationSystems& es, const string& system_name);
void apply_initial_velocity(EquationSystems& es, const string& system_name);
Real
initial_position(const Point& p, const Parameters& /*parameters*/, const string& /*system*/, const string& var_name)
{
//    if (var_name == "x")
//    {
//        return p(0);
//    }
//    else if (var_name == "y")
//    {
//        return p(1);
//    }
//    else if (var_name == "z")
//    {
//        return p(2);
//    }
//    else
//    {
//        return numeric_limits<double>::quiet_NaN();
//    }
    return 0.0;
}
Real
initial_velocity(const Point& /*p*/,
                 const Parameters& /*parameters*/,
                 const string& /*system*/,
                 const string& /*var_name*/)
{
    // could include something more complicated
    return 0.0;
}
Real dirichlet_function(const Real x,
			const Real y,
			const Real z,
			const Real t)
{
  return gx*sin(t);
}
void dirichlet_fcn_wrapper(DenseVector<Real>& output,
			   const Point& p,
			   Real time)
{
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  output(0) = dirichlet_function(x,y,0, time);
  output(1) = dirichlet_function(x,y,0, time);
  return;
}
int
main(int argc, const char* argv[])
{
    LibMeshInit init(argc, argv);
    GetPot input_file("input.2d");
    dim = 2;
    int ps = 64;
    int ps_y = ps;
    FEFamily family = LAGRANGE;
    // physical parameters
    if (input_file.search(1, "rho")) rho     = input_file.next(rho);
    if (input_file.search(1, "eta")) eta     = input_file.next(eta);
    if (input_file.search(1, "mu")) mu       = input_file.next(mu);
    if (input_file.search(1, "kappa")) kappa = input_file.next(kappa);
    if (input_file.search(1, "nu"))
    {
      nu = input_file.next(nu);
      if (nu < 0.0 || nu >= .5)
      {
	nu = .3;
	cout << "Possion ratio supplied was not between 0 and .5. Using a Poisson ratio of mu = " << nu << endl;
      }
      //if nu, the Poisson ratio, specified redefine kappa as follows
      kappa = 2*mu*(1+nu)/(3*(1-2*nu));
      cout << "The bulk modulus was calculated as kappa = " << kappa << endl;
    }
    // traction boundary conditions
    if (input_file.search(1, "T1x")) T1x = input_file.next(T1x);
    if (input_file.search(1, "T1y")) T1y = input_file.next(T1y);
    if (input_file.search(1, "T2x")) T2x = input_file.next(T2x);
    if (input_file.search(1, "T2y")) T2y = input_file.next(T2y);
    if (input_file.search(1, "T3x")) T3x = input_file.next(T3x);
    if (input_file.search(1, "T3y")) T3y = input_file.next(T3y);
    if (input_file.search(1, "T4x")) T4x = input_file.next(T4x);
    if (input_file.search(1, "T4y")) T4y = input_file.next(T4y);
    // body force
    if (input_file.search(1, "fx")) fx = input_file.next(fx);
    if (input_file.search(1, "fy")) fy = input_file.next(fy);

    //Dirichlet force, applied to left hand side
    if (input_file.search(1, "gx")) gx = input_file.next(gx);
    if (input_file.search(1, "gy")) gy = input_file.next(gy);
    if (input_file.search(1, "ps")) ps = input_file.next(ps);
    if (input_file.search(1, "ps_y")) ps_y = input_file.next(ps_y);

    // time-stepping parameters
    Real dt = 1e-4;
    unsigned int n_time_steps = 1e5;
    unsigned int stride = 10;
    if (input_file.search(1, "dt")) dt = input_file.next(dt);
    if (input_file.search(1, "n_time_steps")) n_time_steps = input_file.next(n_time_steps);
    if (input_file.search(1, "stride")) stride = input_file.next(stride);
    std::string order_str                        = "FIRST";
    if (input_file.search(1, "order")) order_str = input_file.next(order_str);
    Order order                                  = Utility::string_to_enum<Order>(order_str);

    std::string elem_type_str                            = "QUAD4";
    if (input_file.search(1, "elem_type")) elem_type_str = input_file.next(elem_type_str);
    ElemType elem_type                                   = Utility::string_to_enum<ElemType>(elem_type_str);
    Mesh mesh(init.comm());
    Real halfwidth = dim > 1 ? .2 : 0.;
    Real halfheight = dim > 2 ? 1. : 0.;
    std::string mesh_name;
    if (input_file.search(1, "mesh_name"))
    {
      mesh_name = input_file.next(mesh_name);
      mesh.read(mesh_name);
    }
    else
    {
    // Build grid mesh
      MeshTools::Generation::build_cube(mesh,
					ps,
					(dim > 1) ? ps_y : 0,
					(dim > 2) ? ps : 0,
					-1., 1.,
					-halfwidth, halfwidth,
					-halfheight,halfheight,
					elem_type);
    }
    mesh.print_info();

    //define an analytic function for the Dirichlet BC
    AnalyticFunction<Real> dirichlet_fcn(dirichlet_fcn_wrapper);
    dirichlet_fcn.init();



    EquationSystems equation_systems(mesh);
    // create systems
    ExplicitSystem& pos_system = equation_systems.add_system<ExplicitSystem>("Position");
    unsigned int x_var = pos_system.add_variable("displacemensx", order, family);
    unsigned int y_var = pos_system.add_variable("displacemensy", order, family);
    pos_system.attach_init_function(apply_initial_position);
    LinearImplicitSystem& vel_system = equation_systems.add_system<LinearImplicitSystem>("Velocity");
    unsigned int u_var = vel_system.add_variable("vx", order, family);
    unsigned int v_var = vel_system.add_variable("vy", order, family);
    vel_system.attach_assemble_function(assemble_equilibrium);
    vel_system.attach_init_function(apply_initial_velocity);
    //define Dirichlet BC's
    std::set<boundary_id_type> boundary_ids;
    //in 2d bdry 0 is the bottom, bdry 1 is the right, bdry 2 is the top, bdry 3 is the left
//    boundary_ids.insert(0);
    //boundary_ids.insert(1);

     if (dim >=2)
     {
//        boundary_ids.insert(2);
       boundary_ids.insert(3);
     }

    std::vector<unsigned int> variables(2);
    variables[0] = u_var;
    variables[1] = v_var;

    DirichletBoundary dirichlet_bc(boundary_ids,
				   variables,
				   &dirichlet_fcn);
    vel_system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);


    // initialize and print info of system
    equation_systems.init();
    pos_system.assemble();
    //jac_system.assemble();
    equation_systems.print_info();
    // write out IC's, write all timesteps to same file
    string exodus_filename = "output.e";
    ExodusII_IO exodus_io(mesh);
    vector<Real> J (n_time_steps);
    Real max_abs_J = 1;
    Real current_J = 1;
    const unsigned int num_elem = mesh.n_elem();

    std::cout << "dt: " << dt << ", dt_cfl: " << 0.2 / ps_y / std::sqrt(kappa/rho) << std::endl;
    // time step size and number
    Real time = 0;
    for (unsigned int time_step = 0; time_step < n_time_steps; time_step++)
    {
        Parameters& parameters = equation_systems.parameters;
        parameters.set<Real>("time") = time;
        parameters.set<Real>("dt") = dt;
        parameters.set<unsigned int>("time step") = time_step;
        // set old solution to the current solution
        vel_system.solve();
        pos_system.solution->add(dt, *vel_system.solution);
        pos_system.solution->localize(*pos_system.current_local_solution);

	//Average determinant is (integral of det)/(volume of domain)
	current_J    = integral_J/num_elem;
	J[time_step] = current_J;
	if (abs(current_J) > abs(max_abs_J)) max_abs_J = current_J;


	if (time_step % stride == 0)
        {
            cout << "current time: " << time << endl;
            exodus_io.write_timestep(exodus_filename, equation_systems, time_step / stride + 1, time);
	    std::cout << "Average Determinant of Deformation Gradient: " << integral_J/num_elem << std::endl;
	    std::cout << "Largest Local Deformation: " << max_int_J << std::endl;
        }

        time += dt;
    }
    std::cout << "The largest determinant in magnitude was " << max_abs_J << std::endl;

    return 0;
}
void assemble_equilibrium(EquationSystems& es, const string& system_name)
{
    libmesh_assert_equal_to(system_name, "Velocity");
    const MeshBase& mesh        = es.get_mesh();
    const unsigned int dim      = mesh.mesh_dimension();
    LinearImplicitSystem& vel_system = es.get_system<LinearImplicitSystem>("Velocity");
    ExplicitSystem& pos_system = es.get_system<ExplicitSystem>("Position");

    //ExplicitSystem& jac_system        = es.get_system<ExplicitSystem>("Jacobian");
    //const unsigned int jac_system_num = jac_system.number();  // identifier number for the system
    //NumericVector<Real>& j_solution   = *jac_system.solution;

    const unsigned int u_var = vel_system.variable_number("vx");
    const unsigned int v_var = vel_system.variable_number("vy");
    // const reference to FE tpe for 1st varible in system
    FEType fe_type = vel_system.variable_type(0);
    // build FE object of type fe_type (AutoPtr is like a pointer)
    AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
    AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
    // use a fifth order Gauss Quadrature rule on interior and boundary
    QGauss qrule(dim, FIFTH);
    QGauss qface(dim - 1, FIFTH);
    fe->attach_quadrature_rule(&qrule);
    fe_face->attach_quadrature_rule(&qface);
    const vector<Point>& qpoint = fe->get_xyz();
    const vector<Real>& JxW = fe->get_JxW();
    const vector<vector<Real> >& phi = fe->get_phi();
    const vector<vector<RealGradient> >& dphi = fe->get_dphi();
    const vector<Point>& qpoint_face = fe_face->get_xyz();
    const vector<Point>& normal_face = fe_face->get_normals();
    const vector<Real>& JxW_face = fe_face->get_JxW();
    const vector<vector<Real> >& phi_face = fe_face->get_phi();
    const DofMap& dof_map = vel_system.get_dof_map();
    // data structures to store element matrix and RHS
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
    DenseSubMatrix<Number> Kuu(Ke), Kuv(Ke), Kvu(Ke), Kvv(Ke);
    DenseSubVector<Number> Fu(Fe), Fv(Fe);
    vector<dof_id_type> dof_indices;
    vector<dof_id_type> dof_indices_u;
    vector<dof_id_type> dof_indices_v;
    const Real dt = es.parameters.get<Real>("dt");
    const Real t = es.parameters.get<Real>("time");
    //the intergal of the elements jacobian determinant
    integral_J = 0.0;
    max_int_J  = 0.0;

    // create iterator for the active elements
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (; el != end_el; ++el)
    {
        // pointer to the element
        const Elem* elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        dof_map.dof_indices(elem, dof_indices_u, u_var);
        dof_map.dof_indices(elem, dof_indices_v, v_var);
        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_u.size();
        // calculates location of quad points and element basis functions for
        // current element
        fe->reinit(elem);
        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);
        // reposition the submatrices
        Kuu.reposition(u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition(u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
        Kvu.reposition(v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition(v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        Fu.reposition(u_var * n_u_dofs, n_u_dofs);
        Fv.reposition(v_var * n_u_dofs, n_v_dofs);
	//element volume
	Real elem_vol = 0;
	//average of determinant J over current element
	Real elem_integral_J = 0;

        // loop over all quad points
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            const Real X = qpoint[qp](0);
            const Real Y = qpoint[qp](1);
            Number u = 0;
            Number v = 0;
            Gradient grad_x = 0;
            Gradient grad_y = 0;
            Gradient grad_z = 0;
            // loop over the test functions to construct u and x from previous time
            // step
            for (unsigned int l = 0; l < phi.size(); l++)
            {
                u += phi[l][qp] * (*vel_system.solution)(dof_indices_u[l]);
                v += phi[l][qp] * (*vel_system.solution)(dof_indices_v[l]);
                grad_x += dphi[l][qp] * (*pos_system.solution)(dof_indices_u[l]);
                grad_y += dphi[l][qp] * (*pos_system.solution)(dof_indices_v[l]);
            }
            RealTensor FF = RealTensor(1.0+grad_x(0),
                                       grad_x(1),
                                       0, // grad_x(2),
                                       grad_y(0),
									   1.0+grad_y(1),
                                       0, // grad_y(2),
                                       0,
                                       0,
                                       1.0); // grad_z(0), grad_z(1), grad_z(2));
            RealTensor PP = evaluate_first_PK_stress(FF);
            RealVectorValue Fb = evaluate_body_force(X, Y, 0, t);

	    //add contributions at each quad point
	    elem_integral_J += JxW[qp]*FF.det();
	    elem_vol        += JxW[qp]*1.0;
            // construct the local RHS
            for (unsigned int i = 0; i < phi.size(); i++)
            {
                // here we include the term Mu^n on the right side, beacause we solve
                // for u^{n+1}
                RealVectorValue PP_dphi = PP * dphi[i][qp];
                Fu(i) += JxW[qp] * ((rho / dt) * u * phi[i][qp] - PP_dphi(0) + Fb(0) * phi[i][qp]);
                Fv(i) += JxW[qp] * ((rho / dt) * v * phi[i][qp] - PP_dphi(1) + Fb(1) * phi[i][qp]);

                // construct the local mass matrix
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    Kuu(i, j) += JxW[qp] * ((rho / dt + eta) * phi[i][qp] * phi[j][qp]);
                    Kvv(i, j) += JxW[qp] * ((rho / dt + eta) * phi[i][qp] * phi[j][qp]);
                    // the u, v, and w coupling
                    Kuv(i, j) = 0.0;
                    Kvu(i, j) = 0.0;
                }
            }
        } // end of quadrature point loop

        //if J = det(FF) == 1, then integral_over_elem_e(J) = vol(elem_e)
        //1/vol(elem_e)*integral_over_elem_e(J) = avg J
        elem_integral_J /= elem_vol;
	integral_J      += elem_integral_J;
	//keep track of the largest local deformation
	if (abs(elem_integral_J) > abs(max_int_J)) max_int_J = elem_integral_J;

        // traction BC section
        {
            // normals of the faces of a square, which determine the BC
            Point normal1, normal2, normal3, normal4;
            normal1 = Point(-1.0, 0.0, 0.0);
            normal2 = Point(1.0, 0.0, 0.0);
            normal3 = Point(0.0, -1.0, 0.0);
            normal4 = Point(0.0, 1.0, 0.0);
            for (unsigned int side = 0; side < elem->n_sides(); side++)
            {
                if (elem->neighbor(side) == NULL)
                {
                    fe_face->reinit(elem, side);
                    for (unsigned int qp = 0; qp < qpoint_face.size(); qp++)
                    {
                        for (unsigned int i = 0; i < phi_face.size(); i++)
                        {
                            // boundary 1
//                             if (normal_face[qp] == normal1)
//                             {
//                                 Fu(i) += JxW_face[qp] * (T1x)*phi_face[i][qp];
//                                 Fv(i) += JxW_face[qp] * (T1y)*phi_face[i][qp];
//                             }
                            //boundary 2
                            if (normal_face[qp] == normal2)
                            {
                                Fu(i) += JxW_face[qp] * (T2x)*phi_face[i][qp];
                                Fv(i) += JxW_face[qp] * (T2y)*phi_face[i][qp];
                            }
                            // boundary 3
                            else if (normal_face[qp] == normal3)
                            {
                                Fu(i) += JxW_face[qp] * (T3x)*phi_face[i][qp];
                                Fv(i) += JxW_face[qp] * (T3y)*phi_face[i][qp];
                            }
                            // boundary 4
                            else if (normal_face[qp] == normal4)
                            {
                                Fu(i) += JxW_face[qp] * (T4x)*phi_face[i][qp];
                                Fv(i) += JxW_face[qp] * (T4y)*phi_face[i][qp];
                            }
                        }
                    }
                }
            }
        } // end of traction BC section
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        // element matrix and RHS are now built for this element. Add to global
        // matrix and RHS
        vel_system.matrix->add_matrix(Ke, dof_indices);
        vel_system.rhs->add_vector(Fe, dof_indices);
    }
}
RealTensor evaluate_first_PK_stress(RealTensor FF)
{
    double J = FF.det();
    RealTensor FF_inv_trans = RealTensor(FF(1, 1), -FF(0, 1), 0, -FF(1, 0), FF(0, 0), 0, 0, 0, 1.0) / J;
    double I1 = FF(0,0)*FF(0,0)
    		  + FF(1,0)*FF(1,0)
    		  + FF(2,0)*FF(2,0)
    		  + FF(0,1)*FF(0,1)
    		  + FF(1,1)*FF(1,1)
    		  + FF(2,1)*FF(2,1)
    		  + FF(0,2)*FF(0,2)
    		  + FF(1,2)*FF(1,2)
    		  + FF(2,2)*FF(2,2);
    return mu * std::pow(J, -2.0/3 ) * (FF - I1 / 3.0 * FF_inv_trans) + kappa * std::log(J) * FF_inv_trans;
//    return mu*FF + (kappa*log(J) - mu)*FF_inv_trans;
}
RealVectorValue evaluate_body_force(const Real X, const Real Y, const Real /*Z*/, const Real t)
{
//  return RealVectorValue(fx * sin(X) * exp(-0.1 * t), fy, 0.0);
  return RealVectorValue(0.0, -1e-3);
}
void apply_initial_position(EquationSystems& es, const string& system_name)
{
    libmesh_assert_equal_to(system_name, "Position");
    ExplicitSystem& system = es.get_system<ExplicitSystem>("Position");
    es.parameters.set<Real>("time") = system.time = 0;
    system.project_solution(initial_position, NULL, es.parameters);
    return;
}
void apply_initial_velocity(EquationSystems& es, const string& system_name)
{
    libmesh_assert_equal_to(system_name, "Velocity");
    LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Velocity");
    es.parameters.set<Real>("time") = system.time = 0;
    system.project_solution(initial_velocity, NULL, es.parameters);
    return;
}
