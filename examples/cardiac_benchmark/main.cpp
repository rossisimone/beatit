/*
 * main.cpp
 *
 *  Created on: Nov 11, 2021
 *      Author: srossi
 */

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/transient_system.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// The definition of a geometric element
#include "libmesh/elem.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

// To impose Dirichlet boundary conditions
#include "libmesh/getpot.h"
#include "Util/IO/io.hpp"
#include "Util/TimeData.hpp"

void dyadic(const libMesh::VectorValue<double> &v1, const libMesh::VectorValue<double> &v2, libMesh::TensorValue<double> &out)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out(i, j) = v1(i) * v2(j);
}

enum class Boundary
{
    ENDO, EPI, BASE
};

void assemble_RHS(libMesh::EquationSystems &es, double tn, double dt, double tau, int step);
void assemble_mass(libMesh::EquationSystems &es, double dt, double rho = 1.0);

struct ActiveTension
{
    static double tau;
    static double tau_n;

    static constexpr double sigma0 = 1e5;
    static constexpr double gamma = 5e-3;
    static constexpr int alpha_min = -30;
    static constexpr int alpha_max = 5;
    static constexpr double t_sys = 0.17-0.17;
    static constexpr double t_dias = 0.484-0.17;
    static double S(double delta);
    static double eval_tau_rhs(double time, double tau_tn);
    static void solve_tau(double time, double dt);
};

double ActiveTension::tau = 0.0;
double ActiveTension::tau_n = 0.0;

int main(int argc, char **argv)
{
    // Initialize libMesh and any dependent libraries, like in example 2.
    libMesh::LibMeshInit init(argc, argv);
    // Read input file
    GetPot data = BeatIt::readInputFile(argc, argv);

    libMesh::ParallelMesh mesh(init.comm());
    libMesh::ExodusII_IO importer(mesh);
    std::string mesh_file = data("mesh", "NONE");
    if (mesh_file != "NONE")
    {
        importer.read(mesh_file);
        mesh.read(mesh_file);
    }
    else
    {
        std::cout << "Meshfile: NONE -> FIX IT!" << std::endl;
        throw std::runtime_error("bad mesh");
    }

    std::vector < std::string > elem_var_names = importer.get_elem_var_names();
    std::cout << "Mesh file contains the following elemental variables: " << std::endl;
    for (auto &&v : elem_var_names)
        std::cout << v << std::endl;


    libMesh::EquationSystems es(mesh);

    libMesh::ExplicitSystem &f_system = es.add_system < libMesh::ExplicitSystem > ("fibers");
    f_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    f_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    libMesh::ExplicitSystem &s_system = es.add_system < libMesh::ExplicitSystem > ("sheets");
    s_system.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_system.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    s_system.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    libMesh::ExplicitSystem &n_system = es.add_system < libMesh::ExplicitSystem > ("xfibers");
    n_system.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_system.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    n_system.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);

    libMesh::TransientExplicitSystem &u_system = es.add_system < libMesh::TransientExplicitSystem > ("displacement");
    u_system.add_variable("ux", libMesh::FIRST);
    u_system.add_variable("uy", libMesh::FIRST);
    u_system.add_variable("uz", libMesh::FIRST);
    u_system.add_vector("u1");
    u_system.add_vector("u2");
    u_system.add_vector("u3");

    libMesh::TransientExplicitSystem &v_system = es.add_system < libMesh::TransientExplicitSystem > ("velocity");
    v_system.add_variable("vx", libMesh::FIRST);
    v_system.add_variable("vy", libMesh::FIRST);
    v_system.add_variable("vz", libMesh::FIRST);
    v_system.add_vector("mass");
    v_system.add_vector("v1");
    v_system.add_vector("v2");
    v_system.add_vector("v3");

    libMesh::TransientExplicitSystem &p_system = es.add_system < libMesh::TransientExplicitSystem > ("pressure");
    p_system.add_variable("p", libMesh::FIRST);
    p_system.add_vector("mass");
    p_system.add_vector("p1");
    p_system.add_vector("p2");
    p_system.add_vector("p3");

    std::cout << "initialize equation systems: " << std::endl;
    es.init();

    std::cout << "Copy elemental variables: " << std::endl;
    importer.copy_elemental_solution(f_system, "fibersx", "fibersx");
    importer.copy_elemental_solution(f_system, "fibersy", "fibersy");
    importer.copy_elemental_solution(f_system, "fibersz", "fibersz");
    importer.copy_elemental_solution(s_system, "sheetsx", "sheetsx");
    importer.copy_elemental_solution(s_system, "sheetsy", "sheetsy");
    importer.copy_elemental_solution(s_system, "sheetsz", "sheetsz");
    importer.copy_elemental_solution(n_system, "xfibersx", "xfibersx");
    importer.copy_elemental_solution(n_system, "xfibersy", "xfibersy");
    importer.copy_elemental_solution(n_system, "xfibersz", "xfibersz");


    BeatIt::TimeData data_time;
    data_time.setup(data);
    data_time.print();

    std::string tau_output = data("tau_output", "tau.out");
    std::ofstream tau_output_file(tau_output);
    tau_output_file << data_time.M_time << ", " << ActiveTension::tau << std::endl;

    double rho = data("rho", 1.0);
    assemble_mass(es, data_time.M_dt, rho);

    libMesh::ExodusII_IO exporter(mesh);
    std::string output_file = data("output", "out.e");
    exporter.write_equation_systems(output_file, es);


    bool solve_mechanics = data("solve_mechanics", true);
    bool solve_FE = data("FE", false);

    int exodus_counter = 0;
    for (; data_time.M_iter <= data_time.M_maxIter, data_time.M_time <= data_time.M_endTime;)
    {
        //std::cout << "\nTime: " << data_time.M_time + data_time.M_dt << std::endl;
        ActiveTension::solve_tau(data_time.M_time, data_time.M_dt);
        // ADVANCE:
//        std::cout << "Update" << std::endl;
        *u_system.older_local_solution = *u_system.old_local_solution;
        *v_system.older_local_solution = *v_system.old_local_solution;
        *p_system.older_local_solution = *p_system.old_local_solution;

        *u_system.old_local_solution = *u_system.current_local_solution;
        *v_system.old_local_solution = *v_system.current_local_solution;
        *p_system.old_local_solution = *p_system.current_local_solution;

        // Solve with Heun's method
        if(solve_mechanics)
        {
            //
            int step = 0;
            assemble_RHS(es, data_time.M_time, data_time.M_dt , ActiveTension::tau_n, step);
            v_system.update();

            u_system.solution->add( data_time.M_dt, *v_system.solution);
            u_system.get_vector(step) = *u_system.solution;
            u_system.update();

            p_system.update();

            if(!solve_FE)
            {
                step = 1;
                assemble_RHS(es, data_time.M_time + data_time.M_dt, data_time.M_dt,  ActiveTension::tau, step);
                v_system.update();

                u_system.solution->add( data_time.M_dt, *v_system.solution);
                u_system.get_vector(step) = *u_system.solution;
                u_system.update();

                p_system.update();

                // Average
                step = 2;
                *v_system.solution += *v_system.old_local_solution;
                *v_system.solution *= 0.5;
                v_system.update();

                *u_system.solution += *u_system.old_local_solution;
                *u_system.solution *= 0.5;
                u_system.update();
            }
        }


        data_time.advance();

        if( data_time.M_iter % data_time.M_saveIter == 0)
        {
            std::cout << "Time: " <<  data_time.M_time << " tau = " << ActiveTension::tau << ", tau_n = " << ActiveTension::tau_n <<   std::endl;
            std::cout << "Export solution" << std::endl;
            exodus_counter ++;
            exporter.write_timestep(output_file, es, exodus_counter, data_time.M_time);
        }
        //std::cout << "Write to file" << std::endl;
        //tau_output_file << data_time.M_time << ", " << ActiveTension::tau << std::endl;

    }
    tau_output_file.close();
    //exporter.write_element_data(es);

    return 0;
}

void ActiveTension::solve_tau(double time, double dt)
{
    // Forward Euler Step 1
    tau_n = tau;
    double fn = eval_tau_rhs(time, tau_n);
    double tau_1 = tau_n + dt * fn;
    // Forward Euler Step 2
    double fnp1 = eval_tau_rhs(time + dt, tau_1);
    double tau_2 = tau_1 + dt * fnp1;
    tau = 0.5 * (tau_1 + tau_2);
}

double ActiveTension::eval_tau_rhs(double time, double tau_tn)
{
    double delta_sys = time - t_sys;
    double Sp = S(delta_sys);
    // Use property of tanh as odd function
    double negative_delta_dias = t_dias - time;
    double Sm = S(negative_delta_dias);
    double f = Sp * Sm;
    double a = alpha_max * f + alpha_min * (1 - f);
    double abs_a = std::abs(a);
    double abs_ap = std::max(a, 0.0);
    return -abs_a * tau_tn + sigma0 * abs_ap;
}

double ActiveTension::S(double delta)
{
    return 0.5 * (1 + tanh(delta / gamma));
}


// u1 = un + dt * fn
// u2 = u1 + dt * f1
// unp1 = 0.5 * ( un + u2 )
//      = 0.5 * ( un + u1 + dt * f1)
//      = 0.5 * ( un + un + dt * fn + dt * f1)
//      = un + 0.5 * dt * ( fn + f1 )

void assemble_RHS(libMesh::EquationSystems &es, double tn, double dt, double tau, int step)
{
    // Get stuff for assembly:
    const libMesh::ExplicitSystem &f_system = es.get_system < libMesh::ExplicitSystem > ("fibers");
    const libMesh::ExplicitSystem &s_system = es.get_system < libMesh::ExplicitSystem > ("sheets");
    const libMesh::TransientExplicitSystem &u_system = es.get_system < libMesh::TransientExplicitSystem > ("displacement");
    libMesh::TransientExplicitSystem& v_system = es.get_system < libMesh::TransientExplicitSystem > ("velocity");
    const libMesh::TransientExplicitSystem &p_system = es.get_system < libMesh::TransientExplicitSystem > ("pressure");
    v_system.rhs->close();
    v_system.rhs->zero();

    int dim = 3;
    const libMesh::DofMap &u_dof_map = u_system.get_dof_map();
    const libMesh::DofMap &p_dof_map = p_system.get_dof_map();
    const libMesh::DofMap &f_dof_map = f_system.get_dof_map();
    libMesh::FEType fe_type = u_dof_map.variable_type(0);
    std::unique_ptr < libMesh::FEBase > fe(libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule(dim, libMesh::FIRST);
    fe->attach_quadrature_rule(&qrule);
    std::unique_ptr < libMesh::FEBase > fe_face(libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qface(dim - 1, libMesh::SECOND);
    fe_face->attach_quadrature_rule(&qface);
    const std::vector<libMesh::Real> &JxW = fe->get_JxW();
    const std::vector<libMesh::Point> &q_point = fe->get_xyz();
    const std::vector<std::vector<libMesh::Real>> &phi = fe->get_phi();
    const std::vector<std::vector<libMesh::RealGradient>> &dphi = fe->get_dphi();
    const std::vector<libMesh::Real> &JxW_face = fe_face->get_JxW();
    const std::vector<libMesh::Point> &q_point_face = fe_face->get_xyz();
    const std::vector<std::vector<libMesh::Real>> &phi_face = fe_face->get_phi();
    const std::vector<libMesh::Point> &N = fe_face->get_normals();

    std::vector < libMesh::dof_id_type > u_dof_indices;
    std::vector < libMesh::dof_id_type > p_dof_indices;
    std::vector < libMesh::dof_id_type > f_dof_indices;
    libMesh::TensorValue<double> I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    libMesh::TensorValue<double> F;
    libMesh::TensorValue<double> P;
    libMesh::TensorValue<double> P_visc;
    libMesh::TensorValue<double> P_pass;
    libMesh::TensorValue<double> P_act;
    libMesh::TensorValue<double> P_vol;
    libMesh::TensorValue<double> FinvT;
    libMesh::TensorValue<double> Grad_u;
    libMesh::TensorValue<double> Fdot;
    libMesh::TensorValue<double> FdotT;
    libMesh::TensorValue<double> Grad_test;
    // 2.0 * f \otimes f0
    libMesh::TensorValue<double> dI4f;
    // 2.0 * s \otimes s0
    libMesh::TensorValue<double> dI4s;
    // f \otimes s0 + s \otimes f0
    libMesh::TensorValue<double> fos0;
    libMesh::TensorValue<double> sof0;
    libMesh::TensorValue<double> dI8fs;

    libMesh::VectorValue<double> f0;
    libMesh::VectorValue<double> s0;
    libMesh::VectorValue<double> f;
    libMesh::VectorValue<double> s;
    libMesh::VectorValue<double> un;
    libMesh::VectorValue<double> vn;
    double pn;

    std::vector<double> f_vec;
    std::vector<double> s_vec;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> p;

    libMesh::DenseVector < libMesh::Number > Fe;

    // Boundary conditions:
    std::map<int, Boundary> boundary_condition_map;
    boundary_condition_map[1] = Boundary::ENDO;
    boundary_condition_map[2] = Boundary::EPI;
    boundary_condition_map[3] = Boundary::BASE;
    //std::cout << "Loop over elem" << std::endl;

    for (const auto &elem : es.get_mesh().active_local_element_ptr_range())
    {
        u_dof_map.dof_indices(elem, u_dof_indices);
        p_dof_map.dof_indices(elem, p_dof_indices);
        f_dof_map.dof_indices(elem, f_dof_indices);
        const unsigned int n_dofs = u_dof_indices.size();
        fe->reinit(elem);

        Fe.resize(n_dofs);

        u_system.solution->get(u_dof_indices, u);
        v_system.solution->get(u_dof_indices, v);
        p_system.solution->get(p_dof_indices, p);

        f_system.solution->get(f_dof_indices, f_vec);
        for (int i = 0; i < 3; ++i) f0(i) = f_vec[i];
        s_system.solution->get(f_dof_indices, s_vec);
        for (int i = 0; i < 3; ++i) s0(i) = s_vec[i];

        auto n_qp = qrule.n_points();
        auto n_phi = phi.size();

        for (unsigned int qp = 0; qp < n_qp; qp++)
        {
            // Eval F
            F *= 0.0;
            FinvT *= 0.0;
            Grad_u *= 0.0;
            Fdot *= 0;
            P_pass *= 0.0;
            P_act *= 0.0;
            P_visc *= 0.0;
            P *= 0.0;
            un *= 0.0;
            vn *= 0.0;
            pn = 0;
            for (int i = 0; i < n_phi; ++i)
            {
                pn += p[i] * phi[i][qp];
                for (int idim = 0; idim < 3; ++idim)
                {
                    un(idim) += u[idim * n_phi + i] * phi[i][qp];
                    vn(idim) += v[idim * n_phi + i] * phi[i][qp];
                    for (int jdim = 0; idim < 3; ++idim)
                    {
                        Grad_u(idim, jdim) += u[idim * n_phi + i] * dphi[i][qp](jdim);
                        // Grad v
                        Fdot(idim, jdim) += v[idim * n_phi + i] * dphi[i][qp](jdim);
                    }
                }
            }
            //std::cout << "Eval F" << std::endl;

            F = I + Grad_u;


            double J = F.det();
            //std::cout << "Eval FinvT, J: " << J << std::endl;
            //F.print();
            FinvT = F.inverse().transpose();


            double mu = 1;
            double J2 = J*J;
            double Jm23 = 1.0 / std::cbrt(J2);
            P_act = 0.5 * tau * dI4f;

            double eta = 1e3;
            P_visc = 0.5 * eta * F * (Fdot.transpose() * F + F.transpose() * Fdot);

            // U = k / 4 (J^2 - 1 - 2 ln J)
            //dU = k/2/J ( J*J -  1)
            double kappa = 1e1;
            double dU = 0.5 * kappa / J * ( J2 - 1.0);

            P_vol = J * dU * FinvT;
            //
            double Pa2Dyne = 10.0;
            double a = 59.0 * Pa2Dyne;
            double af = 18472.0 * Pa2Dyne;
            double as = 2481.0 * Pa2Dyne;
            double afs = 216.0 * Pa2Dyne;
            double b = 8.023;
            double bf = 16.026;
            double bs = 11.12;
            double bfs = 11.436;

            // Isotropic
            double I1 = F.contract(F);
            double dW1bar = 0.5 * a * std::exp( b * (Jm23 * I1 - 3.0) );
            P_pass = dW1bar * Jm23 * F;
            // Fibers
            // eval f \otimes f
            f = F * f0;
            dyadic(2.0 * f, f0, dI4f);
            double I4f = f.contract(f);

            // heavyside function
            double k = 100.0;
            double auxf = (1. + std::exp(-k * (I4f - 1.)));
            double Chi_I4f = 1. / auxf;
            double dChi_i4f =  k * (auxf - 1.0 )  / ( auxf * auxf );
            double expI4f = exp(bf * (I4f - 1.0) * (I4f-1.0) );
            double dW4f = af * Chi_I4f * (I4f - 1.0) * expI4f + 0.5 * af / bf * ( expI4f - 1.0 ) * dChi_i4f;
            P_pass += dW4f * dI4f;

            // Sheets
            // eval s \otimes s
            s = F * s0;
            dyadic(2.0 * s, s0, dI4s);
            double I4s = s.contract(s);
            double auxs = (1. + std::exp(-k * (I4s - 1.)));
            double Chi_I4s = 1. / auxs;
            double dChi_i4s =  k * (auxs - 1.0 )  / ( auxs * auxs );
            double expI4s = exp(bs * (I4s - 1.0) * (I4s-1.0) );
            double dW4s = as * Chi_I4s * (I4s - 1.0) * expI4s + 0.5 * as / bs * ( expI4s - 1.0 ) * dChi_i4s;
            P_pass += dW4s * dI4s;

            // fs
            // eval s \otimes s
            dyadic(f, s0, fos0);
            dyadic(s, f0, sof0);
            dI8fs = fos0 + sof0;
            double I8fs = f.contract(s);
            double dW8fs = afs * I8fs * exp(bfs * I8fs * I8fs);
            P_pass += dW8fs * dI8fs;

            // sigma = J P F^T
            // dev[sigma] = sigma - sigma : I / 3
            // dev[sigma] F^-T/J = sigma F^-T/J - sigma : I / 3 F^-T/J
            // dev[sigma] F^-T/J = P - J (P : F) / 3 F^-T/J
            // dev[sigma] F^-T/J = P - (P : F) / 3 * F^-T
            P = P_pass + P_visc + P_act;
            P -= P.contract(F) / 3.0 * FinvT;
            P += P_vol;

            for (int i = 0; i < n_phi; ++i)
            {
                for (int jdim = 0; jdim < dim; jdim++)
                {
                    Grad_test *= 0.0;
                    Grad_test(jdim, 0) = JxW[qp] * dphi[i][qp](0);
                    Grad_test(jdim, 1) = JxW[qp] * dphi[i][qp](1);
                    Grad_test(jdim, 2) = JxW[qp] * dphi[i][qp](2);
                    // Compute  - \nabla \cdot \sigma + f
                    Fe(i + jdim * n_phi) -= P.contract(Grad_test);
                }
            }
        } // end elemental qp loop
          // apply BC
        for (auto side : elem->side_index_range())
        {
            if (elem->neighbor_ptr(side) == nullptr)
            {
                // sidesets:
                // 1 = endo
                // 2 = epi
                // 3 = base
                auto sideset_id = es.get_mesh().boundary_info->boundary_id(elem, side);
                auto it = boundary_condition_map.find(sideset_id);
                if (it != boundary_condition_map.end())
                {
                    fe_face->reinit(elem, side);
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        for (int i = 0; i < n_phi; ++i)
                        {
                            for (int idim = 0; idim < 3; ++idim)
                            {
                                un(idim) += u[idim * n_phi + i] * phi_face[i][qp];
                                vn(idim) += v[idim * n_phi + i] * phi_face[i][qp];
                            }
                        }
                        switch (it->second)
                        {
                            case Boundary::EPI:
                            {
                                double alpha = 0.0;//1e7;
                                double beta = 0.0;//5e2;
                                for (int idim = 0; idim < dim; ++idim)
                                {
                                    for (unsigned int i = 0; i < n_phi; i++)
                                    {
                                        Fe(i + idim * n_phi) -= JxW_face[qp] * (alpha * un(idim) * N[qp](idim) + beta * vn(idim) * N[qp](idim)) * phi_face[i][qp];
                                    }
                                }
                                break;
                            }
                            case Boundary::BASE:
                            {
                                double alpha = 1e14;// 1e4;
                                double beta = 0.0; //5e2;
                                for (int idim = 0; idim < dim; ++idim)
                                {
                                    for (unsigned int i = 0; i < n_phi; i++)
                                    {
                                        Fe(i + idim * n_phi) -= JxW_face[qp] * (alpha * un(idim) + beta * vn(idim)) * phi_face[i][qp];
                                    }
                                }
                                break;
                            }
                            default:
                            {
                                break;
                            }
                        }
                    }
                } // end loop over face qp
            } // end if side == nullptr
        } // end loop over sides

        v_system.rhs->add_vector(Fe, u_dof_indices);
    }

    //std::cout << "Assemble RHS done" << std::endl;
    v_system.get_vector(step).close();
//    double mass_norm =  v_system.get_vector("mass").l2_norm();
//    std::cout << "mass norm: " << mass_norm << std::endl;
    ( *v_system.rhs ) /= v_system.get_vector("mass");
    ( *v_system.solution ) += ( *v_system.rhs );
    v_system.get_vector(step) = (*v_system.solution);
}

void assemble_mass(libMesh::EquationSystems &es, double dt, double rho)
{
    std::cout << "Assemble Lumped Mass Matrix Vectors" << std::endl;
    // Get stuff for assembly:
    libMesh::TransientExplicitSystem &v_system = es.get_system < libMesh::TransientExplicitSystem > ("velocity");
    libMesh::TransientExplicitSystem &p_system = es.get_system < libMesh::TransientExplicitSystem > ("pressure");

    int dim = 3;
    const libMesh::DofMap &v_dof_map = v_system.get_dof_map();
    const libMesh::DofMap &p_dof_map = p_system.get_dof_map();
    libMesh::FEType fe_type = v_dof_map.variable_type(0);
    std::unique_ptr < libMesh::FEBase > fe(libMesh::FEBase::build(dim, fe_type));
    libMesh::QGauss qrule(dim, libMesh::SECOND);
    fe->attach_quadrature_rule(&qrule);
    const std::vector<libMesh::Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<libMesh::Real>> &phi = fe->get_phi();

    std::vector < libMesh::dof_id_type > v_dof_indices;
    std::vector < libMesh::dof_id_type > p_dof_indices;

    libMesh::DenseVector < libMesh::Number > Fe_v;
    libMesh::DenseVector < libMesh::Number > Fe_p;

    for (const auto &elem : es.get_mesh().active_local_element_ptr_range())
    {
        v_dof_map.dof_indices(elem, v_dof_indices);
        p_dof_map.dof_indices(elem, p_dof_indices);
        const unsigned int v_n_dofs = v_dof_indices.size();
        const unsigned int p_n_dofs = p_dof_indices.size();
        fe->reinit(elem);
        Fe_v.resize(v_n_dofs);
        Fe_p.resize(p_n_dofs);

        auto n_qp = qrule.n_points();
        auto n_phi = phi.size();

        for (unsigned int qp = 0; qp < n_qp; qp++)
        {
            for (unsigned int i = 0; i < n_phi; ++i)
            {
                for (unsigned int j = 0; j < n_phi; ++j)
                {
                    double ml = JxW[qp] * rho / dt * phi[i][qp] * phi[j][qp];
                    Fe_p(i) += ml;
                    Fe_v(i) += ml;
                    Fe_v(i + n_phi) += ml;
                    Fe_v(i + 2 * n_phi) += ml;
                }
            }
        }
        v_system.get_vector("mass").add_vector(Fe_v, v_dof_indices);
        p_system.get_vector("mass").add_vector(Fe_p, p_dof_indices);
    }
    v_system.get_vector("mass").close();
    p_system.get_vector("mass").close();
}
