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
 * \file main.cpp
 *
 * \brief Here we test GetPot
 *
 * \author srossi
 *
 * \version 0.0
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 7, 2016
 *
 */

#include "libmesh/getpot.h"
#include "Util/IO/io.hpp"
// Functions to initialize the library.
#include "libmesh/libmesh.h"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/vector_value.h"

//#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"

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
#include "libmesh/linear_solver.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"

double exact_solution(double mu, double alpha, double x, double time,
		int model);
double exact_derivative(double mu, double alpha, double x, double time,
		int model);
void set_exact_ic(double mu, double alpha, double time,
		libMesh::EquationSystems & es, int model);
void set_exact_solution(double mu, double alpha, double time,
		libMesh::EquationSystems & es, int model);
void eval_l2_err(double& l2_V_err, double& l2_Q_err, double dt, double mu,
		double alpha, double time, libMesh::EquationSystems & es, int model);

void solve(libMesh::EquationSystems & es, int model, bool compute_q = false,
		int order = 2);
void computeIonicCurrent(int model, double V, double Q, double h, double alpha,
		double& I, double& dI, bool is_hyperbolic);

int main(int argc, char ** argv)
{
	// Bring in everything from the libMesh namespace

	using namespace libMesh;
	// Initialize libraries, like in example 2.
	LibMeshInit init(argc, argv, MPI_COMM_WORLD);

	// Use the MeshTools::Generation mesh generator to create a uniform
	// 3D grid
	// We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
	// the number of elements on each side is read from the input file
	// Create a mesh, with dimension to be overridden later, on the
	// default MPI communicator.
	libMesh::Mesh mesh(init.comm());

	GetPot commandLine(argc, argv);
	std::string datafile_name = commandLine.follow("data.beat", 2, "-i",
			"--input");
	GetPot data(datafile_name);
	// allow us to use higher-order approximation.
	int numElementsX = data("mesh/elX", 1600);

	int model = data("model", 0); // 0: cubic; 1: McKean
	bool compute_q = data("q", false); // 0: cubic; 1: McKean
	int order = data("ti_order", 2); // Time integrator order (1 or 2)

	MeshTools::Generation::build_line(mesh, numElementsX, -25.0, 25.0);

	// Print information about the mesh to the screen.
	mesh.print_info();

	libMesh::EquationSystems systems(mesh);

	double mu = data("mu", 0.2);
	double alpha = data("alpha", 0.1);
	int ms = data("ms", 10);
	double dt = data("dt", 0.003125);
	typedef libMesh::Real Real;
	systems.parameters.set<Real>("mu") = mu;
	systems.parameters.set<Real>("alpha") = alpha;
	systems.parameters.set<Real>("dt") = dt;
	systems.parameters.set<Real>("dummy") = 42.0;
	systems.parameters.set<Real>("nobody") = 0.0;

	systems.parameters.set<bool>("test") = true;
	typedef libMesh::Real Real;
	systems.parameters.set<Real>("dummy") = 42.0;
	systems.parameters.set<Real>("nobody") = 0.0;

	libMesh::TransientLinearImplicitSystem& Qesys = systems.add_system<
			libMesh::TransientLinearImplicitSystem>("ExactQSystem");
	unsigned int Qe_var = Qesys.add_variable("Qexact", libMesh::FIRST);
	libMesh::TransientLinearImplicitSystem& Vesys = systems.add_system<
			libMesh::TransientLinearImplicitSystem>("ExactVSystem");
	unsigned int Ve_var = Vesys.add_variable("Vexact", libMesh::FIRST);

	libMesh::TransientLinearImplicitSystem& Vsys = systems.add_system<
			libMesh::TransientLinearImplicitSystem>("VSystem");
	unsigned int V_var = systems.get_system("VSystem").add_variable("V",
			libMesh::FIRST);

	Vsys.add_matrix("ML");
	Vsys.add_matrix("M");
	Vsys.add_matrix("K");

	Vsys.add_vector("k1");
	Vsys.add_vector("k2");
	Vsys.add_vector("Vn");
	Vsys.add_vector("V1");
	Vsys.add_vector("V2");
	Vsys.add_vector("I1");
	Vsys.add_vector("I2");
	Vsys.add_vector("MI1");
	Vsys.add_vector("MI2");
	Vsys.add_vector("MVn");
	Vsys.add_vector("MV1");
	Vsys.add_vector("MV2");
	Vsys.add_vector("KV");
	Vsys.add_vector("KV1");
	Vsys.add_vector("KV2");

	libMesh::TransientLinearImplicitSystem& Qsys = systems.add_system<
			libMesh::TransientLinearImplicitSystem>("QSystem");
	unsigned int Q_var = Qsys.add_variable("Q", libMesh::FIRST);

	Qsys.add_matrix("QML");
	Qsys.add_matrix("QM");
	Qsys.add_matrix("QK");

	Qsys.add_vector("Qn");
	Qsys.add_vector("Q1");
	Qsys.add_vector("Q2");
	Qsys.add_vector("Qk1");
	Qsys.add_vector("Qk2");
	Qsys.add_vector("QVn");
	Qsys.add_vector("QV1");
	Qsys.add_vector("QV2");
	Qsys.add_vector("QI1");
	Qsys.add_vector("QI2");
	Qsys.add_vector("QdI1");
	Qsys.add_vector("QdI2");
	Qsys.add_vector("QMI1");
	Qsys.add_vector("QMI2");
	Qsys.add_vector("QMdI1");
	Qsys.add_vector("QMdI2");
	Qsys.add_vector("QMVn");
	Qsys.add_vector("QMQn");
	Qsys.add_vector("QMQ2");
	Qsys.add_vector("QMV1");
	Qsys.add_vector("QMV2");
	Qsys.add_vector("QKV");
	Qsys.add_vector("QKV1");
	Qsys.add_vector("QKV2");
	Qsys.add_vector("QKQ1");
	Qsys.add_vector("QKQ2");

	systems.init();

	systems.print_info();
	double& time = systems.get_system("VSystem").time;
	systems.get_system("VSystem").time = 0.;

	double l2_V_err = 0.0;
	double l2_Q_err = 0.0;
	set_exact_ic(mu, alpha, time, systems, model);

	set_exact_solution(mu, alpha, time, systems, model);
	std::string err_output = "error_" + std::to_string(numElementsX) + "_mu_"
			+ std::to_string(mu) + "_alpha_" + std::to_string(alpha) + "_model_"
			+ std::to_string(model) + ".txt";
	std::ofstream errfile;
	errfile.open(err_output);
	errfile << "Time  err_V  err_Q\n";

	// TIME LOOP
	std::string exodus_filename = data("exofile", "transient_ex1.e");
	libMesh::ExodusII_IO exo(mesh);

	exo.write_equation_systems(exodus_filename, systems);
	exo.append(true);

	unsigned int t_step = 1;
	exo.write_timestep(exodus_filename, systems, t_step, time);
	for (; systems.get_system("VSystem").time < 1.0;)
	{
		std::cout << "step: " << t_step << std::endl;
		t_step++;
		systems.get_system("VSystem").time += dt;

		solve(systems, model, compute_q, order);

		set_exact_solution(mu, alpha, time, systems, model);
		eval_l2_err(l2_V_err, l2_Q_err, dt, mu, alpha, time, systems, model);
		errfile << std::fixed << std::setprecision(15) << time << " "
				<< l2_V_err << " " << l2_Q_err << std::endl;

		exo.write_timestep(exodus_filename, systems, t_step, time);

	}
	return 0;
}

void solve(libMesh::EquationSystems & es, int model, bool compute_q, int order)
{
	using std::unique_ptr;

	const libMesh::MeshBase & mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	libMesh::TransientLinearImplicitSystem & Vsys = es.get_system<
			libMesh::TransientLinearImplicitSystem>("VSystem");
	const unsigned int V_Var = Vsys.variable_number("V");
	libMesh::TransientLinearImplicitSystem & Qsys = es.get_system<
			libMesh::TransientLinearImplicitSystem>("QSystem");
	const unsigned int Q_Var = Qsys.variable_number("Q");
	const libMesh::DofMap & dof_map = Vsys.get_dof_map();
	libMesh::FEType fe_type = dof_map.variable_type(0);

	double mu = es.parameters.get<double>("mu");
	bool is_hyperbolic = (mu > 0.0) ? true : false;
	double dt = es.parameters.get<double>("dt");
	double time = es.get_system("VSystem").time;
	double alpha = 0.1;

	std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));
	libMesh::QGauss qrule(dim, libMesh::SECOND);
	fe->attach_quadrature_rule(&qrule);

	const std::vector<libMesh::Real> & JxW = fe->get_JxW();
	const std::vector<libMesh::Point> & q_point = fe->get_xyz();
	const std::vector<std::vector<libMesh::Real> > & phi = fe->get_phi();
	const std::vector<std::vector<libMesh::RealGradient> > & dphi =
			fe->get_dphi();

	libMesh::DenseMatrix<libMesh::Number> Ke;
	libMesh::DenseMatrix<libMesh::Number> Me;
	libMesh::DenseMatrix<libMesh::Number> MLe;

	std::vector<libMesh::dof_id_type> dof_indices;

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	// Loop over the elements.

	Vsys.get_matrix("M").zero();
	Vsys.get_matrix("ML").zero();
	Vsys.get_matrix("K").zero();
	double h = 0;
	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		h = elem->hmax();
		dof_map.dof_indices(elem, dof_indices);
		fe->reinit(elem);
		MLe.resize(dof_indices.size(), dof_indices.size());
		Me.resize(dof_indices.size(), dof_indices.size());
		Ke.resize(dof_indices.size(), dof_indices.size());

		for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
		{
			for (unsigned int i = 0; i < phi.size(); i++)
			{
				for (unsigned int j = 0; j < phi.size(); j++)
				{
					// Mass term
					Me(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
					MLe(i, i) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
					// stiffness term
					Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
				}
			}
		}

		Vsys.get_matrix("M").add_matrix(Me, dof_indices);
		Vsys.get_matrix("ML").add_matrix(MLe, dof_indices);
		Vsys.get_matrix("K").add_matrix(Ke, dof_indices);
	}
	Vsys.get_matrix("M").close();
	Vsys.get_matrix("ML").close();
	Vsys.get_matrix("K").close();

	Qsys.get_matrix("QM").zero();
	Qsys.get_matrix("QML").zero();
	Qsys.get_matrix("QK").zero();
	Qsys.get_matrix("QM").close();
	Qsys.get_matrix("QML").close();
	Qsys.get_matrix("QK").close();
	Qsys.get_matrix("QM").add(1.0, Vsys.get_matrix("M"));
	Qsys.get_matrix("QML").add(1.0, Vsys.get_matrix("ML"));
	Qsys.get_matrix("QK").add(1.0, Vsys.get_matrix("K"));

	if (compute_q)
	{

		Qsys.matrix->zero();
		Qsys.matrix->close();

		// First step A x = b
		// Compute matrix: A

		// A = (mu+dt/2) * M + dt^2/4 * K
		if (order == 2)
		{
			if (is_hyperbolic)
			{
				Qsys.matrix->add(mu + 0.5 * dt, Qsys.get_matrix("QML"));
				Qsys.matrix->add(0.25 * dt * dt, Qsys.get_matrix("QK"));
			}
			// A = M + dt/2 * K
			else
			{
				Qsys.matrix->add(1.0, Qsys.get_matrix("QML"));
				Qsys.matrix->add(0.5 * dt, Qsys.get_matrix("QK"));
			}
		}
		else
		{
			if (is_hyperbolic)
			{
				Qsys.matrix->add(mu + dt, Qsys.get_matrix("QML"));
				Qsys.matrix->add(dt * dt, Qsys.get_matrix("QK"));
			}
			// A = M + dt * K
			else
			{
				Qsys.matrix->add(1.0, Qsys.get_matrix("QML"));
				Qsys.matrix->add(dt, Qsys.get_matrix("QK"));
			}
		}

		Qsys.get_vector("QKV1").zero();
		Qsys.get_vector("QKV1").close();
		Qsys.get_vector("QKQ1").zero();
		Qsys.get_vector("QKQ1").close();
		Qsys.get_vector("QMQn").zero();
		Qsys.get_vector("QMQn").close();
		Qsys.get_vector("QVn") = *Vsys.current_local_solution;
		Qsys.get_vector("QV1") = *Vsys.current_local_solution;
		Qsys.get_vector("Qn") = *Qsys.current_local_solution;
		Qsys.get_vector("Q1") = *Qsys.current_local_solution;
		Qsys.get_matrix("QK").vector_mult_add(Qsys.get_vector("QKV1"),
				Qsys.get_vector("QV1"));
		Qsys.get_matrix("QK").vector_mult_add(Qsys.get_vector("QKQ1"),
				Qsys.get_vector("Q1"));
		Qsys.get_matrix("QML").vector_mult_add(Qsys.get_vector("QMQn"),
				Qsys.get_vector("Q1"));

		Qsys.get_vector("QI1").zero();
		Qsys.get_vector("QI1").close();
		Qsys.get_vector("QdI1").zero();
		Qsys.get_vector("QdI1").close();
		libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
		const libMesh::MeshBase::const_node_iterator end_node =
				mesh.local_nodes_end();
		for (; node != end_node; ++node)
		{
			const libMesh::Node * nn = *node;
			dof_map.dof_indices(nn, dof_indices, 0);
			double V1 = (*Vsys.current_local_solution)(dof_indices[0]);
			double Q1 = (*Qsys.current_local_solution)(dof_indices[0]);
			double I1 = 0.0;
			double dI1 = 0.0;
			computeIonicCurrent(model, V1, Q1, h, alpha, I1, dI1,
					is_hyperbolic);
			Qsys.get_vector("QI1").set(dof_indices[0], I1);
			Qsys.get_vector("QdI1").set(dof_indices[0], dI1);
		}

		Qsys.get_vector("QI1").close();
		Qsys.get_vector("QMI1").zero();
		Qsys.get_vector("QMI1").close();
		Qsys.get_matrix("QM").vector_mult_add(Qsys.get_vector("QMI1"),
				Qsys.get_vector("QI1"));
		Qsys.get_vector("QdI1").close();
		Qsys.get_vector("QMdI1").zero();
		Qsys.get_vector("QMdI1").close();
		Qsys.get_matrix("QM").vector_mult_add(Qsys.get_vector("QMdI1"),
				Qsys.get_vector("QdI1"));

		Qsys.get_vector("Qk1").zero();
		Qsys.get_vector("Qk1").close();
		// RHS: ( mu - dt /2 ) * ML * Qn -dt * K * Vn - dt^2 / 4 * K * Qn - dt * M * (In + mu * dIn)
		if (order == 2)
		{
			if (is_hyperbolic)
			{
				Qsys.get_vector("Qk1").add(mu - 0.5 * dt,
						Qsys.get_vector("QMQn"));
				Qsys.get_vector("Qk1").add(-dt, Qsys.get_vector("QKV1"));
				Qsys.get_vector("Qk1").add(-0.25 * dt * dt,
						Qsys.get_vector("QKQ1"));
				Qsys.get_vector("Qk1").add(-mu * dt, Qsys.get_vector("QMdI1"));
				Qsys.get_vector("Qk1").add(-dt, Qsys.get_vector("QMI1"));
			}
			// RHS: ML * Qn  - dt / 2 * K * Qn - dt * M * dIn
			else
			{
				Qsys.get_vector("Qk1").add(1.0, Qsys.get_vector("QMQn"));
				Qsys.get_vector("Qk1").add(-0.5 * dt, Qsys.get_vector("QKQ1"));
				Qsys.get_vector("Qk1").add(-1.0 * dt, Qsys.get_vector("QMdI1"));
			}
		}
		else
		{
			if (is_hyperbolic)
			{
				Qsys.get_vector("Qk1").add(mu, Qsys.get_vector("QMQn"));
				Qsys.get_vector("Qk1").add(-dt, Qsys.get_vector("QKV1"));
				Qsys.get_vector("Qk1").add(-mu * dt, Qsys.get_vector("QMdI1"));
				Qsys.get_vector("Qk1").add(-dt, Qsys.get_vector("QMI1"));
			}
			// RHS: ML * Qn  - dt * M * dIn
			else
			{
				Qsys.get_vector("Qk1").add(1.0, Qsys.get_vector("QMQn"));
				Qsys.get_vector("Qk1").add(-1.0 * dt, Qsys.get_vector("QMdI1"));
			}
		}
		Qsys.rhs->zero();
		Qsys.rhs->close();
		Qsys.rhs->add(1.0, Qsys.get_vector("Qk1"));

		double tol = 1e-12;
		double max_iter = 2000;

		std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);

		Qsys.get_vector("Q2").zero();
		Qsys.get_vector("Q2").close();
		std::unique_ptr < libMesh::LinearSolver<libMesh::Number>
				> linearSolver;
		linearSolver = libMesh::LinearSolver<libMesh::Number>::build(es.comm());
		linearSolver->init();

		rval = linearSolver->solve(*Qsys.matrix, nullptr, Qsys.get_vector("Q2"),
				*Qsys.rhs, tol, max_iter);

		Vsys.get_vector("V2").zero();
		Vsys.get_vector("V2").close();
		if (order == 2)
		{
			Vsys.get_vector("V2").add(1.0, *Vsys.current_local_solution);
			Vsys.get_vector("V2").add(0.5 * dt, *Qsys.current_local_solution);
			Vsys.get_vector("V2").add(0.5 * dt, Qsys.get_vector("Q2"));
		}
		else
		{
			Vsys.get_vector("V2").add(1.0, *Vsys.current_local_solution);
			Vsys.get_vector("V2").add(dt, Qsys.get_vector("Q2"));

			*Vsys.current_local_solution = Vsys.get_vector("V2");
			*Qsys.current_local_solution = Qsys.get_vector("Q2");

			*Qsys.solution = *Qsys.current_local_solution;
			Qsys.update();
			*Qsys.old_local_solution = *Qsys.solution;

			*Vsys.solution = *Vsys.current_local_solution;
			Vsys.update();
			*Vsys.old_local_solution = *Vsys.solution;
			return;

		}
		// SECOND STEP

		Qsys.get_vector("QI2").zero();
		Qsys.get_vector("QI2").close();
		Qsys.get_vector("QdI2").zero();
		Qsys.get_vector("QdI2").close();
		node = mesh.local_nodes_begin();
		for (; node != end_node; ++node)
		{
			const libMesh::Node * nn = *node;
			dof_map.dof_indices(nn, dof_indices, 0);
			double V2 = Vsys.get_vector("V2")(dof_indices[0]);
			double I2 = 0.0;
			double Q2 = Qsys.get_vector("Q2")(dof_indices[0]);
			double dI2 = 0.0;
			computeIonicCurrent(model, V2, Q2, h, alpha, I2, dI2,
					is_hyperbolic);
			Qsys.get_vector("QI2").set(dof_indices[0], I2);
			Qsys.get_vector("QdI2").set(dof_indices[0], dI2);
		}

		Qsys.get_vector("QI2").close();
		Qsys.get_vector("QMI2").zero();
		Qsys.get_vector("QMI2").close();
		Qsys.get_matrix("QM").vector_mult_add(Qsys.get_vector("QMI2"),
				Qsys.get_vector("QI2"));
		Qsys.get_vector("QdI2").close();
		Qsys.get_vector("QMdI2").zero();
		Qsys.get_vector("QMdI2").close();
		Qsys.get_matrix("QM").vector_mult_add(Qsys.get_vector("QMdI2"),
				Qsys.get_vector("QdI2"));

		Qsys.get_vector("QKV2").zero();
		Qsys.get_vector("QKV2").close();
		Qsys.get_matrix("QK").vector_mult_add(Qsys.get_vector("QKV2"),
				Vsys.get_vector("V2"));
		Qsys.get_vector("QKQ2").zero();
		Qsys.get_vector("QKQ2").close();
		Qsys.get_matrix("QK").vector_mult_add(Qsys.get_vector("QKQ2"),
				Qsys.get_vector("Q2"));

		Qsys.get_vector("QMQ2").zero();
		Qsys.get_vector("QMQ2").close();
		Qsys.get_matrix("QML").vector_mult_add(Qsys.get_vector("QMQ2"),
				Qsys.get_vector("Q2"));

		Qsys.get_vector("Qk2").zero();
		Qsys.get_vector("Qk2").close();
		// RHS: ML * Qn - dt / mu * M * ( Qn + Q2 ) / 2 - dt / mu * K * (Vn + V2) / 2 - dt / mu * M * ( In + I2 + mu * ( dIn + dI2 ) ) / 2
		if (is_hyperbolic)
		{
			// ML * Qn
			Qsys.get_vector("Qk2").add(1.0, Qsys.get_vector("QMQn"));
			// - dt * M * ( Qn + Q2 ) / 2
			Qsys.get_vector("Qk2").add(-0.5 * dt / mu, Qsys.get_vector("QMQn"));
			Qsys.get_vector("Qk2").add(-0.5 * dt / mu, Qsys.get_vector("QMQ2"));
			// - dt * K * (Vn + V2) / 2
			Qsys.get_vector("Qk2").add(-0.5 * dt / mu, Qsys.get_vector("QKV1"));
			Qsys.get_vector("Qk2").add(-0.5 * dt / mu, Qsys.get_vector("QKV2"));
			// - dt * M * ( In + I2 + mu * ( dIn + dI2 ) ) / 2
			Qsys.get_vector("Qk2").add(-0.5 * dt / mu, Qsys.get_vector("QMI1"));
			Qsys.get_vector("Qk2").add(-0.5 * dt / mu, Qsys.get_vector("QMI2"));
			Qsys.get_vector("Qk2").add(-0.5 * mu * dt / mu,
					Qsys.get_vector("QMdI1"));
			Qsys.get_vector("Qk2").add(-0.5 * mu * dt / mu,
					Qsys.get_vector("QMdI2"));

		}
		// RHS: ML * Qn - dt * K * ( Qn + Q2 ) / 2 - dt * M * ( dIn + dI2 ) / 2
		else
		{
			Qsys.get_vector("Qk2").add(1.0, Qsys.get_vector("QMQn"));
			Qsys.get_vector("Qk2").add(-0.5 * dt, Qsys.get_vector("QKQ1"));
			Qsys.get_vector("Qk2").add(-0.5 * dt, Qsys.get_vector("QKQ2"));
			Qsys.get_vector("Qk2").add(-0.5 * dt, Qsys.get_vector("QMdI1"));
			Qsys.get_vector("Qk2").add(-0.5 * dt, Qsys.get_vector("QMdI2"));
		}

		Qsys.rhs->zero();
		Qsys.rhs->close();
		Qsys.rhs->add(1.0, Qsys.get_vector("Qk2"));

		std::unique_ptr < libMesh::LinearSolver<libMesh::Number>
				> linearSolver2;
		linearSolver2 = libMesh::LinearSolver<libMesh::Number>::build(
				es.comm());
		linearSolver2->init();
		rval = linearSolver2->solve(Qsys.get_matrix("QML"), nullptr,
				Qsys.get_vector("Q2"), *Qsys.rhs, tol, max_iter);

		Vsys.get_vector("V2").zero();
		Vsys.get_vector("V2").close();
		Vsys.get_vector("V2").add(1.0, *Vsys.current_local_solution);
		Vsys.get_vector("V2").add(0.5 * dt, *Qsys.current_local_solution);
		Vsys.get_vector("V2").add(0.5 * dt, Qsys.get_vector("Q2"));

		*Vsys.current_local_solution = Vsys.get_vector("V2");
		*Qsys.current_local_solution = Qsys.get_vector("Q2");

		*Qsys.solution = *Qsys.current_local_solution;
		Qsys.update();
		*Qsys.old_local_solution = *Qsys.solution;

		*Vsys.solution = *Vsys.current_local_solution;
		Vsys.update();
		*Vsys.old_local_solution = *Vsys.solution;
	}
	else
	{
		Vsys.matrix->zero();
		Vsys.matrix->close();
		Vsys.matrix->add(1.0, Vsys.get_matrix("ML"));
		Vsys.matrix->add(0.5 * dt, Vsys.get_matrix("K"));

		Vsys.get_vector("KV1").zero();
		Vsys.get_vector("KV1").close();
		Vsys.get_vector("Vn") = *Vsys.current_local_solution;
		Vsys.get_vector("V1") = *Vsys.current_local_solution;
		Vsys.get_matrix("K").vector_mult_add(Vsys.get_vector("KV1"),
				Vsys.get_vector("V1"));

		Vsys.get_vector("MVn").zero();
		Vsys.get_vector("MVn").close();
		Vsys.get_matrix("ML").vector_mult_add(Vsys.get_vector("MVn"),
				Vsys.get_vector("Vn"));

		Vsys.get_vector("I1").zero();
		Vsys.get_vector("I1").close();
		libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
		const libMesh::MeshBase::const_node_iterator end_node =
				mesh.local_nodes_end();
		for (; node != end_node; ++node)
		{
			const libMesh::Node * nn = *node;
			dof_map.dof_indices(nn, dof_indices, 0);
			double V1 = (*Vsys.current_local_solution)(dof_indices[0]);
			double I1 = 0.0;
			if (0 == model)
				I1 = V1 * (V1 - alpha) * (V1 - 1.0);
			else
			{
				I1 = V1;
				if (V1 > alpha)
					I1 -= 1.0;
			}
			double dI1 = 0.0;
			double Q1 = 0;
			computeIonicCurrent(model, V1, Q1, h, alpha, I1, dI1,
					is_hyperbolic);
			Vsys.get_vector("I1").set(dof_indices[0], I1);
		}

		Vsys.get_vector("I1").close();
		Vsys.get_vector("MI1").zero();
		Vsys.get_vector("MI1").close();
		Vsys.get_matrix("M").vector_mult_add(Vsys.get_vector("MI1"),
				Vsys.get_vector("I1"));

		Vsys.get_vector("k1").zero();
		Vsys.get_vector("k1").close();
		Vsys.get_vector("k1").add(1.0, Vsys.get_vector("MVn"));
		Vsys.get_vector("k1").add(-0.5 * dt, Vsys.get_vector("KV1"));
		Vsys.get_vector("k1").add(-dt, Vsys.get_vector("MI1"));

		Vsys.rhs->zero();
		Vsys.rhs->close();
		Vsys.rhs->add(1.0, Vsys.get_vector("k1"));

		double tol = 1e-12;
		double max_iter = 2000;

		std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);

		Vsys.get_vector("V2").zero();
		Vsys.get_vector("V2").close();
		std::unique_ptr < libMesh::LinearSolver<libMesh::Number>
				> linearSolver;
		linearSolver = libMesh::LinearSolver<libMesh::Number>::build(es.comm());
		linearSolver->init();

		rval = linearSolver->solve(*Vsys.matrix, nullptr, Vsys.get_vector("V2"),
				*Vsys.rhs, tol, max_iter);
		// SECOND STEP

		Vsys.get_vector("I2").zero();
		Vsys.get_vector("I2").close();
		node = mesh.local_nodes_begin();
		for (; node != end_node; ++node)
		{
			const libMesh::Node * nn = *node;
			dof_map.dof_indices(nn, dof_indices, 0);
			double V2 = Vsys.get_vector("V2")(dof_indices[0]);
			double I2 = 0.0;
			if (0 == model)
				I2 = V2 * (V2 - alpha) * (V2 - 1.0);
			else
			{
				if (V2 > alpha)
					I2 = V2;
				if (V2 > alpha)
					I2 -= 1.0;
			}
			double dI2 = 0.0;
			double Q2 = 0;
			computeIonicCurrent(model, V2, Q2, h, alpha, I2, dI2,
					is_hyperbolic);
			Vsys.get_vector("I2").set(dof_indices[0], I2);
		}

		Vsys.get_vector("I2").close();
		Vsys.get_vector("MI2").zero();
		Vsys.get_vector("MI2").close();
		Vsys.get_matrix("M").vector_mult_add(Vsys.get_vector("MI2"),
				Vsys.get_vector("I2"));

		Vsys.get_vector("KV2").zero();
		Vsys.get_vector("KV2").close();
		Vsys.get_matrix("K").vector_mult_add(Vsys.get_vector("KV2"),
				Vsys.get_vector("V2"));

		Vsys.get_vector("k2").zero();
		Vsys.get_vector("k2").close();
		Vsys.get_vector("k2").add(1.0, Vsys.get_vector("MVn"));
		Vsys.get_vector("k2").add(-0.5 * dt, Vsys.get_vector("KV2"));
		Vsys.get_vector("k2").add(-0.5 * dt, Vsys.get_vector("KV1"));
		Vsys.get_vector("k2").add(-0.5 * dt, Vsys.get_vector("MI1"));
		Vsys.get_vector("k2").add(-0.5 * dt, Vsys.get_vector("MI2"));

		Vsys.rhs->zero();
		Vsys.rhs->close();
		Vsys.rhs->add(1., Vsys.get_vector("k2"));

		std::unique_ptr < libMesh::LinearSolver<libMesh::Number>
				> linearSolver2;
		linearSolver2 = libMesh::LinearSolver<libMesh::Number>::build(
				es.comm());
		linearSolver2->init();
		rval = linearSolver2->solve(Vsys.get_matrix("ML"), nullptr,
				*Vsys.current_local_solution, *Vsys.rhs, tol, max_iter);

		*Vsys.solution = *Vsys.current_local_solution;
		Vsys.update();
		*Vsys.old_local_solution = *Vsys.solution;
	}
}

double exact_solution(double mu, double alpha, double x, double time, int model)
{
	double exact_solution;
	double z0 = 0;
	if (model == 0)
	{
		double c = 0.5 * std::sqrt(2) * (1 - 2 * alpha);
		double z = x - c * time;
		exact_solution = 0.5 * (1.0 - std::tanh(0.5 * z / std::sqrt(2)));
	}
	else
	{
		double c = (1 - 2 * alpha)
				/ std::sqrt(mu + (alpha - alpha * alpha) * (mu - 1) * (mu - 1));
		double gamma = c * c * mu - 1.0;
		double beta = -c * (1 + mu);
		double Delta = beta * beta - 4 * gamma;
		double lp = (-beta + std::sqrt(Delta)) / (2 * gamma);
		double lm = (-beta - std::sqrt(Delta)) / (2 * gamma);
		double z = x - c * time;

		if (z >= z0)
		{
			exact_solution = alpha * std::exp(lp * z);
		}
		else
		{
			exact_solution = 1 + (alpha - 1) * std::exp(lm * z);
		}
	}
	return exact_solution;

}

double exact_derivative(double mu, double alpha, double x, double time,
		int model)
{
	double exact_der;
	double z0 = 0;
	if (model == 0)
	{
		double c = 0.5 * std::sqrt(2) * (1 - 2 * alpha);
		double z = x - c * time;
		double th = std::tanh(0.5 * z / std::sqrt(2));
		exact_der = c * 0.25 / std::sqrt(2.0) * (1.0 - th * th);
	}
	else
	{
		double c = (1 - 2 * alpha)
				/ std::sqrt(mu + (alpha - alpha * alpha) * (mu - 1) * (mu - 1));
		double gamma = c * c * mu - 1.0;
		double beta = -c * (1 + mu);
		double Delta = beta * beta - 4 * gamma;
		double lp = (-beta + std::sqrt(Delta)) / (2 * gamma);
		double lm = (-beta - std::sqrt(Delta)) / (2 * gamma);

		double z = x - c * time;

		if (z >= z0)
		{
			exact_der = -c * lp * alpha * std::exp(lp * z);
		}
		else
		{
			exact_der = -c * lm * (alpha - 1) * std::exp(lm * z);
		}
	}
	return exact_der;

}

void set_exact_ic(double mu, double alpha, double time,
		libMesh::EquationSystems & es, int model)
{
	const libMesh::MeshBase& mesh = es.get_mesh();
	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	typedef libMesh::TransientLinearImplicitSystem MonodomainSystem;
	typedef libMesh::TransientExplicitSystem IonicModelSystem;
	typedef libMesh::ExplicitSystem ParameterSystem;

	MonodomainSystem& Vsys = es.get_system<MonodomainSystem>("VSystem");     //V
	MonodomainSystem& Qsys = es.get_system<MonodomainSystem>("QSystem");     //V

	const libMesh::DofMap& dof_map = Vsys.get_dof_map();

	std::vector<libMesh::dof_id_type> dof_indices;

	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map.dof_indices(nn, dof_indices, 0);
		libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));

		double exsol = exact_solution(mu, alpha, p(0), time, model);
		Vsys.solution->set(dof_indices[0], exsol);        //V

		double exder = exact_derivative(mu, alpha, p(0), time, model);
		Qsys.solution->set(dof_indices[0], exder);        //Q
	}

	Vsys.solution->close();
	Vsys.update();
	Qsys.solution->close();
	Qsys.update();
}

void set_exact_solution(double mu, double alpha, double time,
		libMesh::EquationSystems & es, int model)
{
	const libMesh::MeshBase& mesh = es.get_mesh();
	libMesh::MeshBase::const_node_iterator node = mesh.local_nodes_begin();
	const libMesh::MeshBase::const_node_iterator end_node =
			mesh.local_nodes_end();

	typedef libMesh::TransientLinearImplicitSystem MonodomainSystem;
	typedef libMesh::TransientExplicitSystem IonicModelSystem;
	typedef libMesh::ExplicitSystem ParameterSystem;

	MonodomainSystem& Vsys = es.get_system<MonodomainSystem>("ExactVSystem"); //V
	MonodomainSystem& Qsys = es.get_system<MonodomainSystem>("ExactQSystem"); //V

	const libMesh::DofMap& dof_map = Vsys.get_dof_map();

	std::vector<libMesh::dof_id_type> dof_indices;

	for (; node != end_node; ++node)
	{
		const libMesh::Node * nn = *node;
		dof_map.dof_indices(nn, dof_indices, 0);
		libMesh::Point p((*nn)(0), (*nn)(1), (*nn)(2));

		double exsol = exact_solution(mu, alpha, p(0), time, model);
		Vsys.solution->set(dof_indices[0], exsol);        //V

		double exder = exact_derivative(mu, alpha, p(0), time, model);
		Qsys.solution->set(dof_indices[0], exder);        //Q
	}

	Vsys.solution->close();
	Vsys.update();
	Qsys.solution->close();
	Qsys.update();
}
void eval_l2_err(double& l2_V_err, double& l2_Q_err, double dt, double mu,
		double alpha, double time, libMesh::EquationSystems & es, int model)
{

	typedef libMesh::TransientLinearImplicitSystem MonodomainSystem;
	typedef libMesh::TransientExplicitSystem IonicModelSystem;
	typedef libMesh::ExplicitSystem ParameterSystem;

	const libMesh::MeshBase & mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	MonodomainSystem& Vsys = es.get_system<MonodomainSystem>("VSystem");     //V
	Vsys.update();
	MonodomainSystem& Qsys = es.get_system<MonodomainSystem>("QSystem");     //Q
	Qsys.update();
	const libMesh::DofMap & dof_map_monodomain = Vsys.get_dof_map();
	libMesh::FEType fe_type = dof_map_monodomain.variable_type(0);
	std::unique_ptr<libMesh::FEBase> fe(
			libMesh::FEBase::build(dim, fe_type));
	libMesh::QGauss qrule(dim, libMesh::FOURTH);
	fe->attach_quadrature_rule(&qrule);
	const std::vector<libMesh::Real> & JxW = fe->get_JxW();
	const std::vector<libMesh::Point> & q_point = fe->get_xyz();
	const std::vector<std::vector<libMesh::Real> > & phi = fe->get_phi();

	std::vector<libMesh::dof_id_type> dof_indices;

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();

	double err_l2_V2 = 0.0;
	double err_l2_Q2 = 0.0;

	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		const unsigned int elem_id = elem->id();
		dof_map_monodomain.dof_indices(elem, dof_indices);
		fe->reinit(elem);

		for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
		{

			libMesh::Number Vn = 0;
			libMesh::Number Qn = 0;
			for (unsigned int i = 0; i < phi.size(); ++i)
			{
				Qn += phi[i][qp]
						* (*Qsys.current_local_solution)(dof_indices[i]); //Q
				Vn += phi[i][qp]
						* (*Vsys.current_local_solution)(dof_indices[i]); //V
			}

			double x = q_point[qp](0);
			double exact_V = exact_solution(mu, alpha, x, time, model);
			double exact_Q = exact_derivative(mu, alpha, x, time, model);

			err_l2_V2 += JxW[qp] * (exact_V - Vn) * (exact_V - Vn);
			err_l2_Q2 += JxW[qp] * (exact_Q - Qn) * (exact_Q - Qn);

		}
	}

	double err_l2_V = std::sqrt(err_l2_V2);
	double err_l2_Q = 0.0;
	err_l2_Q = std::sqrt(err_l2_Q2);

	l2_V_err += dt * err_l2_V;
	l2_Q_err += dt * err_l2_Q;
}

void computeIonicCurrent(int model, double V, double Q, double h, double alpha,
		double& I, double& dI, bool is_hyperbolic)
{
//is_hyperbolic = true;
// Cubic
	if (0 == model)
	{
		I = V * (V - alpha) * (V - 1.0);
		dI = (V - alpha) * (V - 1.0) + V * (V - 1.0) + V * (V - alpha);
		dI *= Q;
	}
// Piecewise Linear (McKean)
	else
	{
		double H = 0;
		double dH = 0;
		if (V >= alpha - h && V < alpha)
		{
			// Linear Approximation
			//H = alpha * (V2 / h + 1);
			//dH = alpha / h;
			// Quadratic Approximation
			double A = alpha / h / h;
			double B = 2 * alpha * (h - alpha) / h / h;
			double C = (alpha * (alpha * alpha - 2 * alpha * h + h * h)) / h
					/ h;
			H = A * V * V + B * V + C;
			dH = 2 * A * V + B;
			// No approximation
			if (is_hyperbolic)
			{
				H = 0;
				dH = 0;
			}
		}
		if (V >= alpha && V < alpha + h)
		{
			// Linear Approximation
			//H = ( 1 - alpha ) / h * V2 + alpha;
			//dH = ( 1 - alpha ) / h;
			// Quadratic Approximation
			double A = (alpha - 1) / h / h;
			double B = 2 * (1 - alpha) * (h + alpha) / h / h;
			double C = (alpha
					* (alpha * alpha + 2 * alpha * h - alpha + h * h - 2 * h))
					/ h / h;
			H = A * V * V + B * V + C;
			dH = 2 * A * V + B;
			// No approximation
			if (is_hyperbolic)
			{
				H = 1;
				dH = 0;
			}
		}
		if (V >= alpha + h)
			H = 1;
		I = V - H;
		dI = (1 - dH) * Q;
	}
}
