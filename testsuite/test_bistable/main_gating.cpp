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

/*
 * main_gating.cpp
 *
 *  Created on: Aug 11, 2019
 *      Author: srossi
 */

#include "libmesh/getpot.h"
#include "libmesh/utility.h"
#include "libmesh/string_to_enum.h"

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
#include "libmesh/quadrature_trap.h"

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

#include "Util/SpiritFunction.hpp"

using namespace libMesh;



struct Params
{
	double v0;
	double v1;
	double v2;
	double k;
	double Cm;
	double Chi;
	double dt;
	double sigma_l;
	double sigma_t;

	QBase* qrule;
	Order qrule_order;

	bool lumped;

        void print_info()
        {
             std::cout << "V0: " << v0 << ", V1: " << v1 << ", V2: " << v2 << std::endl;
             std::cout << "k : " << k  << ", Cm: " << Cm << ", Chi: " << Chi << std::endl;
             std::cout << "dt: " << dt << ", sigma_l: " << sigma_l << ", sigma_t: " << sigma_t << std::endl;
             std::cout << "lumped: " << lumped << ", qrule_order: " << qrule_order << std::endl;
             qrule->print_info();
        }
};

void solve(	EquationSystems& es,
		Params& param,
			int t_step,
			int time_integrator_order = 1);

int main(int argc, char ** argv)
{
	// Initialize libMesh
	LibMeshInit init(argc, argv, MPI_COMM_WORLD);

	// Create input file:
	GetPot commandLine(argc, argv);
	std::string datafile_name = commandLine.follow("data.beat", 2, "-i",
			"--input");
	GetPot data(datafile_name);

	// Create Mesh:
	std::cout << "Creating mesh" << std::endl;

	libMesh::Mesh mesh(init.comm());
	int numElementsX = data("elX", 10);
	int numElementsY = data("elY", 0);
	int numElementsZ = data("elZ", 0);
	double minX = data("minX", 0.0);
	double maxX = data("maxX", 1.0);
	double minY = data("minY", 0.0);
	double maxY = data("maxY", 1.0);
	double minZ = data("minZ", 0.0);
	double maxZ = data("maxZ", 1.0);

	std::string elType = data("elType", "EDGE2");
        ElemType elt = libMesh::Utility::string_to_enum<ElemType>(elType);

	MeshTools::Generation::build_cube(mesh,
			                          numElementsX, numElementsY, numElementsZ,
									  minX, maxX,
									  minY, maxY,
									  minZ, maxZ, elt);

	// Print information about the mesh to the screen.
	mesh.print_info();
	const unsigned int dim = mesh.mesh_dimension();


	// Time integrator
    double final_time = data("TF", 1.0);
    double dt = data("dt", 0.125);
    std::string time_integrator = data("time_integrator", "RL");
    int time_integrator_order = data("time_integrator_order", 1);


	// Create system:
	std::cout << "Creating system" << std::endl;

	libMesh::EquationSystems systems(mesh);
	libMesh::TransientLinearImplicitSystem& monodomain = systems.add_system<
			libMesh::TransientLinearImplicitSystem>("V");

    // read order of the basis functions
    std::string porder = data("order", "FIRST");
    Order order = Utility::string_to_enum<Order>(porder);
    // read order of finite element family (for DG)
    std::string fefamily = data("fefamily", "LAGRANGE");
    FEFamily family = libMesh::Utility::string_to_enum<libMesh::FEFamily>(fefamily);

	unsigned int V_var = monodomain.add_variable("V", order, family);
	monodomain.add_vector("Iold");

	monodomain.init();
	monodomain.time  = 0.0;

	monodomain.print_info();

	// set initial condition
	std::cout << "Setting initial condition" << std::endl;
    std::string v_ic = data("ic", "");
    if(v_ic != "")
    {
        BeatIt::SpiritFunction monodomain_ic;
        monodomain_ic.read(v_ic);
        monodomain.project_solution(&monodomain_ic);
    }

	// output
	std::cout << "Creating output" << std::endl;

	std::string exodus_filename = data("exofile", "monodomain.exo");
	libMesh::ExodusII_IO exo(mesh);
	exo.write_equation_systems(exodus_filename, systems);
	exo.append(true);

	// Setup Parameters
	Params param;
	param.v0 = data("v0", -85.0);
	param.v1 = data("v1", -57.6);
	param.v2 = data("v2", 30.0);
	param.k = data("k", 1.4e-3);
	param.Cm = data("Cm", 1.0);
	param.Chi = data("Chi", 1.4e3);
	param.sigma_l = data("sigma_l", 1.0);
	param.sigma_t = data("sigma_t", 1.0);
	param.lumped = data("lumped", true);
	param.dt = dt;

    Order qrule_order = Utility::string_to_enum<Order>(data("qrule_order","SECOND"));
	param.qrule_order = qrule_order;
	std::string qrule_type = data("qrule_type", "QGauss");
	if("QGauss" == qrule_type) 
        {
            param.qrule = new QGauss(dim, qrule_order);
        }
        else if ( "" == qrule_type) 
        {
            param.qrule = new QTrap(dim, qrule_order);
        }

        param.qrule->init(elt);
         // TIME LOOP
	std::cout << "Looping over time" << std::endl;

        param.print_info(); 
	unsigned int t_step = 1;

	for (; monodomain.time < final_time;)
	{
		t_step++;
		monodomain.time += dt;

		monodomain.update();
		*monodomain.older_local_solution = *monodomain.old_local_solution;
		*monodomain.old_local_solution = *monodomain.current_local_solution;


		solve(systems, param, t_step);

		exo.write_timestep(exodus_filename, systems, t_step, monodomain.time);


	}

	delete param.qrule;
	return 0;
}


void solve(	EquationSystems& es,
		Params& param,
			int t_step,
			int time_integrator_order)
{
	using libMesh::UniquePtr;

	bool assemble_matrix = true;
	if(t_step > 2) assemble_matrix = false;



	const libMesh::MeshBase & mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	libMesh::TransientLinearImplicitSystem & Vsys = es.get_system<
			libMesh::TransientLinearImplicitSystem>("V");
	const unsigned int V_Var = Vsys.variable_number("V");

	Vsys.rhs->zero();


	const libMesh::DofMap & dof_map = Vsys.get_dof_map();
	libMesh::FEType fe_type = dof_map.variable_type(0);

	double dt = param.dt;
	double time = Vsys.time;

	UniquePtr<libMesh::FEBase> fe_matrix(libMesh::FEBase::build(dim, fe_type));
	libMesh::QGauss qrule(dim, SECOND);
	fe_matrix->attach_quadrature_rule(&qrule);

	const std::vector<libMesh::Real> & JxW = fe_matrix->get_JxW();
	const std::vector<libMesh::Point> & q_point = fe_matrix->get_xyz();
	const std::vector<std::vector<libMesh::Real> > & phi = fe_matrix->get_phi();
	const std::vector<std::vector<libMesh::RealGradient> > & dphi =
			fe_matrix->get_dphi();


        TensorValue<double> sigma;
        sigma(0,0) = param.sigma_l/param.Chi;
        sigma(1,1) = param.sigma_t/param.Chi;
        sigma(2,2) = param.sigma_t/param.Chi;
	if(assemble_matrix)
	{
                Vsys.matrix->zero();
		libMesh::DenseMatrix<libMesh::Number> Ke;

		std::vector<libMesh::dof_id_type> dof_indices;

		libMesh::MeshBase::const_element_iterator el =
				mesh.active_local_elements_begin();
		const libMesh::MeshBase::const_element_iterator end_el =
				mesh.active_local_elements_end();


		for (; el != end_el; ++el)
		{
			const libMesh::Elem * elem = *el;
			dof_map.dof_indices(elem, dof_indices);
			fe_matrix->reinit(elem);
			Ke.resize(dof_indices.size(), dof_indices.size());

			for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
			{
				for (unsigned int i = 0; i < phi.size(); i++)
				{
					for (unsigned int j = 0; j < phi.size(); j++)
					{
						// Mass term
						double mass = param.Cm * JxW[qp] / dt * (phi[i][qp] * phi[j][qp]);
						if(param.lumped)
							Ke(i, i) += mass;
						else
							Ke(i, j) += mass;
						// stiffness term
						Ke(i, j) += JxW[qp] * ( dphi[i][qp] * (sigma * dphi[j][qp] ) );
					}
				}
			}

			Vsys.matrix->add_matrix(Ke, dof_indices, dof_indices);
		}
		Vsys.matrix->close();
	}

	UniquePtr<libMesh::FEBase> fe_current(libMesh::FEBase::build(dim, fe_type));
	fe_current->attach_quadrature_rule(param.qrule);

	const std::vector<libMesh::Real> & JxW_i = fe_current->get_JxW();
	const std::vector<libMesh::Point> & q_point_i = fe_current->get_xyz();
	const std::vector<std::vector<libMesh::Real> > & phi_i = fe_current->get_phi();

	// RHS
	libMesh::DenseVector<libMesh::Number> Fe;
	std::vector<libMesh::dof_id_type> dof_indices;

	libMesh::MeshBase::const_element_iterator el =
			mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el =
			mesh.active_local_elements_end();


	for (; el != end_el; ++el)
	{
		const libMesh::Elem * elem = *el;
		dof_map.dof_indices(elem, dof_indices);
		fe_matrix->reinit(elem);
		fe_current->reinit(elem);
		Fe.resize(dof_indices.size());

		// VN
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
		{
			for (unsigned int i = 0; i < phi.size(); i++)
			{
				for (unsigned int j = 0; j < phi.size(); j++)
				{
					if(param.lumped)
						Fe(i) +=  param.Cm * JxW[qp] / dt * phi[j][qp] * phi[i][qp] * (*Vsys.old_local_solution)(dof_indices[i]);
					else
						Fe(i) +=  param.Cm * JxW[qp] / dt * phi[j][qp] * phi[i][qp] * (*Vsys.old_local_solution)(dof_indices[j]);
				}
			}
		}
		// VN
		for (unsigned int qp = 0; qp < param.qrule->n_points(); qp++)
		{
			double V_old = 0.0;
			for (unsigned int j = 0; j < phi_i.size(); j++)
			{
				V_old += phi_i[j][qp] * (*Vsys.old_local_solution)(dof_indices[j]);
			}
			double I_ion = param.k * (V_old - param.v0) * (V_old - param.v1) * (V_old - param.v2);
			for (unsigned int i = 0; i < phi_i.size(); i++)
			{
				Fe(i) -= JxW_i[qp] * phi_i[i][qp] * I_ion;
			}
		}
		Vsys.rhs->add_vector(Fe,dof_indices);

	}

	Vsys.rhs->close();


	double tol = 1e-12;
	double max_iter = 2000;

	std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);
	libMesh::UniquePtr < libMesh::LinearSolver<libMesh::Number>
			> linearSolver;
	linearSolver = libMesh::LinearSolver<libMesh::Number>::build(es.comm());
	linearSolver->init();

	rval = linearSolver->solve(*Vsys.matrix, nullptr, *Vsys.solution,
			*Vsys.rhs, tol, max_iter);

}


