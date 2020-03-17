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
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"

#include "Util/IO/io.hpp"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"

#include <random>

int main(int argc, char ** argv)
{
	// Import input file
	GetPot data = BeatIt::readInputFile(argc, argv);

	using namespace libMesh;
	LibMeshInit init(argc, argv, MPI_COMM_WORLD);

	// Create/Read mesh
	Mesh mesh(init.comm());

	std::string mesh_name = data("mesh", "NONE");

	if ("NONE" != mesh_name)
	{
		mesh.read(&mesh_name[0]);
	  int num_refs = data("refs", 0);
	  std::cout <<  "refinments: " << num_refs << std::endl;
	  MeshRefinement(mesh).uniformly_refine(num_refs);
	  mesh.prepare_for_use();

		bool set_regions = data("set_regions", false);
		if(set_regions)
		{}
//			double x1 = data("x1", 0.0416849669806048);
//		double z1 = data("z1", -0.0338434083683408);
//		double x2 = data("x2",0.0670653739716849);
//		double z2 = data("z2",-0.0447536820553136);
//
//		double m1 = (z1 - z2) / (x1 - x2);
//		double q1 = z1 - m1 * z1;
//
//		double x3 = data("x3", 0.00869961251998251);
//		double z3 = data("z3", -0.0560889956394845);
//		double x4 = data("x4",0.0488975409960176);
//		double z4 = data("z4",-0.022876343139214);
//		double m2 = (z3 - z4) / (x3 - x4);
//		double q2 = z3 - m2 * x3;
//
//		MeshBase::element_iterator el = mesh.elements_begin();
//		const MeshBase::element_iterator end_el = mesh.elements_end();
//		auto gen = std::bind(std::uniform_int_distribution<>(0, 1),
//				std::default_random_engine());
//
//		double center = 0.0;
//		for (; el != end_el; ++el)
//		{
//			Elem * elem = *el;
//			const Point cent = elem->centroid();
//			double x = cent(0);
//			double y = cent(2);
//
//			if ( y >= m2 * x + q2)
//			{
//				elem->subdomain_id() = 1;
//			}
//			else
//				elem->subdomain_id() = 2;
//
//			if (y >= m1 * x + q1 && 2 == elem->subdomain_id() )
//			{
//				elem->subdomain_id() = 3;
//			}
//
//			if (y >= m2 * x + q2 && 1 == elem->subdomain_id() )
//			{
//				for (unsigned int side = 0; side < elem->n_sides(); side++)
//				{
//					if (elem->neighbor_ptr(side) == libmesh_nullptr)
//					{
//						elem->subdomain_id() = 4;
//						break;
//					}
//				}
//				if (4 == elem->subdomain_id())
//				{
//					for (unsigned int side = 0; side < elem->n_sides(); side++)
//					{
//						if (elem->neighbor_ptr(side) != libmesh_nullptr)
//						{
//							elem->neighbor_ptr(side)->subdomain_id() = 4;
//						}
//					}
//				}
//			}
//		}
//		el = mesh.elements_begin();
//		for (; el != end_el; ++el)
//		{
//			Elem * elem = *el;
//			if (4 == elem->subdomain_id())
//			{
//				for (unsigned int side = 0; side < elem->n_sides();
//						side++)
//				{
//					if (elem->neighbor_ptr(side) != libmesh_nullptr)
//					{
//						elem->neighbor_ptr(side)->subdomain_id() = 5;
//					}
//				}
//			}
//		}
//		el = mesh.elements_begin();
//		for (; el != end_el; ++el)
//		{
//			Elem * elem = *el;
//			if (5 == elem->subdomain_id())
//			{
//				elem->subdomain_id() = 4;
//			}
//		}
	}
	else
	{
		int nelx = data("nelx", 10);
		int nely = data("nely", 10);
		int nelz = data("nelz", 10);
		int maxx = data("maxx", 1.0);
		int maxy = data("maxy", 1.0);
		int maxz = data("maxz", 1.0);
		int minx = data("minx", 1.0);
		int miny = data("miny", 1.0);
		int minz = data("minz", 1.0);

		auto elType = TET4;
		if (nelz == 0)
			elType = TRI3;
		MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, minx, maxx,
				miny, maxy, minz, maxz, elType);
		MeshBase::element_iterator el = mesh.elements_begin();
		const MeshBase::element_iterator end_el = mesh.elements_end();

		double center = 0.0;
		for (; el != end_el; ++el)
		{
			Elem * elem = *el;
			const Point cent = elem->centroid();
			if (cent(0) <= center && cent(1) <= center && cent(2) <= center)
				elem->subdomain_id() = 1;
			else if (cent(0) > center && cent(1) <= center && cent(2) <= center)
				elem->subdomain_id() = 1;
			else if (cent(0) > center && cent(1) > center && cent(2) <= center)
				elem->subdomain_id() = 5;
			else if (cent(0) <= center && cent(1) > center && cent(2) <= center)
				elem->subdomain_id() = 5;
			else if (cent(0) <= center && cent(1) <= center && cent(2) > center)
				elem->subdomain_id() = 1;
			else if (cent(0) > center && cent(1) <= center && cent(2) > center)
				elem->subdomain_id() = 1;
			else if (cent(0) > center && cent(1) > center && cent(2) > center)
				elem->subdomain_id() = 5;
			else if (cent(0) <= center && cent(1) > center && cent(2) > center)
				elem->subdomain_id() = 5;
		}
		double scale = data("scale", 1.0);
		MeshTools::Modification::scale(mesh, scale, scale, scale);
	}

	mesh.print_info();

	BeatIt::TimeData datatime;
	datatime.setup(data, "");
	datatime.print();

	// Create Equations systems;
	libMesh::EquationSystems es(mesh);
	// TODO: import mesh with data

	// Constructor
	int save_iter = 0;
	int save_iter_ve = 1;
	std::string model = data("model", "monowave");
	std::string system_mass = data(model + "/diffusion_mass", "mass");
	std::string iion_mass = data(model + "/reaction_mass", "lumped_mass");
	bool useMidpointMethod = false;
	int step0 = 0;
	int step1 = 1;

	std::cout << "Create bidomain with bath..." << std::endl;
	BeatIt::ElectroSolver* solver =
			BeatIt::ElectroSolver::ElectroFactory::Create(model, es);
	std::string section = data("section", "monowave");
	std::cout << "Calling setup..." << std::endl;
	solver->setup(data, model);
	std::cout << "Calling init ..." << std::endl;
	es.print_info();
	solver->init(0.0);
	solver->save_parameters();
	std::cout << "Assembling matrices" << std::endl;
	solver->assemble_matrices(datatime.M_dt);
	solver->save_potential(save_iter, 0.0);

	double threshold = data("threshold", -10.0);
    int at_save_iter = data("at_save_iter", 25);

	std::cout << "Time loop starts:" << std::endl;
	for (;
			datatime.M_iter < datatime.M_maxIter
					&& datatime.M_time < datatime.M_endTime;)
	{
		datatime.advance();
		std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter
				<< std::endl;
		solver->advance();
		solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0,
				useMidpointMethod, iion_mass);

		solver->solve_diffusion_step(datatime.M_dt, datatime.M_time,
				useMidpointMethod, iion_mass);

		solver->update_activation_time(datatime.M_time, threshold);

		if (0 == datatime.M_iter % datatime.M_saveIter)
		{
			save_iter++;
			solver->save_potential(save_iter, datatime.M_time);
			//solver->save_exo_timestep(save_iter, datatime.M_time);
		}

        if (0 == datatime.M_iter % (at_save_iter * datatime.M_saveIter ) )
        {
            solver->save_activation_times(save_iter);
        }


	}

	solver->evaluate_conduction_velocity();
	solver->save_conduction_velocity(save_iter);
	solver->save_activation_times(save_iter);

//	// eval cvs
//	{
//		typedef libMesh::ExplicitSystem  ExplicitSystem;
//		ExplicitSystem& cv = es.add_system<ExplicitSystem>("CV");
//		cv.add_variable("cv", MONOMIAL, CONSTANT);
//		cv.init();
//
//		MeshBase::element_iterator el = mesh.elements_begin();
//		const MeshBase::element_iterator end_el = mesh.elements_end();
//		auto gen = std::bind(std::uniform_int_distribution<>(0, 1),
//				std::default_random_engine());
//
//		double center = 0.0;
//		for (; el != end_el; ++el)
//		{
//			const libMesh::Elem * elem = *el;
//			const libMesh::Elem * side0 = elem->side(0);
//			const libMesh::Elem * side1 = elem->side(1);
//			const libMesh::Elem * side2 = elem->side(2);
//
//			double a = side0->volume();
//			double b = side1->volume();
//			double c = side2->volume();
//
//			double theta = acos ( dx^)
//
//		}
//	}

	delete solver;

	//create fibers

	// Time loop
	return 0;
}
