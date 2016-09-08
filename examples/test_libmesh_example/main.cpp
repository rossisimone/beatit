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

// Functions to initialize the library.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/numeric_vector.h"

int main (int argc, char ** argv)
{
    using namespace libMesh;
	LibMeshInit init (argc, argv);
	Mesh mesh1(init.comm());
    Mesh mesh2(init.comm());
	std::cout << "Creating Mesh ..." << std::endl;
	MeshTools::Generation::build_cube (mesh1,
                                       2 , 2 , 2,
                                       0., 1.,
                                       0., 1.,
                                       0., 1.,
                                       TET4);
    std::cout << "Creating Mesh ..." << std::endl;
    MeshTools::Generation::build_cube (mesh2,
                                       2 , 2 , 2,
                                       0., 1.,
                                       0., 1.,
                                       0., 1.,
                                       TET4);

    std::cout << "Creating System 1: " << std::endl;
    EquationSystems systems1(mesh1);
    ExplicitSystem& s1 = systems1.add_system<ExplicitSystem>("s1");
    s1.add_variable("s1x", FIRST);
//    s1.add_variable("s1y", FIRST);
    std::cout << "Initializing s1" << std::endl;
    systems1.init();
    std::cout << "Creating System 2: " << std::endl;
    EquationSystems systems2(mesh1);
    ExplicitSystem& s2 = systems2.add_system<ExplicitSystem>("s2");
    s2.add_variable("s2a", FIRST);
    s2.add_variable("s2b", FIRST); // Commenting this line below it runs ok
    std::cout << "Initializing s2" << std::endl;
    systems2.init();
    std::cout << "Exporting s1" << std::endl;
    std::set<std::string> system_names;
    system_names.insert("s1");
    std::vector<std::string> vars_names;
    vars_names.push_back("s1x");
    ExodusII_IO(systems1.get_mesh()).write_nodal_data("s1.e", *s1.solution, vars_names);
    std::cout << "Exporting s2" << std::endl;
    std::vector<std::string> vars_names2;
    vars_names2.push_back("s2a");
    vars_names2.push_back("s2b");
    ExodusII_IO(systems1.get_mesh()).write_nodal_data("s2.e", *s2.solution, vars_names2);
//    ExodusII_IO(systems2.get_mesh()).write_equation_systems("s2.e", systems2);
//    std::cout << "Exporting s3" << std::endl;
//    ExodusII_IO(systems1.get_mesh()).write_equation_systems("s3.e", systems1);

	return 0;
}






