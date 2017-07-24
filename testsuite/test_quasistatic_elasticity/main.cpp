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
 * main.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */
// Basic include files needed for the mesh functionality.
#include "Elasticity/MixedElasticity.hpp"

#include "libmesh/linear_implicit_system.h"

#include "libmesh/wrapped_functor.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "Util/SpiritFunction.hpp"

#include "libmesh/numeric_vector.h"
#include "libmesh/elem.h"
#include "Util/CTestUtil.hpp"
#include "libmesh/dof_map.h"
#include <iomanip>


enum Formulation { Primal,  Mixed, Incompressible };
enum CL { Linear,  NH, Hexagonal };

double reference_norm(Formulation f, CL cl);

int main (int argc, char ** argv)
{
    // Bring in everything from the libMesh namespace

    using namespace libMesh;
    // Initialize libraries, like in example 2.
    LibMeshInit init (argc, argv, MPI_COMM_WORLD);


    // Use the MeshTools::Generation mesh generator to create a uniform
    // 3D grid
    // We build a linear tetrahedral mesh (TET4) on  [0,2]x[0,0.7]x[0,0.3]
    // the number of elements on each side is read from the input file
    GetPot commandLine ( argc, argv );
    std::string datafile_name = commandLine.follow ( "data.pot", 2, "-i", "--input" );
    GetPot data(datafile_name);
    // allow us to use higher-order approximation.
    // Create a mesh, with dimension to be overridden later, on the
    // default MPI communicator.
    libMesh::Mesh mesh(init.comm());

    int elX = data("elX", 10);
    int elY = data("elY", 2);
    //std::string elTypeName = data("elType", "TRI6");
    //std::map<std::string, ElemType> orderMap;
    //orderMap["TRI6"] = TRI6;
    //auto elType = orderMap.find(elTypeName)->second;
    auto elType = TRI6;
    auto order = FIRST;
    if(elType == TRI6 || elType == QUAD9 ) order = SECOND;
    auto elType2 = TRI3;
    if( elType == QUAD9 ) elType2 = QUAD4;

    std::cout << "Element type chosen:  " << elType << std::endl;
    MeshTools::Generation::build_square ( mesh,
                                          elX, elY,
                                          0., 10.,
                                          0., 2.,
                                          elType );

//    std::string meshfile = data("mesh", "NOMESH");
//    mesh.read (&meshfile[0]);

    libMesh::EquationSystems es(mesh);

    Formulation f;
    std::string formulation = data("elasticity/formulation", "NO_FORMULATION");
    std::cout << "Creating " << formulation << std::endl;

    CL cl;
    std::string cl_name = data("elasticity/materials", "NO_MATERIAL");
    std::cout << "Material " << cl_name << std::endl;

    BeatIt::Elasticity* elas = nullptr;
    if(formulation == "primal")
    {
        elas = new BeatIt::Elasticity(es, "Elasticity");
        f = Primal;
        if(cl_name == "linear") cl = Linear;
        else cl = NH;

    }
    else
    {
        elas = new BeatIt::MixedElasticity(es, "Elasticity");
        double nu = 0.0;
        if(cl_name == "linear")
        {
            cl = Linear;
            nu = data("elasticity/materials/linear/nu", -1.0);
        }
        // Hexagonal coded only for the incompressible mixed formulation
        else if(cl_name == "hexagonal")
        {
            cl = Hexagonal;
            //
            nu = 0.5;
        }
        else
        {
            cl = NH;
            nu = data("elasticity/materials/"+cl_name+"/nu", -1.0);
        }
        if(nu < -0.5) return EXIT_FAILURE;
        else if(nu < 0.5) f = Mixed;
        else f = Incompressible;

    }
    std::cout << "Setting up ... " << std::endl;
    elas->setup(data,"elasticity");
    std::cout << "Initializing output ... " << std::endl;
    elas->init_exo_output("solution.exo");
    std::cout << "Solving ... " << std::endl;
    elas->newton();
    if(Hexagonal == cl) elas->evaluate_nodal_I4f();
    elas->save_exo("solution.exo", 1, 1.0);
    libMesh::LinearImplicitSystem& system  =  elas->M_equationSystems.get_system<libMesh::LinearImplicitSystem>(elas->M_myName);

    auto norm = system.solution->linfty_norm ();
    auto reference_value = reference_norm(f, cl);
    std::cout << std::setprecision(20) << "norm is: " << norm << std::endl;

    delete elas;
    return BeatIt::CTest::check_test(norm, reference_value, 1e-10);
}

double reference_norm(Formulation f, CL cl)
{
    double norm;
    switch(cl)
    {
        case Linear:
        {
            switch(f)
            {
                case Primal:
                {
                    norm = 0.00264652492963676;
                    break;
                }
                case Mixed:
                {
                    norm =  52703.24720521792915;
                    break;
                }
                case Incompressible:
                {
                    norm =  50729.94868658449559;
                    break;
                }

                default:
                {
                    std::cout << "Formulation not set correctly!" << std::endl;
                    std::runtime_error("Formulation not set correctly!");
                    break;
                }
            }
            break;
        }
        case NH:
        {
            switch(f)
            {
                case Primal:
                {
                    norm = 0.18980564561565821391;
                    break;
                }
                case Mixed:
                {
                    norm =  12.66136556159523785;
                    break;
                }
                case Incompressible:
                {
                    norm = 12.66534007214748847;
                    break;
                }

                default:
                {
                    std::cout << "Formulation not set correctly!" << std::endl;
                    std::runtime_error("Formulation not set correctly!");
                    break;
                }
            }
            break;
        }
        case Hexagonal:
        {
            switch(f)
            {
				case Incompressible:
				{
					norm = 717138.49841777875554;
					break;
				}
                default:
                {
                    std::cout << "Formulation not set correctly!" << std::endl;
                    std::runtime_error("Formulation not set correctly!");
                    break;
                }
            }
        	break;
        }
        default:
        {
            std::cout << "Constitutile law not set correctly!" << std::endl;
            std::runtime_error("CL not set correctly!");
            break;
        }
    }

    return norm;
}

