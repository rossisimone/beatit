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
 * \file GenerateFibers.cpp
 *
 * \class GenerateFibers
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Sep 15, 2016
 *
 */


#include "PoissonSolver/Poisson.hpp"
#include "libmesh/mesh.h"

#include "libmesh/numeric_vector.h"
#include "Util/GenerateFibers.hpp"
namespace BeatIt
{

namespace Util
{

const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
generate_gradient_field( libMesh::MeshBase& mesh,
                         const GetPot& data,
                         const std::string& section)
{
    libMesh::Mesh new_mesh( dynamic_cast< libMesh::Mesh&>(mesh) );
    libMesh::EquationSystems es(mesh);
    Poisson p(es);
    p.setup(data, section);
    p.assemble_system();
    p.solve_system();
    p.compute_elemental_solution_gradient();

    libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >
    ptr( libMesh::NumericVector<libMesh::Number>::build( mesh.comm() ) );
    const auto& grad_ptr = p.get_gradient();
    return grad_ptr;
}


void
project_function(std::string& function, libMesh::System& sys)
{
    SpiritFunction func;
    func.read(function );
    sys.project_solution(&func);
}



} /* namespace Util */


}/* namespace BeatIt*/
