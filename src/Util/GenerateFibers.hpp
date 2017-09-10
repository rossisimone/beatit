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
 * \file GenerateFibers.hpp
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

#ifndef SRC_UTIL_GENERATEFIBERS_HPP_
#define SRC_UTIL_GENERATEFIBERS_HPP_

#include <string>

class GetPot;

namespace libMesh
{
class MeshBase;
class System;
}
namespace BeatIt
{

namespace Util
{

const libMesh::UniquePtr< libMesh::NumericVector<libMesh::Number> >&
generate_gradient_field( libMesh::MeshBase& mesh,
                         const GetPot& data,
                         const std::string& section = "rule_based_fibers");


void project_function(std::string& function, libMesh::System& sys);

void normalize(double& x, double& y, double& z,
               double X, double Y, double Z);


template< class Vector >
void normalize( Vector& vec, double X = 1.0, double Y = 0.0, double Z = 0.0 )
{
    auto first = vec.first_local_index();
    auto last = vec.last_local_index();

    for(int i = first; i < last; )
    {
        int j = i;
        double v1_x = vec(i);
        i++;
        double v1_y = vec(i);
        i++;
        double v1_z = vec(i);
        i++;
        normalize(v1_x,v1_y,v1_z,X,Y,Z);
        vec.set(j,v1_x);
        j++;
        vec.set(j,v1_y);
        j++;
        vec.set(j,v1_z);
    }

}


}

}

#endif /* SRC_UTIL_GENERATEFIBERS_HPP_ */
