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
 * \file MonodomainUtil.hpp
 *
 *
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
 * Created on: Aug 15, 2016
 *
 */
#ifndef SRC_ELECTROPHYSIOLOGY_MONODOMAIN_MONODOMAINUTIL_HPP_
#define SRC_ELECTROPHYSIOLOGY_MONODOMAIN_MONODOMAINUTIL_HPP_

#include "libmesh/libmesh_common.h"

namespace libMesh
{
class System;
template <typename T> class NumericVector;

}

namespace BeatIt
{

namespace Electrophysiology
{


bool isFullyActivated(libMesh::System& system, double initiValue = -1.0);

void updateActivationTimes( const libMesh::NumericVector<libMesh::Number>& potential,
		                                                 libMesh::System& system,
														 double time,
														 double threshold = 0.8,
														 double initValue = -1.0);


} // Electrophysiology

} // BeatIt


#endif /* SRC_ELECTROPHYSIOLOGY_MONODOMAIN_MONODOMAINUTIL_HPP_ */
