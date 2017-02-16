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
 * MapsToLibMeshTypes.hpp
 *
 *  Created on: Oct 19, 2016
 *      Author: srossi
 */

#ifndef SRC_UTIL_MAPSTOLIBMESHTYPES_HPP_
#define SRC_UTIL_MAPSTOLIBMESHTYPES_HPP_

#include <string>
#include <map>
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

namespace BeatIt
{

static const std::map<int, libMesh::Order> libmesh_order_map =
{
		{ 0, libMesh::CONSTANT },
		{ 1, libMesh::FIRST },
		{ 2, libMesh::SECOND },
};


static const std::map<std::string, libMesh::FEFamily> libmesh_fefamily_map =
{
		{ "lagrange", libMesh::LAGRANGE },
		{ "monomial", libMesh::MONOMIAL },
        { "dg", libMesh::L2_LAGRANGE },
};

} // BeatIt namespace


#endif /* SRC_UTIL_MAPSTOLIBMESHTYPES_HPP_ */
